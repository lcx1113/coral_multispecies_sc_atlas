#!/usr/bin/env python3
"""Shared helpers for the formal scAtlas orthogroup gene-tree workflow."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
TRUE_VALUES = {"true", "t", "1", "yes", "y"}
FALSE_VALUES = {"false", "f", "0", "no", "n"}


def resolve_path(value: str | Path, project_root: Path = PROJECT_ROOT) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else project_root / path


def load_config(config_path: str | Path) -> tuple[dict, Path]:
    path = resolve_path(config_path)
    with path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError(f"Configuration is not a mapping: {path}")
    return config, path


def input_path(config: Mapping, key: str) -> Path:
    try:
        return resolve_path(config["inputs"][key])
    except KeyError as exc:
        raise KeyError(f"Missing inputs.{key}") from exc


def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def open_text(path: str | Path, mode: str = "rt"):
    """Open plain text or gzip; gzip writes are reproducible (mtime=0)."""
    path = Path(path)
    if path.suffix != ".gz":
        return path.open(mode.replace("t", ""), encoding="utf-8", newline="") if "t" in mode else path.open(mode)
    if "r" in mode:
        return gzip.open(path, mode, encoding="utf-8", newline="")
    ensure_dir(path.parent)
    raw_mode = mode.replace("t", "").replace("b", "") + "b"
    raw = path.open(raw_mode)
    gz = gzip.GzipFile(filename="", mode=raw_mode, fileobj=raw, mtime=0)
    return io.TextIOWrapper(gz, encoding="utf-8", newline="")


def read_tsv(path: str | Path) -> list[dict]:
    with open_text(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: str | Path, rows: Iterable[Mapping], fieldnames: Sequence[str]) -> None:
    path = Path(path)
    ensure_dir(path.parent)
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: normalize_scalar(row.get(key, "")) for key in fieldnames})


def write_json(path: str | Path, value: object) -> None:
    path = Path(path)
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")


def normalize_scalar(value: object) -> object:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    return value


def parse_bool(value: object, field: str = "value") -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in TRUE_VALUES:
        return True
    if normalized in FALSE_VALUES:
        return False
    raise ValueError(f"Invalid boolean for {field}: {value!r}")


def parse_fasta(path: str | Path) -> Iterator[tuple[str, str, str]]:
    """Yield first-token ID, full header, and sequence."""
    seq_id = None
    header = ""
    pieces: list[str] = []
    with open_text(path, "rt") as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, header, "".join(pieces)
                header = line[1:]
                if not header.strip():
                    raise ValueError(f"Empty FASTA header at {path}:{line_number}")
                seq_id = header.split()[0]
                pieces = []
            elif seq_id is None:
                if line.strip():
                    raise ValueError(f"Sequence before first header at {path}:{line_number}")
            else:
                pieces.append(line.strip())
    if seq_id is not None:
        yield seq_id, header, "".join(pieces)


def write_fasta(path: str | Path, records: Iterable[tuple[str, str]], width: int = 80) -> None:
    path = Path(path)
    ensure_dir(path.parent)
    with open_text(path, "wt") as handle:
        for record_id, sequence in records:
            handle.write(f">{record_id}\n")
            for start in range(0, len(sequence), width):
                handle.write(sequence[start : start + width] + "\n")


def fasta_dict(path: str | Path) -> dict[str, str]:
    records: dict[str, str] = {}
    for record_id, _header, sequence in parse_fasta(path):
        if record_id in records:
            raise ValueError(f"Duplicate FASTA ID {record_id!r} in {path}")
        records[record_id] = sequence
    return records


def alignment_stats(path: str | Path) -> dict:
    records = fasta_dict(path)
    lengths = {record_id: len(sequence) for record_id, sequence in records.items()}
    unique_lengths = set(lengths.values())
    if len(unique_lengths) > 1:
        raise ValueError(f"Alignment sequences have unequal lengths in {path}: {sorted(unique_lengths)[:10]}")
    alignment_length = next(iter(unique_lengths), 0)
    non_gap = {
        record_id: sum(character not in {"-", "."} for character in sequence)
        for record_id, sequence in records.items()
    }
    return {
        "n_sequences": len(records),
        "alignment_length": alignment_length,
        "minimum_non_gap": min(non_gap.values(), default=0),
        "maximum_non_gap": max(non_gap.values(), default=0),
        "mean_gap_fraction": (
            sum((alignment_length - count) / alignment_length for count in non_gap.values()) / len(non_gap)
            if records and alignment_length else 0.0
        ),
        "ids": set(records),
        "non_gap_by_id": non_gap,
    }


def run_command(
    command: Sequence[str], stdout_path: str | Path, stderr_path: str | Path,
    timeout: int | None = None,
) -> dict:
    """Run a command and persist its streams, timing, status, and timeout state."""
    stdout_path = Path(stdout_path)
    stderr_path = Path(stderr_path)
    ensure_dir(stdout_path.parent)
    ensure_dir(stderr_path.parent)
    started = utc_now()
    timed_out = False
    exception = ""
    exit_code = 1
    try:
        with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
            completed = subprocess.run(
                [str(value) for value in command], stdout=stdout, stderr=stderr,
                timeout=timeout, check=False, env=os.environ.copy(),
            )
        exit_code = completed.returncode
    except subprocess.TimeoutExpired as exc:
        timed_out = True
        exit_code = 124
        exception = repr(exc)
    except Exception as exc:  # preserve a command record even for launch failures
        exception = repr(exc)
    return {
        "command": [str(value) for value in command],
        "command_text": shlex.join([str(value) for value in command]),
        "start_utc": started,
        "end_utc": utc_now(),
        "exit_code": exit_code,
        "timed_out": timed_out,
        "exception": exception,
        "stdout": str(stdout_path.resolve()),
        "stderr": str(stderr_path.resolve()),
    }


def reusable_command_record(record: Mapping, expected_command: Sequence[str], input_sha256: str) -> bool:
    recorded = [str(value) for value in record.get("command", [])]
    expected = [str(value) for value in expected_command]
    return (
        record.get("exit_code") == 0
        and not record.get("timed_out", False)
        and record.get("input_sha256") == input_sha256
        and recorded == expected
    )


def archive_existing(paths: Iterable[str | Path], archive_root: str | Path) -> str:
    """Move existing paths into one timestamped archive and return its path."""
    present = [Path(path) for path in paths if Path(path).exists()]
    if not present:
        return ""
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    destination = Path(archive_root) / stamp
    counter = 1
    while destination.exists():
        destination = Path(archive_root) / f"{stamp}_{counter}"
        counter += 1
    ensure_dir(destination)
    for path in present:
        shutil.move(str(path), str(destination / path.name))
    return str(destination.resolve())


def transcript_from_protein_id(protein_id: str, species: str) -> str:
    """Frozen Broccoli rule inherited from the pilot."""
    prefix = f"{species}_"
    if not protein_id.startswith(prefix):
        raise ValueError(f"Protein ID {protein_id!r} lacks expected prefix {prefix!r}")
    return re.sub(r"\.p\d+$", "", protein_id[len(prefix) :])


def canonical_gene(raw_gene: str, species: str) -> str:
    fixed = raw_gene.replace("-", "_")
    prefix = f"{species}_"
    return fixed if fixed.startswith(prefix) else prefix + fixed


def parse_gtf_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for match in re.finditer(r'(\S+)\s+"([^"]*)"\s*;?', text):
        attributes.setdefault(match.group(1), match.group(2))
    return attributes


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def provenance_record(
    script: str | Path, config_path: str | Path, argv: Sequence[str], inputs: Sequence[str | Path], outputs: Sequence[str | Path]
) -> dict:
    def record(path_value: str | Path) -> dict:
        path = Path(path_value).resolve()
        return {
            "path": str(path),
            "exists": path.is_file(),
            "sha256": sha256_file(path) if path.is_file() else "",
            "size_bytes": path.stat().st_size if path.is_file() else 0,
        }

    return {
        "created_utc": utc_now(),
        "script": record(script),
        "config": record(config_path),
        "argv": list(argv),
        "inputs": [record(value) for value in inputs],
        "outputs": [record(value) for value in outputs],
    }


def species_manifest(config: Mapping, included_only: bool = False) -> list[dict]:
    rows = read_tsv(input_path(config, "species_inputs"))
    if included_only:
        rows = [row for row in rows if parse_bool(row["include_in_tree"], "include_in_tree")]
    return rows


def active_candidate_sources(config: Mapping) -> list[dict]:
    return [row for row in read_tsv(input_path(config, "candidate_sources")) if parse_bool(row["active"], "active")]


def selected_ogs() -> list[str]:
    path = PROJECT_ROOT / "intermediate" / "02_candidates" / "candidate_ogs.tsv"
    return [row["orthogroup"] for row in read_tsv(path) if parse_bool(row.get("include_for_tree", "TRUE"))]


def og_slug(orthogroup: str) -> str:
    if re.fullmatch(r"OG:\d+", orthogroup):
        return orthogroup.replace(":", "_")
    if re.fullmatch(r"FAM_[A-Z0-9_]+", orthogroup):
        return orthogroup
    raise ValueError(f"Invalid orthogroup/family ID: {orthogroup}")


def decision_step(value: str) -> int:
    match = re.search(r"(\d+)", value or "")
    return int(match.group(1)) if match else 10**9


def fail(message: str, code: int = 1) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(code)
