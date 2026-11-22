#!/usr/bin/env python3
"""Step 06: run frozen MAFFT E-INS-i and verify exact leaf retention."""

from __future__ import annotations

import argparse
import json
import sys

from common import (
    PROJECT_ROOT, alignment_stats, archive_existing, fasta_dict, load_config, og_slug,
    parse_bool, reusable_command_record, run_command, selected_ogs, sha256_file,
    write_json, write_tsv,
)


FIELDS = [
    "orthogroup", "og_slug", "input_sequences", "aligned_sequences", "alignment_length",
    "minimum_non_gap", "maximum_non_gap", "mean_gap_fraction", "leaf_ids_match",
    "exit_code", "timed_out", "input_sha256", "output_sha256", "run_action",
    "archive_path", "command", "status", "error_message",
]


def choose_ogs(requested: list[str] | None) -> list[str]:
    active = selected_ogs()
    result = requested or active
    if len(result) != len(set(result)):
        raise ValueError("Duplicate --og values")
    unknown = sorted(set(result) - set(active))
    if unknown:
        raise ValueError(f"Requested OGs are not active candidates: {unknown}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/analysis_config.yaml")
    parser.add_argument("--og", action="append", help="Process this active OG; repeatable")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    config, _ = load_config(args.config)
    cfg = config["alignment"]
    if cfg.get("algorithm") != "E-INS-i" or not parse_bool(cfg.get("genafpair"), "alignment.genafpair"):
        raise ValueError("Step 06 requires the frozen E-INS-i/--genafpair policy")
    ogs = choose_ogs(args.og)
    root = PROJECT_ROOT / "intermediate" / "06_alignments"
    summaries = []
    for og in ogs:
        slug = og_slug(og)
        input_fasta = PROJECT_ROOT / "intermediate" / "05_sequences" / slug / f"{slug}.clean.faa"
        out_dir = root / slug
        output = out_dir / f"{slug}.mafft.aln.faa"
        stderr = out_dir / f"{slug}.mafft.stderr.log"
        record_path = out_dir / f"{slug}.mafft.command.json"
        row = {key: "" for key in FIELDS}
        row.update({"orthogroup": og, "og_slug": slug, "status": "ERROR"})
        try:
            if not input_fasta.is_file():
                raise FileNotFoundError(input_fasta)
            input_ids = set(fasta_dict(input_fasta))
            input_sha256 = sha256_file(input_fasta)
            command = [
                str(cfg["executable"]), "--genafpair", "--maxiterate", str(cfg["maxiterate"]),
                "--quiet", "--thread", str(config["resources"]["mafft_threads_per_og"]),
                str(input_fasta.resolve()),
            ]
            existing = json.loads(record_path.read_text(encoding="utf-8")) if record_path.is_file() else {}
            reusable = output.is_file() and reusable_command_record(existing, command, input_sha256)
            archive_path = ""
            if reusable and not args.force:
                record = existing
                action = "REUSED_EXISTING"
            else:
                archive_path = archive_existing(
                    [output, stderr, record_path], out_dir / "archive"
                )
                record = run_command(command, output, stderr, timeout=int(cfg["timeout_seconds"]))
                record["input_sha256"] = input_sha256
                record["archive_path"] = archive_path
                write_json(record_path, record)
                action = "EXECUTED"
            stats = alignment_stats(output)
            ids_match = stats["ids"] == input_ids
            status = "PASS" if record["exit_code"] == 0 and ids_match and stats["alignment_length"] > 0 else "ERROR"
            error = "" if status == "PASS" else "MAFFT exit/output/leaf validation failed"
            row.update({
                "input_sequences": len(input_ids), "aligned_sequences": stats["n_sequences"],
                "alignment_length": stats["alignment_length"], "minimum_non_gap": stats["minimum_non_gap"],
                "maximum_non_gap": stats["maximum_non_gap"],
                "mean_gap_fraction": f"{stats['mean_gap_fraction']:.8f}", "leaf_ids_match": ids_match,
                "exit_code": record["exit_code"], "timed_out": record.get("timed_out", False),
                "input_sha256": input_sha256, "output_sha256": sha256_file(output),
                "run_action": action, "archive_path": archive_path,
                "command": record["command_text"], "status": status, "error_message": error,
            })
        except Exception as exc:
            row["error_message"] = str(exc)
        summaries.append(row)
        write_tsv(out_dir / "mafft_summary.tsv", [row], FIELDS)
        write_json(out_dir / "mafft_summary.json", {
            "status": row["status"], "orthogroup": og,
            "method": "MAFFT E-INS-i --genafpair --maxiterate 10000",
        })
        print(f"Step 06 {og}: {row['status']} ({row['aligned_sequences'] or 0} sequences x {row['alignment_length'] or 0} sites)", file=sys.stderr, flush=True)
    n_errors = sum(row["status"] != "PASS" for row in summaries)
    if len(ogs) > 1 or set(ogs) == set(selected_ogs()):
        write_tsv(root / "mafft_summary.tsv", summaries, FIELDS)
        write_json(root / "mafft_summary.json", {
            "status": "PASS" if not n_errors else "FAIL", "n_orthogroups": len(summaries),
            "n_errors": n_errors, "requested_ogs": ogs,
            "method": "MAFFT E-INS-i --genafpair --maxiterate 10000",
        })
    return 1 if n_errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
