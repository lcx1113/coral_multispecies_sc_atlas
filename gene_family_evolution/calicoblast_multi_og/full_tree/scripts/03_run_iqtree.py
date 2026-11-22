#!/usr/bin/env python3
"""Step 08: infer unrooted ML gene trees with IQ-TREE and UFBoot 1000 only."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

from common import (
    PROJECT_ROOT, alignment_stats, archive_existing, load_config, og_slug,
    reusable_command_record, run_command, selected_ogs, sha256_file, write_json, write_tsv,
)


FIELDS = [
    "orthogroup", "og_slug", "n_sequences", "trimmed_sites", "selected_model",
    "ufboot_replicates", "sh_alrt_replicates", "normal_completion",
    "required_outputs_present", "treefile", "contree", "exit_code", "timed_out",
    "input_sha256", "treefile_sha256", "run_action", "archive_path", "executed_command",
    "logged_command", "status", "error_message",
]


def first_match(text: str, patterns: list[str]) -> str:
    for pattern in patterns:
        match = re.search(pattern, text, re.MULTILINE)
        if match:
            return match.group(1).strip()
    return ""


def choose_ogs(requested: list[str] | None) -> list[str]:
    active = selected_ogs()
    result = requested or active
    unknown = sorted(set(result) - set(active))
    if len(result) != len(set(result)) or unknown:
        raise ValueError(f"Invalid or duplicate requested OGs: {unknown}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/analysis_config.yaml")
    parser.add_argument("--og", action="append")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    config, _ = load_config(args.config)
    cfg = config["tree"]
    if int(cfg["ufboot_replicates"]) != 1000 or int(cfg["sh_alrt_replicates"]) != 0:
        raise ValueError("Step 08 requires UFBoot 1000 and no SH-aLRT")
    ogs = choose_ogs(args.og)
    root = PROJECT_ROOT / "intermediate" / "07_iqtree"
    summaries = []
    for og in ogs:
        slug = og_slug(og)
        alignment = PROJECT_ROOT / "intermediate" / "06_alignments" / slug / f"{slug}.clipkit.aln.faa"
        out_dir = root / slug
        out_dir.mkdir(parents=True, exist_ok=True)
        prefix = out_dir / slug
        treefile = Path(str(prefix) + ".treefile")
        contree = Path(str(prefix) + ".contree")
        report = Path(str(prefix) + ".iqtree")
        log = Path(str(prefix) + ".log")
        stdout = out_dir / f"{slug}.iqtree.stdout.log"
        stderr = out_dir / f"{slug}.iqtree.stderr.log"
        command_json = out_dir / f"{slug}.iqtree.command.json"
        row = {key: "" for key in FIELDS}
        row.update({"orthogroup": og, "og_slug": slug, "status": "ERROR"})
        try:
            stats = alignment_stats(alignment)
            if stats["n_sequences"] < int(cfg["minimum_sequences"]) or stats["alignment_length"] < int(cfg["minimum_trimmed_sites"]):
                raise ValueError(f"Alignment below IQ-TREE minimum: {stats['n_sequences']} sequences, {stats['alignment_length']} sites")
            three_tip_display = stats["n_sequences"] == 3
            command = [
                str(cfg["executable"]), "-s", str(alignment.resolve()), "-m", str(cfg["model_finder"]),
                "-mset", str(cfg["model_set"]),
            ]
            if not three_tip_display:
                command.extend(["-bb", str(cfg["ufboot_replicates"]), "-bnni"])
            command.extend([
                "-nt", str(config["resources"]["iqtree_threads_per_og"]), "-pre", str(prefix.resolve()),
                "-nm", str(cfg["max_iterations"]), "-nstop", str(cfg["stop_generations"]),
                "-me", str(cfg["likelihood_epsilon"]), "-seed", str(config["project"]["random_seed"]),
            ])
            input_sha256 = sha256_file(alignment)
            existing = json.loads(command_json.read_text(encoding="utf-8")) if command_json.is_file() else {}
            required = [treefile, report, log] if three_tip_display else [treefile, contree, report, log]
            reusable = all(path.is_file() and path.stat().st_size for path in required) and reusable_command_record(existing, command, input_sha256)
            archive_path = ""
            if reusable and not args.force:
                record, action = existing, "REUSED_EXISTING"
            else:
                archive_path = archive_existing(
                    list(out_dir.glob(f"{slug}.*")), out_dir / "archive"
                )
                record = run_command(command, stdout, stderr, timeout=int(cfg["timeout_seconds"]))
                record.update({"input_sha256": input_sha256, "archive_path": archive_path})
                write_json(command_json, record)
                action = "EXECUTED"
            log_text = log.read_text(encoding="utf-8", errors="replace") if log.is_file() else ""
            report_text = report.read_text(encoding="utf-8", errors="replace") if report.is_file() else ""
            combined = report_text + "\n" + log_text
            completed = "Analysis results written to:" in combined and "Date and Time:" in combined
            selected_model = first_match(combined, [
                r"Best-fit model according to BIC:\s+([^\s]+)", r"Model of substitution:\s+([^\s]+)",
            ])
            logged_command = first_match(log_text, [r"^Command:\s+(.+)$"])
            outputs_present = all(path.is_file() and path.stat().st_size for path in required)
            status = "PASS" if record["exit_code"] == 0 and completed and outputs_present and selected_model else "ERROR"
            row.update({
                "n_sequences": stats["n_sequences"], "trimmed_sites": stats["alignment_length"],
                "selected_model": selected_model,
                "ufboot_replicates": 0 if three_tip_display else cfg["ufboot_replicates"],
                "sh_alrt_replicates": cfg["sh_alrt_replicates"], "normal_completion": completed,
                "required_outputs_present": outputs_present, "treefile": str(treefile.resolve()),
                "contree": str(contree.resolve()), "exit_code": record["exit_code"],
                "timed_out": record.get("timed_out", False), "input_sha256": input_sha256,
                "treefile_sha256": sha256_file(treefile) if treefile.is_file() else "",
                "run_action": action, "archive_path": archive_path,
                "executed_command": record["command_text"], "logged_command": logged_command,
                "status": status, "error_message": "" if status == "PASS" else "IQ-TREE completion/output/model validation failed",
            })
        except Exception as exc:
            row["error_message"] = str(exc)
        summaries.append(row)
        write_tsv(out_dir / "iqtree_summary.tsv", [row], FIELDS)
        write_json(out_dir / "iqtree_summary.json", {
            "status": row["status"], "orthogroup": og,
            "support_policy": "UFBoot 1000 only; three-tip display trees have zero testable splits and no bootstrap",
        })
        print(f"Step 08 {og}: {row['status']} (model={row['selected_model'] or 'NA'})", file=sys.stderr, flush=True)
    n_errors = sum(row["status"] != "PASS" for row in summaries)
    if len(ogs) > 1 or set(ogs) == set(selected_ogs()):
        write_tsv(root / "iqtree_summary.tsv", summaries, FIELDS)
        write_json(root / "iqtree_summary.json", {
            "status": "PASS" if not n_errors else "FAIL", "n_orthogroups": len(summaries),
            "n_errors": n_errors, "requested_ogs": ogs,
            "support_policy": "UFBoot 1000 only; three-tip display trees have zero testable splits and no bootstrap",
        })
    return 1 if n_errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
