#!/usr/bin/env python3
"""Run SC11 all-member combined-family trees on the pruned SC11 species tree."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from common import PROJECT_ROOT, load_config, read_tsv, sha256_file, write_json, write_tsv


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/analysis_config.yaml")
    parser.add_argument("--from-stage", type=int, default=0, choices=range(0, 11))
    parser.add_argument("--to-stage", type=int, default=10, choices=range(0, 11))
    parser.add_argument("--force-generax", action="store_true")
    args = parser.parse_args()
    if args.from_stage > args.to_stage:
        parser.error("--from-stage cannot exceed --to-stage")

    config, config_path = load_config(args.config)
    track_config, track_config_path = load_config("config/tip_tracks.yaml")
    tracks = track_config["multiple_og_tracks"]
    python = str(Path(config["environments"]["gene_tree_python"]).resolve())
    rscript = str(Path(config["environments"]["rscript"]).resolve())
    scope = str(config["multiple_og"]["reconciliation_scope"])
    config_abs = str(config_path.resolve())
    family_ids = [
        row["family_id"] for row in read_tsv(PROJECT_ROOT / config["multiple_og"]["family_definitions"])
        if row.get("active", "").upper() == "TRUE"
    ]
    family_ids = list(dict.fromkeys(family_ids))
    og_args = [item for family_id in family_ids for item in ("--og", family_id)]

    stages: list[tuple[int, str, list[str]]] = [
        (0, "prepare_sc11_external_reference_inputs", [python, str(PROJECT_ROOT / "scripts/prepare_reference_inputs.py"), "--config", config_abs]),
        (1, "mafft_einsi", [python, str(PROJECT_ROOT / "scripts/01_run_mafft.py"), "--config", config_abs, *og_args]),
        (2, "clipkit", [python, str(PROJECT_ROOT / "scripts/02_run_clipkit.py"), "--config", config_abs, *og_args]),
        (3, "iqtree_ufboot", [python, str(PROJECT_ROOT / "scripts/03_run_iqtree.py"), "--config", config_abs, *og_args]),
        (4, "prepare_reconciliation", [python, str(PROJECT_ROOT / "scripts/04_prepare_reconciliation_inputs.py"), "--config", config_abs, "--scope-label", scope, *og_args]),
        (5, "generax_eval", [python, str(PROJECT_ROOT / "scripts/05_run_generax_eval.py"), "--config", config_abs, "--scope-label", scope] + (["--force"] if args.force_generax else [])),
        (6, "parse_reconciliation", [python, str(PROJECT_ROOT / "scripts/06_parse_reconciliation_dls.py"), "--config", config_abs, "--scope-label", scope]),
        (7, "prepare_calicoblast_expression_cdd_tracks", [python, str(PROJECT_ROOT / "scripts/prepare_reference_tracks.py"), "--config", config_abs, "--track-config", str(track_config_path.resolve())]),
        (8, "plot_circular_expression_trees", [rscript, str(PROJECT_ROOT / "scripts/07_plot_family_trees.R"), config_abs, scope, str(track_config_path.resolve())]),
        (9, "plot_expression_domains", [rscript, str(PROJECT_ROOT / "scripts/09_plot_family_expression_domains.R"), config_abs, scope, str(track_config_path.resolve())]),
        (10, "validate_expression_domains", [python, str(PROJECT_ROOT / "scripts/validate_10_tip_tracks.py"), "--config", config_abs, "--track-config", str(track_config_path.resolve())]),
    ]

    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    status_path = log_dir / "pipeline_stage_status.tsv"
    status_rows = []
    if args.from_stage > 0 and status_path.is_file():
        for prior in read_tsv(status_path):
            if int(prior["stage"]) < args.from_stage:
                prior["stage"] = int(prior["stage"])
                prior["exit_code"] = int(prior["exit_code"])
                status_rows.append(prior)
    for stage, label, command in stages:
        if not args.from_stage <= stage <= args.to_stage:
            continue
        started = utc_now()
        log_path = log_dir / f"stage_{stage:02d}_{label}.log"
        print(f"START stage {stage}: {label}", flush=True)
        with log_path.open("w", encoding="utf-8", newline="\n") as log_handle:
            log_handle.write("COMMAND\t" + " ".join(command) + "\n")
            process = subprocess.Popen(
                command, cwd=PROJECT_ROOT, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, text=True, bufsize=1,
            )
            assert process.stdout is not None
            for line in process.stdout:
                sys.stdout.write(line)
                sys.stdout.flush()
                log_handle.write(line)
            exit_code = process.wait()
        row = {
            "stage": stage, "label": label, "started_utc": started,
            "finished_utc": utc_now(), "exit_code": exit_code,
            "status": "PASS" if exit_code == 0 else "FAIL",
            "log": str(log_path.resolve()), "log_sha256": sha256_file(log_path),
            "command": " ".join(command),
        }
        status_rows.append(row)
        write_tsv(status_path, status_rows, list(row))
        print(f"{'PASS' if exit_code == 0 else 'FAIL'} stage {stage}: {label}", flush=True)
        if exit_code != 0:
            write_json(log_dir / "pipeline_run_summary.json", {
                "status": "FAIL", "failed_stage": stage, "failed_label": label,
                "scope_label": scope, "family_ids": family_ids,
                "config": config_abs, "stages": status_rows,
            })
            return exit_code or 1

    final_status = "PASS"
    plot_manifest = PROJECT_ROOT / config["multiple_og"]["plot_output_directory"] / scope / "family_tree_plot_manifest.tsv"
    if args.to_stage >= 8:
        if not plot_manifest.is_file():
            final_status = "FAIL"
        else:
            plot_rows = read_tsv(plot_manifest)
            if len(plot_rows) != len(family_ids) or any(row["status"] != "PASS" for row in plot_rows):
                final_status = "FAIL"
    integrated_manifest = PROJECT_ROOT / tracks["plot_output_directory"] / "family_expression_domain_plot_manifest.tsv"
    if args.to_stage >= 9:
        if not integrated_manifest.is_file():
            final_status = "FAIL"
        else:
            integrated_rows = read_tsv(integrated_manifest)
            if len(integrated_rows) != len(family_ids) or any(row["status"] != "PASS" for row in integrated_rows):
                final_status = "FAIL"
    summary = {
        "status": final_status, "scope_label": scope, "family_ids": family_ids,
        "n_families": len(family_ids), "from_stage": args.from_stage,
        "to_stage": args.to_stage, "config": config_abs,
        "plot_manifest": str(plot_manifest.resolve()) if plot_manifest.is_file() else "",
        "integrated_plot_manifest": str(integrated_manifest.resolve()) if integrated_manifest.is_file() else "",
        "track_config": str(track_config_path.resolve()),
        "stages": status_rows,
    }
    write_json(log_dir / "pipeline_run_summary.json", summary)
    print(f"Combined-family pipeline: {final_status} ({len(family_ids)} families)", flush=True)
    return 0 if final_status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
