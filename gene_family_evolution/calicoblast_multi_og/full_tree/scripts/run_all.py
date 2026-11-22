#!/usr/bin/env python3
"""Run the SC11 + Hvul/Aaur reference inference and display workflows with audits."""
from __future__ import annotations
import argparse
import fcntl
import os
import subprocess
import sys
from pathlib import Path
from common import PROJECT_ROOT, load_config, read_tsv, sha256_file, utc_now, write_json, write_tsv

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from-stage", type=int, default=0, choices=range(12))
    parser.add_argument("--to-stage", type=int, default=11, choices=range(12))
    args = parser.parse_args()
    if args.from_stage > args.to_stage:
        parser.error("from-stage must not exceed to-stage")
    cfg, _ = load_config("config/analysis_config.yaml")
    py = cfg["environments"]["gene_tree_python"]
    scope = cfg["multiple_og"]["reconciliation_scope"]
    families = list(dict.fromkeys(r["family_id"] for r in read_tsv(PROJECT_ROOT / cfg["multiple_og"]["family_definitions"])))
    og_args = [v for f in families for v in ("--og", f)]
    commands = []
    for i in range(11):
        stage_cmds = [[py, "scripts/run_pipeline.py", "--from-stage", str(i), "--to-stage", str(i)]]
        if i == 0:
            stage_cmds += [[py, "scripts/validate_reference_scope.py", "--phase", "inputs"]]
        if i == 5:
            stage_cmds += [[py, "scripts/validate_steps_01_05.py", "--scope-label", scope, *og_args]]
        if i == 6:
            stage_cmds += [[py, "scripts/validate_06_reconciliation_dls.py", "--scope-label", scope]]
        commands.append((i, stage_cmds))
    commands.append((11, [
        [py, "scripts/capture_software_versions.py"],
        [py, "display_control_194/scripts/run_pipeline.py"],
        [py, "scripts/validate_reference_scope.py", "--phase", "final"],
    ]))
    logs = PROJECT_ROOT / "logs"
    logs.mkdir(exist_ok=True)
    lock = (logs / "workflow.lock").open("a")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise SystemExit("Another workflow launcher is running in this directory")
    status_file = logs / "workflow_status.json"
    summary = {"status": "RUNNING", "pid": os.getpid(), "started_utc": utc_now(),
               "from_stage": args.from_stage, "to_stage": args.to_stage,
               "families": families, "stages": []}
    if args.from_stage and status_file.exists():
        import json
        old = json.loads(status_file.read_text())
        summary["previous_stages"] = old.get("stages", [])
        summary["stages"] = [r for r in old.get("stages", []) if r["stage"] < args.from_stage and r["status"] == "PASS"]
    sources = [p for folder in ("config", "scripts", "display_control_194/config", "display_control_194/scripts")
               for p in (PROJECT_ROOT / folder).glob("*") if p.is_file()]
    write_tsv(logs / "run_code_config_manifest.tsv", [
        {"path": str(p), "sha256": sha256_file(p)} for p in sources], ["path", "sha256"])
    write_json(status_file, summary)
    for stage, cmds in commands:
        if not args.from_stage <= stage <= args.to_stage:
            continue
        row = {"stage": stage, "started_utc": utc_now(), "status": "RUNNING", "commands": cmds}
        summary["current_stage"] = stage
        summary["stages"].append(row)
        write_json(status_file, summary)
        for k, cmd in enumerate(cmds):
            print(f"START stage {stage} command {k}: {' '.join(cmd)}", flush=True)
            with (logs / f"workflow_{stage:02d}_{k}.log").open("w") as log:
                result = subprocess.run(cmd, cwd=PROJECT_ROOT, stdout=log, stderr=subprocess.STDOUT)
            if result.returncode:
                row.update(status="FAIL", exit_code=result.returncode, finished_utc=utc_now())
                summary.update(status="FAIL", finished_utc=utc_now())
                write_json(status_file, summary)
                print(f"FAIL stage {stage}: see logs/workflow_{stage:02d}_{k}.log", flush=True)
                return result.returncode
        row.update(status="PASS", exit_code=0, finished_utc=utc_now())
        write_json(status_file, summary)
        print(f"PASS stage {stage}", flush=True)
    import json
    freeze = json.loads((PROJECT_ROOT / "inputs/frozen_members/selected_references/freeze_manifest.json").read_text())
    summary["deferred_families"] = freeze.get("deferred_families", [])
    summary["approved_scope_status"] = "PASS" if args.to_stage == 11 else "PARTIAL"
    summary["deferred_scope_status"] = "USER_CONFIRMED_DEFERRED"
    summary.update(status=("PASS_WITH_DEFERRED_FAMILIES" if summary["deferred_families"] else "PASS") if args.to_stage == 11 else "PARTIAL", finished_utc=utc_now())
    write_json(status_file, summary)
    print(summary["status"], flush=True)
    return 0

if __name__ == "__main__":
    sys.exit(main())
