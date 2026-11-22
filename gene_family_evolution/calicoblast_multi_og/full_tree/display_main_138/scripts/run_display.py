#!/usr/bin/env python3
"""Select representatives, prune/plot and validate; never call inference tools."""
import json,subprocess,sys
import yaml
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1]
now=lambda:datetime.now(timezone.utc).isoformat()
logs=ROOT/"logs";logs.mkdir(exist_ok=True)
config=yaml.safe_load((ROOT.parent/"display_control_194/config/analysis_config.yaml").read_text())
commands=[
 ["select",[sys.executable,"scripts/00_select_reference_representatives.py"]],
 ["plot",[config["visualization"]["rscript"],"scripts/01_prune_and_plot_expression_domains.R",str(ROOT/"config/analysis_config.yaml")]],
 ["validate",[sys.executable,"scripts/02_validate_outputs.py","--config",str(ROOT/"config/analysis_config.yaml")]],
]
status=dict(status="RUNNING",started_utc=now(),stages=[])
save=lambda:(logs/"run_status.json").write_text(json.dumps(status,indent=2)+"\n")
save()
for name,cmd in commands:
 row=dict(stage=name,command=cmd,status="RUNNING",started_utc=now());status["stages"].append(row);save()
 print("START",name,flush=True)
 with (logs/(name+".log")).open("w") as log:p=subprocess.run(cmd,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
 row.update(status="PASS" if not p.returncode else "FAIL",exit_code=p.returncode,finished_utc=now());save()
 if p.returncode:
  status.update(status="FAIL",finished_utc=now());save();raise SystemExit(p.returncode)
 print("PASS",name,flush=True)
status.update(status="PASS",finished_utc=now());save()
