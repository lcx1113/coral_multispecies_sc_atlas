#!/usr/bin/env python3
"""Run the display-pruning workflow from selection through validation."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]


def run(command: list[str]) -> None:
    print("RUN", " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(ROOT / "config" / "analysis_config.yaml"))
    args = parser.parse_args()
    config_path = str(Path(args.config).resolve())
    config = yaml.safe_load(Path(config_path).read_text())
    python = sys.executable
    run([python, str(ROOT / "scripts" / "00_select_tips.py"), "--config", config_path])
    run([config["visualization"]["rscript"], str(ROOT / "scripts" / "01_prune_and_plot_expression_domains.R"), config_path])
    run([python, str(ROOT / "scripts" / "02_validate_outputs.py"), "--config", config_path])
    print("194-tip comparison display pipeline: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

