#!/usr/bin/env python3
"""Select all markers plus one highest-log2FC tip per species in every OG."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def as_float(value: str, default: float) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(ROOT / "config" / "analysis_config.yaml"))
    args = parser.parse_args()
    config_path = Path(args.config).resolve()
    config = yaml.safe_load(config_path.read_text())
    source_root = Path(config["source"]["project"])
    expression_path = source_root / config["source"]["expression"]
    family_path = ROOT / config.get("family_definitions", "config/family_sets.tsv")
    expression = read_tsv(expression_path)
    definitions = [r for r in read_tsv(family_path) if r["active"].upper() == "TRUE"]
    expected_pairs = {(r["family_id"], r["original_orthogroup"]) for r in definitions}

    keyed: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in expression:
        keyed[(row["family_id"], row["original_orthogroup"])].append(row)
    if {key for key in keyed if key[1] != "NA"} != expected_pairs:
        missing = sorted(expected_pairs - set(keyed))
        extra = sorted(set(keyed) - expected_pairs)
        raise SystemExit(f"Family/OG mismatch; missing={missing}; extra={extra}")

    selected: set[tuple[str, str]] = set()
    reason_sets: dict[tuple[str, str], set[str]] = defaultdict(set)
    og_summaries: list[dict[str, object]] = []
    species_summaries: list[dict[str, object]] = []
    for family_id, original_og in sorted(keyed):
        rows = keyed[(family_id, original_og)]
        markers = [r for r in rows if r["positive_marker"].upper() == "TRUE"]
        keep_by_leaf: dict[str, dict[str, str]] = {}
        for row in markers:
            keep_by_leaf[row["leaf_id"]] = row
            reason_sets[(family_id, row["leaf_id"])].add("positive_marker")
        rows_by_species: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            rows_by_species[row["tree_species"]].append(row)
        for species in sorted(rows_by_species):
            species_rows = rows_by_species[species]
            if species in config["selection"].get("retain_all_species_without_expression", []):
                for row in species_rows:
                    if row["expression_observation_status"] != "no_expression_dataset" or row["positive_marker"] != "NA":
                        raise ValueError(f"Unexpected expression for exempt species: {row['leaf_id']}")
                    keep_by_leaf[row["leaf_id"]] = row
                    reason_sets[(family_id, row["leaf_id"])].add("retain_all_external_reference_no_expression_filter")
                species_summaries.append({
                    "family_id": family_id, "original_orthogroup": original_og,
                    "tree_species": species, "n_source_tips": len(species_rows),
                    "n_observed_tips": 0, "n_marker_tips": 0,
                    "top_selection_mode": "retain_all_external_reference_no_expression_filter",
                    "top_leaf_id": "", "top_is_marker": "NA", "top_avg_log2FC": "NA",
                    "top_pct_target": "NA", "top_p_val_adj": "NA",
                    "n_tied_at_highest_log2FC": 0, "n_retained_tips": len(species_rows),
                    "retained_leaf_ids": ";".join(sorted(r["leaf_id"] for r in species_rows)),
                })
                continue
            species_markers = [r for r in species_rows if r["positive_marker"].upper() == "TRUE"]
            observed = [r for r in species_rows
                        if r["expression_observation_status"] == "observed"
                        and r["avg_log2FC"] not in {"", "NA", "NaN", "nan"}]
            tie_count = 0
            if observed:
                ranked = sorted(
                    observed,
                    key=lambda r: (
                        -as_float(r["avg_log2FC"], float("-inf")),
                        -as_float(r["pct_target"], float("-inf")),
                        as_float(r["p_val_adj"], float("inf")),
                        r["leaf_id"],
                    ),
                )
                top = ranked[0]
                best_fc = as_float(top["avg_log2FC"], float("-inf"))
                tie_count = sum(as_float(r["avg_log2FC"], float("-inf")) == best_fc for r in observed)
                top_mode = "species_top_observed_avg_log2FC"
            else:
                top = sorted(species_rows, key=lambda r: r["leaf_id"])[0]
                top_mode = "species_fallback_no_observed_lexical_first"
            keep_by_leaf[top["leaf_id"]] = top
            reason_sets[(family_id, top["leaf_id"])].add(top_mode)
            species_summaries.append({
                "family_id": family_id,
                "original_orthogroup": original_og,
                "tree_species": species,
                "n_source_tips": len(species_rows),
                "n_observed_tips": len(observed),
                "n_marker_tips": len(species_markers),
                "top_selection_mode": top_mode,
                "top_leaf_id": top["leaf_id"],
                "top_is_marker": "TRUE" if top["positive_marker"].upper() == "TRUE" else "FALSE",
                "top_avg_log2FC": top["avg_log2FC"],
                "top_pct_target": top["pct_target"],
                "top_p_val_adj": top["p_val_adj"],
                "n_tied_at_highest_log2FC": tie_count,
                "n_retained_tips": len({r["leaf_id"] for r in species_markers} | {top["leaf_id"]}),
                "retained_leaf_ids": ";".join(sorted({r["leaf_id"] for r in species_markers} | {top["leaf_id"]})),
            })
        keep = list(keep_by_leaf.values())
        for row in keep:
            selected.add((family_id, row["leaf_id"]))
        og_summaries.append({
            "family_id": family_id,
            "original_orthogroup": original_og,
            "n_source_tips": len(rows),
            "n_species": len(rows_by_species),
            "n_observed_tips": sum(r["expression_observation_status"] == "observed" for r in rows),
            "n_marker_tips": len(markers),
            "selection_mode": "all_markers_union_per_species_top_log2FC",
            "n_retained_tips": len(keep),
            "n_marker_retained_tips": sum(r["positive_marker"].upper() == "TRUE" for r in keep),
            "n_top_only_nonmarker_tips": sum(
                r["positive_marker"].upper() != "TRUE" and
                "species_top_observed_avg_log2FC" in reason_sets[(family_id, r["leaf_id"])]
                for r in keep
            ),
            "n_no_observed_fallback_tips": sum(
                "species_fallback_no_observed_lexical_first" in reason_sets[(family_id, r["leaf_id"])]
                for r in keep
            ),
            "retained_leaf_ids": ";".join(sorted(r["leaf_id"] for r in keep)),
        })

    audit: list[dict[str, object]] = []
    for row in expression:
        key = (row["family_id"], row["leaf_id"])
        out = dict(row)
        out["retain_tip"] = "TRUE" if key in selected else "FALSE"
        out["selection_reason"] = ";".join(sorted(reason_sets.get(key, {"not_selected_nonmarker_non_top"})))
        audit.append(out)

    selection_dir = ROOT / config["outputs"]["selection_dir"]
    selection_dir.mkdir(parents=True, exist_ok=True)
    audit_path = selection_dir / "tip_selection_audit.tsv"
    summary_path = selection_dir / "og_selection_summary.tsv"
    species_summary_path = selection_dir / "species_og_selection_summary.tsv"
    write_tsv(audit_path, audit, list(expression[0]) + ["retain_tip", "selection_reason"])
    write_tsv(summary_path, og_summaries, list(og_summaries[0]))
    write_tsv(species_summary_path, species_summaries, list(species_summaries[0]))
    for family_id in sorted({r["family_id"] for r in definitions}):
        family_tips = sorted(leaf for fam, leaf in selected if fam == family_id)
        (selection_dir / f"{family_id}.retained_tips.txt").write_text("".join(f"{x}\n" for x in family_tips))

    source_files = []
    for role, path in [
        ("analysis_config", config_path),
        ("family_definitions", family_path),
        ("source_expression", expression_path),
    ]:
        source_files.append({"role": role, "path": str(path), "sha256": sha256(path), "size_bytes": path.stat().st_size})
    write_tsv(selection_dir / "selection_input_files.tsv", source_files, list(source_files[0]))
    result = {
        "status": "PASS",
        "n_families": len({r["family_id"] for r in definitions}),
        "n_ogs": sum(r["original_orthogroup"] != "NA" for r in og_summaries),
        "n_reference_groups": sum(r["original_orthogroup"] == "NA" for r in og_summaries),
        "n_source_tips": len(expression),
        "n_retained_tips": len(selected),
        "n_marker_retained_tips": sum(r["retain_tip"] == "TRUE" and r["positive_marker"].upper() == "TRUE" for r in audit),
        "n_species_og_groups": len(species_summaries),
        "n_top_only_nonmarker_tips": sum(
            r["retain_tip"] == "TRUE" and r["positive_marker"].upper() != "TRUE" and
            "species_top_observed_avg_log2FC" in r["selection_reason"] for r in audit
        ),
        "n_no_observed_species_fallback_tips": sum(
            r["retain_tip"] == "TRUE" and "species_fallback_no_observed_lexical_first" in r["selection_reason"]
            for r in audit
        ),
        "selection_rule": "SC11: all markers union top observed log2FC per OG/species; Hvul/Aaur: retain all",
    }
    (selection_dir / "tip_selection_summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
