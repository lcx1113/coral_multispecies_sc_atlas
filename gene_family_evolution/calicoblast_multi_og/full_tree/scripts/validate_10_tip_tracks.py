#!/usr/bin/env python3
"""Validate SC11 expression/CDD tracks, circular trees, and domain PDFs."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path

from common import PROJECT_ROOT, load_config, read_tsv, sha256_file, write_json, write_tsv


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/analysis_config.yaml")
    parser.add_argument("--track-config", default="config/tip_tracks.yaml")
    parser.add_argument("--og", action="append", help="Validate an isolated completed-family subset using per-family leaf lookups")
    args = parser.parse_args()
    config, _ = load_config(args.config)
    track_config, _ = load_config(args.track_config)
    tracks = track_config["multiple_og_tracks"]
    track_root = PROJECT_ROOT / tracks["output_directory"]
    plot_root = PROJECT_ROOT / tracks["plot_output_directory"]

    active = ([r for family in args.og for r in read_tsv(PROJECT_ROOT / "intermediate" / "06_alignments" / family / "active_tree_leaf_lookup.tsv")]
              if args.og else read_tsv(PROJECT_ROOT / "intermediate" / "06_alignments" / "active_tree_leaf_lookup.tsv"))
    expression = read_tsv(track_root / "tip_relative_expression.tsv")
    domain_status = read_tsv(track_root / "cdd_sequence_status.tsv")
    domain_hits = read_tsv(track_root / "cdd_domain_hits.tsv")
    coverage = read_tsv(track_root / "tip_track_coverage.tsv")
    scope = str(config["multiple_og"]["reconciliation_scope"])
    nodes = read_tsv(PROJECT_ROOT / "results" / "03_reconciliation" / scope / "reconciliation_nodes.tsv")
    manifest = read_tsv(plot_root / "family_expression_domain_plot_manifest.tsv")
    circular_root = (
        PROJECT_ROOT / config["multiple_og"]["plot_output_directory"] / scope
    )
    circular_manifest = read_tsv(circular_root / "family_tree_plot_manifest.tsv")
    checks = []

    def check(check_id: str, condition: bool, observed: object, expected: object, detail: str = "") -> None:
        checks.append({
            "check_id": check_id, "status": "PASS" if condition else "FAIL",
            "observed": observed, "expected": expected, "detail": detail,
        })

    active_keys = {(row["orthogroup"], row["leaf_id"]) for row in active}
    expression_keys = {(row["family_id"], row["leaf_id"]) for row in expression}
    status_keys = {(row["family_id"], row["leaf_id"]) for row in domain_status}
    node_keys = {
        (row["orthogroup"], row["leaf_id"])
        for row in nodes if row["node_kind"] == "leaf"
    }
    check("active_key_unique", len(active_keys) == len(active), len(active_keys), len(active))
    check("expression_key_unique", len(expression_keys) == len(expression), len(expression_keys), len(expression))
    check("domain_status_key_unique", len(status_keys) == len(domain_status), len(status_keys), len(domain_status))
    check("expression_exact_active_keys", expression_keys == active_keys, len(expression_keys), len(active_keys))
    check("domain_status_exact_active_keys", status_keys == active_keys, len(status_keys), len(active_keys))
    check("reconciliation_leaf_exact_active_keys", node_keys == active_keys, len(node_keys), len(active_keys))

    active_by_key = {(row["orthogroup"], row["leaf_id"]): row for row in active}
    original_og_match = all(
        row["original_orthogroup"] == active_by_key[(row["family_id"], row["leaf_id"])]["original_orthogroup"]
        for row in expression
    )
    check("original_orthogroup_preserved", original_og_match, original_og_match, True)
    check("celltype_calicoblast_only", {row["celltype"] for row in expression} == {"calicoblast"}, ";".join(sorted({row["celltype"] for row in expression})), "calicoblast")

    allowed_states = {"observed", "no_expression_dataset", "gene_not_tested", "celltype_not_present"}
    observed_rows = [row for row in expression if row["expression_observation_status"] == "observed"]
    missing_rows = [row for row in expression if row["expression_observation_status"] != "observed"]
    expression_numeric = True
    for row in observed_rows:
        try:
            values = [float(row[field]) for field in ["avg_log2FC", "pct_target", "pct_background", "p_val", "p_val_adj"]]
            expression_numeric &= 0 <= values[1] <= 1 and 0 <= values[2] <= 1
        except ValueError:
            expression_numeric = False
    check("expression_states_allowed", {row["expression_observation_status"] for row in expression} <= allowed_states, ";".join(sorted({row["expression_observation_status"] for row in expression})), ";".join(sorted(allowed_states)))
    check("observed_expression_numeric", expression_numeric, len(observed_rows), "all observed numeric")
    check("missing_expression_not_zero", all(row["avg_log2FC"] == "NA" for row in missing_rows), sum(row["avg_log2FC"] == "NA" for row in missing_rows), len(missing_rows))

    hits_by_key = Counter((row["family_id"], row["leaf_id"]) for row in domain_hits)
    hit_keys = set(hits_by_key)
    check("domain_hit_keys_subset_active", hit_keys <= active_keys, len(hit_keys - active_keys), 0)
    domain_count_match = all(
        hits_by_key[key] == int(row["n_representative_domains"])
        for key, row in {
            (row["family_id"], row["leaf_id"]): row for row in domain_status
        }.items()
    )
    check("domain_hit_count_matches_status", domain_count_match, domain_count_match, True)
    domain_state_match = all(
        (row["domain_status"] == "domain_hit") == (hits_by_key[(row["family_id"], row["leaf_id"])] > 0)
        for row in domain_status
    )
    check("domain_state_matches_hits", domain_state_match, domain_state_match, True)
    domain_coordinates_valid = all(
        1 <= int(row["domain_from"]) <= int(row["domain_to"]) <= int(row["protein_length"])
        and int(row["overlap_lane"]) >= 1
        for row in domain_hits
    )
    check("domain_coordinates_valid", domain_coordinates_valid, len(domain_hits), "all within protein length")
    check("cdd_version_frozen", {row["cdd_database_version"] for row in domain_status} == {"CDD.v3.21"}, ";".join(sorted({row["cdd_database_version"] for row in domain_status})), "CDD.v3.21")

    coverage_by_family = {row["family_id"]: row for row in coverage}
    manifest_by_family = {row["family_id"]: row for row in manifest}
    circular_by_family = {row["family_id"]: row for row in circular_manifest}
    families = sorted({row["orthogroup"] for row in active})
    check("coverage_family_set", set(coverage_by_family) == set(families), len(coverage_by_family), len(families))
    check("plot_manifest_family_set", set(manifest_by_family) == set(families), len(manifest_by_family), len(families))
    check("circular_manifest_family_set", set(circular_by_family) == set(families), len(circular_by_family), len(families))
    expression_by_family = defaultdict(list)
    status_by_family = defaultdict(list)
    for row in expression:
        expression_by_family[row["family_id"]].append(row)
    for row in domain_status:
        status_by_family[row["family_id"]].append(row)
    for family_id in families:
        states = Counter(row["expression_observation_status"] for row in expression_by_family[family_id])
        cov = coverage_by_family[family_id]
        plot = manifest_by_family[family_id]
        circular = circular_by_family[family_id]
        expected_tips = len(expression_by_family[family_id])
        check(f"{family_id}_coverage_tip_count", int(cov["n_active_tips"]) == expected_tips, cov["n_active_tips"], expected_tips)
        check(f"{family_id}_plot_tip_count", int(plot["n_tips"]) == expected_tips, plot["n_tips"], expected_tips)
        check(f"{family_id}_plot_observed_count", int(plot["n_observed_tips"]) == states["observed"], plot["n_observed_tips"], states["observed"])
        expected_hit_tips = sum(row["domain_status"] == "domain_hit" for row in status_by_family[family_id])
        check(f"{family_id}_plot_domain_tip_count", int(plot["n_domain_hit_tips"]) == expected_hit_tips, plot["n_domain_hit_tips"], expected_hit_tips)
        check(f"{family_id}_circular_tip_count", int(circular["n_tips"]) == expected_tips, circular["n_tips"], expected_tips)
        check(
            f"{family_id}_circular_observed_count",
            int(circular["n_observed_tips"]) == states["observed"],
            circular["n_observed_tips"], states["observed"],
        )
        expected_markers = sum(row["positive_marker"] == "TRUE" for row in expression_by_family[family_id])
        check(
            f"{family_id}_circular_marker_count",
            int(circular["n_marker_tips"]) == expected_markers,
            circular["n_marker_tips"], expected_markers,
        )
        check(
            f"{family_id}_circular_loss_hidden",
            circular["loss_displayed"] == "FALSE",
            circular["loss_displayed"], "FALSE",
        )
        for variant in ["scientific_abbreviation", "original_identifier"]:
            pdf = Path(circular[f"{variant}_pdf"])
            pdf_ok = (
                circular["status"] == "PASS"
                and circular[f"{variant}_ghostscript_status"] == "0"
                and pdf.is_file()
                and pdf.stat().st_size == int(circular[f"{variant}_size_bytes"])
                and sha256_file(pdf) == circular[f"{variant}_sha256"]
            )
            check(
                f"{family_id}_{variant}_circular_pdf",
                pdf_ok, str(pdf), "existing hash-matched Ghostscript PASS PDF",
            )
        pdf = Path(plot["pdf"])
        pdf_ok = (
            plot["status"] == "PASS" and plot["ghostscript_status"] == "0"
            and pdf.is_file() and pdf.stat().st_size == int(plot["size_bytes"])
            and sha256_file(pdf) == plot["sha256"]
        )
        check(f"{family_id}_pdf_manifest", pdf_ok, str(pdf), "existing hash-matched Ghostscript PASS PDF")

    validation_root = track_root / "validation"
    validation_root.mkdir(parents=True, exist_ok=True)
    out_tsv = validation_root / "integrated_tracks_validation.tsv"
    out_json = validation_root / "integrated_tracks_validation_summary.json"
    write_tsv(out_tsv, checks, ["check_id", "status", "observed", "expected", "detail"])
    failures = [row for row in checks if row["status"] != "PASS"]
    write_json(out_json, {
        "status": "PASS" if not failures else "FAIL", "n_checks": len(checks),
        "n_failures": len(failures), "failed_checks": [row["check_id"] for row in failures],
        "n_active_tips": len(active_keys), "n_expression_observed": len(observed_rows),
        "n_cdd_domain_hits": len(domain_hits),
        "n_domain_pdfs": len(manifest),
        "n_circular_pdfs": 2 * len(circular_manifest),
    })
    print(f"Integrated tip-track validation: {'PASS' if not failures else 'FAIL'} ({len(checks)} checks; {len(failures)} failures)")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
