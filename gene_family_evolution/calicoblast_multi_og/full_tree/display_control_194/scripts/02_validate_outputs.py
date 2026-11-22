#!/usr/bin/env python3
"""Validate selection, pruned trees, inherited annotations, tracks and PDFs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from collections import defaultdict
from pathlib import Path

import yaml
from Bio import Phylo


ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(ROOT / "config" / "analysis_config.yaml"))
    args = parser.parse_args()
    config = yaml.safe_load(Path(args.config).read_text())
    source_root = Path(config["source"]["project"])
    selection_dir = ROOT / config["outputs"]["selection_dir"]
    tree_dir = ROOT / config["outputs"]["tree_dir"]
    validation_dir = ROOT / config["outputs"]["validation_dir"]
    validation_dir.mkdir(parents=True, exist_ok=True)

    audit = read_tsv(selection_dir / "tip_selection_audit.tsv")
    og_summary = read_tsv(selection_dir / "og_selection_summary.tsv")
    article_dir = ROOT / config["outputs"]["article_figure_dir"]
    article_manifest = read_tsv(article_dir / "article_final_figure_manifest.tsv")
    article_labels = read_tsv(article_dir / "article_final_tip_labels.tsv")
    mapping = read_tsv(tree_dir / "all_families.pruned_node_annotation_mapping.tsv")
    source_expression = read_tsv(source_root / config["source"]["expression"])
    source_status = read_tsv(source_root / config["source"]["cdd_sequence_status"])
    source_hits = read_tsv(source_root / config["source"]["cdd_domain_hits"])
    definitions = [r for r in read_tsv(ROOT / config.get("family_definitions", "config/family_sets.tsv")) if r["active"].upper() == "TRUE"]
    checks: list[dict[str, object]] = []

    def check(check_id: str, condition: bool, observed: object, expected: object, detail: str = "") -> None:
        checks.append({"check_id": check_id, "status": "PASS" if condition else "FAIL",
                       "observed": observed, "expected": expected, "detail": detail})

    expected_families = {r["family_id"] for r in definitions}
    expected_ogs = {(r["family_id"], r["original_orthogroup"]) for r in definitions}
    selected = [r for r in audit if r["retain_tip"] == "TRUE"]
    selected_keys = {(r["family_id"], r["leaf_id"]) for r in selected}
    check("audit_source_row_count", len(audit) == len(source_expression), len(audit), len(source_expression))
    check("audit_key_unique", len({(r["family_id"], r["leaf_id"]) for r in audit}) == len(audit), len(audit), "all unique")
    check("selection_family_set", {r["family_id"] for r in selected} == expected_families, len({r["family_id"] for r in selected}), len(expected_families))
    check("selection_og_set", {(r["family_id"], r["original_orthogroup"]) for r in selected if r["original_orthogroup"] != "NA"} == expected_ogs, len({(r["family_id"], r["original_orthogroup"]) for r in selected}), len(expected_ogs))

    audit_by_og: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    audit_by_species_og: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in audit:
        audit_by_og[(row["family_id"], row["original_orthogroup"])].append(row)
        audit_by_species_og[(row["family_id"], row["original_orthogroup"], row["tree_species"])].append(row)
    rule_ok = True
    no_observed_group_count = 0
    for key, rows in audit_by_species_og.items():
        markers = {r["leaf_id"] for r in rows if r["positive_marker"] == "TRUE"}
        retained = {r["leaf_id"] for r in rows if r["retain_tip"] == "TRUE"}
        if key[2] in config["selection"].get("retain_all_species_without_expression", []):
            rule_ok &= retained == {r["leaf_id"] for r in rows}
            rule_ok &= all(r["expression_observation_status"] == "no_expression_dataset" and r["positive_marker"] == "NA" for r in rows)
            continue
        observed = [r for r in rows if r["expression_observation_status"] == "observed" and r["avg_log2FC"] not in {"", "NA"}]
        if observed:
            top = sorted(observed, key=lambda r: (-float(r["avg_log2FC"]), -float(r["pct_target"]), float(r["p_val_adj"]), r["leaf_id"]))[0]
        else:
            no_observed_group_count += 1
            top = sorted(rows, key=lambda r: r["leaf_id"])[0]
        rule_ok &= retained == markers | {top["leaf_id"]}
    check("selection_rule_exact", rule_ok, rule_ok, True)
    check("species_og_group_count", len(audit_by_species_og) == sum(int(r["n_species"]) for r in og_summary), len(audit_by_species_og), sum(int(r["n_species"]) for r in og_summary))
    check("no_observed_species_group_count", no_observed_group_count == sum(int(r["n_no_observed_fallback_tips"]) for r in og_summary), no_observed_group_count, sum(int(r["n_no_observed_fallback_tips"]) for r in og_summary))

    article_by_family = {r["family_id"]: r for r in article_manifest}
    check("article_manifest_family_set", set(article_by_family) == expected_families, len(article_by_family), len(expected_families))
    article_label_keys = {(r["family_id"], r["leaf_id"]) for r in article_labels}
    check("article_tip_label_key_unique", len(article_label_keys) == len(article_labels), len(article_label_keys), len(article_labels))
    check("article_tip_label_exact_selected_set", article_label_keys == selected_keys, len(article_label_keys), len(selected_keys))
    article_label_format_ok = all(
        r["article_label"] == (f"{r['tree_species']}|{r['protein_id']}" if r["tree_species"] in {"Hvul", "Aaur"} else f"{r['tree_species']}|{r['Uniprot_name'] or 'NA'}|{r['gene_id'] or 'NA'}")
        for r in article_labels
    )
    check("article_tip_label_format", article_label_format_ok, article_label_format_ok, "species|Uniprot_name|gene_id")
    selected_by_family: dict[str, set[str]] = defaultdict(set)
    for family_id, leaf_id in selected_keys:
        selected_by_family[family_id].add(leaf_id)
    status_keys = {(r["family_id"], r["leaf_id"]) for r in source_status}
    hit_keys = {(r["family_id"], r["leaf_id"]) for r in source_hits}
    selected_missing = [r for r in selected if r["expression_observation_status"] != "observed"]
    check("selected_missing_only_deterministic_fallback_or_external_reference", all("species_fallback_no_observed" in r["selection_reason"] or "retain_all_external_reference_no_expression_filter" in r["selection_reason"] for r in selected_missing), len(selected_missing), "flagged fallback or all external_reference retained")
    check("selected_cdd_status_covered", selected_keys <= status_keys, len(selected_keys - status_keys), 0)
    check("selected_cdd_hits_subset_status", hit_keys <= status_keys, len(hit_keys - status_keys), 0)

    for family_id in sorted(expected_families):
        article_row = article_by_family[family_id]
        tree_path = Path(article_row["tree"])
        tree = Phylo.read(tree_path, "newick")
        tree_tips = {tip.name for tip in tree.get_terminals()}
        check(f"{family_id}_tree_tip_set", tree_tips == selected_by_family[family_id], len(tree_tips), len(selected_by_family[family_id]))
        check(f"{family_id}_manifest_tip_count", int(article_row["n_retained_tips"]) == len(tree_tips), article_row["n_retained_tips"], len(tree_tips))
        check(f"{family_id}_tree_hash", tree_path.is_file() and sha256(tree_path) == article_row["tree_sha256"], tree_path, "existing hash-matched tree")
        article_pdf_path = Path(article_row["pdf"])
        article_gs = subprocess.run([config["visualization"]["ghostscript"], "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=nullpage", str(article_pdf_path)], capture_output=True)
        article_pdf_ok = article_pdf_path.is_file() and sha256(article_pdf_path) == article_row["sha256"] and article_gs.returncode == 0 and int(article_row["ghostscript_status"]) == 0
        check(f"{family_id}_article_pdf", article_pdf_ok, article_pdf_path, "existing hash-matched Ghostscript PASS article PDF")
        og_by_tip = {
            r["leaf_id"]: r["original_orthogroup"]
            for r in selected if r["family_id"] == family_id
        }
        expected_cross_og_nodes = sum(
            len({og_by_tip[t.name] for t in clade.get_terminals()}) > 1
            for clade in tree.get_nonterminals()
        )
        check(
            f"{family_id}_cross_og_node_count",
            int(article_row["n_cross_og_annotated_nodes"]) == expected_cross_og_nodes,
            article_row["n_cross_og_annotated_nodes"], expected_cross_og_nodes,
        )
        check(
            f"{family_id}_pure_og_events_omitted",
            int(article_row["n_pure_og_internal_events_displayed"]) == 0,
            article_row["n_pure_og_internal_events_displayed"], 0,
        )
        check(
            f"{family_id}_article_same_tip_count",
            int(article_row["n_retained_tips"]) == len(tree_tips),
            article_row["n_retained_tips"], len(tree_tips),
        )

    ambiguous = [r for r in mapping if r["mapping_status"] != "mapped_unique_informative_original_node"]
    ambiguous_blank = all(not r["inherited_node_id"] and not r["inherited_event_type"] and not r["inherited_iqtree_ufboot"] for r in ambiguous)
    check("ambiguous_annotations_hidden", ambiguous_blank, len(ambiguous), "all ambiguous annotations blank")
    check("node_mapping_family_set", {r["family_id"] for r in mapping} == expected_families, len({r["family_id"] for r in mapping}), len(expected_families))

    write_tsv(validation_dir / "validation_checks.tsv", checks, ["check_id", "status", "observed", "expected", "detail"])
    failures = [r for r in checks if r["status"] == "FAIL"]
    summary = {
        "status": "PASS" if not failures else "FAIL",
        "n_checks": len(checks),
        "n_failures": len(failures),
        "failed_checks": [r["check_id"] for r in failures],
        "n_families": len(expected_families),
        "n_ogs": len(expected_ogs),
        "n_source_tips": len(audit),
        "n_retained_tips": len(selected),
        "n_marker_tips": sum(r["positive_marker"] == "TRUE" for r in selected),
        "n_species_top_nonmarker_tips": sum(r["positive_marker"] != "TRUE" and "species_top_observed" in r["selection_reason"] for r in selected),
        "n_no_observed_fallback_tips": sum("species_fallback_no_observed" in r["selection_reason"] for r in selected),
        "n_pruned_internal_nodes": len(mapping),
        "n_unique_inherited_internal_nodes": sum(r["mapping_status"] == "mapped_unique_informative_original_node" for r in mapping),
        "n_hidden_ambiguous_or_unmapped_internal_nodes": len(ambiguous),
        "n_article_final_pdfs": len(article_manifest),
        "n_article_cross_og_annotated_nodes": sum(int(r["n_cross_og_annotated_nodes"]) for r in article_manifest),
        "n_article_tip_labels": len(article_labels),
    }
    (validation_dir / "validation_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, sort_keys=True))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
