#!/usr/bin/env python3
"""Independently validate formal Steps 06-10 for one explicit OG scope."""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from pathlib import Path

from Bio import Phylo

from common import (
    PROJECT_ROOT, alignment_stats, fasta_dict, load_config, og_slug, parse_bool,
    read_tsv, selected_ogs, sha256_file, write_json, write_tsv,
)


SMOKE_OGS = ["OG:000304", "OG:005580"]


def tree_tips(path: Path) -> set[str]:
    tree = Phylo.read(path, "newick")
    names = [tip.name for tip in tree.get_terminals()]
    if None in names or len(names) != len(set(names)):
        raise ValueError(f"Missing/duplicate tree tips in {path}")
    return set(names)


def row_map(path: Path) -> dict[str, dict]:
    rows = read_tsv(path)
    if len(rows) != len({row["orthogroup"] for row in rows}):
        raise ValueError(f"Duplicate orthogroup rows in {path}")
    return {row["orthogroup"]: row for row in rows}


def per_og_row_map(root: Path, filename: str, ogs: list[str], fallback: Path) -> dict[str, dict]:
    result = {}
    fallback_rows = None
    for og in ogs:
        path = root / og_slug(og) / filename
        if path.is_file():
            rows = read_tsv(path)
        else:
            if fallback_rows is None:
                fallback_rows = row_map(fallback)
            rows = [fallback_rows[og]] if og in fallback_rows else []
        if len(rows) != 1 or rows[0]["orthogroup"] != og:
            raise ValueError(f"Expected one per-OG summary row for {og}: {path}")
        result[og] = rows[0]
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/analysis_config.yaml")
    parser.add_argument("--og", action="append")
    parser.add_argument("--smoke-test", action="store_true")
    parser.add_argument("--scope-label", default="")
    args = parser.parse_args()
    config, _ = load_config(args.config)
    if args.smoke_test and args.og:
        parser.error("--smoke-test and --og are mutually exclusive")
    ogs = SMOKE_OGS if args.smoke_test else (args.og or selected_ogs())
    scope_label = "smoke_test" if args.smoke_test else args.scope_label
    checks = []

    def check(check_id: str, condition: bool, observed: object, expected: object, detail: str = "") -> None:
        checks.append({
            "check_id": check_id, "status": "PASS" if condition else "FAIL",
            "observed": observed, "expected": expected, "detail": detail,
        })

    expected = set(ogs)
    alignment_root = PROJECT_ROOT / "intermediate" / "06_alignments"
    iqtree_root = PROJECT_ROOT / "intermediate" / "07_iqtree"
    mafft = per_og_row_map(alignment_root, "mafft_summary.tsv", ogs, alignment_root / "mafft_summary.tsv")
    clipkit = per_og_row_map(alignment_root, "clipkit_summary.tsv", ogs, alignment_root / "clipkit_summary.tsv")
    iqtree = per_og_row_map(iqtree_root, "iqtree_summary.tsv", ogs, iqtree_root / "iqtree_summary.tsv")
    check("06_summary_family_set", set(mafft) == expected, ";".join(sorted(mafft)), ";".join(sorted(expected)))
    check("07_summary_family_set", set(clipkit) == expected, ";".join(sorted(clipkit)), ";".join(sorted(expected)))
    check("08_summary_family_set", set(iqtree) == expected, ";".join(sorted(iqtree)), ";".join(sorted(expected)))
    for og in ogs:
        slug = og_slug(og)
        source = PROJECT_ROOT / "intermediate" / "05_sequences" / slug / f"{slug}.clean.faa"
        mafft_path = PROJECT_ROOT / "intermediate" / "06_alignments" / slug / f"{slug}.mafft.aln.faa"
        clip_path = PROJECT_ROOT / "intermediate" / "06_alignments" / slug / f"{slug}.clipkit.aln.faa"
        raw_clip_path = PROJECT_ROOT / "intermediate" / "06_alignments" / slug / f"{slug}.clipkit.raw.aln.faa"
        active_lookup_path = PROJECT_ROOT / "intermediate" / "06_alignments" / slug / "active_tree_leaf_lookup.tsv"
        removed_path = PROJECT_ROOT / "intermediate" / "06_alignments" / slug / "clipkit_removed_sequences.tsv"
        tree_path = PROJECT_ROOT / "intermediate" / "07_iqtree" / slug / f"{slug}.treefile"
        source_ids = set(fasta_dict(source))
        mafft_stats = alignment_stats(mafft_path)
        clip_stats = alignment_stats(clip_path)
        raw_clip_stats = alignment_stats(raw_clip_path)
        active_ids = {row["leaf_id"] for row in read_tsv(active_lookup_path)}
        removed_ids = {row["leaf_id"] for row in read_tsv(removed_path)}
        check(f"06_{slug}_status", mafft[og]["status"] == "PASS", mafft[og]["status"], "PASS")
        check(f"06_{slug}_leaves", mafft_stats["ids"] == source_ids, len(mafft_stats["ids"]), len(source_ids))
        check(f"07_{slug}_status", clipkit[og]["status"] == "PASS", clipkit[og]["status"], "PASS")
        check(f"07_{slug}_raw_leaves", raw_clip_stats["ids"] == source_ids, len(raw_clip_stats["ids"]), len(source_ids))
        check(f"07_{slug}_active_leaves", clip_stats["ids"] == active_ids == source_ids - removed_ids and not (active_ids & removed_ids), len(active_ids), len(source_ids) - len(removed_ids))
        check(f"07_{slug}_removed_all_gap", all(raw_clip_stats["non_gap_by_id"][leaf] == 0 for leaf in removed_ids), len(removed_ids), "all non_gap=0")
        check(f"07_{slug}_sites", 0 < clip_stats["alignment_length"] <= mafft_stats["alignment_length"], clip_stats["alignment_length"], f"1..{mafft_stats['alignment_length']}")
        check(f"08_{slug}_status", iqtree[og]["status"] == "PASS", iqtree[og]["status"], "PASS")
        check(f"08_{slug}_tree_leaves", tree_tips(tree_path) == active_ids, len(tree_tips(tree_path)), len(active_ids))
        check(f"08_{slug}_model", bool(iqtree[og]["selected_model"]), iqtree[og]["selected_model"], "nonempty")

    input_base = PROJECT_ROOT / "intermediate" / "08_reconciliation_inputs"
    input_root = input_base / scope_label if scope_label else input_base
    family_manifest = read_tsv(input_root / "generax" / "family_manifest.tsv")
    check("09_family_set", {row["orthogroup"] for row in family_manifest} == expected, len(family_manifest), len(expected))
    split_rows = read_tsv(input_root / "iqtree_split_support.tsv")
    split_counts = Counter(row["orthogroup"] for row in split_rows)
    for og in ogs:
        n_leaves = int(iqtree[og]["n_sequences"])
        check(f"09_{og_slug(og)}_split_count", split_counts[og] == n_leaves - 3, split_counts[og], n_leaves - 3)
        og_splits = [row for row in split_rows if row["orthogroup"] == og]
        support_ok = all(
            (row["ufboot_status"] == "reported" and 0 <= float(row["ufboot"]) <= 100)
            or (row["ufboot_status"] == "missing_support" and row["ufboot"] == "NA")
            for row in og_splits
        )
        check(f"09_{og_slug(og)}_support_contract", support_ok, "valid" if og_splits else "empty", "reported numeric or NA+missing_support")
        check(f"09_{og_slug(og)}_branch_length", all(float(row["branch_length"]) >= 0 for row in og_splits), "valid" if og_splits else "empty", "all >=0")
    species_tree_version = str(config["project"]["species_tree_version"])
    species_tree = input_root / "common" / f"{species_tree_version}.generax.nwk"
    species = Phylo.read(species_tree, "newick")
    expected_names = {
        row["canonical_node_name"] for row in read_tsv(PROJECT_ROOT / "intermediate" / "01_species_tree" / "species_tree_node_contract.resolved.tsv")
        if row["node_kind"] == "internal"
    }
    observed_names = {clade.name for clade in species.get_nonterminals()}
    check("09_species_internal_names", observed_names == expected_names, len(observed_names), len(expected_names))
    expected_root = config["reconciliation"]["species_tree_root"]
    check("09_species_root_name", species.root.name == expected_root, species.root.name, expected_root)
    file_manifest = read_tsv(input_root / "reconciliation_input_files.tsv")
    bad_hashes = [row["path"] for row in file_manifest if not Path(row["path"]).is_file() or sha256_file(row["path"]) != row["sha256"]]
    check("09_input_manifest_hashes", not bad_hashes, len(bad_hashes), 0, ";".join(bad_hashes[:3]))

    result_base = PROJECT_ROOT / "results" / "02_generax"
    result_root = result_base / scope_label if scope_label else result_base
    run_rows = read_tsv(result_root / "generax_run_summary.tsv")
    family_rows = read_tsv(result_root / "generax_family_summary.tsv")
    check("10_single_run", len(run_rows) == 1, len(run_rows), 1)
    check("10_run_status", len(run_rows) == 1 and run_rows[0]["status"] == "PASS", run_rows[0]["status"] if run_rows else "missing", "PASS")
    check("10_family_set", {row["orthogroup"] for row in family_rows} == expected, len(family_rows), len(expected))
    check("10_family_status", all(row["status"] == "PASS" for row in family_rows), sum(row["status"] != "PASS" for row in family_rows), 0)
    check("10_topology_fixed", all(parse_bool(row["unrooted_split_set_identical_to_iqtree"]) for row in family_rows), sum(not parse_bool(row["unrooted_split_set_identical_to_iqtree"]) for row in family_rows), 0)
    check("10_no_transfer", all(int(row["transfer_events"]) == 0 and int(row["transfer_loss_events"]) == 0 for row in family_rows), "all zero" if family_rows else "empty", "all zero")
    check("10_eval_ll_not_evidence", len(run_rows) == 1 and not parse_bool(run_rows[0]["likelihood_evidence"]), run_rows[0]["likelihood_evidence"] if run_rows else "missing", "FALSE")

    suffix = f".{scope_label}" if scope_label else ""
    output_root = PROJECT_ROOT / "results" / "00_input_audit"
    out_tsv = output_root / f"steps_06_10_validation{suffix}.tsv"
    out_json = output_root / f"steps_06_10_validation{suffix}.json"
    write_tsv(out_tsv, checks, ["check_id", "status", "observed", "expected", "detail"])
    failures = [row for row in checks if row["status"] != "PASS"]
    write_json(out_json, {
        "status": "PASS" if not failures else "FAIL", "scope_label": scope_label,
        "requested_ogs": ogs, "n_checks": len(checks), "n_failures": len(failures),
        "failed_checks": [row["check_id"] for row in failures],
    })
    print(f"Steps 06-10 validation: {'PASS' if not failures else 'FAIL'} ({len(checks)} checks; {len(failures)} failures)")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
