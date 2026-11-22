#!/usr/bin/env python3
"""Step 09: prepare GeneRax EVAL inputs and canonical IQ-TREE split evidence."""

from __future__ import annotations

import argparse
import hashlib
import re
import shutil
from collections import defaultdict
from pathlib import Path

from Bio import Phylo

from common import (
    PROJECT_ROOT, alignment_stats, ensure_dir, fasta_dict, load_config, og_slug,
    parse_bool, read_tsv, selected_ogs, sha256_file, write_json, write_tsv,
)


TOOL_FORBIDDEN_LABEL_CHARS = set(" \t\r\n,():;[]'")
SPLIT_FIELDS = [
    "orthogroup", "split_id", "side_a_hash", "side_b_hash", "n_side_a", "n_side_b",
    "side_a_leaves", "side_b_leaves", "ufboot", "ufboot_status", "branch_length", "source_tree",
    "source_tree_sha256",
]


def tree_tips(tree) -> list[str]:
    names = [leaf.name for leaf in tree.get_terminals()]
    if None in names or len(names) != len(set(names)):
        raise ValueError("Tree has missing or duplicate tip labels")
    return names


def leaf_hash(names: tuple[str, ...]) -> str:
    return hashlib.sha256("\n".join(names).encode("utf-8")).hexdigest()


def canonical_split(descendants: set[str], all_leaves: set[str]) -> tuple[tuple[str, ...], tuple[str, ...]]:
    complement = all_leaves - descendants
    if not descendants or not complement or descendants & complement or descendants | complement != all_leaves:
        raise ValueError("Invalid bipartition")
    sides = [tuple(sorted(descendants)), tuple(sorted(complement))]
    sides.sort(key=lambda side: (len(side), side))
    return sides[0], sides[1]


def split_identifier(side_a: tuple[str, ...], side_b: tuple[str, ...]) -> str:
    return hashlib.sha256(("\n".join(side_a) + "\n||\n" + "\n".join(side_b)).encode("utf-8")).hexdigest()


def iqtree_split_rows(path: Path, orthogroup: str, support_config: dict) -> list[dict]:
    tree = Phylo.read(path, "newick")
    leaves = set(tree_tips(tree))
    if len(tree.root.clades) != 3:
        raise ValueError(f"Expected IQ-TREE unrooted display root degree 3: {path}")
    if any(len(clade.clades) not in (0, 2) for clade in tree.find_clades() if clade is not tree.root):
        raise ValueError(f"Expected strictly binary unrooted topology: {path}")
    source_sha = sha256_file(path)
    result, seen = [], set()
    for clade in tree.get_nonterminals():
        if clade is tree.root:
            continue
        side_a, side_b = canonical_split({tip.name for tip in clade.get_terminals()}, leaves)
        split_id = split_identifier(side_a, side_b)
        if split_id in seen:
            raise ValueError(f"Duplicate internal split in {path}")
        seen.add(split_id)
        if clade.confidence is None:
            ufboot = str(support_config["missing_ufboot_value"])
            ufboot_status = str(support_config["missing_ufboot_status"])
        else:
            if not 0 <= float(clade.confidence) <= 100:
                raise ValueError(f"Invalid UFBoot on IQ-TREE internal split in {path}")
            ufboot = f"{float(clade.confidence):.12g}"
            ufboot_status = str(support_config["reported_ufboot_status"])
        if clade.branch_length is None or clade.branch_length < 0:
            raise ValueError(f"Missing/invalid branch length on IQ-TREE internal split in {path}")
        result.append({
            "orthogroup": orthogroup, "split_id": split_id,
            "side_a_hash": leaf_hash(side_a), "side_b_hash": leaf_hash(side_b),
            "n_side_a": len(side_a), "n_side_b": len(side_b),
            "side_a_leaves": ";".join(side_a), "side_b_leaves": ";".join(side_b),
            "ufboot": ufboot, "ufboot_status": ufboot_status,
            "branch_length": f"{clade.branch_length:.12g}",
            "source_tree": str(path.resolve()), "source_tree_sha256": source_sha,
        })
    expected = len(leaves) - 3
    if len(result) != expected:
        raise ValueError(f"Expected {expected} internal unrooted splits; observed {len(result)} in {path}")
    return sorted(result, key=lambda row: row["split_id"])


def rooted_clades(tree) -> set[frozenset[str]]:
    all_tips = set(tree_tips(tree))
    return {
        frozenset(tip.name for tip in clade.get_terminals())
        for clade in tree.get_nonterminals()
        if clade is not tree.root and 1 < len(clade.get_terminals()) < len(all_tips)
    }


def validate_canonical_species_tree(path: Path, config: dict) -> tuple[set[str], list[str]]:
    tree = Phylo.read(path, "newick")
    tips = set(tree_tips(tree))
    labels = [clade.name for clade in tree.get_nonterminals()]
    expected_labels = {
        row["canonical_node_name"] for row in read_tsv(PROJECT_ROOT / "intermediate" / "01_species_tree" / "species_tree_node_contract.resolved.tsv")
        if row["node_kind"] == "internal"
    }
    if len(tips) != int(config["preflight"]["expected_tree_species"]):
        raise ValueError("Canonical species-tree tip count differs from frozen expectation")
    if len(labels) != int(config["preflight"]["expected_internal_nodes"]) or set(labels) != expected_labels:
        raise ValueError("Canonical species-tree internal labels differ from the configured tree contract")
    if tree.root.name != config["reconciliation"]["species_tree_root"] or len(tree.root.clades) != 2:
        raise ValueError("Canonical species tree root is not the frozen Cnidaria bifurcation")
    return tips, labels


def write_space_map(path: Path, rows: list[dict]) -> None:
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(f"{row['leaf_id']} {row['species']}\n")


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
    parser.add_argument("--scope-label", default="", help="Optional isolated subdirectory, e.g. smoke_test")
    args = parser.parse_args()
    config, config_path = load_config(args.config)
    if args.scope_label and not re.fullmatch(r"[A-Za-z0-9_.-]+", args.scope_label):
        raise ValueError("Unsafe --scope-label")
    support = config["support"]
    expected_support = {
        "export_iqtree_split_support": True,
        "split_definition": "canonical_unrooted_bipartition",
        "split_hash_algorithm": "sha256",
        "split_leaf_order": "lexical_sorted_stable_leaf_ids",
        "collapse_tree_active_in_reconciliation": False,
        "missing_ufboot_value": "NA",
        "missing_ufboot_status": "missing_support",
        "reported_ufboot_status": "reported",
        "require_branch_length_when_ufboot_missing": True,
    }
    for key, expected in expected_support.items():
        observed = parse_bool(support[key], f"support.{key}") if isinstance(expected, bool) else support[key]
        if observed != expected:
            raise ValueError(f"support.{key} must be {expected!r}")
    reconciliation = config["reconciliation"]
    if reconciliation["method"] != "GeneRax_EVAL_only" or reconciliation["gene_tree_variant"] != "uncollapsed_only":
        raise ValueError("Step 09 requires uncollapsed GeneRax EVAL-only inputs")
    ogs = choose_ogs(args.og)
    base = PROJECT_ROOT / "intermediate" / "08_reconciliation_inputs"
    root = base / args.scope_label if args.scope_label else base
    common_dir = root / "common"
    species_tree_version = str(config["project"]["species_tree_version"])
    source_species_tree = PROJECT_ROOT / "intermediate" / "01_species_tree" / f"{species_tree_version}.nwk"
    species_tips, species_labels = validate_canonical_species_tree(source_species_tree, config)
    species_tree = common_dir / f"{species_tree_version}.generax.nwk"
    ensure_dir(common_dir)
    shutil.copyfile(source_species_tree, species_tree)
    copied = Phylo.read(species_tree, "newick")
    if set(tree_tips(copied)) != species_tips or rooted_clades(copied) != rooted_clades(Phylo.read(source_species_tree, "newick")):
        raise ValueError("Species-tree topology changed during GeneRax input copy")

    summaries, split_rows, blocks, file_manifest = [], [], [], []
    for og in ogs:
        slug = og_slug(og)
        source_tree = PROJECT_ROOT / "intermediate" / "07_iqtree" / slug / f"{slug}.treefile"
        alignment = PROJECT_ROOT / "intermediate" / "06_alignments" / slug / f"{slug}.clipkit.aln.faa"
        og_dir = root / slug
        tree_copy = og_dir / "trees" / f"{slug}.ml_unrooted.nwk"
        ensure_dir(tree_copy.parent)
        shutil.copyfile(source_tree, tree_copy)
        active_lookup = PROJECT_ROOT / "intermediate" / "06_alignments" / slug / "active_tree_leaf_lookup.tsv"
        rows = sorted(read_tsv(active_lookup), key=lambda row: row["leaf_id"])
        tree_ids = set(tree_tips(Phylo.read(tree_copy, "newick")))
        alignment_ids = set(fasta_dict(alignment))
        lookup_ids = {row["leaf_id"] for row in rows}
        if tree_ids != alignment_ids or tree_ids != lookup_ids:
            raise ValueError(f"Leaf mismatch before reconciliation input preparation for {og}")
        if any(set(row["leaf_id"]) & TOOL_FORBIDDEN_LABEL_CHARS for row in rows):
            raise ValueError(f"Gene labels contain GeneRax-forbidden characters for {og}")
        mapped_species = {row["species"] for row in rows}
        if not mapped_species <= species_tips:
            raise ValueError(f"Mapped species absent from canonical species tree for {og}")
        current_splits = iqtree_split_rows(tree_copy, og, support)
        split_rows.extend(current_splits)
        human_map = og_dir / "mappings" / f"{slug}.leaf_species.tsv"
        generax_map = og_dir / "mappings" / f"{slug}.generax.map"
        human_rows = [{
            "orthogroup": og, "leaf_id": row["leaf_id"], "species": row["species"],
            "protein_id": row["protein_id"], "transcript_id": row["transcript_id"],
            "canonical_gene_id": row["canonical_gene_id"], "is_orphan": row["is_orphan"],
            "is_multi_og": row["is_multi_og"],
        } for row in rows]
        write_tsv(human_map, human_rows, [
            "orthogroup", "leaf_id", "species", "protein_id", "transcript_id",
            "canonical_gene_id", "is_orphan", "is_multi_og",
        ])
        write_space_map(generax_map, rows)
        per_og_iq_summary = PROJECT_ROOT / "intermediate" / "07_iqtree" / slug / "iqtree_summary.tsv"
        if per_og_iq_summary.is_file():
            iq_rows = read_tsv(per_og_iq_summary)
        else:
            iq_rows = [
                row for row in read_tsv(PROJECT_ROOT / "intermediate" / "07_iqtree" / "iqtree_summary.tsv")
                if row["orthogroup"] == og
            ]
        iq_row = iq_rows[0] if len(iq_rows) == 1 else None
        if not iq_row or iq_row["status"] != "PASS" or not iq_row["selected_model"]:
            raise ValueError(f"Missing valid IQ-TREE summary/model for {og}")
        block = {
            "orthogroup": og, "family": slug, "starting_gene_tree": str(tree_copy.resolve()),
            "alignment": str(alignment.resolve()), "mapping": str(generax_map.resolve()),
            "subst_model": iq_row["selected_model"],
        }
        blocks.append(block)
        for role, path in [
            ("source_iqtree_treefile", source_tree), ("primary_uncollapsed_tree", tree_copy),
            ("trimmed_alignment", alignment), ("human_leaf_species_map", human_map),
            ("generax_map", generax_map), ("active_tree_leaf_lookup", active_lookup),
        ]:
            file_manifest.append({"orthogroup": og, "file_role": role, "path": str(path.resolve()), "sha256": sha256_file(path)})
        stats = alignment_stats(alignment)
        reported_support = [float(row["ufboot"]) for row in current_splits if row["ufboot_status"] == "reported"]
        summaries.append({
            "orthogroup": og, "family": slug, "n_leaves": len(tree_ids),
            "n_species": len(mapped_species), "alignment_sites": stats["alignment_length"],
            "iqtree_model": iq_row["selected_model"],
            "n_internal_unrooted_splits": len(current_splits),
            "n_splits_with_support": len(reported_support),
            "n_splits_missing_support": len(current_splits) - len(reported_support),
            "minimum_reported_ufboot": f"{min(reported_support):.12g}" if reported_support else "NA",
            "primary_tree": str(tree_copy.resolve()), "generax_mapping": str(generax_map.resolve()),
            "species_tree_version": config["project"]["species_tree_version"],
            "analysis_role": str(config["reconciliation"].get("analysis_role", "GENERAX_EVAL_PRIMARY")), "status": "PASS",
        })

    split_path = root / "iqtree_split_support.tsv"
    write_tsv(split_path, split_rows, SPLIT_FIELDS)
    generax_dir = root / "generax"
    ensure_dir(generax_dir)
    families_path = generax_dir / "families.txt"
    with families_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("[FAMILIES]\n")
        for block in blocks:
            handle.write(f"- {block['family']}\n")
            handle.write(f"starting_gene_tree = {block['starting_gene_tree']}\n")
            handle.write(f"alignment = {block['alignment']}\n")
            handle.write(f"mapping = {block['mapping']}\n")
            handle.write(f"subst_model = {block['subst_model']}\n")
    family_manifest = generax_dir / "family_manifest.tsv"
    write_tsv(family_manifest, blocks, ["orthogroup", "family", "starting_gene_tree", "alignment", "mapping", "subst_model"])
    for role, path in [
        ("canonical_species_tree", species_tree), ("iqtree_split_support", split_path),
        ("generax_families", families_path), ("generax_family_manifest", family_manifest),
    ]:
        file_manifest.append({"orthogroup": "ALL", "file_role": role, "path": str(path.resolve()), "sha256": sha256_file(path)})
    write_tsv(root / "reconciliation_input_files.tsv", file_manifest, ["orthogroup", "file_role", "path", "sha256"])
    write_tsv(root / "reconciliation_input_summary.tsv", summaries, [
        "orthogroup", "family", "n_leaves", "n_species", "alignment_sites", "iqtree_model",
        "n_internal_unrooted_splits", "n_splits_with_support", "n_splits_missing_support",
        "minimum_reported_ufboot", "primary_tree", "generax_mapping",
        "species_tree_version", "analysis_role", "status",
    ])
    write_json(root / "reconciliation_input_summary.json", {
        "status": "PASS", "scope_label": args.scope_label, "requested_ogs": ogs,
        "n_orthogroups": len(summaries), "n_iqtree_internal_splits": len(split_rows),
        "species_tree": str(species_tree.resolve()), "species_tree_sha256": sha256_file(species_tree),
        "n_species": len(species_tips), "n_named_internal_species_nodes": len(species_labels),
        "primary_policy": "GeneRax 1.2.2 EVAL on uncollapsed IQ-TREE unrooted topology",
        "support_policy": "canonical unrooted bipartition SHA256 with UFBoot and branch length",
        "config": str(config_path.resolve()), "config_sha256": sha256_file(config_path),
    })
    print(f"Step 09 PASS: {len(blocks)} families; {len(split_rows)} canonical IQ-TREE splits; scope={args.scope_label or 'production'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
