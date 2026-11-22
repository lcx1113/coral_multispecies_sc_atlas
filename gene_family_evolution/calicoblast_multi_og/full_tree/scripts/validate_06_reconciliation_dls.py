#!/usr/bin/env python3
"""Independently validate Step 11 node/event/DLS tables and IQ-TREE evidence joins."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from Bio import Phylo

from common import PROJECT_ROOT, load_config, read_tsv, sha256_file, write_json, write_tsv


SUPPORT_SCOPE = "iqtree_unrooted_bipartition_not_reconciliation_event"


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def leafset_hash(leaves: tuple[str, ...]) -> str:
    return sha256_text("\n".join(leaves))


def stable_node_id(orthogroup: str, leaves: tuple[str, ...]) -> str:
    return "GN_" + sha256_text(orthogroup + "\n" + "\n".join(leaves))[:20]


def stable_event_id(orthogroup: str, context: str) -> str:
    return "EV_" + sha256_text(orthogroup + "\n" + context)[:20]


def split_species_locations(path: Path) -> dict[str, tuple[bool, tuple[str, ...]]]:
    tree = Phylo.read(path, "newick")
    return {
        clade.name: (clade.is_terminal(), tuple(sorted(tip.name for tip in clade.get_terminals())))
        for clade in tree.find_clades(order="preorder")
    }


def xml_path_key(value: str) -> tuple[int, ...]:
    if value == "root":
        return ()
    return tuple(int(part) for part in value.removeprefix("root.").split("."))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/analysis_config.yaml")
    parser.add_argument("--scope-label", default="")
    args = parser.parse_args()
    config, _ = load_config(args.config)
    output_base = PROJECT_ROOT / "results" / "03_reconciliation"
    output_root = output_base / args.scope_label if args.scope_label else output_base
    nodes = read_tsv(output_root / "reconciliation_nodes.tsv")
    events = read_tsv(output_root / "reconciliation_events.tsv")
    dl_rows = read_tsv(output_root / "dls_events.tsv")
    loss_edges = read_tsv(output_root / "loss_edges.tsv")
    summaries = read_tsv(output_root / "reconciliation_summary.tsv")
    inputs = read_tsv(output_root / "step11_input_files.tsv")
    outputs = read_tsv(output_root / "step11_output_files.tsv")
    generax_base = PROJECT_ROOT / "results" / "02_generax"
    generax_root = generax_base / args.scope_label if args.scope_label else generax_base
    step10 = read_tsv(generax_root / "generax_family_summary.tsv")
    input_base = PROJECT_ROOT / "intermediate" / "08_reconciliation_inputs"
    input_root = input_base / args.scope_label if args.scope_label else input_base
    splits = read_tsv(input_root / "iqtree_split_support.tsv")
    species_tree_version = str(config["project"]["species_tree_version"])
    species_tree = input_root / "common" / f"{species_tree_version}.generax.nwk"
    species_locations = split_species_locations(species_tree)
    expected_ogs = [row["orthogroup"] for row in step10 if row["status"] == "PASS"]
    checks = []

    def check(orthogroup: str, name: str, passed: bool, details: str) -> None:
        checks.append({
            "orthogroup": orthogroup,
            "check": name,
            "status": "PASS" if passed else "ERROR",
            "details": details,
        })

    check("ALL", "config_eval_only", config["reconciliation"].get("method") == "GeneRax_EVAL_only" and config["reconciliation"].get("strategy") == "EVAL", "GeneRax EVAL only")
    check("ALL", "exact_family_set", {row["orthogroup"] for row in summaries} == set(expected_ogs) and len(summaries) == len(expected_ogs), f"observed={len(summaries)} expected={len(expected_ogs)}")
    check("ALL", "global_node_ids_unique", len({row["node_id"] for row in nodes}) == len(nodes), f"n={len(nodes)}")
    check("ALL", "global_event_ids_unique", len({row["event_id"] for row in events}) == len(events), f"n={len(events)}")
    check("ALL", "all_outputs_uncompressed_tsv", all(not str(path).endswith(".gz") for path in output_root.glob("*.tsv")), "human-review reconciliation tables")

    input_hashes_ok = all(Path(row["path"]).is_file() and sha256_file(row["path"]) == row["sha256"] for row in inputs)
    output_hashes_ok = all(Path(row["path"]).is_file() and sha256_file(row["path"]) == row["sha256"] and Path(row["path"]).stat().st_size == int(row["size_bytes"]) for row in outputs)
    check("ALL", "input_fingerprints_exact", input_hashes_ok, f"n={len(inputs)}")
    check("ALL", "output_fingerprints_exact", output_hashes_ok, f"n={len(outputs)}")
    run_summary_path = output_root / "step11_run_summary.json"
    run_summary = json.loads(run_summary_path.read_text()) if run_summary_path.is_file() else {}
    check("ALL", "run_summary_pass", run_summary.get("status") == "PASS" and run_summary.get("parser_version") == "1.1", str(run_summary.get("status")))
    check("ALL", "run_summary_counts", run_summary.get("n_gene_nodes") == len(nodes) and run_summary.get("n_events") == len(events) and run_summary.get("n_duplication_loss_rows") == len(dl_rows), f"nodes={len(nodes)} events={len(events)} dl={len(dl_rows)}")

    step10_by_og = {row["orthogroup"]: row for row in step10 if row["status"] == "PASS"}
    summary_by_og = {row["orthogroup"]: row for row in summaries}
    split_by_og: defaultdict[str, dict[str, dict]] = defaultdict(dict)
    for row in splits:
        split_by_og[row["orthogroup"]][row["split_id"]] = row
    nodes_by_og: defaultdict[str, list[dict]] = defaultdict(list)
    events_by_og: defaultdict[str, list[dict]] = defaultdict(list)
    dl_by_og: defaultdict[str, list[dict]] = defaultdict(list)
    for row in nodes:
        nodes_by_og[row["orthogroup"]].append(row)
    for row in events:
        events_by_og[row["orthogroup"]].append(row)
    for row in dl_rows:
        dl_by_og[row["orthogroup"]].append(row)

    for orthogroup in expected_ogs:
        og_nodes = nodes_by_og[orthogroup]
        og_events = events_by_og[orthogroup]
        node_lookup = {row["node_id"]: row for row in og_nodes}
        event_lookup = {row["event_id"]: row for row in og_events}
        roots = [row for row in og_nodes if row["is_root"] == "TRUE"]
        leaves = [row for row in og_nodes if row["node_kind"] == "leaf"]
        internals = [row for row in og_nodes if row["node_kind"] == "internal"]
        n_leaves = len(leaves)
        check(orthogroup, "rooted_binary_node_counts", len(roots) == 1 and len(internals) == n_leaves - 1 and len(og_nodes) == 2 * n_leaves - 1, f"root={len(roots)} internal={len(internals)} leaf={n_leaves}")

        stable_ids_ok = True
        graph_ok = True
        descendants_ok = True
        location_ok = True
        for row in og_nodes:
            descendant_leaves = tuple(row["descendant_leaves"].split(";")) if row["descendant_leaves"] else ()
            stable_ids_ok &= row["node_id"] == stable_node_id(orthogroup, descendant_leaves) and row["descendant_leaf_hash"] == leafset_hash(descendant_leaves)
            expected_location = species_locations.get(row["species_location"])
            location_ok &= expected_location is not None
            if expected_location is not None:
                location_ok &= row["species_location_is_tip"] == ("TRUE" if expected_location[0] else "FALSE")
                location_ok &= row["species_location_descendants"] == ";".join(expected_location[1])
                location_ok &= int(row["n_species_location_descendants"]) == len(expected_location[1])
            children = row["child_node_ids"].split(";") if row["child_node_ids"] else []
            graph_ok &= all(child in node_lookup and node_lookup[child]["parent_node_id"] == row["node_id"] for child in children)
            graph_ok &= int(row["n_children"]) == len(children)
            if row["node_kind"] == "leaf":
                graph_ok &= len(children) == 0 and row["leaf_id"] == descendant_leaves[0]
            else:
                graph_ok &= len(children) == 2
                union = tuple(sorted({leaf for child in children for leaf in node_lookup[child]["descendant_leaves"].split(";")}))
                descendants_ok &= union == descendant_leaves
        check(orthogroup, "stable_node_ids_and_leaf_hashes", stable_ids_ok, f"n={len(og_nodes)}")
        check(orthogroup, "parent_child_graph_exact", graph_ok and roots[0]["parent_node_id"] == "", f"n={len(og_nodes)}")
        check(orthogroup, "descendant_sets_reconstruct", descendants_ok, "internal set = union(children)")
        check(orthogroup, "species_locations_exact", location_ok, f"species-tree nodes={len(species_locations)}")

        leaf_map_path = input_root / step10_by_og[orthogroup]["family"] / "mappings" / f"{step10_by_og[orthogroup]['family']}.leaf_species.tsv"
        leaf_metadata = {row["leaf_id"]: row for row in read_tsv(leaf_map_path)}
        leaf_metadata_ok = all(
            row["leaf_id"] in leaf_metadata
            and row["leaf_species"] == leaf_metadata[row["leaf_id"]]["species"]
            and row["protein_id"] == leaf_metadata[row["leaf_id"]]["protein_id"]
            and row["canonical_gene_id"] == leaf_metadata[row["leaf_id"]]["canonical_gene_id"]
            for row in leaves
        ) and {row["leaf_id"] for row in leaves} == set(leaf_metadata)
        check(orthogroup, "leaf_metadata_exact", leaf_metadata_ok, f"n={len(leaves)}")

        mapped_nodes = [row for row in internals if row["split_id"]]
        multiplicity = Counter(row["split_id"] for row in mapped_nodes)
        root_edge_nodes = [row for row in mapped_nodes if row["split_mapping_status"] == "matched_root_edge_pair"]
        trivial_root_edge_nodes = [
            row for row in internals
            if row["split_mapping_status"] == "not_applicable_trivial_root_edge"
        ]
        support_ok = True
        for row in og_nodes:
            if row["node_kind"] == "leaf":
                support_ok &= row["split_mapping_status"] == "not_applicable_leaf" and not row["split_id"] and row["iqtree_ufboot"] == "NA" and row["ufboot_status"] == "not_applicable"
            elif row["is_root"] == "TRUE":
                support_ok &= row["split_mapping_status"] == "not_applicable_root" and not row["split_id"] and row["iqtree_ufboot"] == "NA" and row["ufboot_status"] == "not_applicable"
            elif row["split_mapping_status"] == "not_applicable_trivial_root_edge":
                support_ok &= (
                    not row["split_id"]
                    and row["iqtree_ufboot"] == "NA"
                    and row["ufboot_status"] == "not_applicable"
                    and row["support_scope"] == "not_applicable_trivial_root_edge"
                    and row["parent_node_id"] == roots[0]["node_id"]
                    and int(row["n_descendant_leaves"]) == n_leaves - 1
                )
            else:
                support = split_by_og[orthogroup].get(row["split_id"])
                support_ok &= support is not None and row["support_scope"] == SUPPORT_SCOPE and row["support_is_event_confidence"] == "FALSE"
                if support is not None:
                    support_ok &= row["iqtree_ufboot"] == support["ufboot"] and row["ufboot_status"] == support["ufboot_status"]
                    support_ok &= (
                        (row["ufboot_status"] == "reported" and 0 <= float(row["iqtree_ufboot"]) <= 100)
                        or (row["ufboot_status"] == "missing_support" and row["iqtree_ufboot"] == "NA")
                    )
                    support_ok &= float(row["iqtree_split_branch_length"]) == float(support["branch_length"])
                    support_ok &= row["iqtree_source_tree_sha256"] == support["source_tree_sha256"]
        internal_root_edge_ok = (
            len(root_edge_nodes) == 2
            and not trivial_root_edge_nodes
            and len({row["split_id"] for row in root_edge_nodes}) == 1
            and all(row["parent_node_id"] == roots[0]["node_id"] for row in root_edge_nodes)
            and Counter(multiplicity.values())[2] == 1
            and Counter(multiplicity.values())[1] == len(split_by_og[orthogroup]) - 1
            and set(multiplicity) == set(split_by_og[orthogroup])
        )
        root_children = roots[0]["child_node_ids"].split(";")
        terminal_root_edge_ok = (
            not root_edge_nodes
            and len(trivial_root_edge_nodes) == 1
            and sum(node_lookup[child]["node_kind"] == "leaf" for child in root_children) == 1
            and Counter(multiplicity.values())[2] == 0
            and Counter(multiplicity.values())[1] == len(split_by_og[orthogroup])
            and set(multiplicity) == set(split_by_og[orthogroup])
        )
        check(orthogroup, "iqtree_support_values_exact", support_ok, f"mapped_nodes={len(mapped_nodes)}")
        check(
            orthogroup, "root_edge_support_mapping_explicit",
            internal_root_edge_ok or terminal_root_edge_ok,
            f"distinct={len(multiplicity)} pair_nodes={len(root_edge_nodes)} trivial_root_nodes={len(trivial_root_edge_nodes)}",
        )

        primary_events = [row for row in og_events if row["is_primary_gene_node_event"] == "TRUE"]
        primary_by_node = {row["node_id"]: row for row in primary_events}
        primary_ok = len(primary_events) == len(og_nodes) and set(primary_by_node) == set(node_lookup)
        primary_ids_ok = True
        for node_id, node in node_lookup.items():
            event = primary_by_node.get(node_id)
            if event is None:
                primary_ok = False
                continue
            primary_ok &= event["event_type"] == node["event_type"] and event["species_location"] == node["species_location"]
            primary_ids_ok &= event["event_id"] == stable_event_id(orthogroup, f"PRIMARY|{node_id}|{event['event_type']}")
        check(orthogroup, "one_primary_xml_event_per_gene_node", primary_ok and primary_ids_ok, f"events={len(primary_events)} nodes={len(og_nodes)}")

        auxiliary = [row for row in og_events if row["event_type"] in {"SL", "L"}]
        auxiliary_ok = True
        for event in auxiliary:
            auxiliary_ok &= event["is_auxiliary_reconciliation_event"] == "TRUE" and event["node_id"] == ""
            auxiliary_ok &= event["edge_parent_node_id"] in node_lookup and event["edge_child_node_id"] in node_lookup
            if event["edge_child_node_id"] in node_lookup and event["event_type"] == "SL":
                auxiliary_ok &= event["descendant_extant_leaves"] == node_lookup[event["edge_child_node_id"]]["descendant_leaves"]
            if event["event_type"] == "L":
                sl = event_lookup.get(event["associated_sl_event_id"])
                auxiliary_ok &= sl is not None and sl["event_type"] == "SL"
                if sl is not None:
                    auxiliary_ok &= event["edge_parent_node_id"] == sl["edge_parent_node_id"] and event["edge_child_node_id"] == sl["edge_child_node_id"]
                auxiliary_ok &= event["n_descendant_extant_leaves"] == "0"
        check(orthogroup, "auxiliary_sl_loss_edge_mapping", auxiliary_ok, f"n={len(auxiliary)}")

        auxiliary_ids_ok = True
        sl_groups: defaultdict[tuple[str, str, str], list[dict]] = defaultdict(list)
        for event in auxiliary:
            if event["event_type"] == "SL":
                sl_groups[(event["edge_parent_node_id"], event["edge_child_node_id"], event["species_location"])].append(event)
        for group, values in sl_groups.items():
            for ordinal, event in enumerate(sorted(values, key=lambda row: xml_path_key(row["xml_path"])), 1):
                expected_id = stable_event_id(orthogroup, "SL|" + "|".join(group) + f"|{ordinal}")
                auxiliary_ids_ok &= event["event_id"] == expected_id
        losses_by_sl: defaultdict[str, list[dict]] = defaultdict(list)
        for event in auxiliary:
            if event["event_type"] == "L":
                losses_by_sl[event["associated_sl_event_id"]].append(event)
        for sl_event_id, values in losses_by_sl.items():
            for ordinal, event in enumerate(sorted(values, key=lambda row: xml_path_key(row["xml_path"])), 1):
                expected_id = stable_event_id(orthogroup, f"L|{sl_event_id}|{event['species_location']}|{ordinal}")
                auxiliary_ids_ok &= event["event_id"] == expected_id
        check(orthogroup, "stable_auxiliary_event_ids", auxiliary_ids_ok, f"SL={sum(len(v) for v in sl_groups.values())} L={sum(len(v) for v in losses_by_sl.values())}")

        event_counts = Counter(row["event_type"] for row in og_events)
        step10_row = step10_by_og[orthogroup]
        expected_counts = {
            "S": int(step10_row["speciation_events"]),
            "SL": int(step10_row["speciation_loss_events"]),
            "D": int(step10_row["duplication_events"]),
            "L": int(step10_row["speciation_loss_events"]) + int(step10_row["standalone_loss_events"]),
            "Leaf": int(step10_row["leaf_events"]),
        }
        check(orthogroup, "event_counts_match_step10", all(event_counts[key] == value for key, value in expected_counts.items()), str(dict(event_counts)))

        expected_dl_ids = {row["event_id"] for row in og_events if row["event_type"] in {"D", "L"}}
        observed_dl_ids = {row["event_id"] for row in dl_by_og[orthogroup]}
        dl_ok = expected_dl_ids == observed_dl_ids and len(observed_dl_ids) == len(dl_by_og[orthogroup])
        for row in dl_by_og[orthogroup]:
            source = event_lookup[row["event_id"]]
            evidence_id = source["node_id"] if source["event_type"] == "D" else source["edge_child_node_id"]
            evidence = node_lookup[evidence_id]
            dl_ok &= row["evidence_node_id"] == evidence_id and row["split_id"] == evidence["split_id"]
            dl_ok &= row["iqtree_ufboot"] == evidence["iqtree_ufboot"] and row["ufboot_status"] == evidence["ufboot_status"] and row["support_is_event_confidence"] == "FALSE"
        check(orthogroup, "duplication_loss_subset_exact", dl_ok, f"n={len(observed_dl_ids)}")

        summary = summary_by_og[orthogroup]
        summary_ok = (
            int(summary["n_leaves"]) == n_leaves
            and int(summary["n_internal_nodes"]) == len(internals)
            and int(summary["n_gene_nodes"]) == len(og_nodes)
            and int(summary["n_speciation_events"]) == event_counts["S"]
            and int(summary["n_duplication_events"]) == event_counts["D"]
            and int(summary["n_speciation_loss_events"]) == event_counts["SL"]
            and int(summary["n_loss_events"]) == event_counts["L"]
            and summary["root_node_id"] == roots[0]["node_id"]
            and summary["root_species_location"] == roots[0]["species_location"]
            and summary["status"] == "PASS"
        )
        check(orthogroup, "family_summary_exact", summary_ok, str(summary))

    validation_dir = output_root / "validation"
    loss_event_ids = {row["event_id"] for row in events if row["event_type"] == "L"}
    check("ALL", "loss_edge_table_exact", {row["event_id"] for row in loss_edges} == loss_event_ids and len(loss_edges) == len(loss_event_ids), f"n={len(loss_edges)}")
    validation_path = validation_dir / "step11_validation.tsv"
    write_tsv(validation_path, checks, ["orthogroup", "check", "status", "details"])
    n_errors = sum(row["status"] == "ERROR" for row in checks)
    write_json(
        validation_dir / "step11_validation.summary.json",
        {"status": "PASS" if n_errors == 0 else "FAIL", "n_checks": len(checks), "n_errors": n_errors},
    )
    print(f"Step 11 validation: {'PASS' if n_errors == 0 else 'FAIL'} ({len(checks)} checks, {n_errors} errors)")
    return 1 if n_errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
