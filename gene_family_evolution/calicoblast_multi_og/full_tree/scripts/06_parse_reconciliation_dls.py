#!/usr/bin/env python3
"""Step 11: parse GeneRax EVAL into stable gene-node, event, and DLS tables."""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter, defaultdict
from pathlib import Path
from xml.etree import ElementTree

from Bio import Phylo

from common import (
    PROJECT_ROOT,
    ensure_dir,
    load_config,
    read_tsv,
    sha256_file,
    write_json,
    write_tsv,
)


PARSER_VERSION = "1.1"
METHOD = "GeneRax_1.2.2_EVAL_UndatedDL"
NODE_ID_SCHEME = "sha256(orthogroup+rooted_descendant_leaf_set),first20"
EVENT_ID_SCHEME = "sha256(orthogroup+event_context),first20"
SUPPORT_SCOPE = "iqtree_unrooted_bipartition_not_reconciliation_event"

NODE_FIELDS = [
    "orthogroup", "family", "node_id", "node_kind", "is_root", "parent_node_id",
    "child_node_ids", "n_children", "event_type", "species_location",
    "species_location_is_tip", "n_species_location_descendants",
    "species_location_descendants", "n_descendant_leaves", "descendant_leaf_hash",
    "descendant_leaves", "n_descendant_species", "descendant_species", "leaf_id",
    "leaf_species", "protein_id", "canonical_gene_id", "reconciled_branch_length",
    "nhx_B", "nhx_embedded_label", "split_id", "split_mapping_status",
    "split_mapping_multiplicity", "iqtree_ufboot", "ufboot_status", "iqtree_split_branch_length",
    "support_scope", "support_is_event_confidence", "iqtree_source_tree",
    "iqtree_source_tree_sha256", "reconciled_nhx", "reconciled_xml",
    "reconciliation_method", "parser_version",
]

EVENT_FIELDS = [
    "orthogroup", "family", "event_id", "event_type", "raw_xml_event",
    "is_primary_gene_node_event", "is_auxiliary_reconciliation_event", "node_id",
    "edge_parent_node_id", "edge_child_node_id", "associated_sl_event_id",
    "species_location", "species_location_is_tip", "n_species_location_descendants",
    "species_location_descendants", "xml_clade_name", "xml_depth", "xml_path",
    "n_descendant_extant_leaves", "descendant_extant_leaf_hash",
    "descendant_extant_leaves", "source_xml", "source_xml_sha256",
    "reconciliation_method", "parser_version",
]

DL_FIELDS = [
    "orthogroup", "family", "event_id", "event_type", "node_id",
    "edge_parent_node_id", "edge_child_node_id", "associated_sl_event_id",
    "species_location", "species_location_is_tip", "n_species_location_descendants",
    "species_location_descendants", "evidence_node_id", "split_id",
    "split_mapping_status", "split_mapping_multiplicity", "iqtree_ufboot", "ufboot_status",
    "iqtree_split_branch_length", "support_scope", "support_is_event_confidence",
    "source_xml", "reconciliation_method",
]

SUMMARY_FIELDS = [
    "orthogroup", "family", "n_leaves", "n_internal_nodes", "n_gene_nodes",
    "n_speciation_events", "n_duplication_events", "n_speciation_loss_events",
    "n_loss_events", "n_support_mapped_nodes", "n_distinct_mapped_splits",
    "n_unique_split_nodes", "n_root_edge_pair_nodes", "root_node_id",
    "root_species_location", "root_edge_split_id", "status",
]

LOSS_EDGE_FIELDS = [
    "orthogroup", "family", "event_id", "associated_sl_event_id",
    "edge_parent_node_id", "edge_child_node_id", "species_location",
    "species_location_descendants", "source_xml", "reconciliation_method",
]


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def leafset_hash(leaves: tuple[str, ...]) -> str:
    return sha256_text("\n".join(leaves))


def stable_node_id(orthogroup: str, leaves: tuple[str, ...]) -> str:
    return "GN_" + sha256_text(orthogroup + "\n" + "\n".join(leaves))[:20]


def stable_event_id(orthogroup: str, context: str) -> str:
    return "EV_" + sha256_text(orthogroup + "\n" + context)[:20]


def parse_nhx_comment(comment: str | None) -> dict[str, str]:
    if not comment or not comment.startswith("&&NHX:"):
        raise ValueError(f"Missing GeneRax NHX annotation: {comment!r}")
    result = {}
    for field in comment[len("&&NHX:") :].split(":"):
        if "=" not in field:
            raise ValueError(f"Malformed NHX field {field!r}")
        key, value = field.split("=", 1)
        result[key] = value
    for key in ["S", "D", "H", "B"]:
        if key not in result:
            raise ValueError(f"Missing NHX {key} in {comment!r}")
    return result


def canonical_split(descendants: set[str], all_leaves: set[str]) -> tuple[tuple[str, ...], tuple[str, ...]]:
    sides = [tuple(sorted(descendants)), tuple(sorted(all_leaves - descendants))]
    sides.sort(key=lambda side: (len(side), side))
    return sides[0], sides[1]


def split_identifier(side_a: tuple[str, ...], side_b: tuple[str, ...]) -> str:
    return sha256_text("\n".join(side_a) + "\n||\n" + "\n".join(side_b))


def xml_tag(element) -> str:
    return element.tag.rsplit("}", 1)[-1]


def direct_children(element, name: str) -> list:
    return [child for child in list(element) if xml_tag(child) == name]


def direct_text(element, name: str) -> str:
    matches = direct_children(element, name)
    return (matches[0].text or "").strip() if matches else ""


def build_species_locations(species_tree: Path) -> dict[str, dict]:
    tree = Phylo.read(species_tree, "newick")
    result = {}
    for clade in tree.find_clades(order="preorder"):
        if not clade.name:
            raise ValueError(f"Unnamed species-tree node in {species_tree}")
        if clade.name in result:
            raise ValueError(f"Duplicate species-tree node name {clade.name}")
        descendants = tuple(sorted(tip.name for tip in clade.get_terminals()))
        result[clade.name] = {
            "is_tip": clade.is_terminal(),
            "descendants": descendants,
        }
    return result


def location_fields(location: str, species_locations: dict[str, dict]) -> dict:
    if location not in species_locations:
        raise ValueError(f"GeneRax speciesLocation absent from species tree: {location}")
    metadata = species_locations[location]
    return {
        "species_location_is_tip": metadata["is_tip"],
        "n_species_location_descendants": len(metadata["descendants"]),
        "species_location_descendants": ";".join(metadata["descendants"]),
    }


def parse_gene_nodes(
    orthogroup: str,
    family: str,
    nhx_path: Path,
    xml_path: Path,
    leaf_metadata: dict[str, dict],
    split_lookup: dict[str, dict],
    species_locations: dict[str, dict],
) -> tuple[list[dict], dict[tuple[str, ...], dict]]:
    tree = Phylo.read(nhx_path, "newick")
    all_leaves = {tip.name for tip in tree.get_terminals()}
    if None in all_leaves or len(all_leaves) != len(tree.get_terminals()):
        raise ValueError(f"Missing or duplicate leaves in {nhx_path}")
    if set(leaf_metadata) != all_leaves:
        raise ValueError(f"Leaf lookup differs from reconciled NHX for {orthogroup}")
    expected_root_id = stable_node_id(orthogroup, tuple(sorted(all_leaves)))

    temporary = []

    def walk(clade, parent_id: str, is_root: bool) -> str:
        leaves = tuple(sorted(tip.name for tip in clade.get_terminals()))
        node_id = stable_node_id(orthogroup, leaves)
        annotations = parse_nhx_comment(clade.comment)
        if annotations["H"] != "N":
            raise ValueError(f"Unexpected transfer annotation at {node_id}")
        if clade.is_terminal():
            event_type = "Leaf"
            if annotations["D"] != "N":
                raise ValueError(f"Leaf marked as duplication at {node_id}")
        else:
            event_type = "D" if annotations["D"] == "Y" else "S"
            if annotations["D"] not in {"Y", "N"}:
                raise ValueError(f"Unexpected NHX D value at {node_id}: {annotations['D']}")
        child_ids = [walk(child, node_id, False) for child in clade.clades]
        descendant_species = tuple(sorted({leaf_metadata[leaf]["species"] for leaf in leaves}))
        leaf_id = leaves[0] if clade.is_terminal() else ""
        metadata = leaf_metadata.get(leaf_id, {})
        label = clade.name if clade.name is not None else clade.confidence
        row = {
            "orthogroup": orthogroup,
            "family": family,
            "node_id": node_id,
            "node_kind": "leaf" if clade.is_terminal() else "internal",
            "is_root": is_root,
            "parent_node_id": parent_id,
            "child_node_ids": ";".join(child_ids),
            "n_children": len(child_ids),
            "event_type": event_type,
            "species_location": annotations["S"],
            **location_fields(annotations["S"], species_locations),
            "n_descendant_leaves": len(leaves),
            "descendant_leaf_hash": leafset_hash(leaves),
            "descendant_leaves": ";".join(leaves),
            "n_descendant_species": len(descendant_species),
            "descendant_species": ";".join(descendant_species),
            "leaf_id": leaf_id,
            "leaf_species": metadata.get("species", ""),
            "protein_id": metadata.get("protein_id", ""),
            "canonical_gene_id": metadata.get("canonical_gene_id", ""),
            "reconciled_branch_length": "" if clade.branch_length is None else f"{clade.branch_length:.12g}",
            "nhx_B": annotations["B"],
            "nhx_embedded_label": "" if label is None else str(label),
            "split_id": "",
            "split_mapping_status": "not_applicable_leaf" if clade.is_terminal() else ("not_applicable_root" if is_root else "pending"),
            "split_mapping_multiplicity": 0,
            "iqtree_ufboot": "NA",
            "ufboot_status": "not_applicable" if clade.is_terminal() or is_root else "missing_support",
            "iqtree_split_branch_length": "",
            "support_scope": "not_applicable" if clade.is_terminal() or is_root else SUPPORT_SCOPE,
            "support_is_event_confidence": False,
            "iqtree_source_tree": "",
            "iqtree_source_tree_sha256": "",
            "reconciled_nhx": str(nhx_path.resolve()),
            "reconciled_xml": str(xml_path.resolve()),
            "reconciliation_method": METHOD,
            "parser_version": PARSER_VERSION,
        }
        if not clade.is_terminal() and not is_root:
            side_a, side_b = canonical_split(set(leaves), all_leaves)
            if len(side_a) == 1:
                if parent_id != expected_root_id or len(leaves) != len(all_leaves) - 1:
                    raise ValueError(
                        f"Unexpected non-root trivial split for {orthogroup} node {node_id}"
                    )
                row.update(
                    {
                        "split_mapping_status": "not_applicable_trivial_root_edge",
                        "ufboot_status": "not_applicable",
                        "support_scope": "not_applicable_trivial_root_edge",
                    }
                )
            else:
                split_id = split_identifier(side_a, side_b)
                if split_id not in split_lookup:
                    raise ValueError(f"No Step 09 split match for {orthogroup} node {node_id}")
                support = split_lookup[split_id]
                row.update(
                    {
                        "split_id": split_id,
                        "iqtree_ufboot": support["ufboot"],
                        "ufboot_status": support["ufboot_status"],
                        "iqtree_split_branch_length": support["branch_length"],
                        "iqtree_source_tree": support["source_tree"],
                        "iqtree_source_tree_sha256": support["source_tree_sha256"],
                    }
                )
        temporary.append(row)
        return node_id

    root_id = walk(tree.root, "", True)
    if len({row["node_id"] for row in temporary}) != len(temporary):
        raise ValueError(f"Stable node ID collision in {orthogroup}")
    split_multiplicity = Counter(row["split_id"] for row in temporary if row["split_id"])
    for row in temporary:
        if not row["split_id"]:
            continue
        multiplicity = split_multiplicity[row["split_id"]]
        row["split_mapping_multiplicity"] = multiplicity
        if multiplicity == 1:
            row["split_mapping_status"] = "matched_unique"
        elif multiplicity == 2 and row["parent_node_id"] == root_id:
            row["split_mapping_status"] = "matched_root_edge_pair"
        else:
            raise ValueError(f"Unexpected split multiplicity {multiplicity} at {row['node_id']}")
    temporary.sort(key=lambda row: (row["orthogroup"], row["node_id"]))
    by_descendants = {
        tuple(row["descendant_leaves"].split(";")): row for row in temporary
    }
    return temporary, by_descendants


def parse_xml_events(
    orthogroup: str,
    family: str,
    xml_path: Path,
    nodes_by_descendants: dict[tuple[str, ...], dict],
    species_locations: dict[str, dict],
) -> list[dict]:
    root = ElementTree.parse(xml_path).getroot()
    rec_gene_trees = [element for element in root.iter() if xml_tag(element) == "recGeneTree"]
    if len(rec_gene_trees) != 1:
        raise ValueError(f"Expected one recGeneTree in {xml_path}")
    phylogenies = direct_children(rec_gene_trees[0], "phylogeny")
    gene_roots = direct_children(phylogenies[0], "clade") if len(phylogenies) == 1 else []
    if len(gene_roots) != 1:
        raise ValueError(f"Expected one reconciled gene root in {xml_path}")

    records: list[dict] = []
    by_path: dict[tuple[int, ...], dict] = {}

    def walk(clade, path: tuple[int, ...]) -> dict:
        event_containers = direct_children(clade, "eventsRec")
        event_elements = list(event_containers[0]) if len(event_containers) == 1 else []
        if len(event_elements) != 1:
            raise ValueError(f"Expected one XML event at path {path} in {xml_path}")
        raw_event = xml_tag(event_elements[0])
        if raw_event not in {"speciation", "duplication", "loss", "leaf"}:
            raise ValueError(f"Unsupported XML event {raw_event} in {xml_path}")
        location = event_elements[0].attrib.get("speciesLocation", "")
        child_records = [walk(child, path + (index,)) for index, child in enumerate(direct_children(clade, "clade"))]
        if raw_event == "leaf":
            leaves = (direct_text(clade, "name"),)
        elif raw_event == "loss":
            leaves = ()
        else:
            leaves = tuple(sorted({leaf for child in child_records for leaf in child["leaves"]}))
        n_gene_children = sum(bool(child["leaves"]) for child in child_records)
        if raw_event == "speciation":
            if n_gene_children >= 2:
                event_type = "S"
            elif n_gene_children == 1 and any(child["raw_xml_event"] == "loss" for child in child_records):
                event_type = "SL"
            else:
                raise ValueError(f"Unclassifiable XML speciation at path {path}")
        else:
            event_type = {"duplication": "D", "loss": "L", "leaf": "Leaf"}[raw_event]
        record = {
            "path_tuple": path,
            "parent_path": path[:-1] if path else None,
            "raw_xml_event": raw_event,
            "event_type": event_type,
            "species_location": location,
            "xml_clade_name": direct_text(clade, "name"),
            "leaves": leaves,
        }
        records.append(record)
        by_path[path] = record
        return record

    walk(gene_roots[0], ())
    primary_types = {"S", "D", "Leaf"}
    for record in records:
        if record["event_type"] in primary_types:
            node = nodes_by_descendants.get(record["leaves"])
            if node is None:
                raise ValueError(f"XML primary event lacks NHX node at {record['path_tuple']}")
            if node["event_type"] != record["event_type"] or node["species_location"] != record["species_location"]:
                raise ValueError(f"NHX/XML event mismatch at {node['node_id']}")
            record["node_id"] = node["node_id"]
        else:
            record["node_id"] = ""

    def nearest_primary_ancestor(record: dict) -> str:
        path = record["parent_path"]
        while path is not None:
            ancestor = by_path[path]
            if ancestor["event_type"] in primary_types:
                return ancestor["node_id"]
            path = ancestor["parent_path"]
        return ""

    sl_groups: defaultdict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for record in records:
        if record["event_type"] != "SL":
            continue
        child_node = nodes_by_descendants.get(record["leaves"])
        if child_node is None:
            raise ValueError(f"SL path lacks downstream NHX node at {record['path_tuple']}")
        record["edge_parent_node_id"] = nearest_primary_ancestor(record)
        record["edge_child_node_id"] = child_node["node_id"]
        if not record["edge_parent_node_id"]:
            raise ValueError(f"SL path lacks upstream NHX node at {record['path_tuple']}")
        sl_groups[(record["edge_parent_node_id"], record["edge_child_node_id"], record["species_location"])].append(record)
    for group, values in sl_groups.items():
        for ordinal, record in enumerate(sorted(values, key=lambda item: item["path_tuple"]), 1):
            context = "SL|" + "|".join(group) + f"|{ordinal}"
            record["event_id"] = stable_event_id(orthogroup, context)

    for record in records:
        event_type = record["event_type"]
        if event_type in primary_types:
            record["event_id"] = stable_event_id(orthogroup, f"PRIMARY|{record['node_id']}|{event_type}")
            record["edge_parent_node_id"] = ""
            record["edge_child_node_id"] = ""
            record["associated_sl_event_id"] = ""
        elif event_type == "SL":
            record["associated_sl_event_id"] = ""
        else:
            parent = by_path.get(record["parent_path"])
            if parent is None or parent["event_type"] != "SL":
                raise ValueError(f"Loss is not directly attached to SL at {record['path_tuple']}")
            record["edge_parent_node_id"] = parent["edge_parent_node_id"]
            record["edge_child_node_id"] = parent["edge_child_node_id"]
            record["associated_sl_event_id"] = parent["event_id"]
            sibling_losses = sorted(
                value["path_tuple"] for value in records
                if value["event_type"] == "L" and value["parent_path"] == record["parent_path"]
            )
            ordinal = sibling_losses.index(record["path_tuple"]) + 1
            record["event_id"] = stable_event_id(
                orthogroup,
                f"L|{parent['event_id']}|{record['species_location']}|{ordinal}",
            )

    if len({record["event_id"] for record in records}) != len(records):
        raise ValueError(f"Stable event ID collision in {orthogroup}")
    source_hash = sha256_file(xml_path)
    output = []
    for record in sorted(records, key=lambda item: item["path_tuple"]):
        location = record["species_location"]
        row = {
            "orthogroup": orthogroup,
            "family": family,
            "event_id": record["event_id"],
            "event_type": record["event_type"],
            "raw_xml_event": record["raw_xml_event"],
            "is_primary_gene_node_event": record["event_type"] in primary_types,
            "is_auxiliary_reconciliation_event": record["event_type"] in {"SL", "L"},
            "node_id": record["node_id"],
            "edge_parent_node_id": record.get("edge_parent_node_id", ""),
            "edge_child_node_id": record.get("edge_child_node_id", ""),
            "associated_sl_event_id": record.get("associated_sl_event_id", ""),
            "species_location": location,
            **location_fields(location, species_locations),
            "xml_clade_name": record["xml_clade_name"],
            "xml_depth": len(record["path_tuple"]),
            "xml_path": "root" if not record["path_tuple"] else "root." + ".".join(map(str, record["path_tuple"])),
            "n_descendant_extant_leaves": len(record["leaves"]),
            "descendant_extant_leaf_hash": leafset_hash(record["leaves"]),
            "descendant_extant_leaves": ";".join(record["leaves"]),
            "source_xml": str(xml_path.resolve()),
            "source_xml_sha256": source_hash,
            "reconciliation_method": METHOD,
            "parser_version": PARSER_VERSION,
        }
        output.append(row)
    return output


def make_duplication_loss_rows(events: list[dict], node_lookup: dict[str, dict]) -> list[dict]:
    rows = []
    for event in events:
        if event["event_type"] not in {"D", "L"}:
            continue
        evidence_node_id = event["node_id"] if event["event_type"] == "D" else event["edge_child_node_id"]
        evidence = node_lookup[evidence_node_id]
        rows.append(
            {
                "orthogroup": event["orthogroup"],
                "family": event["family"],
                "event_id": event["event_id"],
                "event_type": "duplication" if event["event_type"] == "D" else "loss",
                "node_id": event["node_id"],
                "edge_parent_node_id": event["edge_parent_node_id"],
                "edge_child_node_id": event["edge_child_node_id"],
                "associated_sl_event_id": event["associated_sl_event_id"],
                "species_location": event["species_location"],
                "species_location_is_tip": event["species_location_is_tip"],
                "n_species_location_descendants": event["n_species_location_descendants"],
                "species_location_descendants": event["species_location_descendants"],
                "evidence_node_id": evidence_node_id,
                "split_id": evidence["split_id"],
                "split_mapping_status": evidence["split_mapping_status"],
                "split_mapping_multiplicity": evidence["split_mapping_multiplicity"],
                "iqtree_ufboot": evidence["iqtree_ufboot"],
                "ufboot_status": evidence["ufboot_status"],
                "iqtree_split_branch_length": evidence["iqtree_split_branch_length"],
                "support_scope": evidence["support_scope"],
                "support_is_event_confidence": False,
                "source_xml": event["source_xml"],
                "reconciliation_method": METHOD,
            }
        )
    return sorted(rows, key=lambda row: (row["orthogroup"], row["event_id"]))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/analysis_config.yaml")
    parser.add_argument("--scope-label", default="", help="Must match Step 09/10 scope")
    args = parser.parse_args()
    config, config_path = load_config(args.config)
    reconciliation = config["reconciliation"]
    expected_config = {
        "method": "GeneRax_EVAL_only",
        "version": "1.2.2",
        "rec_model": "UndatedDL",
        "strategy": "EVAL",
        "parser_node_id_scheme": "rooted_descendant_leaf_set_sha256_first20",
        "parser_event_id_scheme": "event_context_sha256_first20",
        "parser_support_source": "step09_canonical_unrooted_split_table",
    }
    for key, expected in expected_config.items():
        if str(reconciliation.get(key)) != expected:
            raise ValueError(f"Step 11 requires reconciliation.{key}={expected!r}")

    input_base = PROJECT_ROOT / "intermediate" / "08_reconciliation_inputs"
    input_root = input_base / args.scope_label if args.scope_label else input_base
    generax_base = PROJECT_ROOT / "results" / "02_generax"
    generax_root = generax_base / args.scope_label if args.scope_label else generax_base
    family_summary_path = generax_root / "generax_family_summary.tsv"
    family_rows = read_tsv(family_summary_path)
    if not family_rows or any(row["status"] != "PASS" for row in family_rows):
        raise ValueError("Step 10 family summary is empty or contains failures")
    expected_ogs = [row["orthogroup"] for row in family_rows]
    family_by_og = {row["orthogroup"]: row for row in family_rows}
    if len(family_by_og) != len(expected_ogs):
        raise ValueError("Duplicate families in Step 10 family summary")

    split_path = input_root / "iqtree_split_support.tsv"
    split_rows = read_tsv(split_path)
    splits_by_og: defaultdict[str, dict[str, dict]] = defaultdict(dict)
    for row in split_rows:
        if row["split_id"] in splits_by_og[row["orthogroup"]]:
            raise ValueError(f"Duplicate Step 09 split {row['split_id']}")
        splits_by_og[row["orthogroup"]][row["split_id"]] = row

    species_tree_version = str(config["project"]["species_tree_version"])
    species_tree = input_root / "common" / f"{species_tree_version}.generax.nwk"
    species_locations = build_species_locations(species_tree)
    all_nodes: list[dict] = []
    all_events: list[dict] = []
    summaries: list[dict] = []
    input_files = [
        {"orthogroup": "ALL", "role": "config", "path": str(config_path.resolve()), "sha256": sha256_file(config_path)},
        {"orthogroup": "ALL", "role": "step10_family_summary", "path": str(family_summary_path.resolve()), "sha256": sha256_file(family_summary_path)},
        {"orthogroup": "ALL", "role": "step09_split_support", "path": str(split_path.resolve()), "sha256": sha256_file(split_path)},
        {"orthogroup": "ALL", "role": "species_tree", "path": str(species_tree.resolve()), "sha256": sha256_file(species_tree)},
    ]

    for orthogroup in expected_ogs:
        family_row = family_by_og[orthogroup]
        family = family_row["family"]
        nhx_path = Path(family_row["reconciled_nhx"])
        xml_path = Path(family_row["reconciled_xml"])
        event_counts_path = Path(family_row["event_counts"])
        leaf_map_path = input_root / family / "mappings" / f"{family}.leaf_species.tsv"
        for role, path in [
            ("reconciled_nhx", nhx_path), ("reconciled_xml", xml_path),
            ("generax_event_counts", event_counts_path), ("leaf_metadata", leaf_map_path),
        ]:
            if not path.is_file():
                raise FileNotFoundError(path)
            input_files.append({"orthogroup": orthogroup, "role": role, "path": str(path.resolve()), "sha256": sha256_file(path)})
        leaf_rows = read_tsv(leaf_map_path)
        leaf_metadata = {row["leaf_id"]: row for row in leaf_rows}
        if len(leaf_metadata) != len(leaf_rows):
            raise ValueError(f"Duplicate leaf metadata for {orthogroup}")
        nodes, nodes_by_descendants = parse_gene_nodes(
            orthogroup, family, nhx_path, xml_path, leaf_metadata,
            splits_by_og[orthogroup], species_locations,
        )
        events = parse_xml_events(
            orthogroup, family, xml_path, nodes_by_descendants, species_locations,
        )
        event_counts = Counter(event["event_type"] for event in events)
        expected_counts = {
            "S": int(family_row["speciation_events"]),
            "SL": int(family_row["speciation_loss_events"]),
            "D": int(family_row["duplication_events"]),
            "L": int(family_row["speciation_loss_events"]) + int(family_row["standalone_loss_events"]),
            "Leaf": int(family_row["leaf_events"]),
        }
        if any(event_counts[key] != value for key, value in expected_counts.items()):
            raise ValueError(f"Step 10/XML event count mismatch for {orthogroup}: {event_counts}")
        root_nodes = [node for node in nodes if node["is_root"]]
        root_edge_nodes = [node for node in nodes if node["split_mapping_status"] == "matched_root_edge_pair"]
        trivial_root_edge_nodes = [
            node for node in nodes
            if node["split_mapping_status"] == "not_applicable_trivial_root_edge"
        ]
        root_edge_splits = {node["split_id"] for node in root_edge_nodes}
        internal_root_edge = len(root_edge_nodes) == 2 and len(root_edge_splits) == 1 and not trivial_root_edge_nodes
        terminal_root_edge = not root_edge_nodes and len(trivial_root_edge_nodes) == 1
        if len(root_nodes) != 1 or not (internal_root_edge or terminal_root_edge):
            raise ValueError(f"Unexpected rooted-edge mapping for {orthogroup}")
        summaries.append(
            {
                "orthogroup": orthogroup,
                "family": family,
                "n_leaves": sum(node["node_kind"] == "leaf" for node in nodes),
                "n_internal_nodes": sum(node["node_kind"] == "internal" for node in nodes),
                "n_gene_nodes": len(nodes),
                "n_speciation_events": event_counts["S"],
                "n_duplication_events": event_counts["D"],
                "n_speciation_loss_events": event_counts["SL"],
                "n_loss_events": event_counts["L"],
                "n_support_mapped_nodes": sum(bool(node["split_id"]) for node in nodes),
                "n_distinct_mapped_splits": len({node["split_id"] for node in nodes if node["split_id"]}),
                "n_unique_split_nodes": sum(node["split_mapping_status"] == "matched_unique" for node in nodes),
                "n_root_edge_pair_nodes": len(root_edge_nodes),
                "root_node_id": root_nodes[0]["node_id"],
                "root_species_location": root_nodes[0]["species_location"],
                "root_edge_split_id": next(iter(root_edge_splits)) if root_edge_splits else "",
                "status": "PASS",
            }
        )
        all_nodes.extend(nodes)
        all_events.extend(events)

    all_nodes.sort(key=lambda row: (row["orthogroup"], row["node_id"]))
    all_events.sort(key=lambda row: (row["orthogroup"], row["xml_path"]))
    node_lookup = {row["node_id"]: row for row in all_nodes}
    if len(node_lookup) != len(all_nodes):
        raise ValueError("Node IDs are not globally unique")
    dl_rows = make_duplication_loss_rows(all_events, node_lookup)

    output_base = PROJECT_ROOT / "results" / "03_reconciliation"
    output_root = ensure_dir(output_base / args.scope_label if args.scope_label else output_base)
    nodes_path = output_root / "reconciliation_nodes.tsv"
    events_path = output_root / "reconciliation_events.tsv"
    dl_path = output_root / "dls_events.tsv"
    loss_edges_path = output_root / "loss_edges.tsv"
    summary_path = output_root / "reconciliation_summary.tsv"
    inputs_path = output_root / "step11_input_files.tsv"
    sc11 = set(config["sc11_sampling"]["single_cell_species"])
    for row in all_nodes + all_events + dl_rows:
        taxa = set(row.get("species_location_descendants", "").split(";")) - {""}
        row["event_scope"] = ("within_Anthozoa" if taxa and taxa <= sc11 else
                              "reference_Medusozoa" if taxa and taxa <= {"Hvul", "Aaur"} else
                              "above_Anthozoa_within_sampled_Cnidaria")
    write_tsv(nodes_path, all_nodes, NODE_FIELDS + ["event_scope"])
    write_tsv(events_path, all_events, EVENT_FIELDS + ["event_scope"])
    write_tsv(dl_path, dl_rows, DL_FIELDS + ["event_scope"])
    loss_edges = [
        {key: event.get(key, "") for key in LOSS_EDGE_FIELDS}
        for event in all_events if event["event_type"] == "L"
    ]
    write_tsv(loss_edges_path, loss_edges, LOSS_EDGE_FIELDS)
    write_tsv(summary_path, summaries, SUMMARY_FIELDS)
    write_tsv(inputs_path, input_files, ["orthogroup", "role", "path", "sha256"])
    output_files = []
    for role, path in [
        ("reconciliation_nodes", nodes_path), ("reconciliation_events", events_path),
        ("dls_events", dl_path), ("loss_edges", loss_edges_path),
        ("reconciliation_summary", summary_path), ("step11_input_files", inputs_path),
    ]:
        output_files.append({"role": role, "path": str(path.resolve()), "sha256": sha256_file(path), "size_bytes": path.stat().st_size})
    outputs_path = output_root / "step11_output_files.tsv"
    write_tsv(outputs_path, output_files, ["role", "path", "sha256", "size_bytes"])
    write_json(
        output_root / "step11_run_summary.json",
        {
            "status": "PASS",
            "parser_version": PARSER_VERSION,
            "node_id_scheme": NODE_ID_SCHEME,
            "event_id_scheme": EVENT_ID_SCHEME,
            "support_scope": SUPPORT_SCOPE,
            "scope_label": args.scope_label,
            "n_families": len(summaries),
            "n_gene_nodes": len(all_nodes),
            "n_events": len(all_events),
            "n_duplication_loss_rows": len(dl_rows),
            "outputs": output_files,
        },
    )
    print(
        f"Step 11 parsing: PASS ({len(summaries)} families, {len(all_nodes)} gene nodes, "
        f"{len(all_events)} XML events)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
