#!/usr/bin/env python3
"""Step 10: run frozen GeneRax 1.2.2 EVAL-only reconciliation and validate outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from xml.etree import ElementTree

from Bio import Phylo

from common import (
    PROJECT_ROOT, archive_existing, load_config, parse_bool, read_tsv,
    reusable_command_record, run_command, sha256_file, write_json, write_tsv,
)


EVENT_KEYS = {"S", "SL", "D", "T", "TL", "L", "Leaf"}


def tree_tip_set(path: Path) -> set[str]:
    tree = Phylo.read(path, "newick")
    names = [tip.name for tip in tree.get_terminals()]
    if None in names or len(names) != len(set(names)):
        raise ValueError(f"Missing or duplicate tree tips in {path}")
    return set(names)


def canonical_split(descendants: set[str], leaves: set[str]) -> tuple[tuple[str, ...], tuple[str, ...]]:
    complement = leaves - descendants
    if not descendants or not complement:
        raise ValueError("Invalid bipartition")
    sides = [tuple(sorted(descendants)), tuple(sorted(complement))]
    sides.sort(key=lambda side: (len(side), side))
    return sides[0], sides[1]


def split_identifier(side_a: tuple[str, ...], side_b: tuple[str, ...]) -> str:
    return hashlib.sha256(("\n".join(side_a) + "\n||\n" + "\n".join(side_b)).encode("utf-8")).hexdigest()


def unrooted_split_ids(path: Path) -> set[str]:
    tree = Phylo.read(path, "newick")
    leaves = tree_tip_set(path)
    result = set()
    for clade in tree.find_clades(order="preorder"):
        descendants = {tip.name for tip in clade.get_terminals()}
        complement = leaves - descendants
        if len(descendants) < 2 or len(complement) < 2:
            continue
        result.add(split_identifier(*canonical_split(descendants, leaves)))
    return result


def rooted_clades(path: Path) -> set[frozenset[str]]:
    tree = Phylo.read(path, "newick")
    leaves = tree_tip_set(path)
    return {
        frozenset(tip.name for tip in clade.get_terminals())
        for clade in tree.get_nonterminals()
        if clade is not tree.root and 1 < len(clade.get_terminals()) < len(leaves)
    }


def composite_sha256(paths: list[Path], values: list[str]) -> str:
    digest = hashlib.sha256()
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
        digest.update(str(path.resolve()).encode())
        digest.update(b"\0")
        digest.update(sha256_file(path).encode())
        digest.update(b"\0")
    for value in values:
        digest.update(value.encode())
        digest.update(b"\0")
    return digest.hexdigest()


def parse_event_counts(path: Path) -> dict[str, int]:
    result = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            key, value = line.split(":", 1)
            result[key] = int(value)
    if set(result) != EVENT_KEYS or any(value < 0 for value in result.values()):
        raise ValueError(f"Invalid GeneRax event-count file: {path}")
    return result


def xml_event_counts(path: Path) -> Counter:
    counts = Counter()
    for element in ElementTree.parse(path).getroot().iter():
        tag = element.tag.rsplit("}", 1)[-1]
        if tag in {"speciation", "duplication", "loss", "leaf", "transferBack", "branchingOut"}:
            counts[tag] += 1
    return counts


def parse_stats(path: Path) -> dict[str, float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    result = {}
    for key in ["JointLL", "LibpllLL", "RecLL"]:
        match = re.search(rf"{key}:\s*([-+0-9.eE]+)", text)
        if not match:
            raise ValueError(f"Missing {key} in {path}")
        result[key] = float(match.group(1))
    return result


def expected_family_outputs(prefix: Path, family: str) -> dict[str, Path]:
    root = prefix / "reconciliations"
    return {
        "event_counts": root / f"{family}_eventCounts.txt",
        "species_event_counts": root / f"{family}_speciesEventCounts.txt",
        "reconciled_nhx": root / f"{family}_reconciliated.nhx",
        "reconciled_xml": root / f"{family}_reconciliated.xml",
        "orthogroups": root / f"{family}_orthogroups.txt",
        "orthogroups_all": root / f"{family}_orthogroups_all.txt",
        "transfers": root / f"{family}_transfers.txt",
    }


def require_config(config: dict) -> dict:
    cfg = config["reconciliation"]
    expected = {
        "method": "GeneRax_EVAL_only", "version": "1.2.2", "rec_model": "UndatedDL",
        "strategy": "EVAL", "gene_tree_source": "iqtree_treefile",
        "gene_tree_variant": "uncollapsed_only", "substitution_model_source": "iqtree_selected_model",
        "support_threshold": "not_applicable_fixed_topology_eval",
    }
    for key, value in expected.items():
        if str(cfg.get(key)) != str(value):
            raise ValueError(f"reconciliation.{key} must be {value!r}")
    for key, value in {
        "unrooted_gene_tree": True, "treerecs_active": False, "allow_horizontal_transfer": False,
        "per_family_rates": False, "eval_likelihood_is_evidence": False,
    }.items():
        if parse_bool(cfg.get(key), f"reconciliation.{key}") != value:
            raise ValueError(f"reconciliation.{key} must be {value}")
    if float(cfg["default_duplication_rate"]) != 1.0 or float(cfg["default_loss_rate"]) != 1.0:
        raise ValueError("Frozen GeneRax EVAL default duplication/loss rates must be 1.0/1.0")
    for path_key, hash_key in [
        ("generax_executable", "generax_binary_sha256"), ("generax_mpiexec", "mpiexec_sha256"),
    ]:
        path = Path(config["environments"][path_key])
        if not path.is_file() or sha256_file(path) != cfg[hash_key]:
            raise ValueError(f"Frozen executable/hash check failed for {path_key}: {path}")
    ranks = int(config["resources"]["generax_mpi_ranks"])
    if ranks != 16 or ranks > int(config["resources"]["maximum_total_cpu"]):
        raise ValueError("GeneRax must use 16 ranks within the global CPU budget")
    return cfg


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/analysis_config.yaml")
    parser.add_argument("--scope-label", default="", help="Must match Step 09 scope")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.scope_label and not re.fullmatch(r"[A-Za-z0-9_.-]+", args.scope_label):
        raise ValueError("Unsafe --scope-label")
    config, config_path = load_config(args.config)
    cfg = require_config(config)
    input_base = PROJECT_ROOT / "intermediate" / "08_reconciliation_inputs"
    input_root = input_base / args.scope_label if args.scope_label else input_base
    family_manifest_path = input_root / "generax" / "family_manifest.tsv"
    families_path = input_root / "generax" / "families.txt"
    species_tree_version = str(config["project"]["species_tree_version"])
    species_tree = input_root / "common" / f"{species_tree_version}.generax.nwk"
    split_path = input_root / "iqtree_split_support.tsv"
    manifest = read_tsv(family_manifest_path)
    if not manifest or len({row["family"] for row in manifest}) != len(manifest):
        raise ValueError("Empty or duplicate GeneRax family manifest")
    expected_splits = defaultdict(set)
    for row in manifest:
        expected_splits[row["orthogroup"]] = set()
    for row in read_tsv(split_path):
        expected_splits[row["orthogroup"]].add(row["split_id"])
    if set(expected_splits) != {row["orthogroup"] for row in manifest}:
        raise ValueError("IQ-TREE split table family set differs from GeneRax manifest")

    result_base = PROJECT_ROOT / "results" / "02_generax"
    result_root = result_base / args.scope_label if args.scope_label else result_base
    native_root = result_root / "native"
    native_root.mkdir(parents=True, exist_ok=True)
    prefix = native_root / "generax_eval_v1_2_2"
    stdout = native_root / "generax_eval_v1_2_2.stdout.log"
    stderr = native_root / "generax_eval_v1_2_2.stderr.log"
    command_json = native_root / "generax_eval_v1_2_2.command.json"
    ranks = int(config["resources"]["generax_mpi_ranks"])
    command = [
        str(Path(config["environments"]["generax_mpiexec"]).resolve()), "-np", str(ranks),
        str(Path(config["environments"]["generax_executable"]).resolve()),
        "--families", str(families_path.resolve()), "--species-tree", str(species_tree.resolve()),
        "--rec-model", str(cfg["rec_model"]), "--strategy", "EVAL", "--unrooted-gene-tree",
        "--prefix", str(prefix.resolve()), "--seed", str(cfg["seed"]),
    ]
    input_files = [families_path, family_manifest_path, species_tree, split_path, config_path]
    for row in manifest:
        input_files.extend(Path(row[key]) for key in ["starting_gene_tree", "alignment", "mapping"])
    input_sha = composite_sha256(input_files, [
        "GeneRax-1.2.2", "EVAL", str(ranks), str(cfg["rec_model"]), str(cfg["seed"]),
        "unrooted", "no-per-family-rates", str(cfg["generax_binary_sha256"]), str(cfg["mpiexec_sha256"]),
    ])
    family_outputs = {row["family"]: expected_family_outputs(prefix, row["family"]) for row in manifest}
    run_paths = {
        "stats": prefix / "stats.txt", "starting_species_tree": prefix / "starting_species_tree.newick",
        "inferred_species_tree": prefix / "inferred_species_tree.newick",
        "per_species_coverage": prefix / "perSpeciesCoverage.txt", "generax_log": prefix / "generax.log",
    }
    complete = all(
        path.is_file() and (path.stat().st_size > 0 or name == "transfers")
        for paths in family_outputs.values() for name, path in paths.items()
    ) and all(path.is_file() and path.stat().st_size > 0 for path in run_paths.values())
    existing = json.loads(command_json.read_text(encoding="utf-8")) if command_json.is_file() else {}
    archive_path = ""
    if complete and reusable_command_record(existing, command, input_sha) and not args.force:
        record, action = existing, "REUSED_EXISTING"
    else:
        archive_path = archive_existing([prefix, stdout, stderr, command_json], result_root / "archive")
        record = run_command(command, stdout, stderr, timeout=int(cfg["timeout_seconds"]))
        record.update({
            "input_sha256": input_sha, "scope_label": args.scope_label,
            "strategy": "EVAL", "generax_version": "1.2.2", "archive_path": archive_path,
        })
        write_json(command_json, record)
        action = "EXECUTED"
    run_status = "PASS" if record.get("exit_code") == 0 else "ERROR"
    run_error = "" if run_status == "PASS" else f"GeneRax exit={record.get('exit_code')}: {record.get('exception', '')}"

    family_rows, n_errors = [], 0
    for row in manifest:
        family, og = row["family"], row["orthogroup"]
        paths = family_outputs[family]
        values = {key: 0 for key in EVENT_KEYS}
        xml_counts = Counter()
        status, error = "PASS", ""
        n_input = n_output = output_split_count = 0
        topology_same = False
        try:
            if run_status != "PASS":
                raise ValueError(run_error)
            missing = [name for name, path in paths.items() if not path.is_file() or (path.stat().st_size == 0 and name != "transfers")]
            if missing:
                raise ValueError(f"Missing/empty GeneRax outputs: {missing}")
            input_tips = tree_tip_set(Path(row["starting_gene_tree"]))
            output_tips = tree_tip_set(paths["reconciled_nhx"])
            n_input, n_output = len(input_tips), len(output_tips)
            if input_tips != output_tips:
                raise ValueError("Reconciled NHX leaf set differs from IQ-TREE input")
            output_splits = unrooted_split_ids(paths["reconciled_nhx"])
            output_split_count = len(output_splits)
            topology_same = output_splits == expected_splits[og]
            if not topology_same:
                raise ValueError("GeneRax EVAL changed the IQ-TREE unrooted split set")
            values = parse_event_counts(paths["event_counts"])
            if values["Leaf"] != n_input or values["T"] != 0 or values["TL"] != 0:
                raise ValueError("GeneRax event counts violate leaf/UndatedDL expectations")
            xml_counts = xml_event_counts(paths["reconciled_xml"])
            if xml_counts["duplication"] != values["D"] or xml_counts["leaf"] != values["Leaf"]:
                raise ValueError("RecPhyloXML duplication/leaf counts disagree with eventCounts")
            if xml_counts["speciation"] != values["S"] + values["SL"] or xml_counts["loss"] != values["SL"] + values["L"]:
                raise ValueError("RecPhyloXML speciation/loss counts disagree with eventCounts")
        except Exception as exc:
            status, error = "ERROR", str(exc)
            n_errors += 1
        family_rows.append({
            "orthogroup": og, "family": family, "strategy": "EVAL",
            "n_input_leaves": n_input, "n_output_leaves": n_output,
            "n_iqtree_internal_splits": len(expected_splits[og]),
            "n_reconciled_internal_splits": output_split_count,
            "unrooted_split_set_identical_to_iqtree": topology_same,
            "speciation_events": values["S"], "speciation_loss_events": values["SL"],
            "duplication_events": values["D"], "transfer_events": values["T"],
            "transfer_loss_events": values["TL"], "standalone_loss_events": values["L"],
            "leaf_events": values["Leaf"], "xml_speciation_tags": xml_counts["speciation"],
            "xml_duplication_tags": xml_counts["duplication"], "xml_loss_tags": xml_counts["loss"],
            "xml_leaf_tags": xml_counts["leaf"],
            **{name: str(path.resolve()) for name, path in paths.items()},
            **{f"{name}_sha256": sha256_file(path) if path.is_file() else "" for name, path in paths.items()},
            "status": status, "error_message": error,
        })
        print(f"Step 10 {family}: {status}; splits={len(expected_splits[og])}->{output_split_count}; D={values['D']}; SL={values['SL']}")

    stats_values = {"JointLL": "", "LibpllLL": "", "RecLL": ""}
    if run_status == "PASS":
        try:
            missing = [name for name, path in run_paths.items() if not path.is_file() or path.stat().st_size == 0]
            if missing:
                raise ValueError(f"Missing/empty GeneRax run outputs: {missing}")
            stats_values = parse_stats(run_paths["stats"])
            if any(value != 0 for value in stats_values.values()):
                raise ValueError("GeneRax 1.2.2 EVAL likelihood fields were expected to be zero")
            if rooted_clades(run_paths["starting_species_tree"]) != rooted_clades(species_tree):
                raise ValueError("GeneRax starting species tree differs from Step 09 species tree")
            if rooted_clades(run_paths["inferred_species_tree"]) != rooted_clades(species_tree):
                raise ValueError("GeneRax inferred species tree differs from fixed species tree")
        except Exception as exc:
            run_status, run_error = "ERROR", str(exc)
            n_errors += 1
    run_row = {
        "scope_label": args.scope_label, "strategy": "EVAL", "generax_version": "1.2.2",
        "generax_binary": str(Path(config["environments"]["generax_executable"]).resolve()),
        "generax_binary_sha256": cfg["generax_binary_sha256"],
        "mpiexec": str(Path(config["environments"]["generax_mpiexec"]).resolve()),
        "mpiexec_sha256": cfg["mpiexec_sha256"], "mpi_ranks": ranks,
        "rec_model": cfg["rec_model"], "unrooted_gene_tree": True, "per_family_rates": False,
        "default_duplication_rate": cfg["default_duplication_rate"],
        "default_loss_rate": cfg["default_loss_rate"], "likelihood_evidence": False,
        "joint_ll": stats_values["JointLL"], "libpll_ll": stats_values["LibpllLL"],
        "rec_ll": stats_values["RecLL"], "seed": cfg["seed"],
        "families_file": str(families_path.resolve()), "species_tree": str(species_tree.resolve()),
        "iqtree_split_support": str(split_path.resolve()), "prefix": str(prefix.resolve()),
        "stats": str(run_paths["stats"].resolve()), "command_record": str(command_json.resolve()),
        "input_sha256": input_sha, "exit_code": record.get("exit_code", ""),
        "timed_out": record.get("timed_out", False), "run_action": action,
        "archive_path": archive_path, "status": run_status, "error_message": run_error,
    }
    write_tsv(result_root / "generax_run_summary.tsv", [run_row], list(run_row))
    write_tsv(result_root / "generax_family_summary.tsv", family_rows, list(family_rows[0]))
    final_pass = n_errors == 0 and run_status == "PASS"
    write_json(result_root / "generax_run_summary.json", {
        "status": "PASS" if final_pass else "FAIL", "scope_label": args.scope_label,
        "strategy": "EVAL", "generax_version": "1.2.2", "n_families": len(manifest),
        "n_errors": n_errors, "native_output_root": str(native_root.resolve()),
        "config": str(config_path.resolve()), "config_sha256": sha256_file(config_path),
        "likelihood_interpretation": "EVAL likelihood fields are zero and are not evidence",
    })
    return 0 if final_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
