#!/usr/bin/env python3
"""Step 07: run ClipKIT, then remove and report only all-gap sequences."""

from __future__ import annotations

import argparse
import json
import sys

from common import (
    PROJECT_ROOT, alignment_stats, archive_existing, fasta_dict, load_config, og_slug,
    parse_fasta, read_tsv, reusable_command_record, run_command, selected_ogs,
    sha256_file, write_fasta, write_json, write_tsv,
)


SUMMARY_FIELDS = [
    "orthogroup", "og_slug", "n_sequences_before", "n_sequences_clipkit_raw",
    "n_removed_all_gap", "n_sequences_active", "sites_before", "sites_after",
    "retained_site_fraction", "minimum_non_gap_active", "maximum_non_gap_active",
    "mean_gap_fraction_active", "raw_leaf_ids_match_mafft", "active_leaf_contract_valid",
    "exit_code", "timed_out", "input_sha256", "raw_output_sha256",
    "active_output_sha256", "run_action", "archive_path", "command", "status", "error_message",
]
PER_SEQUENCE_FIELDS = [
    "orthogroup", "leaf_id", "species", "protein_id", "transcript_id",
    "canonical_gene_id", "trimmed_non_gap_aa", "trimmed_gap_fraction",
    "all_gap_after_clipkit", "included_after_trimming", "removal_reason",
]
REMOVAL_FIELDS = [
    "orthogroup", "leaf_id", "species", "protein_id", "transcript_id",
    "canonical_gene_id", "trimmed_non_gap_aa", "reason", "mafft_sha256",
    "clipkit_raw_sha256",
]


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
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    config, _ = load_config(args.config)
    cfg = config["trimming"]
    if cfg.get("remove_all_gap_after_clipkit") is not True or cfg.get("all_gap_definition") != "non_gap_aa_eq_0":
        raise ValueError("Step 07 requires the frozen all-gap-removal contract")
    ogs = choose_ogs(args.og)
    root = PROJECT_ROOT / "intermediate" / "06_alignments"
    lookup_rows = read_tsv(PROJECT_ROOT / "intermediate" / "04_candidate_members" / "tree_leaf_lookup.tsv")
    lookup = {(row["orthogroup"], row["leaf_id"]): row for row in lookup_rows}
    summaries, per_sequence, removals, active_lookup_all = [], [], [], []

    for og in ogs:
        slug, out_dir = og_slug(og), root / og_slug(og)
        input_alignment = out_dir / f"{slug}.mafft.aln.faa"
        raw_output = out_dir / f"{slug}.clipkit.raw.aln.faa"
        active_output = out_dir / f"{slug}.clipkit.aln.faa"
        stdout, stderr = out_dir / f"{slug}.clipkit.stdout.log", out_dir / f"{slug}.clipkit.stderr.log"
        record_path = out_dir / f"{slug}.clipkit.command.json"
        row = {key: "" for key in SUMMARY_FIELDS}
        row.update({"orthogroup": og, "og_slug": slug, "status": "ERROR"})
        og_per_sequence, og_removals, active_rows = [], [], []
        try:
            before = alignment_stats(input_alignment)
            input_sha256 = sha256_file(input_alignment)
            command = [
                str(cfg["executable"]), str(input_alignment.resolve()), "-m", str(cfg["mode"]),
                "-g", str(cfg["gap_threshold"]), "-o", str(raw_output.resolve()), "-l",
            ]
            existing = json.loads(record_path.read_text(encoding="utf-8")) if record_path.is_file() else {}
            reusable = raw_output.is_file() and reusable_command_record(existing, command, input_sha256)
            archive_path = ""
            if reusable and not args.force:
                record, action = existing, "REUSED_EXISTING"
            else:
                archive_path = archive_existing([raw_output, active_output, stdout, stderr, record_path], out_dir / "archive")
                record = run_command(command, stdout, stderr, timeout=int(cfg["timeout_seconds"]))
                record.update({"input_sha256": input_sha256, "archive_path": archive_path})
                write_json(record_path, record)
                action = "EXECUTED"
            raw = alignment_stats(raw_output)
            raw_hash = sha256_file(raw_output)
            raw_ids_match = before["ids"] == raw["ids"]
            raw_records = list(parse_fasta(raw_output))
            removed_ids = {leaf_id for leaf_id, count in raw["non_gap_by_id"].items() if count == 0}
            active_records = [(leaf_id, sequence) for leaf_id, _header, sequence in raw_records if leaf_id not in removed_ids]
            write_fasta(active_output, active_records)
            active = alignment_stats(active_output)
            expected_active_ids = raw["ids"] - removed_ids
            active_contract = active["ids"] == expected_active_ids and not (active["ids"] & removed_ids)
            retained = raw["alignment_length"] / before["alignment_length"] if before["alignment_length"] else 0.0
            minimum_sequences = int(config["tree"]["minimum_sequences"])
            status = "PASS" if (
                record["exit_code"] == 0 and raw_ids_match and active_contract
                and active["n_sequences"] >= minimum_sequences and active["alignment_length"] > 0
                and retained >= float(cfg["minimum_retained_site_fraction"])
            ) else "ERROR"
            for leaf_id, _header, sequence in raw_records:
                source = lookup.get((og, leaf_id))
                if source is None:
                    raise ValueError(f"No Step 04 lookup for {og}/{leaf_id}")
                non_gap = raw["non_gap_by_id"][leaf_id]
                gap_fraction = (raw["alignment_length"] - non_gap) / raw["alignment_length"] if raw["alignment_length"] else 1.0
                removed = leaf_id in removed_ids
                seq_row = {
                    "orthogroup": og, "leaf_id": leaf_id, "species": source["species"],
                    "protein_id": source["protein_id"], "transcript_id": source["transcript_id"],
                    "canonical_gene_id": source["canonical_gene_id"], "trimmed_non_gap_aa": non_gap,
                    "trimmed_gap_fraction": f"{gap_fraction:.8f}", "all_gap_after_clipkit": removed,
                    "included_after_trimming": not removed,
                    "removal_reason": "all_gap_after_clipkit" if removed else "",
                }
                og_per_sequence.append(seq_row)
                if removed:
                    og_removals.append({
                        "orthogroup": og, "leaf_id": leaf_id, "species": source["species"],
                        "protein_id": source["protein_id"], "transcript_id": source["transcript_id"],
                        "canonical_gene_id": source["canonical_gene_id"], "trimmed_non_gap_aa": non_gap,
                        "reason": "all_gap_after_clipkit", "mafft_sha256": input_sha256,
                        "clipkit_raw_sha256": raw_hash,
                    })
                else:
                    active_rows.append(source)
            row.update({
                "n_sequences_before": before["n_sequences"], "n_sequences_clipkit_raw": raw["n_sequences"],
                "n_removed_all_gap": len(removed_ids), "n_sequences_active": active["n_sequences"],
                "sites_before": before["alignment_length"], "sites_after": active["alignment_length"],
                "retained_site_fraction": f"{retained:.8f}", "minimum_non_gap_active": active["minimum_non_gap"],
                "maximum_non_gap_active": active["maximum_non_gap"], "mean_gap_fraction_active": f"{active['mean_gap_fraction']:.8f}",
                "raw_leaf_ids_match_mafft": raw_ids_match, "active_leaf_contract_valid": active_contract,
                "exit_code": record["exit_code"], "timed_out": record.get("timed_out", False),
                "input_sha256": input_sha256, "raw_output_sha256": raw_hash,
                "active_output_sha256": sha256_file(active_output), "run_action": action,
                "archive_path": archive_path, "command": record["command_text"], "status": status,
                "error_message": "" if status == "PASS" else "ClipKIT/all-gap active-leaf validation failed",
            })
        except Exception as exc:
            row["error_message"] = str(exc)
        summaries.append(row)
        per_sequence.extend(og_per_sequence)
        removals.extend(og_removals)
        active_lookup_all.extend(active_rows)
        write_tsv(out_dir / "clipkit_summary.tsv", [row], SUMMARY_FIELDS)
        write_tsv(out_dir / "clipkit_per_sequence_qc.tsv", og_per_sequence, PER_SEQUENCE_FIELDS)
        write_tsv(out_dir / "clipkit_removed_sequences.tsv", og_removals, REMOVAL_FIELDS)
        lookup_fields = list(lookup_rows[0]) if lookup_rows else []
        write_tsv(out_dir / "active_tree_leaf_lookup.tsv", sorted(active_rows, key=lambda value: value["leaf_id"]), lookup_fields)
        write_json(out_dir / "clipkit_summary.json", {
            "status": row["status"], "orthogroup": og,
            "member_policy": "remove_only_all_gap_after_clipkit_and_report",
            "n_removed_all_gap": len(og_removals),
        })
        print(f"Step 07 {og}: {row['status']} ({row['n_sequences_before'] or 0} -> {row['n_sequences_active'] or 0} active tips)", file=sys.stderr, flush=True)

    n_errors = sum(row["status"] != "PASS" for row in summaries)
    if len(ogs) > 1 or set(ogs) == set(selected_ogs()):
        write_tsv(root / "clipkit_summary.tsv", summaries, SUMMARY_FIELDS)
        write_tsv(root / "clipkit_per_sequence_qc.tsv", per_sequence, PER_SEQUENCE_FIELDS)
        write_tsv(root / "clipkit_removed_sequences.tsv", removals, REMOVAL_FIELDS)
        lookup_fields = list(lookup_rows[0]) if lookup_rows else []
        write_tsv(root / "active_tree_leaf_lookup.tsv", sorted(active_lookup_all, key=lambda value: (value["orthogroup"], value["leaf_id"])), lookup_fields)
        write_json(root / "clipkit_summary.json", {
            "status": "PASS" if not n_errors else "FAIL", "n_orthogroups": len(summaries),
            "n_errors": n_errors, "n_removed_all_gap": len(removals), "requested_ogs": ogs,
            "member_policy": "remove_only_all_gap_after_clipkit_and_report",
        })
    return 1 if n_errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
