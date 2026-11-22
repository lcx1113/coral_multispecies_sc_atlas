#!/usr/bin/env python3
"""Reduce reference display only, using maximal pure-reference clades in the full tree."""
from pathlib import Path
import sys,json,hashlib
from collections import Counter
import yaml
from Bio import Phylo
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT.parent
sys.path.insert(0,str(SOURCE/"scripts"))
from common import read_tsv,write_tsv,write_json,sha256_file,fasta_dict,utc_now
REF={"Hvul","Aaur"}
AA=set("ACDEFGHIKLMNPQRSTVWY")

def main():
    baseline=SOURCE/"display_control_194/results/01_tip_selection"
    original=read_tsv(baseline/"tip_selection_audit.tsv")
    cfg=yaml.safe_load((SOURCE/"display_control_194/config/analysis_config.yaml").read_text())
    cfg["family_definitions"]="config/family_sets.tsv"
    cfg["selection"]["reference_display_policy"]="one_per_species_per_maximal_pure_reference_clade_in_full_tree"
    cfg["selection"]["retain_all_species_without_expression"]=[]
    cfg["selection"]["reference_species"]=["Hvul","Aaur"]
    cfg["selection"]["coral_selection"]="exact_copy_of_parent_display_selection"
    cfg["selection"]["representative_ranking"]=["trimmed_alignment_canonical_AA_fraction_desc","mean_patristic_distance_same_species_same_clade_asc","protein_id_lexical_asc"]
    (ROOT/"config/analysis_config.yaml").write_text(yaml.safe_dump(cfg,sort_keys=False))
    nodes=read_tsv(SOURCE/cfg["source"]["reconciliation_nodes"])
    definitions=read_tsv(ROOT/"config/family_sets.tsv")
    families=list(dict.fromkeys(r["family_id"] for r in definitions))
    proofs=[];clades=[];summary=[]
    snapshot_paths={SOURCE/"display_control_194/config/analysis_config.yaml",SOURCE/cfg["source"]["reconciliation_nodes"]}
    snapshot_paths.update(p for p in (SOURCE/"display_control_194/results").rglob("*") if p.is_file())
    for k in ["expression","cdd_domain_hits","cdd_sequence_status","reconciliation_events","reconciliation_summary"]:
        snapshot_paths.add(SOURCE/cfg["source"][k])
    frozen=read_tsv(SOURCE/"inputs/frozen_members/selected_references/reference_lookup.tsv")
    freeze_keys={(r["family_id"],r["leaf_id"]) for r in frozen if r["complete_core_combination_pass"]=="TRUE"}
    selected=set()
    for f in families:
        ns=[r for r in nodes if r["orthogroup"]==f]
        nm={r["node_id"]:r for r in ns}
        ls={r["leaf_id"]:r for r in ns if r["node_kind"]=="leaf"}
        refs={k for k,r in ls.items() if r["leaf_species"] in REF}
        desc=lambda r:set(r["descendant_leaves"].split(";"))
        groups=[r for r in ns if desc(r)<=refs and
            (r["parent_node_id"] not in nm or not desc(nm[r["parent_node_id"]])<=refs)]
        assert set().union(*(desc(g) for g in groups))==refs
        assert sum(len(desc(g)) for g in groups)==len(refs)
        path=SOURCE/Path(ns[0]["reconciled_nhx"]);snapshot_paths.add(path)
        tree=Phylo.read(path,"newick")
        assert {t.name for t in tree.get_terminals()}==set(ls)
        alnpath=SOURCE/f"intermediate/06_alignments/{f}/{f}.clipkit.aln.faa"
        aln=fasta_dict(alnpath);snapshot_paths.add(alnpath)
        lengths={len(s) for s in aln.values()};assert len(lengths)==1
        sites=next(iter(lengths))
        baseline_f=[r for r in original if r["family_id"]==f]
        keepcoral=[r for r in baseline_f if r["tree_species"] not in REF and r["retain_tip"]=="TRUE"]
        for r in keepcoral:selected.add((f,r["leaf_id"]))
        nrefs=0
        for group in sorted(groups,key=lambda g:g["node_id"]):
            leaves=desc(group)
            reps=[]
            for sp in sorted({ls[k]["leaf_species"] for k in leaves}):
                ids=sorted(k for k in leaves if ls[k]["leaf_species"]==sp)
                scored=[]
                for leaf in ids:
                    assert (f,leaf) in freeze_keys
                    seq=aln[leaf].upper()
                    standard=sum(c in AA for c in seq)
                    dist=sum(tree.distance(leaf,other) for other in ids if other!=leaf)/(len(ids)-1) if len(ids)>1 else 0.0
                    assert dist>=0
                    scored.append(dict(family_id=f,reference_clade_node_id=group["node_id"],species=sp,
                        leaf_id=leaf,protein_id=ls[leaf]["protein_id"],n_same_species_candidates=len(ids),
                        alignment_sites=sites,n_standard_amino_acids=standard,
                        n_gap=sum(c in "-." for c in seq),n_other_characters=sum(c not in AA and c not in "-." for c in seq),
                        standard_amino_acid_fraction=standard/sites,mean_patristic_distance_same_species_clade=dist,
                        unflagged_core_requirement_pass=True))
                ranked=sorted(scored,key=lambda r:(-r["n_standard_amino_acids"],r["mean_patristic_distance_same_species_clade"],r["protein_id"]))
                rep=ranked[0]["leaf_id"];selected.add((f,rep));reps.append(rep);nrefs+=1
                for rank,row in enumerate(ranked,1):
                    row.update(selection_rank=rank,is_representative=rank==1,representative_leaf_id=rep,
                        represented_members=";".join(ids),n_represented_members=len(ids),
                        reason="highest_alignment_valid_AA_then_minimum_mean_distance_then_protein_ID")
                    proofs.append(row)
            clades.append(dict(family_id=f,reference_clade_node_id=group["node_id"],
                clade_ufboot=group["iqtree_ufboot"],ufboot_status=group["ufboot_status"],
                grouping_basis="maximal_pure_reference_subtree_in_full_reconciled_tree_not_verified_subfamily",
                full_clade_size=len(leaves),n_Hvul=sum(ls[k]["leaf_species"]=="Hvul" for k in leaves),
                n_Aaur=sum(ls[k]["leaf_species"]=="Aaur" for k in leaves),n_representatives=len(reps),
                representative_leaf_ids=";".join(sorted(reps)),full_clade_members=";".join(sorted(leaves))))
        summary.append(dict(family_id=f,n_full_tree_tips=len(ls),n_pure_reference_clades=len(groups),
            n_reference_before=len(refs),n_reference_after=nrefs,n_reference_hidden=len(refs)-nrefs,
            n_SC11_display_unchanged=len(keepcoral),n_display_before=sum(r["retain_tip"]=="TRUE" for r in baseline_f),
            n_display_after=len(keepcoral)+nrefs))
    pm={(r["family_id"],r["leaf_id"]):r for r in proofs}
    audit=[]
    for row in original:
        r=row.copy();k=(r["family_id"],r["leaf_id"])
        if r["tree_species"] in REF:
            assert k in pm and r["retain_tip"]=="TRUE"
            r["retain_tip"]="TRUE" if k in selected else "FALSE"
            r["selection_reason"]="external_reference_representative_no_expression_filter" if k in selected else "hidden_reference_redundant_within_full_tree_clade_species"
            for field in ["reference_clade_node_id","representative_leaf_id","n_represented_members"]:
                r[field]=pm[k][field]
        else:
            r.update(reference_clade_node_id="",representative_leaf_id="",n_represented_members="")
            assert (r["retain_tip"]=="TRUE")== (k in selected)
        audit.append(r)
    out=ROOT/cfg["outputs"]["selection_dir"]
    write_tsv(out/"tip_selection_audit.tsv",audit,list(audit[0]))
    write_tsv(out/"reference_representative_audit.tsv",proofs,list(proofs[0]))
    write_tsv(out/"reference_clades.tsv",clades,list(clades[0]))
    write_tsv(out/"family_pruning_summary.tsv",summary,list(summary[0]))
    for family in families:
        (out/f"{family}.retained_tips.txt").write_text("".join(x+"\n" for f,x in sorted(selected) if f==family))
    og=read_tsv(baseline/"og_selection_summary.tsv")
    species=read_tsv(baseline/"species_og_selection_summary.tsv")
    for rows in [og,species]:
        for r in rows:
            if r["original_orthogroup"]!="NA":continue
            subset=[a for a in audit if a["family_id"]==r["family_id"] and a["original_orthogroup"]=="NA"
                    and ("tree_species" not in r or a["tree_species"]==r["tree_species"])]
            kept=[a for a in subset if a["retain_tip"]=="TRUE"]
            r["n_retained_tips"]=len(kept)
            r["retained_leaf_ids"]=";".join(sorted(a["leaf_id"] for a in kept))
            field="top_selection_mode" if "tree_species" in r else "selection_mode"
            r[field]="one_reference_per_full_tree_clade_species"
    write_tsv(out/"og_selection_summary.tsv",og,list(og[0]))
    write_tsv(out/"species_og_selection_summary.tsv",species,list(species[0]))
    snap=[dict(path=str(p),sha256=sha256_file(p)) for p in sorted(snapshot_paths)]
    oldsnap=out/"parent_preservation_manifest.tsv"
    if oldsnap.exists():assert read_tsv(oldsnap)==snap,"Parent files changed since first selection"
    else:write_tsv(oldsnap,snap,["path","sha256"])
    write_json(out/"tip_selection_summary.json",dict(status="PASS",n_source_tips=len(audit),
        n_retained_tips=len(selected),n_reference_before=len(proofs),n_reference_after=sum(r["is_representative"] for r in proofs),
        n_coral_retained=sum(r["retain_tip"]=="TRUE" and r["tree_species"] not in REF for r in audit),
        selection_rule=cfg["selection"],created_utc=utc_now()))
    print(summary,flush=True)

if __name__=="__main__":main()
