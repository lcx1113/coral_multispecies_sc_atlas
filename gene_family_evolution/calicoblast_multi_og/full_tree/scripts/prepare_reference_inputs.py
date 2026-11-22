#!/usr/bin/env python3
"""Construct named SC11 + two-reference species tree and frozen full-family inputs."""
import argparse
import hashlib
from collections import defaultdict
from pathlib import Path
from Bio import Phylo
from Bio.Phylo.Newick import Tree, Clade
from common import PROJECT_ROOT as ROOT, load_config, read_tsv, write_tsv, write_json, write_fasta, fasta_dict, sha256_file

def table(path,rows,fields=None):
    write_tsv(path,rows,fields or list(rows[0]))

def species_tree(cfg):
    sc=set(cfg["sc11_sampling"]["single_cell_species"])
    source=Path(cfg["sc11_sampling"]["source_named_species_tree"])
    tree=Phylo.read(source,"newick")
    for tip in list(tree.get_terminals()):
        if tip.name not in sc: tree.prune(tip)
    assert tree.root.name=="Anthozoa"
    core={(c.name,tuple(sorted(t.name for t in c.get_terminals()))) for c in tree.get_nonterminals()}
    tree=Tree(root=Clade(name="Cnidaria",clades=[Clade(name="Medusozoa",clades=[Clade(name="Hvul"),Clade(name="Aaur")]),tree.root]),rooted=True)
    for node in tree.find_clades(): node.branch_length=None
    assert len(tree.get_terminals())==13 and len(tree.get_nonterminals())==12
    assert all(len(n.clades)==2 for n in tree.get_nonterminals())
    assert len({n.name for n in tree.find_clades()})==25
    focal=next(n for n in tree.find_clades() if n.name=="Anthozoa")
    assert core=={(c.name,tuple(sorted(t.name for t in c.get_terminals()))) for c in focal.find_clades() if not c.is_terminal()}
    out=ROOT/"intermediate/01_species_tree"; out.mkdir(parents=True,exist_ok=True)
    path=out/(cfg["project"]["species_tree_version"]+".nwk")
    Phylo.write(tree,path,"newick",plain=True)
    parents={id(ch):p for p in tree.find_clades() for ch in p.clades}
    rows=[];edges=[]
    for n in tree.find_clades():
        tips=sorted(t.name for t in n.get_terminals()); parent=parents.get(id(n))
        rows.append(dict(species_tree_version=cfg["project"]["species_tree_version"],node_kind="tip" if n.is_terminal() else "internal",
            canonical_node_name=n.name,stable_node_id=n.name,source_internal_label=n.name,
            parent_canonical_node_name=parent.name if parent else "",child_canonical_node_names=";".join(ch.name for ch in n.clades),
            branch_length="",is_root=n is tree.root,n_descendant_species=len(tips),descendant_species=";".join(tips),
            descendant_set_sha256=hashlib.sha256("\n".join(tips).encode()).hexdigest(),
            name_sources="original22_SC11_pruning" if n.name in {x[0] for x in core}|sc else "user_Metazoa_backbone",
            event_scope="within_Anthozoa" if set(tips)<=sc else "reference_Medusozoa" if set(tips)<={"Hvul","Aaur"} else "above_Anthozoa_within_sampled_Cnidaria"))
        if parent: edges.append(dict(parent=parent.name,child=n.name,child_kind="tip" if n.is_terminal() else "internal",branch_length=""))
    table(out/"species_tree_node_contract.resolved.tsv",rows);table(out/"species_tree_edges.tsv",edges)
    write_json(out/"species_tree_build.summary.json",dict(status="PASS",n_species=13,n_internal_nodes=12,n_edges=24,root="Cnidaria",focal_ancestor="Anthozoa",source_tree=str(source),source_tree_sha256=sha256_file(source),output_tree_sha256=sha256_file(path),branch_lengths="omitted"))

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--config",default="config/analysis_config.yaml");a=ap.parse_args()
    cfg,_=load_config(a.config);species_tree(cfg)
    search=ROOT/"inputs/frozen_members"
    frozen=search/"selected_references/freeze_manifest.json"
    import json
    manifest=json.loads(frozen.read_text())
    assert manifest["status"]=="PASS"
    assert manifest["rules_sha256"]==sha256_file(ROOT/"config/reference_selection_rules.yaml")
    decisions=search/"reference_member_decisions.tsv"
    assert sha256_file(decisions)==manifest["decisions_sha256"]
    refs=read_tsv(search/"selected_references/reference_lookup.tsv")
    assert manifest["selected_fasta_sha256"]==sha256_file(search/"selected_references/references.faa")
    refseq=fasta_dict(search/"selected_references/references.faa")
    query=read_tsv(search/"query_lookup.tsv")
    defs=[r for r in read_tsv(ROOT/cfg["multiple_og"]["family_definitions"]) if r["active"]=="TRUE"]
    families=sorted({r["family_id"] for r in defs})
    assert set(families)<=set(manifest["approved_families"])
    members=[];summaries=[]
    for family in families:
        rows=[r.copy() for r in query if r["family_id"]==family]
        seq=fasta_dict(search/f"queries/{family}.faa")
        template={k:"NA" for k in rows[0]}
        for ref in [r for r in refs if r["family_id"]==family]:
            row=template.copy();row.update(ref)
            row.update(orthogroup=family,original_orthogroup="NA",member_role="external_reference",
                family_name=rows[0]["family_name"],og_color="#BDBDBD",og_order=99,
                canonical_gene_id="NA",transcript_id="NA",raw_gene_id="NA",broccoli_gene="NA",
                gene_name="NA",orthogroup_name="NA",is_orphan="NA",is_multi_og="NA",
                mapping_method="external_FASTA_first_token_no_gene_assignment",
                resolution_status="EXTERNAL_PROTEIN_RESOLVED",object_species_alias=ref["species"],
                fasta_header=ref["protein_id"],gene_id_matches_broccoli="NA")
            rows.append(row);seq[ref["leaf_id"]]=refseq[ref["leaf_id"]]
        assert len(seq)==len(rows)==len({r["leaf_id"] for r in rows})
        out=ROOT/f"intermediate/05_sequences/{family}";out.mkdir(parents=True,exist_ok=True)
        fasta=out/f"{family}.clean.faa";write_fasta(fasta,sorted(seq.items()))
        fields=list(dict.fromkeys(k for r in rows for k in r))
        table(out/"family_members.tsv",rows,fields);members+=rows
        summaries.append(dict(family_id=family,family_name=rows[0]["family_name"],
            original_orthogroups=";".join(d["original_orthogroup"] for d in defs if d["family_id"]==family),
            n_members=len(rows),n_species=len({r["species"] for r in rows}),combined_fasta=str(fasta),
            combined_fasta_sha256=sha256_file(fasta),status="PASS"))
    fields=list(dict.fromkeys(k for r in members for k in r))
    for name in ["tree_leaf_lookup","family_members"]:table(ROOT/f"intermediate/04_candidate_members/{name}.tsv",members,fields)
    table(ROOT/"intermediate/02_candidates/candidate_ogs.tsv",[dict(orthogroup=f,include_for_tree=True) for f in families])
    table(ROOT/"intermediate/02_candidates/family_summary.tsv",summaries)
    write_json(ROOT/"results/00_input_audit/family_input_summary.json",dict(status="PASS",n_families=len(families),n_members=len(members),n_species_tree_tips=13,n_species_tree_internal_nodes=12,root="Cnidaria",focal="Anthozoa"))
    print("PASS reference family inputs",len(members),"members",flush=True)

if __name__=="__main__": main()
