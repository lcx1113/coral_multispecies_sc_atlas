#!/usr/bin/env python3
"""Validate frozen inputs and current tree/track outputs; export full membership and Newick."""
import argparse
import json
import re
from collections import Counter
from pathlib import Path
from Bio import Phylo
from common import PROJECT_ROOT as ROOT, load_config,read_tsv,write_tsv,write_json,sha256_file,fasta_dict

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--phase",choices=["inputs","final"],required=True);a=ap.parse_args()
    cfg,_=load_config("config/analysis_config.yaml");scope=cfg["multiple_og"]["reconciliation_scope"];sc=set(cfg["sc11_sampling"]["single_cell_species"])
    m=read_tsv(ROOT/"intermediate/04_candidate_members/tree_leaf_lookup.tsv");fams={r["family_id"] for r in m}
    q=[r for r in read_tsv(ROOT/"inputs/frozen_members/query_lookup.tsv") if r["family_id"] in fams]
    key=lambda r:(r["family_id"],r["leaf_id"])
    checks=[]
    def check(name,ok,detail=""):
        checks.append(dict(check_id=name,status="PASS" if ok else "FAIL",detail=str(detail)))
    check("input_keys_unique",len({key(r) for r in m})==len(m))
    run_defs=read_tsv(ROOT/cfg["multiple_og"]["family_definitions"])
    check("runtime_family_definitions_match_inputs",{r["family_id"] for r in run_defs}==fams)
    check("runtime_original_OG_count",len({r["original_orthogroup"] for r in run_defs})==cfg["preflight"]["expected_unique_ogs"])
    check("SC11_full_member_set_unchanged",{key(r) for r in m if r["species"] in sc}=={key(r) for r in q})
    check("no_Nvec_members",all(r["species"]!="Nvec" for r in m))
    ext=[r for r in m if r["member_role"]=="external_reference"]
    check("reference_identity_not_original_OG",all(r["original_orthogroup"]=="NA" and r["canonical_gene_id"]=="NA" and r["species"] in {"Hvul","Aaur"} for r in ext))
    frozen=ROOT/"inputs/frozen_members/selected_references/freeze_manifest.json"
    freeze=json.loads(frozen.read_text())
    check("frozen_selection_hash",freeze["decisions_sha256"]==sha256_file(ROOT/"inputs/frozen_members/reference_member_decisions.tsv"))
    check("frozen_selection_rules_hash",freeze["rules_sha256"]==sha256_file(ROOT/"config/reference_selection_rules.yaml"))
    decisions=read_tsv(ROOT/"inputs/frozen_members/reference_member_decisions.tsv")
    includes=[r for r in decisions if r["decision"]=="INCLUDE" and r["family_id"] in fams]
    check("selected_references_equal_frozen_inclusions",{key(r) for r in ext}=={key(r) for r in includes})
    rules,_=load_config("config/reference_selection_rules.yaml")
    check("strict_reference_rules_satisfied",all(r["dual70_pass"]=="TRUE" and r["complete_core_combination_pass"]=="TRUE"
        and float(r["best_qualifying_query_coverage"])>=rules["query_coverage_min"]
        and float(r["best_qualifying_subject_coverage"])>=rules["subject_coverage_min"] for r in includes))
    original_strict=read_tsv(ROOT/"inputs/frozen_members/selected_reference_members.tsv")
    check("reference_members_exact_frozen_strict90",{key(r) for r in ext}=={key(r) for r in original_strict} and len(ext)==90)
    core_hits=read_tsv(ROOT/"inputs/frozen_members/domain_completeness_audit.tsv")
    def complete_core(r):
        names={h["domain_name"] for h in core_hits if h["family_id"]==r["family_id"] and h["protein_id"]==r["protein_id"] and h["incomplete"]=="-"}
        rule=rules["domain_rules"][r["family_id"]]
        return bool(names&set(rule["any_of"])) and set(rule.get("all_of",[]))<=names and sum(n.startswith("C2") and n.endswith("_Ferlin") for n in names)>=rule.get("specific_C2_Ferlin_count_min",0)
    check("strict_core_combination_rederived_from_domain_hits",all(complete_core(r) for r in includes))
    check("strict_family_input_counts",{f:sum(r["family_id"]==f for r in m) for f in fams}=={"FAM_CA":95,"FAM_FERLIN":27,"FAM_MEP1B":215,"FAM_SLC4":39})
    check("frozen_selected_fasta_hash",freeze["selected_fasta_sha256"]==sha256_file(ROOT/"inputs/frozen_members/selected_references/references.faa"))
    st=Phylo.read(ROOT/f"intermediate/01_species_tree/{cfg['project']['species_tree_version']}.nwk","newick")
    check("species_tree13",{n.name for n in st.get_terminals()}==sc|{"Hvul","Aaur"})
    check("species_tree12_internals",len(st.get_nonterminals())==12)
    check("root_focal_decoupled",st.root.name=="Cnidaria" and st.common_ancestor(list(sc)).name=="Anthozoa")
    for f in sorted(fams):
        seq=fasta_dict(ROOT/f"intermediate/05_sequences/{f}/{f}.clean.faa")
        old=fasta_dict(ROOT/f"inputs/frozen_members/queries/{f}.faa")
        check(f+"_SC11_sequences_identical",all(seq.get(k)==v for k,v in old.items()))
    if a.phase=="final":
        track=ROOT/f"results/05_tip_tracks/{scope}"
        expr=read_tsv(track/"tip_relative_expression.tsv");active=read_tsv(ROOT/"intermediate/06_alignments/active_tree_leaf_lookup.tsv")
        ak={key(r) for r in active};check("tracks_equal_active_leaves",{key(r) for r in expr}==ak)
        original={key(r):r for r in read_tsv(ROOT/"inputs/merged_tracks/tip_relative_expression.tsv")}
        fields=["avg_log2FC","pct_target","pct_background","p_val","p_val_adj","umifrac","positive_marker","gene_id","Uniprot_name","expression_observation_status"]
        check("SC11_expression_unchanged",all(all(r[k]==original[key(r)][k] for k in fields) for r in expr if r["tree_species"] in sc))
        for name in ["cdd_sequence_status","cdd_domain_hits"]:
            before=read_tsv(ROOT/f"inputs/merged_tracks/{name}.tsv")
            after=read_tsv(track/(name+".tsv"))
            before=[r for r in before if key(r) in ak and r["species"] in sc]
            after=[r for r in after if r["species"] in sc]
            columns=list(before[0])
            canon=lambda rows:sorted(tuple(r.get(k,"") for k in columns) for r in rows)
            check("SC11_"+name+"_unchanged",canon(before)==canon(after))
        for name in ["cdd_sequence_status","cdd_domain_hits"]:
            before=[r for r in read_tsv(ROOT/f"inputs/merged_tracks/{name}.tsv") if key(r) in ak and r["species"] in {"Hvul","Aaur"}]
            after=[r for r in read_tsv(track/(name+".tsv")) if r["species"] in {"Hvul","Aaur"}]
            columns=list(before[0])
            canon=lambda rows:sorted(tuple(r.get(k,"") for k in columns) for r in rows)
            check("reference_"+name+"_unchanged",canon(before)==canon(after))
        check("external_expression_missing_not_zero",all(r["expression_observation_status"]=="no_expression_dataset" and all(r[k]=="NA" for k in fields[:7]) for r in expr if r["tree_species"] in {"Hvul","Aaur"}))
        audit=read_tsv(ROOT/"display_control_194/results/01_tip_selection/tip_selection_audit.tsv")
        chosen={key(r) for r in audit if r["retain_tip"]=="TRUE"}
        external_active={key(r) for r in active if r["member_role"]=="external_reference"}
        check("all_active_reference_tips_displayed",external_active<=chosen)
        summary=[]
        seqs={f:fasta_dict(ROOT/f"intermediate/05_sequences/{f}/{f}.clean.faa") for f in fams}
        em={key(r):r for r in read_tsv(ROOT/"inputs/merged_tracks/tip_relative_expression.tsv")}
        am={key(r):r for r in audit}
        exports=[]
        for r in m:
            k=key(r);out={**r,**em[k],**am.get(k,{})}
            out.update(in_active_tree=k in ak,in_final_tree=k in chosen,
                selection_reason=am[k]["selection_reason"] if k in am else "all_gap_after_clipkit",
                protein_sequence=seqs[r["family_id"]][r["leaf_id"]])
            exports.append(out)
        write_tsv(ROOT/"results/tree_member.tsv",exports,list(dict.fromkeys(k for r in exports for k in r)))
        nodes=read_tsv(ROOT/f"results/03_reconciliation/{scope}/reconciliation_nodes.tsv")
        event=read_tsv(ROOT/f"results/03_reconciliation/{scope}/reconciliation_events.tsv")
        valid={n.name for n in st.find_clades()}
        check("event_species_locations_valid",all(r["species_location"] in valid for r in nodes+event))
        check("event_scopes_present",all(r["event_scope"] in {"within_Anthozoa","reference_Medusozoa","above_Anthozoa_within_sampled_Cnidaria"} for r in nodes+event))
        plots={r["family_id"]:r for r in read_tsv(ROOT/"display_control_194/results/05_article_figures/article_final_figure_manifest.tsv")}
        structures=[]
        for f in sorted(fams):
            full=Path(next(r["reconciled_nhx"] for r in nodes if r["orthogroup"]==f));newick=re.sub(r"\[[^\]]*\]","",full.read_text().strip())
            path=ROOT/f"results/trees/{f}.full_tree.nwk";path.parent.mkdir(parents=True,exist_ok=True);path.write_text(newick+"\n")
            ft=Phylo.read(path,"newick");pt=Phylo.read(plots[f]["tree"],"newick")
            retained={leaf for family,leaf in chosen if family==f}
            check(f+"_full_tree_active_set",{n.name for n in ft.get_terminals()}=={r["leaf_id"] for r in active if r["family_id"]==f})
            check(f+"_display_tree_retained_set",{n.name for n in pt.get_terminals()}==retained)
            projected={frozenset(t.name for t in n.get_terminals())&retained for n in ft.get_nonterminals()}
            check(f+"_display_exact_topology_projection",{s for s in projected if len(s)>=2}=={frozenset(t.name for t in n.get_terminals()) for n in pt.get_nonterminals()})
            structures.append(dict(family_id=f,full_tree_path=str(path),full_reconciled_nhx=str(full),full_tree_newick=newick,display_tree_path=plots[f]["tree"],display_tree_newick=Path(plots[f]["tree"]).read_text().strip()))
            fm=[r for r in m if r["family_id"]==f];fa=[r for r in audit if r["family_id"]==f and r["retain_tip"]=="TRUE"]
            summary.append(dict(family_id=f,full_tips=len(fm),SC11_full_tips=sum(r["species"] in sc for r in fm),
                Hvul_full_tips=sum(r["species"]=="Hvul" for r in fm),Aaur_full_tips=sum(r["species"]=="Aaur" for r in fm),
                active_tips=sum(r["family_id"]==f for r in active),final_tips=len(fa),final_marker_tips=sum(r["positive_marker"]=="TRUE" for r in fa)))
        write_tsv(ROOT/"results/tree_structure.tsv",structures,list(structures[0]))
        write_tsv(ROOT/"results/family_summary.tsv",summary,list(summary[0]))
        rscript=(ROOT/"display_control_194/scripts/01_prune_and_plot_expression_domains.R").read_text()
        check("article_nonmarker_no_outline",'shape = 21, colour = "transparent", stroke = 0' in rscript)
        code=[p for folder in ["scripts","config","display_control_194/scripts","display_control_194/config"] for p in (ROOT/folder).glob("*") if p.is_file()]
        write_tsv(ROOT/"logs/final_code_config_manifest.tsv",[dict(path=str(p),sha256=sha256_file(p)) for p in code],["path","sha256"])
    out=ROOT/"results/validation";write_tsv(out/(a.phase+"_checks.tsv"),checks,["check_id","status","detail"])
    failed=[r["check_id"] for r in checks if r["status"]!="PASS"]
    result=dict(status="FAIL" if failed else "PASS",phase=a.phase,n_checks=len(checks),failed_checks=failed,completed_families=sorted(fams))
    write_json(out/(a.phase+"_summary.json"),result);print(result,flush=True)
    return bool(failed)
if __name__=="__main__":raise SystemExit(main())
