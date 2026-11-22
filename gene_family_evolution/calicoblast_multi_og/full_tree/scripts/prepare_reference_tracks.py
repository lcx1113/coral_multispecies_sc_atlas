#!/usr/bin/env python3
"""Map the supplied expression statistics and precomputed CDD domains onto active leaves."""
import argparse
from collections import Counter,defaultdict
from common import PROJECT_ROOT as ROOT, load_config,read_tsv,write_tsv,write_json

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--config",default="config/analysis_config.yaml");ap.add_argument("--track-config",default="config/tip_tracks.yaml")
    ap.add_argument("--og",action="append")
    a=ap.parse_args();cfg,_=load_config(a.config);tc,_=load_config(a.track_config);tc=tc["multiple_og_tracks"]
    active=([r for f in a.og for r in read_tsv(ROOT/f"intermediate/06_alignments/{f}/active_tree_leaf_lookup.tsv")]
            if a.og else read_tsv(ROOT/"intermediate/06_alignments/active_tree_leaf_lookup.tsv"))
    keys={(r["family_id"],r["leaf_id"]) for r in active};assert len(keys)==len(active)
    out=ROOT/tc["output_directory"]; outputs={}
    for name in ["tip_relative_expression","cdd_domain_hits","cdd_sequence_status"]:
        src=read_tsv(ROOT/f"inputs/merged_tracks/{name}.tsv")
        rows=[r for r in src if (r["family_id"],r["leaf_id"]) in keys]
        if name!="cdd_domain_hits":
            assert {(r["family_id"],r["leaf_id"]) for r in rows}==keys and len(rows)==len(keys)
        write_tsv(out/(name+".tsv"),rows,list(src[0]));outputs[name]=rows
    summary=[]
    for f in sorted({r["family_id"] for r in active}):
        e=[r for r in outputs["tip_relative_expression"] if r["family_id"]==f]
        s=[r for r in outputs["cdd_sequence_status"] if r["family_id"]==f]
        states=Counter(r["expression_observation_status"] for r in e)
        summary.append(dict(family_id=f,n_active_tips=len(e),n_observed=states["observed"],
            n_no_expression_dataset=states["no_expression_dataset"],n_gene_not_tested=states["gene_not_tested"],
            n_celltype_not_present=states["celltype_not_present"],n_other_expression_state=0,
            n_domain_hit_tips=sum(r["domain_status"]=="domain_hit" for r in s),
            n_no_domain_hit_tips=sum(r["domain_status"]=="no_domain_hit" for r in s),
            n_cdd_domain_hits=sum(r["family_id"]==f for r in outputs["cdd_domain_hits"]),status="PASS"))
    write_tsv(out/"tip_track_coverage.tsv",summary,list(summary[0]))
    write_json(out/"tip_track_run_summary.json",dict(status="PASS",n_families=len(summary),n_active_tips=len(active),
        expression_recomputed=False,SC11_cdd_recomputed=False,external_cdd_newly_computed=False,missing_expression_is_zero=False))
    print("PASS combined tracks",len(active))
if __name__=="__main__":main()
