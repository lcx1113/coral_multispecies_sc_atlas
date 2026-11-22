from pathlib import Path
import csv,json
from label_rules import simplify,expand,truncate
R=Path(__file__).resolve().parents[1]
def read(p):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def put(p,rr):
 with Path(p).open('w') as f:
  w=csv.DictWriter(f,list(rr[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
rr={}
for direction in ['up','down']:
 for r in read(R/f'results/plot_data/majority_all_celltypes/{direction}/row_annotations.tsv'):rr[r['orthogroup']]=r
cfg=json.loads((R/'config/analysis_config.json').read_text())
ann={x['orthogroup']:x for x in read(Path(cfg['annotation_root'])/'orthogroup_final_annotations.csv')}
res=[]
for og,r in sorted(rr.items()):
 raw=ann[og]['first_annotation'];assert raw==r['first_annotation']
 short=simplify(raw);assert expand(short)==raw
 show=truncate(short)
 res.append(dict(orthogroup=og,first_annotation=raw,display_annotation=short,plot_annotation=show,TF_display_label=og+' | '+show,abbreviated=raw!=short,truncated=len(short)>68,original_length=len(raw),display_length=len(short),plot_length=len(show)))
put(R/'results/annotation_display_mapping.tsv',res)
fractions={x['orthogroup']:x for x in read(R/'results/TF_OG.member_fraction_all_candidates.tsv')}
put(R/'results/TF_OG.plotted_member_fraction_and_first_annotation.tsv',[dict(fractions[og],first_annotation=ann[og]['first_annotation']) for og in sorted(rr)])
mem=read(R/'results/tf_og_member_evidence.tsv');put(R/'results/TF_OG.plotted_member_evidence.tsv',[r for r in mem if r['orthogroup'] in rr])
print('ANNOTATIONS',len(rr),'OGs',sum(x['abbreviated'] for x in res),'shortened',sum(x['truncated'] for x in res),'truncated')
