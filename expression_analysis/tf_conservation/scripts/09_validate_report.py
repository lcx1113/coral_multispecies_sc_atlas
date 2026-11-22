from pathlib import Path
import csv,gzip,json,hashlib,subprocess,math
from label_rules import simplify,expand,truncate
R=Path(__file__).resolve().parents[1]
def read(p):
 p=Path(p);op=gzip.open if p.suffix=='.gz' else open
 with op(str(p),'rt') as f:return list(csv.DictReader(f,delimiter='\t'))
def put(p,rr):
 rr=list(rr)
 with Path(p).open('w') as f:
  w=csv.DictWriter(f,list(rr[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rr)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for c in iter(lambda:f.read(8388608),b''):h.update(c)
 return h.hexdigest()
def num(s):return float('nan') if s in ('NA','') else float(s)
def same(a,b):
 if a==b:return True
 try:return (math.isnan(num(a)) and math.isnan(num(b))) or math.isclose(num(a),num(b),rel_tol=1e-13,abs_tol=1e-13)
 except ValueError:return False
checks=[]
def check(name,ok,detail):
 checks.append(dict(check=name,status='PASS' if ok else 'FAIL',detail=detail));assert ok,(name,detail)
cfg=json.loads((R/'config/analysis_config.json').read_text())
check('statistical_validation',all(x['status']=='PASS' for x in read(R/'results/statistical_validation.tsv')),'native marker votes, fixed denominators and member annotations')
units=read(R/'config/original_combination_manifest.tsv');cols=read(R/'results/plot_columns.tsv');thresholds=read(R/'config/comparison_thresholds.tsv')
check('excluded_types_and_fixed_denominators',len(units)==164 and len(cols)==157 and len(thresholds)==19 and not ({x['comparison_id'] for x in thresholds}|{x['cell_type'] for x in units})&set(cfg['excluded_analysis_targets']),'both epidermis_like targets excluded from statistics, votes and plots; other denominators unchanged')
check('filter_ct_and_unknown_policy',all(x['column_id'] not in cfg['filter_ct'] for x in cols) and 'unknown26' not in {x['cell_type'] for x in units} and 'Xesp|unknown5' not in {x['column_id'] for x in cols},'unknown5 recomputed but hidden; unknown26 background only')
annotations=read(R/'results/annotation_display_mapping.tsv')
check('annotation_shortening_roundtrip',all(expand(x['display_annotation'])==x['first_annotation'] and truncate(x['display_annotation'])==x['plot_annotation'] for x in annotations),'original first_annotation preserved with lossless numeric-prefix shortening before truncation')
fractions=read(R/'results/TF_OG.plotted_member_fraction_and_first_annotation.tsv')
check('TF80_gene_fraction',all(5*int(x['n_tf_genes'])>=4*int(x['n_unique_species_gene']) for x in fractions),'unique species x gene denominator, >=80%')
figs=read(R/'results/figure_manifest.tsv');keys=read(R/'config/plotting_manifest.tsv');pdfs=[];pointcounts=[];cells=0
for k in keys:
 key=k['figure_key'];p=R/'results/plot_data'/key;tp=R/'results/transposed_plot_data'/key
 rr=read(p/'row_annotations.tsv');cc=read(p/'column_annotations.tsv')
 if not rr:continue
 for name in [k['matrix_name'],'marker_state_matrix.tsv','plot_value_source_matrix.tsv']:
  native=read(p/name);trans=read(tp/name)
  assert len(trans)==len(cc) and len(native)==len(rr)
  assert [x['orthogroup'] for x in native]==[x['orthogroup'] for x in rr]
  assert [x['column_id'] for x in trans]==[x['column_id'] for x in cc]
  assert all(same(native[j][c['column_id']],trans[i][r['orthogroup']]) for i,c in enumerate(cc) for j,r in enumerate(rr)),key
  cells+=len(rr)*len(cc)
 native=read(p/k['matrix_name']);st=read(p/'marker_state_matrix.tsv')
 points=read(R/'results/marker_display'/key/'rendered_points.tsv');mask=read(R/'results/marker_display'/key/'marker_display_mask.tsv')
 expected=set()
 for j,r in enumerate(rr):
  for i,c in enumerate(cc):
   yes=st[j][c['column_id']]=='up' and num(native[j][c['column_id']])>0
   assert (mask[i][r['orthogroup']]=='TRUE')==yes
   if yes:expected.add((c['column_id'],r['orthogroup']))
 assert len(points)==len(expected) and {(x['column_id'],x['orthogroup']) for x in points}==expected
 assert all(x['marker_state']=='up' and num(x['plot_log2FC'])>0 for x in points)
 g=read(R/'results/layout_annotations'/key/'render_geometry.tsv')[0]
 assert float(g['row_gap_mm'])==float(g['column_gap_mm'])==1.2
 assert float(g['marker_size_mm'])==(1.0 if k['overview']=='TRUE' else 1.2)
 assert g['annotation_cell_borders']=='TRUE' and g['annotation_border_color']=='#333333' and float(g['annotation_border_lwd_pt'])==float(g['border_lwd_pt'])==0.45 and g['annotation_border_lty']=='solid'
 assert all(float(g[x])%360==270 for x in ['column_names_rot','column_title_rot','conservation_text_rot','conservation_name_rot'])
 assert g['bottom_OG_text_just']=='left' and g['top_group_title_horizontal_just']=='1'
 assert g['actual_row_order_is_parent_column_order']=='TRUE' and g['actual_column_order_is_parent_row_order']=='TRUE'
 assert g['conservation_name_displayed']=='TRUE' and g['internal_grid_removed']=='TRUE'
 assert all(x in g['top_annotation_names'].split(';') for x in ['conservation','Cell_type'])
 assert all(g[x]=='Cell_type' for x in ['bottom_annotation_names','left_annotation_names','right_annotation_names'])
 pointcounts.append(dict(figure_key=key,n_points=len(points),n_down_points=0,n_original_up_down_points=g['parent_dot_count']))
for f in figs:
 p=Path(f['figure']);q=subprocess.run(['/usr/bin/gs','-q','-dNOPAUSE','-dBATCH','-sDEVICE=nullpage',str(p)],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
 pdfs.append(dict(figure=str(p),size_bytes=p.stat().st_size,status='PASS' if q.returncode==0 else 'FAIL',error=q.stderr.decode()[:1000]))
put(R/'results/pdf_validation.tsv',pdfs);put(R/'results/marker_display_counts.tsv',pointcounts)
check('PDFs_readable',len(figs)==28 and all(x['status']=='PASS' for x in pdfs),'28 PDFs; 26 with data and 2 empty-result panels')
check('exact_transposition',True,str(cells)+' values/states/source entries checked with original axes ordered')
check('only_up_marker_points',True,str(sum(x['n_points'] for x in pointcounts))+' displayed points validated, no downregulated points')
check('layout_and_actual_draw_order',True,'four-sided strips, equal 1.2 mm gaps, 270-degree axis text, conservation k/N and original sorting rules')
val=read(R/'results/all_species.member_contribution_validation.tsv')
check('native_member_reconstruction',len(val)==11 and all(float(x['raw_UMI_reconstruction_max_error'])==0 and float(x['normalized_mean_reconstruction_max_error'])<1e-8 for x in val),'all 11 species; RNA-to-OG UMI sums exact')
versions=read(R/'results/version_summary.tsv');bydir={x['direction']:x for x in versions}
# Current exclusion scope must be absent from every active data set and output path.
excluded=set(cfg['excluded_analysis_targets'])
for name in ['TF_OG.original_species_celltype_details.tsv.gz','TF_OG.comparison_original_combination_details.tsv.gz','TF_OG.comparison_species_votes.tsv.gz','TF_OG.comparison_conservation.tsv.gz','TF_OG.selected_associations.with_member_support.tsv','all_OG.comparison_conservation.tsv.gz','all_species.OG_TF_member_contribution.tsv.gz']:
 rows=read(R/'results'/name)
 assert not any(row.get('cell_type') in excluded or row.get('comparison_id') in excluded for row in rows),name
assert not any(any(ct in str(p.relative_to(R)) for ct in excluded) for p in (R/'figures').rglob('*.pdf'))
check('excluded_targets_absent_from_active_outputs',True,'no excluded target in statistics, votes, selected associations, member data or figures')
put(R/'results/output_validation.tsv',checks)
summary=dict(status='PASS',n_checks=len(checks),n_PDFs=len(figs),n_nonempty_PDFs=sum(int(x['n_OGs'])>0 for x in figs),n_up_OGs=int(bydir['up']['n_OGs']),n_down_OGs=int(bydir['down']['n_OGs']),n_unique_plotted_OGs=len(annotations),n_plot_rows=len(cols),n_original_combinations=len(units),n_comparisons=len(thresholds),n_TF80_candidates=len({x['orthogroup'] for x in read(R/'results/TF_OG.comparison_conservation.tsv.gz')}),n_down_marker_points=0,statistics_species=cfg['species_order'],reuse_comparisons=cfg['reuse_comparisons'])
(R/'results/final_validation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
manifest=[]
for p in sorted(R.rglob('*')):
 if p.is_file() and p.name!='output_file_manifest.tsv' and '__pycache__' not in p.parts:
  manifest.append(dict(relative_path=str(p.relative_to(R)),size_bytes=p.stat().st_size,sha256=sha(p)))
put(R/'results/output_file_manifest.tsv',manifest)
print(json.dumps(summary),flush=True)
