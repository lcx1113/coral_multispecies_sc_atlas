from pathlib import Path
import subprocess,sys,json
R=Path(__file__).resolve().parents[1];RS='Rscript'
(R/'logs').mkdir(exist_ok=True)
cfg=json.loads((R/'config/analysis_config.json').read_text())
def run(script,*args):
 cmd=([RS,'--vanilla'] if script.endswith('.R') else [sys.executable])+[str(R/'scripts'/script)]+list(args)
 suffix='_'+args[1] if script=='01_statistics.R' else ''
 with (R/'logs'/('runner_'+script+suffix+'.log')).open('w') as f:subprocess.run(cmd,cwd=R,stdout=f,stderr=subprocess.STDOUT,check=True)
 print('COMPLETE',script,flush=True)
for species in cfg['species_order']:run('01_statistics.R',str(R),species)
for s,args in [('02_select.R',(str(R),)),('03_plot_data.R',(str(R),)),('04_labels.py',()),('run_members.py',()),('07_plot.R',(str(R),)),('08_validate_statistics.R',(str(R),)),('09_validate_report.py',())]:run(s,*args)
