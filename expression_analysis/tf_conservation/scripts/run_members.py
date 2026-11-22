from pathlib import Path
import subprocess,json,concurrent.futures
R=Path(__file__).resolve().parents[1];cfg=json.loads((R/'config/analysis_config.json').read_text())
RS='Rscript'
def run(s):
 with (R/f'logs/05_members_{s}.log').open('w') as f:
  subprocess.run([RS,'--vanilla',str(R/'scripts/05_member_contribution.R'),str(R),s],stdout=f,stderr=subprocess.STDOUT,check=True)
 print('MEMBER COMPLETE',s,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(run,cfg['species_order']))
with (R/'logs/06_finalize.log').open('w') as f:subprocess.run([RS,'--vanilla',str(R/'scripts/06_finalize.R'),str(R)],stdout=f,stderr=subprocess.STDOUT,check=True)
print('ALL MEMBER EVIDENCE COMPLETE',flush=True)
