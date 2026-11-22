#!/usr/bin/env python3
"""Record this run's installed package versions and executable fingerprints."""
import json
import subprocess
from pathlib import Path
from common import PROJECT_ROOT, load_config, sha256_file, write_tsv

cfg, _ = load_config("config/analysis_config.yaml")
rows = []
for env_key in ("gene_tree_prefix","generax_prefix","r_prefix"):
    prefix = Path(cfg["environments"][env_key])
    for file in sorted((prefix / "conda-meta").glob("*.json")):
        rec = json.loads(file.read_text())
        if rec["name"] in {"mafft","clipkit","iqtree","generax","openmpi","python","biopython","r-base",
                           "r-ggplot2","bioconductor-ggtree","r-ape","r-ggnewscale","r-scales","r-digest","r-yaml"}:
            rows.append(dict(environment=str(prefix),software=rec["name"],version=rec["version"],
                             build=rec.get("build",""),evidence=str(file),sha256=sha256_file(file),
                             evidence_type="conda_metadata_not_necessarily_loaded_R_library",exit_code=""))
for name, executable in [
    ("MAFFT",cfg["alignment"]["executable"]),("ClipKIT",cfg["trimming"]["executable"]),
    ("IQ-TREE",cfg["tree"]["executable"]),("GeneRax",cfg["environments"]["generax_executable"]),
    ("OpenMPI",cfg["environments"]["generax_mpiexec"]),
]:
    path = Path(executable)
    version = subprocess.run([str(path),"--version"],capture_output=True,text=True,timeout=20)
    log = PROJECT_ROOT / "logs" / (name.lower().replace("-","_") + ".version.txt")
    log.write_text(version.stdout + version.stderr)
    rows.append(dict(environment=str(path.parent.parent),software=name,version="see native version log",
                     build="",evidence=str(log),sha256=sha256_file(path),
                     evidence_type="native_command",exit_code=version.returncode))
r_expr = 'for(p in c("ggtree","ggplot2","ape","ggnewscale","scales","digest","yaml")) cat(paste("ACTUAL_PACKAGE",p,as.character(packageVersion(p)),find.package(p),sep="\\t"),"\\n",sep="")'
loaded = subprocess.run([cfg["environments"]["rscript"],"-e",r_expr],capture_output=True,text=True,timeout=60)
r_log = PROJECT_ROOT / "logs/r_loaded_packages.txt"
r_log.write_text(loaded.stdout + loaded.stderr)
for line in loaded.stdout.splitlines():
    if line.startswith("ACTUAL_PACKAGE\t"):
        _, name, version, location = line.split("\t")
        rows.append(dict(environment=cfg["environments"]["r_prefix"],software=name,version=version,
                         build="",evidence=location,sha256=sha256_file(r_log),
                         evidence_type="R_packageVersion_and_find.package",exit_code=loaded.returncode))
write_tsv(PROJECT_ROOT / "logs/software_versions.tsv",rows,
          ["environment","software","version","build","evidence","sha256","evidence_type","exit_code"])
print("Software version capture PASS",flush=True)
