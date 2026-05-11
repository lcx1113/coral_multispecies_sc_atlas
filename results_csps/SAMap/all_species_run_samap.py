#!/usr/bin/env python
# coding: utf-8

import os
os.environ["OMP_NUM_THREADS"] = "32"
# In[1]:


from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd


# ## Loading in raw data

# In[2]:


fn1  = './h5ad/dat.Adig.seurat_final.h5ad'
fn2  = './h5ad/dat.Ahya.seurat_final.h5ad'
fn3  = './h5ad/dat.Amur.seurat_final.h5ad'
fn4  = './h5ad/dat.Fspp.seurat_final.h5ad'
fn5  = './h5ad/dat.Gfas.seurat_final.h5ad'
fn6  = './h5ad/dat.Hcoe.seurat_final.h5ad'
fn7  = './h5ad/dat.Pcru.seurat_final.h5ad'
fn8  = './h5ad/dat.Pdam.seurat_final.h5ad'
fn9  = './h5ad/dat.Plut.seurat_final.h5ad'
fn10 = './h5ad/dat.Spis.seurat_final.h5ad'
fn11 = './h5ad/dat.Xesp.seurat_final.h5ad'


# In[3]:


filenames = {
    'ad': fn1,
    'ah': fn2,
    'am': fn3,
    'fs': fn4,
    'gf': fn5,
    'hc': fn6,
    'pc': fn7,
    'pd': fn8,
    'pl': fn9,
    'sp': fn10,
    'xe': fn11
}


# In[ ]:


sm = SAMAP(
        filenames,
        f_maps = '/share/home/Hu_lab_JX/workplace/samap/20260429/extend_longest_maps/',
        save_processed=True ,
        keys = {
            'ad':'cell_type',
            'ah':'cell_type',
            'am':'cell_type',
            'fs':'cell_type',
            'gf':'cell_type',
            'hc':'cell_type',
            'pc':'cell_type',
            'pd':'cell_type',
            'pl':'cell_type',
            'sp':'cell_type',
            'xe':'cell_type'
        },
    )
sm.run(neigh_from_keys = {
        'ad':True,
        'ah':True,
        'am':True,
        'fs':True,
        'gf':True,
        'hc':True,
        'pc':True,
        'pd':True,
        'pl':True,
        'sp':True,
        'xe':True
})


# In[ ]:


sm.run(pairwise=True)
samap = sm.samap # SAM object with three species stitched together


# In[ ]:


sm.scatter()


# In[ ]:


import pickle
with open('./sm_all_genome_species.pkl', 'wb') as file:
    pickle.dump(sm, file)
print("对象已保存到本地文件。")


# Post-SAMap analysis

# In[ ]:


import pickle


with open('./sm_all_genome_species.pkl', 'rb') as file:
    sm = pickle.load(file)

print("对象已成功加载。")


# In[ ]:


import scanpy as sc
import matplotlib.pyplot as plt
import networkx as nx
import umap
from scipy.sparse import csr_matrix
import seaborn as sns



# In[ ]:


samap = sm.samap
adata = samap.adata
adata


# In[ ]:


sc.pl.umap(adata, color=['species'], title='UMAP projection colored by species', size=10)

import scanpy as sc
import matplotlib.pyplot as plt


sc.pl.umap(
    adata, 
    color=['species'], 
    title='Combined Species UMAP', 
    size=5,            
    alpha=0.6,          
    frameon=False, 
    show=False,
    save='_combined_all.pdf' 
)


import os
import scanpy as sc
import matplotlib.pyplot as plt


output_dir = 'species_individual_plots'


if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"创建文件夹: {output_dir}")


if not isinstance(adata.obs['species'].dtype, pd.CategoricalDtype):
    adata.obs['species'] = adata.obs['species'].astype('category')

all_species = adata.obs['species'].cat.categories


for spec in all_species:
    print(f"正在导出物种: {spec} ...")
    
   
    fig, ax = plt.subplots(figsize=(7, 6))
    
    
    sc.pl.umap(
        adata, 
        color=['species'], 
        groups=[spec],     
        title=f'UMAP - Species: {spec}',
        size=15,          
        frameon=False,
        ax=ax,
        show=False         
    )
    

    save_path = os.path.join(output_dir, f'umap_{spec}.pdf')
    
    
    plt.savefig(save_path, format='pdf', bbox_inches='tight')
    
    
    plt.close(fig)

print("\n所有物种图片已导出完成！")


# ## Calculating cell type mapping scores

# sankey_plot

# In[ ]:


keys = {
            'ad':'cell_type',
            'ah':'cell_type',
            'am':'cell_type',
            'fs':'cell_type',
            'gf':'cell_type',
            'hc':'cell_type',
            'pc':'cell_type',
            'pd':'cell_type',
            'pl':'cell_type',
            'sp':'cell_type',
            'xe':'cell_type'
        }
D,MappingTable = get_mapping_scores(sm,keys,n_top = 0)

D.head()
MappingTable.head()


# In[ ]:


MappingTable.index


# In[ ]:


def export_mapping_table_for_sankey(mapping_table, output_file='mapping_table.csv'):
    """
    
    
    Parameters:
    mapping_table: pandas.DataFrame - SAMap生成的mapping table
    output_file: str - 输出CSV文件名
    """
    
    if not isinstance(mapping_table, pd.DataFrame):
        raise ValueError("mapping_table must be a pandas DataFrame")
    
    
    mapping_table.to_csv(output_file, index=True, header=True)
    print(f"Mapping table exported to {output_file}")
    print(f"Table shape: {mapping_table.shape}")
    print(f"Sample of the data:")
    print(mapping_table.iloc[:5, :5])
    
    return output_file


export_mapping_table_for_sankey(MappingTable, 'all_genome_species_samap_mapping_table.csv')








