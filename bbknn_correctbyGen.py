# pip install -i https://pypi.tuna.tsinghua.edu.cn/simple datatable  (use Tsinghua mirror)

# ##############
import scanpy as sc
import os
import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.pyplot as plt
import scrublet as scr
import scanpy.external as sce

# Set working directory
import os
os.chdir("/work")
adata = sc.read_h5ad('input.h5ad')  # 4245220 × 18566
adata.obs['orig.ident'] = ['_'.join(x.split('_')[:6]) for x in adata.obs.index]
adata.obs['paper'] = ['_'.join(x.split('_')[:1]) for x in adata.obs.index]


# %% 1) Quality control
sc.pp.filter_cells(adata, min_genes = 200)                           
sc.pp.filter_cells(adata, max_genes = 6000)                          
sc.pp.filter_genes(adata, min_cells=3)                               
adata.var["mito"] = adata.var_names.str.startswith("MT-")            
sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)    
adata = adata[adata.obs['pct_counts_mito'] < 10, :]                  

# %% 2) Normalization and highly variable gene (HVG) selection
adata.layers["counts"] = adata.X.copy()                              
sc.pp.normalize_total(adata, target_sum=1e4)                       
sc.pp.log1p(adata)                                                   
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key="paper")
# Visualize HVGs
sc.pl.highly_variable_genes(adata)
plt.savefig('highly_variable_genes.png')
adata.raw = adata

# %% 3) PCA
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mito"])        
sc.pp.scale(adata, max_value=10)                                     
# PCA
sc.tl.pca(adata, svd_solver="arpack")
sc.pl.pca_variance_ratio(adata, log=True)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)                 
plt.savefig('pca_variance_ratio_log_scale.png')                      

# %% 4) Batch correction
# Run BBKNN 
sc.external.pp.bbknn(adata, batch_key='paper')

# %% 5) Clustering (Leiden)
sc.tl.leiden(adata, resolution=1)

# PAGA
sc.tl.paga(adata)
sc.pl.paga(adata, plot=False)
sc.tl.umap(adata, init_pos='paga')                                 
sc.pl.paga(adata, plot=True)
plt.savefig('paga_plot.png')


# %% 6) UMAP embedding and coloring by Leiden clusters
sc.tl.umap(adata)
sc.tl.rank_genes_groups(adata, 'leiden')
sc.pl.umap(adata, color='leiden')                                    
# sc.tl.umap(adata, min_dist=0.5, spread=0.1, copy=False)
plt.figure(figsize=(20, 15))
sc.pl.umap(adata, color='leiden')
plt.savefig("umap.png", bbox_inches='tight')

# --------- Marker visualization ---------
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
plt.savefig("marker.png")

top100_markers = pd.DataFrame()
# Iterate over each cluster
for group in adata.uns['rank_genes_groups']['names'].dtype.names:
    # Extract markers for each cluster
    group_markers = pd.DataFrame({
        'gene': adata.uns['rank_genes_groups']['names'][group],
        'scores': adata.uns['rank_genes_groups']['scores'][group],
        'logfoldchanges': adata.uns['rank_genes_groups']['logfoldchanges'][group],
        'pvals': adata.uns['rank_genes_groups']['pvals'][group],
        'pvals_adj': adata.uns['rank_genes_groups']['pvals_adj'][group]
    })
    # Keep upregulated genes only: logFC > 0.25, adjusted p < 0.05, positive score
    group_markers_filtered = group_markers[
        (group_markers['logfoldchanges'] > 0.25) &
        (group_markers['pvals_adj'] < 0.05) &
        (group_markers['scores'] > 0)
    ].copy()
    # Take top 100 by score
    group_markers_filtered = group_markers_filtered.nlargest(n=100, columns='scores')
    # Add cluster label
    group_markers_filtered['cluster'] = group
    # Append to global table
    top100_markers = pd.concat([top100_markers, group_markers_filtered], ignore_index=True)

# Reorder columns as follows:
top100_markers = top100_markers[['gene', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj', 'cluster']]
# Save to CSV
top100_markers.to_csv("top100_markers.csv", index=False)

# %% Save processed data to 'result.h5ad'
adata.write('result.h5ad')
adata.obs.to_csv("meta.csv")


# # %% Plot markers on UMAP
# sc.pl.umap(adata, color=['CD3D', 'CD3E', 'CD4','CD8A','CD8B','LYZ','CST3','ANK3','ZEB2','CD79A',"CD79B",'MS4A1','JCHAIN'])
# plt.savefig("Umap_marker.png")

# marker_genes = ['CD3D', 'CD3E', 'CD4','CD8A','CD8B','LYZ','CST3','ANK3','ZEB2','CD79A',"CD79B",'MS4A1','JCHAIN']
# sc.pl.dotplot(adata, marker_genes, groupby='leiden', dendrogram=True, show=False)
# plt.savefig("dotplot.jpg", dpi=300, bbox_inches='tight')

# # ############ label = TRUE #############
# plt.figure(figsize=(20, 15))
# sc.pl.umap(adata, color='leiden', legend_loc='on data', show=False)
# plt.savefig("umap2.png", bbox_inches='tight')

# plt.figure(figsize=(20, 15))
# sc.pl.umap(adata, color='celltype', legend_loc='on data', show=False)
# plt.savefig("umap_celltype.png", bbox_inches='tight')

# plt.figure(figsize=(20, 15))
# sc.pl.umap(adata, color='anno', legend_loc='on data', show=False)
# plt.savefig("umap_anno.jpg", bbox_inches='tight')

# plt.figure(figsize=(20, 15))
# sc.pl.umap(
#     adata,
#     color='anno',
#     legend_loc=None,        # Do not draw legend
#     title='',               # No title
#     show=False,
#     frameon=False           # Turn off axis frame
# )
# # Also make sure ticks/labels are hidden
# plt.axis('off')             # Hide axes completely
# plt.savefig("umap_anno_clean.jpg", dpi=300, bbox_inches='tight', pad_inches=0)


