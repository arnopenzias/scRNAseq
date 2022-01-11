#import numpy as np
import pandas as pd
import scanpy as sc
#import matplotlib
#import leidenalg

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, facecolor='white')

results_file = 'write/zikv.h5ad'  # the file that will store the analysis results

adata = sc.read_10x_mtx(
    'Zikv ctrl/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)    

adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`

adata

sc.pl.highest_expr_genes(adata, n_top=20) #ver os genes mais expressos, nota-se o ruido dos mt e malat

#filtro básico

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
#aqui vc ve mais fácil os treshholds dos genes pra vc cortar eles fora

#filtrar os outliers cortando o objeto
adata = adata[adata.obs.n_genes_by_counts < 7000, :]
adata = adata[adata.obs.pct_counts_mt < 24, :]

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
#normalizar e logaritimizar pra poder comparar as células

#ver os genes mais expressos depois desse pré processamento
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

#You can get back an AnnData of the object in .raw by calling .raw.to_adata().

adata.raw = adata

#The result of the previous highly-variable-genes detection is stored as an annotation in .var.highly_variable and auto-detected by PCA and hence, sc.pp.neighbors and subsequent manifold/graph tools

#Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance.

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

#Scale each gene to unit variance. Clip values exceeding standard deviation 10.

sc.pp.scale(adata, max_value=10)

#fazer o pca

sc.tl.pca(adata, svd_solver='arpack')

sc.pl.pca_variance_ratio(adata, log=True)

adata.write(results_file)
adata

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
#Let us compute the neighborhood graph of cells using the PCA representation of the data matrix. You might simply use default values here. For the sake of reproducing Seurat’s results, let’s take the following values.

sc.tl.umap(adata)

sc.pl.umap(adata, color=['Dnm1l', 'Fis1', 'Mff'])

sc.pl.umap(adata, color=['Mfn1', 'Mfn2', 'Opa1'])


sc.tools.leiden(adata, key_added="clusters")


sc.pl.umap(adata, color=['clusters'])
adata.write(results_file)

#finding marker genes and stats
sc.tl.rank_genes_groups(adata, 'clusters', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
adata.write(results_file)


adata = sc.read(results_file)
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)



result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(5)

#comparar genes contra os grupos
adata = sc.read(results_file)
sc.pl.violin(adata, ['Dnm1l', 'Fis1', 'Mff'], groupby='clusters')


######## daqui pra frente deu ruim !!!!

#marcar os clusters
#criar os genes pra marcar
marker_genes = ['Id3', 'Clu', 'Rpl32', 'Egfr', 'Dlx2', 'Nes', 'Sox2,' ,'Gfap', 'Rbfox3',  'Cdk4', 'Cdk1', 'Tipin', 'Cenpa', 'Pcna', 'Mki67', 'Top2a', 'Mcm5']
new_cluster_names = [
    'Id3', 'Clu', 'Rpl32', 'Egfr', 'Dlx2', 'Nes', 'Sox2', 'Sox2', 'Gfap', 'Rbfox3' 'NSCs',
    'B', 'CD8 T',
    'NK', 'FCGR3A Monocytes',
    'Dendritic', 'Megakaryocytes']
adata.rename_categories('leiden', new_cluster_names)























