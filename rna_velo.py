#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 18:50:17 2021

@author: Patrick
"""
#importing...

import os
import scvelo as scv

#setting WD
cwd = os.getcwd()
print("Current working directory: {0}".format(cwd))
os.chdir('/home/patrick/analysis_playground/velocyto/control_A/velocyto')

#Working

adata = scv.read("control_A_BM.h5ad")
adata

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)

scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

scv.tl.velocity(adata)

scv.tl.velocity_graph(adata)

scv.pl.velocity_embedding_stream(adata, basis="umap", color="seurat_clusters", dpi=600)

scv.pl.velocity_embedding(adata, basis="umap", color="seurat_clusters", arrow_length=4, arrow_size=4, dpi=600)

scv.tl.recover_dynamics(adata)

scv.tl.latent_time(adata)

scv.pl.scatter(adata, color="latent_time", color_map="gnuplot", dpi=600)

top_genes = adata.var["fit_likelihood"].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby="latent_time", col_color="seurat_clusters", n_convolve=100, figsize=(24, 12), col_cluster=(True))
