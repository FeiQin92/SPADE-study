import warnings
warnings.filterwarnings('ignore')

import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import math
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import SpaGCN as spg
#In order to read in image data, we need to install some package. Here we recommend package "opencv"
#inatll opencv in python
#!pip3 install opencv-python
import cv2
import os
import time
#import rpy2.robjects as robjects


spg.__version__

#from scanpy import read_10x_h5

def paste0(*args, sep=""):
    # Join the elements in args using the specified separator (sep)
    return sep.join(map(str, args))


spatial=pd.read_csv(paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/MOB/info_MOB_Perm1000.txt"),sep="\t",header=None, na_filter=False, index_col=0)

for Perm_i in range(0, 1): 
    print(paste0("Perm_i ############ ", Perm_i))
	
    adata=pd.read_csv(paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/MOB/count_MOB.spa.txt"),sep=" ",header=0, na_filter=False, index_col=0)
    adata = ad.AnnData(adata)
    adata=ad.AnnData.transpose(adata)

    adata.obs["x1"]=spatial[2*Perm_i+1]
    adata.obs["x2"]=spatial[2*Perm_i+2]
    adata.obs["x_array"]=adata.obs["x1"]
    adata.obs["y_array"]=adata.obs["x2"]

    img=cv2.imread("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/MOB/HE_Rep11_MOB.png")
    s=1
    b=2
    adj=spg.calculate_adj_matrix(x=adata.obs["x_array"],y=adata.obs["y_array"], x_pixel=adata.obs["x_array"], y_pixel=adata.obs["y_array"], image=img, beta=b, alpha=s, histology=False)

    #adata.var_names_make_unique()
    #spg.prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros
    #spg.prefilter_specialgenes(adata)
    #Normalize and take log for UMI
    #sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)

    p=0.5
    #Find the l value given p
    l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

    #For this toy data, we set the number of clusters=7 since this tissue has 7 layers
    n_clusters=2
    #Set seed
    r_seed=t_seed=n_seed=100
    #Search for suitable resolution
    res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)
    #res=0.7

    clf=spg.SpaGCN()
    clf.set_l(l)
    #Set seed
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)

    #Run
    clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
    y_pred, prob=clf.predict()
    adata.obs["pred"]= y_pred
    adata.obs["pred"]=adata.obs["pred"].astype('category')
    #Do cluster refinement(optional)
    #shape="hexagon" for Visium data, "square" for ST data.

    refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj, shape="hexagon")
    adata.obs["refined_pred"]=refined_pred
    adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
    #Save results

    #Use domain 0 as an example
    target=0
    #Set filtering criterials
    min_in_group_fraction=0.8
    min_in_out_group_ratio=1
    min_fold_change=1.5
    #Search radius such that each spot in the target domain has approximately 10 neighbors on average
    #adj_2d=spg.calculate_adj_matrix(x=adata.obs["x_array"], y=adata.obs["y_array"], histology=False)
    start, end= np.quantile(adj[adj!=0],q=0.001), np.quantile(adj[adj!=0],q=0.1)
    r=spg.search_radius(target_cluster=target, cell_id=adata.obs.index.tolist(), x=adata.obs["x_array"], y=adata.obs["y_array"],  pred=adata.obs["pred"].tolist(), start=start, end=end, num_min=10, num_max=14,  max_run=100)
    #Detect neighboring domains

    nbr_domians=spg.find_neighbor_clusters(target_cluster=target,
                                   cell_id=adata.obs.index.tolist(),
                                   x=adata.obs["x_array"].tolist(),
                                   y=adata.obs["y_array"].tolist(),
                                   pred=adata.obs["pred"].tolist(),
                                   radius=r,
                                   ratio=1/2)

    nbr_domians=nbr_domians[0:3]
    de_genes_info=spg.rank_genes_groups(input_adata=adata,
                                target_cluster=target,
                                nbr_list=nbr_domians,
                                label_col="pred",
                                adj_nbr=True,
                                log=True)

    #adata_list.append(de_genes_info)
                  
    # Convert the list of adata objects to an R list
    #r_adata_list = robjects.ListVector({f"de_genes_info_{i}": adata_item for i, adata_item in enumerate(adata_list)})

    outpath=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/MOB/SpaGCN/SpaGCN_MOB_Perm_", Perm_i, ".txt")
    de_genes_info.to_csv(outpath, sep='\t', index=False)

    # Save the list of adata objects as a .RData file
    # robjects.r['save'](r_adata_list, file=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/MOB/FC", FC, "/SpaGCN/SpaGCN_MOB_FC",FC,".RData"))


print("################  Done SpaGCN analysis #####################")



