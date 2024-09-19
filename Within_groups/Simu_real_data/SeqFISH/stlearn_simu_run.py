#conda activate stlearn
# pip install -U stlearn


import stlearn as st
import pandas as pd
# import scanpy as sc
# Ingore all warnings
import warnings
warnings.filterwarnings("ignore")


Perm_i=0

def paste0(*args, sep=""):
    # Join the elements in args using the specified separator (sep)
    return sep.join(map(str, args))

for FC in ["1.5", "2", "2.5", "3", "3.5"]:

    print(paste0("FC", FC))

    Count=pd.read_csv(paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/Simu_real_data/SeqFISH/count_SeqFISH_FC", FC, ".spa.txt"),sep=" ",header=0, na_filter=False, index_col=0)
    Count=Count.transpose()

    spatial_all=pd.read_csv(paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/Simu_real_data/SeqFISH/info_SeqFISH_200_FC", FC, "_Perm1000.txt"),sep="\t",header=None, na_filter=False, index_col=0)
    spatial_all=pd.DataFrame(spatial_all)

    for Perm_i in range(0, 1001): 
        print(paste0("Perm_i ############ ", Perm_i))

        spatial = spatial_all.iloc[:, (2*Perm_i):(2*Perm_i+2)]

        spatial.columns = ['imagecol', 'imagerow']

        adata = st.create_stlearn(count=Count, spatial=spatial, library_id="Sample_test", scale=1, background_color="white")

        # Preprocessing data
        st.pp.filter_genes(adata,min_cells=3)
        st.pp.normalize_total(adata)
        st.pp.log1p(adata)
        st.pp.scale(adata)

        # Run PCA
        st.em.run_pca(adata,n_comps=10,random_state=0)
        #st.pp.neighbors(adata,n_neighbors=5,use_rep='X_pca',random_state=0)
        #st.tl.clustering.louvain(adata,random_state=0)
        st.tl.clustering.kmeans(adata, n_clusters=2, use_data="X_pca", key_added="X_pca_kmeans")
        obsm=pd.DataFrame(adata.obs)

        outpath=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/Simu_real_data/SeqFISH/FC", FC, "/stlearn/stlearn_clusters_SeqFISH_FC", FC, "_Perm_", Perm_i, ".txt")
        #outpath=paste0("/home/qinf2/GuoshuaiCai/SPADE/stlearn_clusters_", shape, "_FC", FC, "_Perm_", Perm_i, ".txt")

        obsm.to_csv(outpath)

print("################  Done stlearn analysis #####################")


