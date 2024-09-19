import pandas as pd
import NaiveDE
import SpatialDE

shape="hotspot"
Cellprop="0.2"
FC="1.4"
Perm_i=0

def paste0(*args, sep=""):
    # Join the elements in args using the specified separator (sep)
    return sep.join(map(str, args))

counts=pd.read_csv(paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/MOB/count_MOB.spa.txt"),sep=" ",header=0, na_filter=False, index_col=0)
counts=counts.transpose()

spatial_all=pd.read_csv(paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/MOB/info_MOB_Perm1000.txt"),sep="\t",header=None, na_filter=False, index_col=0)
spatial_all=pd.DataFrame(spatial_all)

for Perm_i in range(400, 601): 
    print(paste0("Perm_i ############ ", Perm_i))

    sample_info = spatial_all.iloc[:, (2*Perm_i):(2*Perm_i+2)]
    sample_info.columns = ['x', 'y']
    sample_info = sample_info[['x', 'y']].copy()
    sample_info.loc[:, "total_counts"] =  counts.sum(axis=1)

    norm_expr = NaiveDE.stabilize(counts.T).T
    resid_expr = NaiveDE.regress_out(sample_info, norm_expr.T, 'np.log(total_counts)').T
    sample_resid_expr = resid_expr.sample(n=100, axis=1, random_state=1)
    X = sample_info[['x', 'y']]
    X = pd.DataFrame.to_numpy(X)

    results = SpatialDE.run(X, resid_expr)

    outpath=paste0("/gpfs/gsfs12/users/qinf2/GuoshuaiCai/SPADE/real_data_single_group/MOB/SpatialDE/SpatialDE_MOB_Perm_", Perm_i, ".txt")

    #outpath=paste0("/home/qinf2/GuoshuaiCai/SPADE/SpatialDE_", shape, "_FC", FC, "_Perm_", Perm_i, ".txt")

    results.to_csv(outpath)

print("################  Done spatialDE analysis #####################")


