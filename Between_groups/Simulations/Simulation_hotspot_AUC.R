
library("dplyr")
library("gstat")
library("sp")
library("spacetime")
library("STRbook")
library("tidyr")

library(spatstat)
library(Matrix)
library(MASS)
library(FRK)
source("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\Code\\utilities.R")


setwd("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\data\\SPARK\\raw_data")

MOB <- t(read.table("Rep11_MOB_count_matrix-1.tsv"))
MOB1 <- MOB[,apply(MOB, 2, sum) >= 10]
MOB2 <- MOB1[(apply(MOB1>0, 1, mean) >= 0.3) ,]
set.seed(2022)
MOB3 <- MOB2[, sample(1:ncol(MOB2), 200)]

library(splatter)
# library(SCRIP)
params <- splatEstimate(MOB3)

sim_trend <-  splatSimulate(params=params, batchCells =400)
sim_trend_RC <- counts(sim_trend)
# save(sim_trend_RC, file="E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\simu.data.400Cells.RData")
sim_trend_RC <- get(load("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\simu.data.400Cells.RData"))
sim_trend_RC1 <- sim_trend_RC[apply(sim_trend_RC>1, 1, mean) >= 0.5,]
dim(sim_trend_RC1)

numGenes = 200
lowSignal = 25
highSignal = 25

numCell = 200

for (cellProp in c(0.1, 0.2, 0.3)){
  for (ieffect in 0.5*c(3:8)){
    # itheta = 10
    set.seed(31)
    sim_trend_RC400 <- sim_trend_RC1[sample(1:nrow(sim_trend_RC1), 200), sample(1:ncol(sim_trend_RC1), 400)]
    
    set.seed(31)
    pp <- rpoispp(0.00020,  win=owin(c(1,1000),c(1,1000)))
    low_expr    <- c(1)
    high_expr   <- c(2)
    
    pp = add_markdist_hotspot(pp,low_marks=low_expr, high_marks=high_expr, nMarkers=(lowSignal+highSignal), ncell=numCell, cell_proportion=cellProp)
    
    expr_mat_C <- matrix(NA,nrow=pp$n,ncol=numGenes)
    expr_mat_T <- matrix(NA,nrow=pp$n,ncol=numGenes)
    for(igene in 1:numGenes){
      
      back_ground_cells   <- which(pp$markformat[,igene]==1)
      spiked_in_cells     <- which(pp$markformat[,igene]==2)
      
      sim_trend_RC200_1 <- sim_trend_RC400[,1:200]
      a <- sim_trend_RC200_1[igene,]
      if(igene <= highSignal){
        a[spiked_in_cells]  <- sim_trend_RC200_1[igene, spiked_in_cells]*ieffect
      }else if(igene <=lowSignal+highSignal) {
        a[spiked_in_cells]  <- sim_trend_RC200_1[igene, spiked_in_cells]/ieffect
      }
    
      sim_trend_RC200_2 <- sim_trend_RC400[,201:400]
      b <- sim_trend_RC200_2[igene,]
      set.seed(igene)
      spiked_in_cells_R <- sample(1:pp$n, length(spiked_in_cells))
      if(igene <= highSignal){
        b[spiked_in_cells_R] <- sim_trend_RC200_2[igene, spiked_in_cells_R]*ieffect
      }else if(igene <=lowSignal+highSignal) {
        b[spiked_in_cells_R] <- sim_trend_RC200_2[igene, spiked_in_cells_R]/ieffect
      }
      
      expr_mat_T[,igene] <- a
      expr_mat_C[,igene] <- b
    }
    
    colnames(expr_mat_T) <- colnames(expr_mat_C)  <- paste0("gene",1:numGenes)
    info                <- cbind.data.frame(x=pp$x,y=pp$y)
    info <- round(info)
    
    rownames(info)      <- rownames(expr_mat_T) <- rownames(expr_mat_C) <- paste0("C",1:nrow(expr_mat_T))
    count_T=t(expr_mat_T)
    count_C=t(expr_mat_C)
    
    write.table(pp$markformat,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspot\\Cell_prop_",cellProp,"\\count_hotspot_200_FC_hr_indexT", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
    
    write.table(count_T,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspot\\Cell_prop_",cellProp,"\\count_hotspot_200_FC_hr_T", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
    write.table(count_C,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspot\\Cell_prop_",cellProp,"\\count_hotspot_200_FC_hr_C", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
    
    write.table(info,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspot\\Cell_prop_",cellProp,"\\info_hotspot_200_FC_hr", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
  }
}






### exchange both length scale theta and gamma values
#########################################################
#########################################################
cellProp=0.2
FC=2.0
for (cellProp in c(0.1, 0.2, 0.3)){
  for (FC in 0.5*c(3:8)){
    print(paste0("cellProp=", cellProp, "; FC=", FC))
    # the first 100/300 genes are marker genes
    readcounts_T <- as.matrix(read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspot\\Cell_prop_",cellProp,"\\count_hotspot_200_FC_hr_T", FC, ".txt")))
    readcounts_C <- as.matrix(read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspot\\Cell_prop_",cellProp,"\\count_hotspot_200_FC_hr_C", FC, ".txt")))
    
    info <- read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspot\\Cell_prop_",cellProp,"\\info_hotspot_200_FC_hr", FC, ".txt"))
    rownames(info) <- colnames(readcounts_T) <- colnames(readcounts_C) <- paste0("C",1:ncol(readcounts_T))
    rownames(readcounts_T) <- rownames(readcounts_C) <-  paste0("G",1:nrow(readcounts_T))
    ED <- as.matrix(dist(info))
    lrang <- SPARK::ComputeGaussianPL(ED, compute_distance=FALSE)
    
    library(SPADE)
    regdata_T <- SPADE_norm(readcounts=as.matrix(readcounts_T), info=info)
    regdata_C <- SPADE_norm(readcounts=as.matrix(readcounts_C), info=info)
    res <- SPADE_DE(regdata_T, regdata_C, info, info, mode="Shape&Strength")
   
    save(res, file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspot\\Cell_prop_",cellProp,"\\SPDE_estimate_hr_",cellProp,"_", FC, "EX_ShapeStrength.RData" ))
    
  }
}



AUC.all <- NULL
for (cellProp in c(0.1, 0.2, 0.3)){
  for (FC in 0.5*c(3:8)){
    print(paste0("cellProp=", cellProp, "; FC=", FC))
    # the first 50/200 genes are marker genes
    truth=c(rep(0, 50), rep(1, 150))
    res_DE <- get(load(paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspot\\Cell_prop_",cellProp,"\\Other_DE_methods_",cellProp,"_", FC, ".RData" )))
    AUC.other <- AUC.gen(res=res_DE, truth=truth, nDE=50)
    
    SPDE_estimate_hr <- get(load(paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\hotspot\\Cell_prop_",cellProp,"\\SPDE_estimate_hr_",cellProp,"_", FC, "EX_ShapeStrength.RData" )))
 
    pred <- prediction(SPDE_estimate_hr$Adjust.Pvalue, truth)
    auc_ROCR <- performance(pred, measure = "auc")
    auc_ROCR_SPDE_ex <- auc_ROCR@y.values[[1]]
    auc_ROCR_SPDE_ex
    
    AUC.all1 <- rbind(AUC.other, data.frame(AUC=auc_ROCR_SPDE_ex, Methods="SPADE"))
    
    AUC.all1$CellProp <- cellProp
    AUC.all1$FC <- FC
    AUC.all <- rbind(AUC.all, AUC.all1)
  }
}







