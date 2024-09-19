
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

sim_trend <-  splatSimulate(params=params)
sim_trend_RC <- counts(sim_trend)
# save(sim_trend_RC, file="E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\hotspot\\simu.data.400Cells.RData")
sim_trend_RC <- get(load("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\hotspot\\simu.data.400Cells.RData"))
sim_trend_RC1 <- sim_trend_RC[apply(sim_trend_RC>1, 1, mean) >= 0.5,]
dim(sim_trend_RC1)

numGenes = 200
lowSignal = 25
highSignal = 25

numCell = 200
## Streak pattern with different directions and same coordinates
for (cellProp in c(0.1, 0.2, 0.3)){
  for (ieffect in c(1.5, 2.0, 2.5, 3.0, 3.5, 10)){
    # itheta = 10
    set.seed(31)
    sim_trend_RC400 <- sim_trend_RC1[sample(1:nrow(sim_trend_RC1), 400),]
    
    set.seed(31)
    ppT <- rpoispp(0.00020,  win=owin(c(1,1000),c(1,1000)))
    
    ppC <- ppT
    ppC$x <- ppT$y
    ppC$y <- ppT$x
    
    low_expr    <- c(1)
    high_expr   <- c(2)
    
    set.seed(31)
    pp_T = add_markdist_streak(ppT,low_marks=low_expr, high_marks=high_expr, nMarkers=(lowSignal+highSignal), ncell=numGenes, cell_proportion=cellProp)
    pp_C = add_markdist_streak_h(ppT,low_marks=low_expr, high_marks=high_expr, nMarkers=(lowSignal+highSignal), ncell=numGenes, cell_proportion=cellProp)
    
    expr_mat_C <- matrix(NA,nrow=pp_C$n,ncol=numGenes)
    expr_mat_T <- matrix(NA,nrow=pp_T$n,ncol=numGenes)
    for(igene in 1:numGenes){
      
      spiked_in_cells_T     <- which(pp_T$markformat[,igene]==2)
      spiked_in_cells_C     <- which(pp_C$markformat[,igene]==2)
      
      a <- sim_trend_RC400[igene,]
      if(igene <= highSignal){
        a[spiked_in_cells_T]  <- sim_trend_RC400[igene,spiked_in_cells_T]*ieffect
      }else if(igene <=lowSignal+highSignal) {
        a[spiked_in_cells_T]  <- sim_trend_RC400[igene,spiked_in_cells_T]/ieffect
      }
      
      b <- sim_trend_RC400[igene,]
      if(igene <= highSignal){
        # set.seed(igene+2000)
        b[spiked_in_cells_C] <- sim_trend_RC400[igene,spiked_in_cells_C]*ieffect
      }else if(igene <=lowSignal+highSignal) {
        # set.seed(igene+2000)
        b[spiked_in_cells_C] <- sim_trend_RC400[igene,spiked_in_cells_C]/ieffect
      }
      
      expr_mat_T[,igene] <- a
      expr_mat_C[,igene] <- b
      

    }
    
    colnames(expr_mat_T) <- colnames(expr_mat_C)  <- paste0("gene",1:numGenes)
    info_T                <- cbind.data.frame(x=ppT$x,y=ppT$y)
    info_T <- round(info_T)
    
    info_C                <- cbind.data.frame(x=ppC$x,y=ppC$y)
    info_C <- round(info_C)
    
    rownames(info_T)      <- rownames(expr_mat_T) <- rownames(expr_mat_C) <- paste0("C",1:nrow(expr_mat_T))
    rownames(info_C)      <- rownames(expr_mat_T) <- rownames(expr_mat_C) <- paste0("C",1:nrow(expr_mat_T))
    
    count_T=t(expr_mat_T)
    count_C=t(expr_mat_C)
    
    write.table(pp_T$markformat,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streak\\Cell_prop_",cellProp,"\\count_streak_200_FC_hr_indexT", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
    write.table(pp_C$markformat,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streak\\Cell_prop_",cellProp,"\\count_streak_200_FC_hr_indexC", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
    
    write.table(count_T,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streak\\Cell_prop_",cellProp,"\\count_streak_200_FC_hr_T", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
    write.table(count_C,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streak\\Cell_prop_",cellProp,"\\count_streak_200_FC_hr_C", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
    
    write.table(info_T,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streak\\Cell_prop_",cellProp,"\\info_streak_200_FC_hr_T", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
    write.table(info_T,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streak\\Cell_prop_",cellProp,"\\info_streak_200_FC_hr_C", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
    }
}





## Streak pattern with the same directions and same coordinates
for (cellProp in c(0.1, 0.2, 0.3)){
  for (ieffect in c(1.5, 2.0, 2.5, 3.0, 3.5, 10)){
    # itheta = 10
    set.seed(31)
    sim_trend_RC400 <- sim_trend_RC1[sample(1:nrow(sim_trend_RC1), 400),]

    set.seed(31)
    pp <- rpoispp(0.00020,  win=owin(c(1,1000),c(1,1000)))
    low_expr    <- c(1)
    high_expr   <- c(2)

    set.seed(31)
    pp_T = add_markdist_streak2(pp,low_marks=low_expr, high_marks=high_expr, nMarkers=(lowSignal+highSignal), 
                                ncell=numGenes, cell_proportion=cellProp, seed=31)

    set.seed(69)
    pp_C = add_markdist_streak2(pp,low_marks=low_expr, high_marks=high_expr, nMarkers=(lowSignal+highSignal), 
                                ncell=numGenes, cell_proportion=cellProp, seed=69)

    expr_mat_C <- matrix(NA,nrow=pp_C$n,ncol=numGenes)
    expr_mat_T <- matrix(NA,nrow=pp_T$n,ncol=numGenes)
    for(igene in 1:numGenes){

      spiked_in_cells_T     <- which(pp_T$markformat[,igene]==2)
      spiked_in_cells_C     <- which(pp_C$markformat[,igene]==2)

      a <- sim_trend_RC400[igene,]
      if(igene <= highSignal){
        a[spiked_in_cells_T]  <- sim_trend_RC400[igene,spiked_in_cells_T]*ieffect
      }else if(igene <=lowSignal+highSignal) {
        a[spiked_in_cells_T]  <- sim_trend_RC400[igene,spiked_in_cells_T]/ieffect
      }

      b <- sim_trend_RC400[igene,]
      if(igene <= highSignal){
        # set.seed(igene+2000)
        b[spiked_in_cells_C] <- sim_trend_RC400[igene,spiked_in_cells_C]*ieffect
      }else if(igene <=lowSignal+highSignal) {
        # set.seed(igene+2000)
        b[spiked_in_cells_C] <- sim_trend_RC400[igene,spiked_in_cells_C]/ieffect
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

    write.table(pp_T$markformat,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streak\\Cell_prop_",cellProp,"\\count_streak_200_FC_hr_samedir_indexT", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
    write.table(pp_C$markformat,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streak\\Cell_prop_",cellProp,"\\count_streak_200_FC_hr_samedir_indexC", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
    
    write.table(count_T,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streak\\Cell_prop_",cellProp,"\\count_streak_200_FC_hr_samedir_T", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
    write.table(count_C,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streak\\Cell_prop_",cellProp,"\\count_streak_200_FC_hr_samedir_C", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")

    write.table(info,file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streak\\Cell_prop_",cellProp,"\\info_streak_200_FC_hr", as.character(ieffect), ".txt"),row.names = F,col.names = F,quote = F,sep="\t")
  }
}




library(SPADE)
### exchange both length scale theta and gamma values
#########################################################
#########################################################
#########################################################
cellProp=0.3
FC=3.5
for (cellProp in c(0.1, 0.2, 0.3)){
  for (FC in c(1.5, 2.0, 2.5, 3.0, 3.5)){
    print(paste0("cellProp=", cellProp, "; FC=", FC))
 
    readcounts_T <- as.matrix(read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streak\\Cell_prop_",cellProp,"\\count_streak_200_FC_hr_T", FC, ".txt")))
    readcounts_C <- as.matrix(read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streak\\Cell_prop_",cellProp,"\\count_streak_200_FC_hr_C", FC, ".txt")))
    
    info_T <- read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streak\\Cell_prop_",cellProp,"\\info_streak_200_FC_hr_T", FC, ".txt"))
    rownames(info_T) <- colnames(readcounts_T) <- colnames(readcounts_C) <- paste0("C",1:ncol(readcounts_T))
    
    info_C <- read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streak\\Cell_prop_",cellProp,"\\info_streak_200_FC_hr_C", FC, ".txt"))
    rownames(info_C) <- colnames(readcounts_T) <- colnames(readcounts_C) <- paste0("C",1:ncol(readcounts_T))
    
    rownames(readcounts_T) <- rownames(readcounts_C) <-  paste0("G",1:nrow(readcounts_T))
    ED <- as.matrix(dist(info_T))
    lrang <- SPARK::ComputeGaussianPL(ED, compute_distance=FALSE)
    
    info_C <- info_T
    info_C$V1 <- info_T$V2
    info_C$V2 <- info_T$V1

    library(spatialDE)
    regdata_T <- SPADE_norm(readcounts=as.matrix(readcounts_T), info=info_T)
    regdata_C <- SPADE_norm(readcounts=as.matrix(readcounts_C), info=info_T)
    res <- SPADE_DE(readcounts1=regdata_T, readcounts2=regdata_C, location1=info_T, location2=info_T, mode="Shape&Strength")

    save(res, file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streakEX\\Cell_prop_",cellProp,"\\SPDE_estimate_hr_",cellProp,"_", FC, "EX.diffDir_ShapeStrength.RData" ))
    
  }
}





cellProp=0.2
FC=2.0
for (cellProp in c(0.1, 0.2, 0.3)){
  for (FC in c(1.5, 2.0, 2.5, 3.0, 3.5)){
    print(paste0("cellProp=", cellProp, "; FC=", FC))
  
    readcounts_T <- as.matrix(read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streak\\Cell_prop_",cellProp,"\\count_streak_200_FC_hr_sameDir_T", FC, ".txt")))
    readcounts_C <- as.matrix(read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streak\\Cell_prop_",cellProp,"\\count_streak_200_FC_hr_sameDir_C", FC, ".txt")))
    
    info <- read.table(file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streak\\Cell_prop_",cellProp,"\\info_streak_200_FC_hr", FC, ".txt"))
    rownames(info) <- colnames(readcounts_T) <- colnames(readcounts_C) <- paste0("C",1:ncol(readcounts_T))
    rownames(readcounts_T) <- rownames(readcounts_C) <-  paste0("G",1:nrow(readcounts_T))
    ED <- as.matrix(dist(info))
    lrang <- SPARK::ComputeGaussianPL(ED, compute_distance=FALSE)
    
    library(spatialDE)
    regdata_T <- SPADE_norm(readcounts=as.matrix(readcounts_T), info=info)
    regdata_C <- SPADE_norm(readcounts=as.matrix(readcounts_C), info=info)
    res <- SPADE_DE(regdata_T, regdata_C, info, info, mode="Shape&Strength")

    save(res, file=paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streakEX\\Cell_prop_",cellProp,"\\SPDE_estimate_hr_",cellProp,"_", FC, "EX.sameDir_ShapeStrength.RData" ))
    
  }
}






FPR.all <- NULL
for (cellProp in c(0.1, 0.2, 0.3)){
  for (FC in c(1.5, 2.0, 2.5, 3.0, 3.5)){
    print(paste0("cellProp=", cellProp, "; FC=", FC))
    # the first 50/200 genes are marker genes
    truth=c(rep(0, 50), rep(1, 150))
    
    SPDE_estimate_hr <- get(load(paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streakEX\\Cell_prop_",cellProp,"\\SPDE_estimate_hr_",cellProp,"_", FC, "EX.sameDir_ShapeStrength.RData" )))
    FPR_sameDir <- sum(SPDE_estimate_hr$Adjust.Pvalue < 0.05)/length(truth)
    
    SPDE_estimate_hr <- get(load(paste0("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\simulation\\DEanalysis\\streakEX\\Cell_prop_",cellProp,"\\SPDE_estimate_hr_",cellProp,"_", FC, "EX.diffDir_ShapeStrength.RData" )))
    FPR_diffDir <- sum(SPDE_estimate_hr$Adjust.Pvalue < 0.05)/length(truth)
    
    FPR.all1 <- rbind(data.frame(FPR=c(FPR_sameDir, FPR_diffDir), Methods=c("Same_dir", "Diff_dir")))
    
    FPR.all1$CellProp <- cellProp
    FPR.all1$FC <- FC
    FPR.all <- rbind(FPR.all, FPR.all1)
  }
}




