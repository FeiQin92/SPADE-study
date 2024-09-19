# library(curl)
# for (filename in c(
#              "Adult.rds", "Control_Juv.rds", "Meta.rds",
#                   "Stage44.rds", "Stage54.rds", "Stage57.rds")){
#   url <- paste0("https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000056/RDS/", filename)
#   destfile <- paste0("E:/Stereo-seq/", filename)
#   curl_download(url =url ,destfile=destfile,
#                 quiet=FALSE, mode="wb")
#   }


# url <- ("https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000056/RDS/5DPI_1.rds")
# destfile <- paste0("E:/Stereo-seq/", "5DPI_1.rds")
# curl_download(url =url ,destfile=destfile,
#               quiet=FALSE, mode="wb")



library(Seurat)
DPI_2D_1 <- readRDS("E:\\Stereo-seq\\2DPI_1.rds")
DPI_2D_1_data <- data.matrix(DPI_2D_1@assays$Spatial@counts)
DPI_2D_1_info <- DPI_2D_1@images$Injury_2DPI_rep1_SS200000147BL_D5@coordinates[,c(3, 2)]
dim(DPI_2D_1_data)
dim(DPI_2D_1_info)

DPI_2D_1_data1 <- DPI_2D_1_data[apply(DPI_2D_1_data>0, 1, mean)>0.2,]
rm(DPI_2D_1)
rm(DPI_2D_1_data)
DPI_2D_1_data2 <- DPI_2D_1_data1[,apply(DPI_2D_1_data1>0, 2, sum)>1000]
DPI_2D_1_info1 <- -DPI_2D_1_info
DPI_2D_1_info2 <- DPI_2D_1_info1[apply(DPI_2D_1_data1>0, 2, sum)>1000,]



DPI_5D_1 <- readRDS("E:\\Stereo-seq\\5DPI_1.rds")
DPI_5D_1_data <- data.matrix(DPI_5D_1@assays$Spatial@counts)
DPI_5D_1_info <- DPI_5D_1@images$Injury_5DPI_rep1_SS200000147BL_D2@coordinates[, c(2, 3)]

dim(DPI_5D_1_data)
dim(DPI_5D_1_info)

DPI_5D_1_data1 <- DPI_5D_1_data[apply(DPI_5D_1_data>0, 1, mean)>0.2,]
rm(DPI_5D_1)
rm(DPI_5D_1_data)
DPI_5D_1_data2 <- DPI_5D_1_data1[,apply(DPI_5D_1_data1>0, 2, sum)>1000]
DPI_5D_1_info1 <- DPI_5D_1_info
DPI_5D_1_info1$row <- -DPI_5D_1_info$row
DPI_5D_1_info2 <- DPI_5D_1_info1[apply(DPI_5D_1_data1>0, 2, sum)>1000,]

dim(DPI_5D_1_data2)
dim(DPI_5D_1_info2)


tables5<- as.matrix(readxl::read_excel("E:\\Stereo-seq\\science.abp9444_tables_s1_to_s7\\science.abp9444_table_s5.xlsx"))
tables <- readxl::read_excel("E:\\Stereo-seq\\science.abp9444_tables_s1_to_s7\\science.abp9444_table_s2.xlsx")
tables1 <- as.matrix(tables[tables$`Slice Name` %in% c("2DPI-1", "10DPI-1"), 11:30])
genelist <- NULL
for (i in c(1:ncol(tables1))){
  genelist <- c(genelist, as.character(tables1[,i]))
}


## Extract overlap genes between two groups
over_gene <- intersect(rownames(DPI_2D_1_data2), rownames(DPI_5D_1_data2))

DPI_2D_1_data <- DPI_2D_1_data2[over_gene, ]
DPI_5D_1_data <- DPI_5D_1_data2[over_gene, ]
dim(DPI_2D_1_data)
dim(DPI_5D_1_data)

setwd("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\data\\Stereo-seq")

# save(DPI_2D_1_data, file="DPI_2D_1_data.RData")
# save(DPI_2D_1_info2, file="DPI_2D_1_info.RData")

# save(DPI_5D_1_data, file="DPI_5D_1_data.RData")
# save(DPI_5D_1_info2, file="DPI_5D_1_info.RData")


load("DPI_2D_1_data.RData")
load("DPI_5D_1_data.RData")
load("DPI_2D_1_info.RData")
load("DPI_5D_1_info.RData")




### Exchange two parameters including tao_1 and theta ######
#####################################
library(SPADE)
DPI_2D <- SPADE_norm(readcounts=DPI_2D_1_data, info=DPI_2D_1_info2)
DPI_5D <- SPADE_norm(readcounts=DPI_5D_1_data, info=DPI_5D_1_info2)
res <- SPADE_DE(DPI_2D, DPI_5D, DPI_2D_1_info2, DPI_5D_1_info2, mode="Shape&Strength")




## Makers from the original study, can be downloaded from orginal paper
setwd("E:\\DB\\Dropbox\\Qinfei\\Research\\Spatial scRNA\\references\\science.abp9444_tables_s1_to_s7")
markers <- readxl::read_xlsx("science.abp9444_table_s2.xlsx")
markers <- markers[markers$`Slice Name` %in% c("2DPI-1","5DPI-1"),]
marker.list <- NULL
for (i in 1:nrow(markers)){
  marker.list <- c(marker.list, as.character(markers[i,11:30]))
}
marker.list <- unique(marker.list)


