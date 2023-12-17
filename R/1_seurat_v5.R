library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(scCustomize)
library(SingleCellExperiment)
library(scDblFinder)

pth <- file.path("data","Gl261_Immune_CAR")
fls <- list.files(pth)

lst <- map(fls, function(x) {
  map(list.files(file.path(pth,x,"outs","per_sample_outs")), function(i) {
    tryCatch({
    .df <- Read10X_h5(file.path(pth,x,"outs","per_sample_outs",i,"count","sample_filtered_feature_bc_matrix.h5"))[["Gene Expression"]] %>% 
      CreateSeuratObject() %>% 
      AddMetaData(ifelse(grepl("CD45_2",i),"CD45_2","panCD45"), "sort_gate") %>% 
      AddMetaData(x, "group") %>% 
      AddMetaData(gsub("_.*","",i), "sample") %>% 
      AddMetaData(PercentageFeatureSet(., pattern = "^mt-"),"percent.mt")
    #.df[["multiplex"]] <- CreateAssayObject(Read10X_h5(file.path(pth,x,"outs","per_sample_outs",i,"count","sample_filtered_feature_bc_matrix.h5"))[[2]])
    sce <- SingleCellExperiment(assays=list(counts=.df@assays$RNA$counts))
    sce <- scDblFinder(sce)
    .df <- .df[,sce@colData$scDblFinder.class == "singlet"]
    .df
    
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  })
})

names(lst) <- fls

lst <- map(lst, Merge_Seurat_List)
all <- Merge_Seurat_List(lst, add.cell.ids=fls)

## normalize data
all <- all %>%
  subset(all, subset = all$percent.mt < 10) %>% 
  NormalizeData(variable.features.n = 10000) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters()

## Plot the CAR expression
FeaturePlot(all, "CAR")
FeaturePlot(all, "RQR8")

## sctransform
all <- all %>% 
  SCTransform(variable.features.n = 10000) %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters()

## plot car cells
FeaturePlot(all[,all[["RNA"]]@counts["CAR",]>0], "CAR")
FeaturePlot(all, "RQR8")

#SaveH5Seurat(all, filename = file.path("data","seurat_object_all.H5Seurat"))
SaveH5Seurat(all, filename = file.path("data","seurat_v5_object_all_sctransform_norm.H5Seurat"))
