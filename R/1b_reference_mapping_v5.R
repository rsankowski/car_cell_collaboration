library(SingleCellExperiment)
library(scDblFinder)
library(tidyverse)
library(data.table)
library(assertthat)
library(SeuratDisk)
library(Seurat)

if (!file.exists(file.path("data", "v5_reference_mapped_all.RData"))) {
  ## Antunes data is downloaded from url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163120
## 10x v3 chemistry data
antunes_v3 <- read.csv(file.path("data","Antunes","GSE163120_Citeseq_Mouse.GBM.KO4_WT4.filtered.RNA.feature.bc.matrix.csv.gz"))
rownames(antunes_v3) <- antunes_v3$X
antunes_v3 <- antunes_v3[, grepl("2$", colnames(antunes_v3))]

antunes_adt <- data.table::fread(file.path("data","Antunes","GSE163120_Citeseq_Mouse.GBM.KO4_WT4.filtered.ADT.feature.bc.matrix.csv.gz")) %>% 
  as.data.frame()
rownames(antunes_adt) <- antunes_adt$V1
antunes_adt <- antunes_adt[, grepl("2$", colnames(antunes_adt))]
colnames(antunes_adt) <- gsub("-",".", colnames(antunes_adt))

#load metadata
metadata_v3 <- fread(file.path("data","Antunes","GSE163120_annot.Citeseq_Mouse.GBM.KO4_WT4.csv.gz")) %>% 
  as.data.frame()
metadata_v3$cell <- gsub("-",".",metadata_v3$cell)
rownames(metadata_v3) <- metadata_v3$cell
metadata_v3 <- metadata_v3[grepl("2$", rownames(metadata_v3)),]

assert_that(identical(rownames(metadata_v3),colnames(antunes_v3)))
assert_that(identical(colnames(antunes_v3), colnames(antunes_adt)))

## load new data
all <- LoadH5Seurat(file.path("data","seurat_object_all_sctransform_norm.H5Seurat"))

if (!file.exists(file.path("data", "Antunes","antunes_v3_multimodal.h5seurat"))) {## Multimodal mapping of the v3 RNA ADT data from Antunes et al
  ## url: https://satijalab.org/seurat/articles/multimodal_vignette.html
  antunes_v3 <- CreateSeuratObject(counts = antunes_v3) %>% 
    AddMetaData(metadata = metadata_v3) %>% 
    AddMetaData(metadata = "Antunes_et_al_10X_v3", col.name = "dataset")
  adt_assay <- CreateAssayObject(counts = antunes_adt)
  
  antunes_v3[["ADT"]] <- adt_assay
  
  ## subset for immune cells only
  antunes_v3 <- subset(antunes_v3, cells=colnames(antunes_v3)[antunes_v3$cluster != "Oligo"])
  
  ## exclude doublets
  sce <- SingleCellExperiment(assays=list(counts=antunes_v3@assays$RNA$counts))
  sce <- scDblFinder(sce)
  antunes_v3 <- antunes_v3[,sce@colData$scDblFinder.class == "singlet"]
  
  ##
  antunes_v3 <- antunes_v3 %>% SCTransform() %>% RunPCA()
  
  ## analyze adt
  DefaultAssay(antunes_v3) <- 'ADT'
  VariableFeatures(antunes_v3) <- rownames(antunes_v3[["ADT"]])
  antunes_v3 <- NormalizeData(antunes_v3, normalization.method = 'CLR', margin = 2) %>% 
    ScaleData() %>% RunPCA(reduction.name = 'apca')
  
  ## find multimodal neighbors
  antunes_v3 <- FindMultiModalNeighbors(
    antunes_v3, reduction.list = list("pca", "apca"), 
    dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
  )
  
  ## run UMAP
  antunes_v3 <- RunUMAP(antunes_v3, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", return.model=TRUE) # the last argument is from url: https://github.com/satijalab/seurat/issues/3615
  antunes_v3 <- FindClusters(antunes_v3, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)
  
  ## plot data
  DimPlot(antunes_v3, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
  DefaultAssay(antunes_v3) <- 'SCT'
  FeaturePlot(antunes_v3, "Hexb")
  
  ## runs SPCA
  antunes_v3 <- antunes_v3 %>% RunSPCA(assay = 'SCT', graph = 'wsnn')
  
  ## save
  save(antunes_v3, file = file.path("data", "Antunes", "antunes_v3_multimodal.RData"))
  #SaveH5Seurat( antunes_v3, file.path("data", "Antunes", "antunes_v3_multimodal.h5seurat"), overwrite = T)
} else {
  load(file.path("data", "Antunes","antunes_v3_multimodal.RData"))
}

  ## run reference mapping
  ## url: https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
  ## Computing a cached neighbor index
  antunes_v3 <- FindNeighbors(
    object = antunes_v3,
    reduction = "spca",
    graph = 'wsnn',
    dims = 1:50,
    graph.name = "spca.annoy.neighbors", 
    k.param = 50,
    cache.index = TRUE,
    return.neighbor = TRUE,
    l2.norm = TRUE
  )

  anchors <- FindTransferAnchors(
        reference = antunes_v3,
        query = all,
        normalization.method = "SCT",
        reference.reduction = "spca",
        reference.neighbors = "spca.annoy.neighbors", 
        dims = 1:50)
   
  all <- MapQuery(
      anchorset = anchors, 
      query = all,
      reference = antunes_v3, 
      refdata = list(
        cluster = "cluster", 
        predicted_ADT = "ADT"),
      reference.reduction = "spca",
      reduction.model = "wnn.umap")
  
all <- all %>% 
  RunUMAP(reduction = 'ref.spca', dims = 1:50) %>%  
  FindNeighbors(reduction = 'ref.spca', dims = 1:50) %>% 
  FindClusters(reduction = 'ref.spca', dims = 1:50) 

DimPlot(all, group.by="predicted.cluster", label=T)

## calculate cell cycle scores
SaveH5Seurat( all, file.path("data", "v5_reference_mapped_all.H5seurat")) #, overwrite = T
save( all, file = file.path("data", "v5_reference_mapped_all.RData")) 
} else {
  load(file.path("data", "v5_reference_mapped_all.RData"))
}

## filter for mt content and rerun UMAP
all <- all %>% 
  AddMetaData(PercentageFeatureSet(., pattern = "^mt-"),"percent.mt")

all <- subset(all, subset = percent.mt < 5 & nFeature_RNA < 5000)

all <- all %>% 
  RunUMAP(reduction = 'ref.spca', dims = 1:50) %>%  
  FindNeighbors(reduction = 'ref.spca', dims = 1:50) %>% 
  FindClusters(reduction = 'ref.spca', dims = 1:50) 

save( all, file = file.path("data", "v5_reference_mapped_all_filtered.RData"))
