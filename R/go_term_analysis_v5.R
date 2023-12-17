library(clusterProfiler)
library(tidyverse)
library(assertthat)
library(SeuratDisk)
library(Seurat)
library(tidyquant)
library(extrafont)

load(file.path("data", "v5_reference_mapped_all_filtered.RData"))
all$group[all$group == "MF"] <- "BMDM"
all$group <- factor(all$group, levels = c("PBS","T","NK","BMDM"))

##subset data for macrophages
macs <- subset(all, idents = as.character(c(0,10,6,13,8,3,21,26, ## TAMs
                                            2,7,18,11,5, ## Tr. Mono
                                            20,19
)))

## go analysis
Idents(macs) <- macs$group
if (!file.exists(file.path("data","cluster_markers_macs_v5.RData"))) {
  macs_markers <- FindAllMarkers(macs, logfc.threshold = .2,recorrect_umi=FALSE) ## https://github.com/satijalab/seurat/issues/6427
  save(macs_markers, file = file.path(file.path("data","cluster_markers_macs_v5.RData")))
  write.csv(macs_markers, file.path("data","cluster_markers_macs_v5.csv"))
} else {
  load(file.path("data","cluster_markers_macs_v5.RData"))
}

## subset data
macs_markers <- macs_markers %>% 
  filter(p_val_adj < .05 & avg_log2FC > 0)

## plot

genes <- bitr(unique(macs_markers$gene), fromType = "SYMBOL",
              toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
              OrgDb = 'org.Mm.eg.db')

colnames(genes)[1] <- "gene"

macs_markers <- macs_markers %>% 
  left_join(genes, relationship = "many-to-many") %>% 
  distinct(cluster, gene, .keep_all = T) %>%
  group_by(cluster) %>% 
  top_n(100, wt=avg_log2FC) %>% 
  na.omit()

macs_markers$cluster <- factor(macs_markers$cluster, levels = c("PBS","T","NK","BMDM"))

go_terms_bp <- compareCluster(ENTREZID ~ cluster , 
                              data=macs_markers, 
                              fun = "enrichGO",
                              OrgDb='org.Mm.eg.db',
                              ont="BP"
)
dotplot(go_terms_bp,
        showCategory = 6,
        size = NULL,
        split = NULL,
        font.size = 12,
        title = "",
        strip_width = 6) +
  theme(text = element_text(family = "Arial",size = 12),
        axis.text.y = element_text(family = "Arial",size = 12),
        axis.text.x = element_text(family = "Arial",size = 12,angle = 45, vjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
                          legend.title = element_text(face = "bold"))

ggsave(file.path("plots","others","macs_GO_term_BP_dotplot_v5.pdf"), device = cairo_pdf, height = 162, width = 125, units = "mm")

## t cells
tcells <- subset(all, idents = as.character(c(27,31,12))) 

## go analysis
Idents(tcells) <- tcells$group
if (!file.exists(file.path("data","cluster_markers_tcells_v5.RData"))) {
  tcells_markers <- FindAllMarkers(tcells, logfc.threshold = .2,recorrect_umi=FALSE) ## https://github.com/satijalab/seurat/issues/6427
  save(tcells_markers, file = file.path(file.path("data","cluster_markers_tcells_v5.RData")))
  write.csv(tcells_markers, file.path("data","cluster_markers_tcells_v5.csv"))
} else {
  load(file.path("data","cluster_markers_tcells_v5.RData"))
}

## subset markers
tcells_markers <- tcells_markers %>% 
  filter(p_val_adj < .05 & avg_log2FC > 0)

## plot
genes <- bitr(unique(tcells_markers$gene), fromType = "SYMBOL",
              toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
              OrgDb = 'org.Mm.eg.db')

colnames(genes)[1] <- "gene"

tcells_markers <- tcells_markers %>% 
  left_join(genes, relationship = "many-to-many") %>% 
  distinct(cluster, gene, .keep_all = T) %>%
  group_by(cluster) %>% 
  top_n(100, wt=avg_log2FC) %>% 
  na.omit()

tcells_markers$cluster <- factor(tcells_markers$cluster, levels = c("PBS","T","NK","BMDM"))

go_terms_bp <- compareCluster(ENTREZID ~ cluster , 
                              data=tcells_markers, 
                              fun = "enrichGO",
                              OrgDb='org.Mm.eg.db',
                              ont="BP"
)
dotplot(go_terms_bp,
        showCategory = 6,
        size = NULL,
        split = NULL,
        font.size = 12,
        title = "",
        strip_width = 6) +
  theme(text = element_text(family = "Arial",size = 12),
        axis.text.y = element_text(family = "Arial",size = 12),
        axis.text.x = element_text(family = "Arial",size = 12,angle = 45, vjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(face = "bold"))

ggsave(file.path("plots","others","tcells_GO_term_BP_dotplot_v5.pdf"), device = cairo_pdf, height = 162, width = 125, units = "mm")

## nk 
nk <- subset(all, idents = "9")

## go analysis
Idents(nk) <- nk$group
if (!file.exists(file.path("data","cluster_markers_nk_v5.RData"))) {
  nk_markers <- FindAllMarkers(nk, logfc.threshold = .2,recorrect_umi=FALSE) ## https://github.com/satijalab/seurat/issues/6427
  save(nk_markers, file = file.path(file.path("data","cluster_markers_nk_v5.RData")))
  write.csv(nk_markers, file.path("data","cluster_markers_nk_v5.csv"))
} else {
  load(file.path("data","cluster_markers_nk_v5.RData"))
}

## subset markers
nk_markers <- nk_markers %>% 
  filter(p_val_adj < .05 & avg_log2FC > 0)

## plot
genes <- bitr(unique(nk_markers$gene), fromType = "SYMBOL",
              toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
              OrgDb = 'org.Mm.eg.db')

colnames(genes)[1] <- "gene"

nk_markers <- nk_markers %>% 
  left_join(genes, relationship = "many-to-many") %>% 
  distinct(cluster, gene, .keep_all = T) %>%
  group_by(cluster) %>% 
  top_n(100, wt=avg_log2FC) %>% 
  na.omit()

nk_markers$cluster <- factor(nk_markers$cluster, levels = c("PBS","T","NK","BMDM"))

go_terms_bp <- compareCluster(ENTREZID ~ cluster , 
                              data=nk_markers, 
                              fun = "enrichGO",
                              OrgDb='org.Mm.eg.db',
                              ont="BP"
)
dotplot(go_terms_bp,
        showCategory = 6,
        size = NULL,
        split = NULL,
        font.size = 12,
        title = "",
        strip_width = 6) +
  theme(axis.text = element_text(family = "Arial",size = 12),
        text = element_text(family = "Arial",size = 12),
        axis.text.y = element_text(family = "Arial",size = 12),
        axis.text.x = element_text(family = "Arial",size = 12,angle = 45, vjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(face = "bold"))

ggsave(file.path("plots","others","nk_GO_term_BP_dotplot_v5.pdf"), device = cairo_pdf, height = 150, width = 120, units = "mm")

## other cell types
## granulocytes
gran <- subset(all, idents = c("2","9","17","14"))
#gran <- JoinLayers(gran)

## go analysis
Idents(gran) <- gran$group
if (!file.exists(file.path("data","cluster_markers_gran_v5.RData"))) {
  gran_markers <- FindAllMarkers(gran, logfc.threshold = .2,recorrect_umi=FALSE) ## https://github.com/satijalab/seurat/issues/6427
  save(gran_markers, file = file.path(file.path("data","cluster_markers_gran_v5.RData")))
  write.csv(gran_markers, file.path("data","cluster_markers_gran_v5.csv"))
} else {
  load(file.path("data","cluster_markers_gran_v5.RData"))
}

## subset markers
gran_markers <- gran_markers %>% 
  filter(p_val_adj < .05 & avg_log2FC > 0)

## plot
genes <- bitr(unique(gran_markers$gene), fromType = "SYMBOL",
              toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
              OrgDb = 'org.Mm.eg.db')

colnames(genes)[1] <- "gene"

gran_markers <- gran_markers %>% 
  left_join(genes, relationship = "many-to-many") %>% 
  distinct(cluster, gene, .keep_all = T) %>%
  group_by(cluster) %>% 
  top_n(100, wt=avg_log2FC) %>% 
  na.omit()

gran_markers$cluster <- factor(gran_markers$cluster, levels = c("PBS","T","MF","NK"))

go_terms_bp <- compareCluster(ENTREZID ~ cluster , 
                              data=gran_markers, 
                              fun = "enrichGO",
                              OrgDb='org.Mm.eg.db',
                              ont="BP"
)
dotplot(go_terms_bp)  
ggsave(file.path("plots","others","gran_GO_term_BP_dotplot_v5.pdf"), height = 7, width = 5)

go_terms_mf <- compareCluster(ENTREZID ~ cluster, 
                              data=gran_markers, 
                              fun = "enrichGO",
                              OrgDb='org.Mm.eg.db',
                              ont="MF"
)
dotplot(go_terms_mf)                         
ggsave(file.path("plots","others","gran_GO_term_MF_dotplot_v5.pdf"), height = 7, width = 6)

## DCs
dcs <- subset(all, idents = c("23","15","25","31","19","16"))
#dcs <- JoinLayers(dcs)

## go analysis
Idents(dcs) <- dcs$group
if (!file.exists(file.path("data","cluster_markers_dcs_v5.RData"))) {
  dcs_markers <- FindAllMarkers(dcs, logfc.threshold = .2,recorrect_umi=FALSE) ## https://github.com/satijalab/seurat/issues/6427
  save(dcs_markers, file = file.path(file.path("data","cluster_markers_dcs_v5.RData")))
  write.csv(dcs_markers, file.path("data","cluster_markers_dcs_v5.csv"))
} else {
  load(file.path("data","cluster_markers_dcs_v5.RData"))
}

## subset markers
dcs_markers <- dcs_markers %>% 
  filter(p_val_adj < .05 & avg_log2FC > 0)

## plot
genes <- bitr(unique(dcs_markers$gene), fromType = "SYMBOL",
              toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
              OrgDb = 'org.Mm.eg.db')

colnames(genes)[1] <- "gene"

dcs_markers <- dcs_markers %>% 
  left_join(genes, relationship = "many-to-many") %>% 
  distinct(cluster, gene, .keep_all = T) %>%
  group_by(cluster) %>% 
  top_n(100, wt=avg_log2FC) %>% 
  na.omit()

dcs_markers$cluster <- factor(dcs_markers$cluster, levels = c("PBS","T","MF","NK"))

go_terms_bp <- compareCluster(ENTREZID ~ cluster , 
                              data=dcs_markers, 
                              fun = "enrichGO",
                              OrgDb='org.Mm.eg.db',
                              ont="BP"
)
dotplot(go_terms_bp)  
ggsave(file.path("plots","others","dcs_GO_term_BP_dotplot_v5.pdf"), height = 7, width = 5)

go_terms_mf <- compareCluster(ENTREZID ~ cluster, 
                              data=dcs_markers, 
                              fun = "enrichGO",
                              OrgDb='org.Mm.eg.db',
                              ont="MF"
)
dotplot(go_terms_mf)                         
ggsave(file.path("plots","others","dcs_GO_term_MF_dotplot_v5.pdf"), height = 7, width = 5)

## all cells
Idents(all) <- all$group
if (!file.exists(file.path("data","cluster_markers_all_cells_by_group_v5.RData"))) {
  all_markers_by_gr <- FindAllMarkers(all, logfc.threshold = .2,recorrect_umi=FALSE) ## https://github.com/satijalab/seurat/issues/6427
  save(all_markers_by_gr, file = file.path(file.path("data","cluster_markers_all_cells_by_group_v5.RData")))
  write.csv(all_markers_by_gr, file.path("data","cluster_markers_all_cells_by_group_v5.csv"))
} else {
  load(file.path("data","cluster_markers_all_cells_by_group_v5.RData"))
}

## subset markers
all_markers_by_gr <- all_markers_by_gr %>% 
  filter(p_val_adj < .05 & avg_log2FC > 0)

## plot
genes <- bitr(unique(all_markers_by_gr$gene), fromType = "SYMBOL",
              toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
              OrgDb = 'org.Mm.eg.db')

colnames(genes)[1] <- "gene"

all_markers_by_gr <- all_markers_by_gr %>% 
  left_join(genes, relationship = "many-to-many") %>% 
  distinct(cluster, gene, .keep_all = T) %>%
  na.omit()

all_markers_by_gr$cluster <- factor(all_markers_by_gr$cluster, levels = c("PBS","T","MF","NK"))

go_terms_bp <- compareCluster(ENTREZID ~ cluster , 
                              data=all_markers_by_gr, 
                              fun = "enrichGO",
                              OrgDb='org.Mm.eg.db',
                              ont="BP"
)
dotplot(go_terms_bp)  
ggsave(file.path("plots","others","all_cells_by_group_GO_term_BP_dotplot_v5.pdf"), height = 7, width = 5)

go_terms_mf <- compareCluster(ENTREZID ~ cluster, 
                              data=all_markers_by_gr, 
                              fun = "enrichGO",
                              OrgDb='org.Mm.eg.db',
                              ont="MF"
)
dotplot(go_terms_mf)                         
ggsave(file.path("plots","others","all_cells_by_group_GO_term_MF_dotplot_v5.pdf"), height = 7, width = 6)
