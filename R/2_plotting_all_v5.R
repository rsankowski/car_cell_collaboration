library(readxl)
library(UCell)
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(tidyquant)
library(ggpointdensity)
library(viridis)
library(ggrastr)
library(extrafont)
library(scCustomize)

source(file.path("R","functions.R"))

#all <- LoadH5Seurat(file.path("data","v5_reference_mapped_all.H5Seurat"))
load(file.path("data", "v5_reference_mapped_all_filtered.RData"))

if (F) {
  ## reorder clusters
  order_clusters <- data.frame(seurat_clusters= all$seurat_clusters, row.names = rownames(all[[]])) %>%
    bind_cols(as.data.frame(t(all[["SCT"]]$scale.data))) %>%
    group_by(seurat_clusters) %>%
    summarize_all(.funs=mean) %>%
    as.data.frame()
  
  rownames(order_clusters) <- order_clusters$seurat_clusters
  order_clusters <- order_clusters$seurat_clusters[hclust(dist(order_clusters[,-1]))$order]
  
  #reorder the clusters
  levels(all) <- rev(order_clusters)
  all$seurat_clusters <- factor(all$seurat_clusters, levels = levels(all))
}

## order groups
all$group[all$group == "MF"] <- "BMDM"
all$group <- factor(all$group, levels = c("PBS","T","NK","BMDM"))

## QC Metrics
VlnPlot(all, c("nCount_RNA","nFeature_RNA","percent.mt"),stack = T, flip = T, split.by = "group")

order_clusters3 <- as.character(c(0,10,6,13,8,3,21,26, ## TAMs
                                  2,7,18,11,5, ## Tr. Mono
                                  20,19, ## Mono
                                  16,22,23,24,30,17, ## DCs
                                  28, ## MG
                                  27,31,12, ## T cells 
                                  9, # NK cells
                                  29, ## B cells
                                  25,4,1,15, ## granulocytes
                                  32, ## Mast cells
                                  14 ## low quality cells
)) #rev(order_clusters[!order_clusters %in% c("6","29","30")])
levels(all) <- order_clusters3

walk(unique(all$group), function(x) {
  p <- VlnPlot(all[,all$group == x], 
                       c("nCount_RNA","nFeature_RNA","percent.mt"),
               stack = T,
               flip = T,
               cols = unname(palette_light())) +
    theme(text = element_text(family = "Arial",size = 18),
          axis.text.x = element_text(family = "Arial",size = 18, angle = 45),
          axis.text.y = element_text(family = "Arial",size = 18),
          line = element_line(linewidth = .1),
          rect = element_rect(linewidth = .1)) +
    NoLegend()
  print(p)
  
  ggsave(file.path("plots","others",paste0(x,"_qc_violin_plt.pdf")), height = 192, width = 180, units = "mm", device = cairo_pdf)
})

## plot clusters
DimPlot(all, label = T) +
  scale_color_manual(values = c(colors_many,colors_fig)[-c(14,21)]) +
  theme_void()
ggsave(file.path("plots","umap","clusters_all_v5.pdf"))

## no legend
DimPlot(all, label = T) +
  #scale_color_tq() +
  scale_color_manual(values = c(colors_many,colors_fig)[-c(14,21)]) +
  theme_void() +
  NoLegend()

ggsave(file.path("plots","umap","clusters_all_no_legend_v5.pdf"))

## no legend
DimPlot(all, label = T, raster = T) +
  #scale_color_tq() +
  scale_color_manual(values = c(colors_many,colors_fig)[-c(14,21)]) +
  theme_void() +
  NoLegend()

ggsave(file.path("plots","umap","clusters_all_no_legend_raster_v5.pdf"))

## plot sort_gate
DimPlot(all, label = T, group.by = "sort_gate") +
  scale_color_brewer(palette = "Set1") +
  theme_void()
ggsave(file.path("plots","umap","sort_gate_all_v5.pdf"))

## plot group
DimPlot(all, group.by = "group") +
  scale_color_brewer(palette = "Set2") +
  theme_void() +
  facet_wrap(~ group, ncol=2)
ggsave(file.path("plots","umap","group_all_v5.pdf"))

## plot predicted.cluster
pr_cl <- DimPlot(all, label = T, group.by = "predicted.cluster", pt.size = .1) +
  #scale_color_brewer(palette = "Set2") +
  scale_color_manual(values = colors_fig) +
  theme_void() +
  NoLegend()

rasterise(pr_cl)
ggsave(file.path("plots","umap","predicted.cluster_all_v5.pdf"))

## hypergeometric testing
hyper_tests <- all[[]] %>% 
  hyper_test_n(var1 = "seurat_clusters",var2 = "group") 

hyper_tests %>% 
  mutate('-log adj.\nP value' = -log(padj),
         Cluster=factor(Cluster, levels=order_clusters3)) %>% 
  filter(Significance!="n.s.",
         Cluster != "14") %>% 
  ggplot(aes(factor(enrichment_var, levels=c("PBS","T","NK","BMDM")), Cluster, size=`-log adj.\nP value`)) +
  geom_point() +
  theme_linedraw() +
  labs(y="Cluster", x="Treatment") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7),
        axis.text = element_text(family = "Arial",size = 6),
        text = element_text(family = "Arial",size = 6),
        axis.title = element_text(size = 6, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "right", 
        legend.justification = "bottom") +
  scale_size_continuous(range = c(.1,3)) ## https://stackoverflow.com/questions/20155400/defining-minimum-point-size-in-ggplot2-geom-point

ggsave(file.path("plots","others","hyper_test_groups_adj_p_values_v5.pdf"), device = cairo_pdf, height = 61, width = 38, units = "mm")

## 2d density plots of the conditions
walk(unique(all$group), function(x) {
  plt <- data.frame(all@reductions$umap@cell.embeddings, "group"=factor(all$group, levels=c("PBS","T","NK","BMDM"))) %>%
    filter(group == x)  %>%
    ggplot(aes(umap_1,umap_2, z=(group))) +
    #stat_summary_hex(bins=200, fun=mean) +
    #facet_wrap(~group) +
    geom_pointdensity(adjust = .5) +
    theme_void() +
    scale_color_viridis(option = "A", direction = -1) 

  plt2 <- rasterise(plt)
  print(plt2)
  ggsave(file.path("plots","umap",paste0("v5_cell_density_",x,".pdf")))

})
  sigs <- read_excel(file.path("data","cell_signatures.xlsx"), skip = 1)
sigs <- as.list(sigs)
sigs <- map(sigs, na.omit)
  ## plot cell signatures
signatures <- list(
  "DAM1"= c("Apoe", "Ctsb","Lyz2","Fth1","Ctsd","Spp1","Clec7a","Itgax","Lpl","Cts7","Timp2","H2-D1","B2m","Ctsb","Tyrobp","Ctsl","Cd63","Csf1","Cd9","Ctsa","Ank", "Ctsz","Axl","Serpine2","Ccl6","Cd68","Cadm1","Cd52","Gusb","Hif1a"),
  "DAM2"= c("Apoe", "Ctsb","Lyz2","Fth1","Ctsd","Spp1","Clec7a","Itgax","Lpl","Cts7","Lgals3","Gpnmb","Igf1","Mif","Cd33","Abca1","Srgap2","Timp2","H2-D1","B2m","Ctsb","Tyrobp","Ctsl","Cd63","Csf1","Cd9","Ctsa","Ank", "Ctsz","Axl","Serpine2","Ccl6","Cd68","Cadm1","Cd52","Gusb","Hif1a","Lilrb4a","Trem2","Serinc3"))

signatures <- c(sigs, signatures)

## add signatures
sig_scores <- ScoreSignatures_UCell(all[["SCT"]]$counts, features=signatures)
all <- AddMetaData(all, as.data.frame(sig_scores))
walk(colnames(all[[]])[16:45], function(x) {
  plt <- FeaturePlot(all, x, raster = T, raster.dpi = c(300,400)) +
    theme_void()
  print(plt)
  ggsave(file.path("plots","umap",paste0(x,"_signature.pdf")))
})


