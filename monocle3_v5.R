library(ggpubr)
library(dendextend)
library(ComplexHeatmap)
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(SeuratWrappers)
library(RColorBrewer)
library(circlize)
library(UCell)
library(clusterProfiler)
library(monocle3)

source(file.path("R","functions.R"))

load(file.path("data", "v5_reference_mapped_all_filtered.RData"))
if (!file.exists( file.path("data","monocle3_object_macrophages.RData"))) {
  cds <- as.cell_data_set(all) %>% 
  cluster_cells()

  plot_cells(cds, color_cells_by = "seurat_clusters")
  
  ## learn the trajectory
  cds <- learn_graph(cds)
  plot_cells(cds,
             color_cells_by = "predicted.cluster",
             label_groups_by_cluster=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE)
  
  ## calculate pseudotime
  cds <- order_cells(cds)
  
  plt <- plot_cells(cds,
             color_cells_by = "pseudotime",
             label_cell_groups=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE,
             graph_label_size=1.5) +
    theme_void()
  
  rasterise(plt)
  
  ggsave(file.path("plots","umap","pseudotime_macrophages.pdf"))
  
  ## choose cells 
  cds_subset <- choose_cells(cds)
  save(cds, file = file.path("data","monocle3_object_macrophages.RData"))
} else {
  load( file.path("data","monocle3_object_macrophages.RData"))
}

## differential gene expression 
## From url: https://github.com/cole-trapnell-lab/monocle-release/issues/295
subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
genes <- row.names(subset(subset_pr_test_res, q_value == 0 & morans_I > 0.25))
genes

pt.matrix <- exprs(cds_subset)[match(genes,rownames(rowData(cds_subset))),order(pseudotime(cds_subset))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds_subset, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes;

#Ward.D2 Hierarchical Clustering
hthc <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

print(hthc)

pdf(file.path("plots","heatmaps","pseudotime_macrophages.pdf"), height = 10, width = 7)
print(hthc)
dev.off()

## genes from heatmap
dend <- row_dend(hthc) 

dend1 <- color_branches(dend, k = 6, groupLabels = TRUE)
plot(dend1)

dend2 <- as.hclust(dend1)
dend3 <- cutree(dend2,6)

module_genes <- data.frame(gene=dend2$labels[dend2$order],
                 cluster=dend3[dend2$order])
lst <- split(module_genes, module_genes$cluster)

names(lst) <- paste0("module",1:6)

## calculate module scores
DefaultAssay(all) <- "RNA"
all <- JoinLayers(all)
module_scores <- ScoreSignatures_UCell(all[["RNA"]]$counts, features=lst)
all <- AddMetaData(all, as.data.frame(module_scores))

## plot pseudotime composition
a <- pseudotime(cds_subset)
ps_time <- data.frame("pseudotime"=a, condition=cds_subset@colData$group) %>% 
  arrange(pseudotime) %>% 
  mutate(p_bin=ceiling(1:length(pseudotime)/500))

ps_time %>% 
  ggplot(aes(x=p_bin, fill=condition)) +
  geom_bar(position = "fill")

## plot modules
walk(colnames(all[[]])[46:51], function(x) {
  plt <- FeaturePlot(all, x) +
    theme_void()
  print(plt)
  
  ggsave(file.path("plots","umap",paste0(x,"_monocle_module_expr.pdf")))
})

DimPlot(all, label=T)

cell_ids <- list("Module_6" = data.frame(cell=rownames(cds_subset@colData)[cds_subset@colData$seurat_clusters %in% c("0","10","18","11")],
                                         group=cds_subset@colData$group[cds_subset@colData$seurat_clusters %in% c("0","10","18","11")]),
            "Module_5" = data.frame(cell=rownames(cds_subset@colData)[cds_subset@colData$seurat_clusters %in% c("0","10")],
                                    group=cds_subset@colData$group[cds_subset@colData$seurat_clusters %in% c("0","10")]), 
            "Module_4" = data.frame(cell=rownames(cds_subset@colData)[cds_subset@colData$seurat_clusters %in% c("8","3","26")],
                                    group=cds_subset@colData$group[cds_subset@colData$seurat_clusters %in% c("8","3","26")]),
            "Module_3" = data.frame(cell=rownames(cds_subset@colData)[cds_subset@colData$seurat_clusters %in% c("2","20")],
                                      group=cds_subset@colData$group[cds_subset@colData$seurat_clusters %in% c("2","20")]),
            "Module_2" = data.frame(cell=rownames(cds_subset@colData)[cds_subset@colData$seurat_clusters %in% c("0","10","8","3","21","26")],
                                    group=cds_subset@colData$group[cds_subset@colData$seurat_clusters %in% c("0","10","8","3","21","26")]),
            "Module_1" = data.frame(cell=rownames(cds_subset@colData)[cds_subset@colData$seurat_clusters %in% c("2","5","11","18","7")],
                                    group=cds_subset@colData$group[cds_subset@colData$seurat_clusters %in% c("2","5","11","18","7")]))

## bar plot relative abundance
cell_ids %>% 
  bind_rows(.id="module") %>% 
  group_by(module, group) %>% 
  summarise(count=n()) %>% 
  mutate(ratio=count/sum(count),
         #module=factor(module, levels=names(cell_ids)),
         group=factor(group, levels=c("PBS","T","NK","BMDM"))) %>% 
  ggplot(aes(x=module, y=ratio,fill=group)) +
  geom_col(position = "fill") +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  theme_pubclean() +
  theme(text = element_text(size=20)) +
  scale_y_continuous(breaks = c(0,.5,1)) +
  scale_y_reverse()

ggsave(file.path("plots","others","monocle_ratios_per_module.pdf"))

test_hyper <- cell_ids %>% 
  bind_rows(.id="module") %>% 
  hyper_test_n(var1 = "module", var2 = "group") %>% 
  write.csv(file.path("data","hyper_test_group_monocle_module.csv"))

## GO term analysis modules
dend <- row_dend(hthc) 
dend2 <- as.hclust(dend)
dend3 <- cutree(dend2,6)

module_genes <- data.frame(gene=dend2$labels[dend2$order],
                           cluster=dend3[dend2$order])

## plot

genes <- bitr(unique(module_genes$gene), fromType = "SYMBOL",
              toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
              OrgDb = 'org.Mm.eg.db')

colnames(genes)[1] <- "gene"

module_genes <- module_genes %>% 
  left_join(genes, relationship = "many-to-many") %>% 
  distinct(cluster, gene, .keep_all = T) %>%
  na.omit()

write.csv(module_genes, file.path("data","monocle_TAMs_module_genes.csv"))
#module_genes$cluster <- factor(module_genes$cluster, levels = c("PBS","T","MF","NK"))

go_terms_bp <- compareCluster(ENTREZID ~ cluster , 
                              data=module_genes, 
                              fun = "enrichGO",
                              OrgDb='org.Mm.eg.db',
                              ont="BP"
)
dotplot(go_terms_bp)  
ggsave(file.path("plots","others","monocle_modules_GO_term_BP_dotplot_v5.pdf"), height = 7, width = 5)

go_terms_mf <- compareCluster(ENTREZID ~ cluster, 
                              data=module_genes, 
                              fun = "enrichGO",
                              OrgDb='org.Mm.eg.db',
                              ont="MF"
)
dotplot(go_terms_mf)                     
ggsave(file.path("plots","others","monocle_modules_GO_term_MF_dotplot_v5.pdf"), height = 7, width = 5)

