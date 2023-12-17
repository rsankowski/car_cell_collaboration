library(RaceID)
library(SeuratObject)
library(tidyverse)
library(assertthat)
library(FateID)
library(readxl)
library(RColorBrewer)
library(ggpubr)
library(tidyquant)
library(ggpubr)
library(Seurat)
library(SeuratDisk)

load(file.path("data", "v5_reference_mapped_all_filtered.RData"))
source(file.path("R","functions.R"))

DefaultAssay(all) <- "RNA"
all <- SeuratObject::JoinLayers(all)

#extract metadata
metadata <- all@meta.data
metadata$seurat_clusters <- factor(metadata$seurat_clusters, levels = levels(all))

date <- Sys.Date()

set.seed(79106)

if (!file.exists(file.path("data","sc_all_v5.RData"))) {
  
  sc <- SCseq(all[["RNA"]]$counts)
  
  #filter data
  sc <- filterdata(sc, 
                   mintotal=100,
                   minnumber = 1,
                   knn=10,
                   minexpr = 1)
  
  # 2.Run Seurat with filtering such that the same cells are retained
  assert_that(length(colnames(all)) == length(colnames(sc@ndata)))
  
  # 3.Re-initialize RaceID output with Seurat data:
  
  part <- as.numeric(as.character(all@meta.data$seurat_clusters))
  d <- as.matrix( dist(all@reductions$pca@cell.embeddings) )
  tsne <- as.data.frame(all@reductions$umap@cell.embeddings)
  names(part) <- colnames(sc@ndata)
  
  n <- colnames(sc@ndata)
  part <- part[n]
  
  # partition
  sc@cpart <- sc@cluster$kpart <- part
  # distances
  sc@distances <- d[n,n]
  # tsne
  sc@tsne <- tsne[n,]
  rm(d, tsne)
  
  sc@medoids <- compmedoids(sc, sc@cpart)
  
  #reorder the clusters
  #idx <- which(order_clusters %in% unique(as.character(sc@cpart)))
  #sc@cpart <- factor(sc@cpart, levels = order_clusters)
  
  sc@fcol <- rep_len(unname(palette_light()), length.out = length(unique(sc@cpart)))#c(colors_pat, colors_many, colors_fig)[-2][idx]
  
  save(sc, file=file.path("data","sc_all_v5.RData"))
  
} else {
  load(file.path("data","sc_all_v5.RData"))
}

sc@fcol <- colors_many

#StemID
if (!file.exists(file.path("data","ltr_all.RData"))){ #data/ltr-larger-clusters.RData
  ltr <- Ltree(sc)
  
  #convert clusters in integers
  ltr@sc@cpart <- as.numeric(as.character(ltr@sc@cpart)) +1
  names(ltr@sc@cpart) <- colnames(ltr@sc@ndata)
  
  ltr <- compentropy(ltr)
  ltr <- projcells(ltr,nmode=TRUE,fr=FALSE) #400
  ltr <- projback(ltr,pdishuf=100)
  ltr <- lineagegraph(ltr)
  ltr <- comppvalue(ltr,pthr=0.2)
  
  save(ltr, file = file.path("data","ltr_all.RData"))
} else {
  load(file.path("data","ltr_all.RData"))
}

#plot entropy umap
.ent <- data.frame(ltr@sc@tsne, Entropy=ltr@entropy, Cluster=ltr@sc@cpart-1)

plot_continuous(param = "Entropy")
ggsave(file.path("plots","umap","10x_all","entropy.pdf"), width = 11.2, height = 8.11, useDingbats = F)

## define colors
my_colors <- c(colors_many, colors,unname(palette_green()))[-c(14,21,28,37)]
names(my_colors) <- levels(all)

.ent %>% 
  filter(Cluster %in% c(20,7)) %>% 
  mutate(Cluster=factor(Cluster,levels=c("7","20")),
         delta_entr = Entropy -min(Entropy)) %>% 
  #group_by(Cluster)  %>% 
  #summarise(mean_entr=mean(Entropy)) %>%
  ggplot(aes(x=Cluster, y=delta_entr, fill=Cluster,color=Cluster)) +
  geom_beeswarm(size=2,pch=19,alpha=.5)+
  stat_summary(fun = mean, geom="crossbar") +
  scale_color_manual(values = my_colors[c("7","20")]) +
  scale_fill_manual(values = my_colors[c("7","20")]) +
  scale_y_continuous(breaks = c(0,.1,.2)) +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 20)) +
  stat_compare_means() +
  labs(y="delta Entropy") +
  expand_limits(y=.23)

ggsave(file.path(file.path("plots","others","monocyte_engraftment_entropy.pdf")), height = 3, width = 6)




x <- compscore(ltr,scthr=0.2)
plotgraph(ltr,showCells=FALSE)

pdf(file.path("plots","umap","lineage-graph-control-all-cells.pdf"), width = 11.2, height = 8.11, useDingbats = F)
plotgraph(ltr,showCells=FALSE)
dev.off()

#x <- compscore(ltr,scthr=0.9)
plotdistanceratio(ltr)
plotspantree(ltr)
plotprojections(ltr)

#lineage tree for moDCs
monocytes <- c(10,8,2,1)
mg <- c(17,5,4)

#pseudotemporal monocytes
n <- cellsfromtree(ltr,monocytes)
x <- getfdata(ltr@sc)

fs  <- filterset(x,n=n$f, minexpr = 1, minnumber = 10)

if (!file.exists(file.path("data","s1d-ctrl-monocytes.Robj"))) {
  s1d <- getsom(fs,nb=1000,alpha=.5)
  save(s1d, file = file.path("data","s1d-ctrl-monocytes.Robj"))
} else {
  load(file.path("data","s1d-ctrl-monocytes.Robj"))
}

ps  <- procsom(s1d,corthr=.8,minsom=5)
y    <- ltr@sc@cpart[n$f]
fcol <- sc@fcol
plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)

pdf(file.path("plots","heatmaps","monocytes-trajectory-heatmap.pdf"), width = 8.57, height = 5.79, useDingbats = F)
plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
dev.off()

pdf(file.path("plots","heatmaps","monocytes-outline-all-trajectory-heatmap.pdf"), width = 8.57, height = 5.79, useDingbats = F)
plotheatmap2(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
dev.off()

png(file.path("plots","heatmaps","map-monocytes-trajectory-heatmap.png"))
image(t(as.matrix(ps$all.z)), col = rev(colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100)), axes = FALSE, ylim = c(-0.02,1))
dev.off()


#export node genes
modules <- data.frame('Node' = NA, 'Genes' = NA)
for (i in 1: max(ps$nodes)) {
  gene_names <- names(ps$nodes)[ps$nodes == i]
  gene_names <- gsub('_.*', '', gene_names)
  modules2 <- data.frame('Node' = NA, 'Genes' = gene_names)
  modules2$Node <- rep(as.character(i), nrow(modules2))
  modules <- rbind(na.omit(modules), modules2)
}

write_csv(modules, file.path("data","nodes-stemid-monocytes-vector-monocytes.csv"))

## plot samples along trajectory
traj <- all[[]][names(ps$nodes.e),] %>% 
  mutate(trajectory=1:nrow(.),
         bin=ceiling(trajectory/50)) %>% 
  group_by(bin,group, sample) %>% 
  summarise(count=n()) %>% 
  group_by(sample) %>% 
  mutate(sample_size=sum(count)) %>% 
  mutate(norm_count=count/sample_size)

## sanity check that all samples add up to 1
traj %>% group_by(sample) %>% summarise(sum(norm_count))

traj %>% 
  ggplot(aes(x=bin, y=norm_count, color=group)) +
    #geom_point() +
    geom_smooth(size=2) +
    theme_bw() +
  facet_wrap(~ group) +
  stat_cor(method = "spearman",
           label.y = .0225) +
  labs(title = "spearman correlation")
  
ggsave(file.path("plots","others","monocyte_trajectory_normalized_cell_ratio.pdf")) 

traj %>% 
  ggplot(aes(x=bin, y=norm_count, color=group)) +
  #geom_point() +
  geom_smooth(size=2) +
  theme_bw() +
  facet_wrap(~ group) +
  stat_cor(method = "pearson",
           label.y = .0225) +
  labs(title = "pearson correlation")

ggsave(file.path("plots","others","monocyte_trajectory_normalized_cell_ratio_pearson.pdf")) 

#pseudotemporal mg
n <- cellsfromtree(ltr,mg)
x <- getfdata(ltr@sc)

fs  <- filterset(x,n=n$f, minexpr = 1, minnumber = 10)

if (!file.exists(file.path("data","s1d-ctrl-mg.Robj"))) {
  s1d <- getsom(fs,nb=1000,alpha=.5)
  save(s1d, file = file.path("data","s1d-ctrl-mg.Robj"))
} else {
  load(file.path("data","s1d-ctrl-mg.Robj"))
}

ps  <- procsom(s1d,corthr=.8,minsom=5)
y    <- ltr@sc@cpart[n$f]
fcol <- sc@fcol
plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)

pdf(file.path("plots","heatmaps","mg-trajectory-heatmap.pdf"), width = 8.57, height = 5.79, useDingbats = F)
plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
dev.off()

pdf(file.path("plots","heatmaps","mg-outline-all-trajectory-heatmap.pdf"), width = 8.57, height = 5.79, useDingbats = F)
plotheatmap2(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
dev.off()

png(file.path("plots","heatmaps","map-mg-trajectory-heatmap.png"))
image(t(as.matrix(ps$all.z)), col = rev(colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100)), axes = FALSE, ylim = c(-0.02,1))
dev.off()


#export node genes
modules <- data.frame('Node' = NA, 'Genes' = NA)
for (i in 1: max(ps$nodes)) {
  gene_names <- names(ps$nodes)[ps$nodes == i]
  gene_names <- gsub('_.*', '', gene_names)
  modules2 <- data.frame('Node' = NA, 'Genes' = gene_names)
  modules2$Node <- rep(as.character(i), nrow(modules2))
  modules <- rbind(na.omit(modules), modules2)
}

write_csv(modules, file.path("data","nodes-stemid-mg-vector-mg.csv"))

## plot samples along trajectory
traj <- all[[]][names(ps$nodes.e),] %>% 
  mutate(trajectory=1:nrow(.),
         bin=ceiling(trajectory/50)) %>% 
  group_by(bin,group, sample) %>% 
  summarise(count=n()) %>% 
  group_by(sample) %>% 
  mutate(sample_size=sum(count)) %>% 
  mutate(norm_count=count/sample_size)

## sanity check that all samples add up to 1
traj %>% group_by(sample) %>% summarise(sum(norm_count))

traj %>% 
  ggplot(aes(x=bin, y=norm_count, color=group)) +
  #geom_point() +
  geom_smooth(size=2) +
  theme_bw() +
  facet_wrap(~ group) +
  stat_cor(method = "spearman",
           label.y = .0225) +
  labs(title = "spearman correlation")

ggsave(file.path("plots","others","mg_trajectory_normalized_cell_ratio_spearman.pdf")) 

traj %>% 
  ggplot(aes(x=bin, y=norm_count, color=group)) +
  #geom_point() +
  geom_smooth(size=2) +
  theme_bw() +
  facet_wrap(~ group) +
  stat_cor(method = "pearson",
           label.y = .0225) +
  labs(title = "pearson correlation")

ggsave(file.path("plots","others","mg_trajectory_normalized_cell_ratio_pearson.pdf")) 
