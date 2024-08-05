rm(list = ls())
source("objects.R")
soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")
aoi <- readRDS("sync_all/sth/atac/merge_v1/ao.Rds")
wd <- paste0("publication/sth/pub3.1a_c20_UMAP")
dir.create(wd)

#>>>>>>>> We first plot UMAPs with c20 highlighted.
pl <- plot.panel.list(soi, 
                      panel.list = "seurat_clusters", 
                      binarize.panel = T, binarize.items = "20",
                      order = F, 
                      assay = "RNA", invisible = T, 
                      plot.out = paste0(wd, "/UMAP/highlight_c20.pdf"),
                      split.by = "tissue", use.split.as.title = T, split.order = c("BM", "SP"),
                      raster = T, label = T, label.size = 1.8,
                      publication = T, 
                      italic.title = F,
                      sub.width = 1.2, sub.height = 1.2)

pl <- plot.panel.list(aoi, ao.embedding = "UMAP",
                      panel.list = "seurat_clusters", 
                      binarize.panel = T, binarize.items = "20",
                      order = F, 
                      assay = "RNA", invisible = T, 
                      plot.out = paste0(wd, "/UMAP/highlight_c20_ATAC.pdf"),
                      split.by = "tissue", use.split.as.title = T, split.order = c("BM", "SP"),
                      raster = T, label = T, label.size = 1.8,
                      publication = T, 
                      italic.title = F,
                      sub.width = 1.2, sub.height = 1.2)

pl <- plot.panel.list(soi, sample = SAMPLES.SP, hide.legend = T,
                      panel.list = "seurat_clusters", 
                      binarize.panel = T, binarize.items = "20",
                      order = F, 
                      assay = "RNA", invisible = T, 
                      plot.out = paste0(wd, "/UMAP/highlight_c20_SP.pdf"),
                      # split.by = "tissue", use.split.as.title = T, split.order = c("BM", "SP"),
                      raster = T, label = T, label.size = 1.8,
                      publication = T, 
                      italic.title = F,
                      sub.width = 1.5, sub.height = 1.5)

# Now we print a zoomed-in version:
soi.zoom <- soi[, soi$seurat_clusters %in% c("1", "9", "16", "20") & 
                  soi$sample %in% SAMPLES.SP]

cells.zoom <- soi.zoom@reductions$umap@cell.embeddings %>% as.data.frame() %>% 
  filter(UMAP_1 > 0, UMAP_2 > -1) %>% rownames()
  
soi.zoom <- soi.zoom[, cells.zoom]

pl <- plot.panel.list(soi.zoom,
                      hide.legend = T,
                      panel.list = "seurat_clusters", 
                      binarize.panel = T, binarize.items = "20",
                      order = F, 
                      assay = "RNA", invisible = T, 
                      plot.out = paste0(wd, "/UMAP/highlight_c20_SP_zoom.pdf"),
                      # split.by = "tissue", use.split.as.title = T, split.order = c("BM", "SP"),
                      raster = T, label = T, label.size = 2.5,
                      publication = T, 
                      italic.title = F,
                      sub.width = 1.5, sub.height = 1.5)

# Plot marker genes: Hlf and Myct1
markers <- c("Hlf", "Myct1", "Ly6a", "Selp")
lapply(c(T, F), function(bOrder) {
  plot.panel.list(soi, sample = SAMPLES.SP,
                  panel.list = markers,
                  order = bOrder, assay = "RNA", invisible = T, 
                  plot.out = paste0(wd, "/UMAP/marker_genes_", bOrder, ".pdf"),
                  raster = T, label = F, publication = T, ymax = NULL, ymin = 0.1, italic.title = T,
                  sub.width = 1.2, sub.height = 1.2, n.col = 4)
  
})

# Plot ATAC-seq UMAP, highlighting clusters 1 or 20
aoi.rm <- aoi[!aoi$seurat_clusters %in% c("22", "19", "13", "21"),]
# aoi.rm$seurat_clusters[aoi.rm$seurat_clusters %in% c("1", "20")] <- "1/20"

pl <- plot.panel.list(aoi.rm, ao.embedding = "UMAP",
                      sample = SAMPLES.SP, hide.legend = T,
                      panel.list = "seurat_clusters", 
                      binarize.panel = T, binarize.items = c("1", "20"),
                      order = T, 
                      assay = "RNA", invisible = T, 
                      plot.out = paste0(wd, "/UMAP/highlight_ATAC_c1_c20_SP.pdf"),
                      # split.by = "tissue", use.split.as.title = T, split.order = c("BM", "SP"),
                      raster = T, label = T, label.size = 1.8,
                      publication = T, 
                      italic.title = F,
                      sub.width = 1.5, sub.height = 1.5)

# Plot BM c20 vs spleen c20, ATAC:
pl <- plot.panel.list(aoi.rm, ao.embedding = "UMAP",
                      hide.legend = T,
                      panel.list = "seurat_clusters", 
                      binarize.panel = T, binarize.items = c("20"),
                      order = T, 
                      assay = "RNA", invisible = T, 
                      plot.out = paste0(wd, "/UMAP/highlight_ATAC.pdf"),
                      split.by = "tissue", use.split.as.title = T, split.order = c("BM", "SP"),
                      raster = T, label = T, label.size = 1.8,
                      publication = T, 
                      italic.title = F,
                      sub.width = 1.5, sub.height = 1.5)

# Plot BM vs spleen c20, RNA:
soi$c20 <- soi$seurat_clusters %>% as.character()
soi$c20[soi$c20 == "20" & soi$tissue == "BM"] <- "20 (BM)"
soi$c20[soi$c20 == "20" & soi$tissue == "SP"] <- "20 (SP)"

n.clusters <- soi$c20 %>% unique() %>% length()
color.map <- rep("gray75", n.clusters)
names(color.map) <- soi$c20 %>% unique() %>% gtools::mixedsort()
color.map["20 (BM)"] <- "orangered"
color.map["20 (SP)"] <- "yellowgreen"

soi$c20 <- factor(soi$c20)

p <- umap.fanc(soi, group.by = "c20", reverse.point.order = T,
               cols = color.map,
               show.legends = F, axis.type = "nothing", 
               pt.size = 0.3, label.size = 2.5,
               label.groups = T,
               plot.out = paste0(wd, "/UMAP/RNA_BM_vs_SP.pdf"),
               height = 1.5, width = 1.5)

soi.zoom2 <- soi[, soi$seurat_clusters %in% c("19", "20", "6")]
cells.zoom2 <- soi.zoom2@reductions$umap@cell.embeddings %>% as.data.frame() %>% 
  filter(UMAP_1 > 6, UMAP_2 > -3) %>% rownames()
soi.zoom2 <- soi.zoom2[, cells.zoom2]

DimPlot(soi.zoom2)

p <- umap.fanc(soi.zoom2, group.by = "c20", reverse.point.order = T,
               cols = color.map,
               show.legends = F, axis.type = "nothing", 
               pt.size = 0.8, label.size = 2.5,
               label.groups = T,
               plot.out = paste0(wd, "/UMAP/RNA_BM_vs_SP_zoom.pdf"),
               height = 1.5, width = 1.5)
