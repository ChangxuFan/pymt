source("objects.R")
wd <- "publication/sth/pub2.1_DC_markers/"
dir.create(wd)

soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")
so1 <- soi[, soi$sample == SAMPLES[1]]

clusters <- c("6", "19")

so1 <- so1[, so1$seurat_clusters %in% clusters]

color.map <- c("#71B000", "#F265E7")
names(color.map) <- clusters

p <- umap.fanc(so1, group.by = "seurat_clusters", 
               show.legends = F, axis.type = "nothing", 
               pt.size = 0.1, label.size = 2.5,
               label.groups = F,
               cols = color.map,
               remove.outlier = c(0, 0, 0, 2), pt.shape = 19,
               plot.out = paste0(wd, "/DC_UMAP.pdf"), width = 0.5, height = 0.5)

dc.markers <- c("Csf1r", "Flt3", "Irf8", "Zeb2", "Id2",
                "H2-Aa", "Batf3", "Sirpa","Siglech")

outliers <- lapply(1:2, function(i) {
  names(which.min(so1@reductions$umap@cell.embeddings[, i]))
}) %>% unlist()

plot.panel.list(panel.list = dc.markers, obj = so1[, !colnames(so1) %in% outliers],
                order = F, assay = "RNA", 
                plot.out = paste0(wd, "/DC_markers.pdf"),
                raster = T, auto.adjust.raster.pt.size = F,
                label = F, label.size = 1.8,
                publication = T, invisible = T, 
                n.col = 5, pdf.by.row = T,
                sub.width = 0.85, sub.height = 0.85,
                pt.size = 5, 
                return.list = F, threads = 1)
