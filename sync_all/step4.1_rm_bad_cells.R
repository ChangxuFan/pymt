source("objects.R")

aoi <- readRDS("sync_all/sth/atac/merge_v1/ao.Rds")
soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")

# get rid of low quality cells from BM_tumor_ba3rep1:

p <- seurat.plot.archr(ao = aoi, ao.embedding = "UMAP", ao.name = "Clusters",
                       shuffle = T)
# bk <- p
aot <- aoi[aoi$Clusters != "C12",]

p <- seurat.plot.archr(ao = aot, ao.embedding = "UMAP", ao.name = "Clusters",
                       shuffle = T)
polygon.df <- data.frame(x = c(0, 3, 3, 0), y = c(-1, 0.5, 2.8, 2.8))

p + geom_polygon(data = polygon.df, mapping = aes(x = x, y = y), 
                 fill = NA, color = "red")
bad.cells <- select.cells.by.polygon(embed.df = aot@embeddings$UMAP$df, polygon.df = polygon.df)

aot <- aot[! aot$cellNames %in% bad.cells, ]
p2 <- seurat.plot.archr(ao = aot, ao.embedding = "UMAP", ao.name = "Clusters",
                        shuffle = T)
p2 <- p2 + theme(aspect.ratio = 1)
ggsave("sync_all/sth/atac/merge_v1/Plots/rmBadCells_Clusters.png", p2, dpi = 100, 
       width = 5, height = 5)

aoi <- aot
soi <- soi[, aoi$cellNames]

p3 <- dim.plot.simple(so = soi, plot.out = "sync_all/sth/rna/int_rpca20/plots/umap_cluster_rmBadCells.png", 
                      group.by = "seurat_clusters", label = T)
system(paste0("mv ", "sync_all/sth/rna/int_rpca20/soi.Rds ", 
              "sync_all/sth/rna/int_rpca20/soi_preRmBadCells_2022-12-16.Rds"))
system(paste0("mv ", "sync_all/sth/atac/merge_v1/ao.Rds ",
              "sync_all/sth/atac/merge_v1/ao_preRmBadCells_2022-12-16.Rds"))

saveRDS(aoi, "sync_all/sth/atac/merge_v1/ao.Rds")
saveRDS(soi, "sync_all/sth/rna/int_rpca20/soi.Rds", compress = F)
