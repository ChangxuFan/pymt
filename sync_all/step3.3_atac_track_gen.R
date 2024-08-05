source("objects.R")
work.dir <- "sync_all/sth/atac/merge_v1/"

aoi <- readRDS(paste0(work.dir, "/ao.Rds"))

lapply(c("seurat_clusters", "Clusters"), function(cluster) {
  getGroupBW.fanc(ArchRProj = aoi, groupBy = cluster, split.by = "Sample",
                  out.dir = paste0(work.dir, "/bw/", cluster, "/"), tileSize = 50)
})
