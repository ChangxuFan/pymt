source("objects.R")
wd <- paste0("sync_all/sth/step5.1_intersect_ATAC_RNA_clustering/")
dir.create(wd)

aoi <- readRDS("sync_all/sth/atac/merge_v1/ao.Rds")
soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")
soi <- seurat.add.archr.meta(so = soi, ao = aoi, metas = "Clusters", overwrite = T)


df <- soi@meta.data[, c("seurat_clusters", "Clusters")]
colnames(df) <- c("sc", "c")

mat <- table(df) %>% as.matrix()
write.csv(mat, paste0(wd, "/confusion_mat_raw.csv"))

# write an algorithm to find the correspondance between RNA and ATAC based clustering

soi <- filter.by.cluster.correspondence(so = soi, cluster.ident = "seurat_clusters", filter.by = "Clusters", 
                                        new.cluster.ident = "scc", step.pct = 0.5, max.steps = 3,
                                        out.dir = wd) 

saveRDS(soi, "sync_all/sth/rna/int_rpca20/soi.Rds")
aoi <- archr.add.seurat.m(ao = aoi, so = soi, metas = "scc", as.factor = T)
saveRDS(aoi, "sync_all/sth/atac/merge_v1/ao.Rds")

#### plot a confusion matrix:
mat <- read.csv(paste0(wd, "/confusion_mat_raw.csv"), row.names = 1) %>% as.matrix()
colnames(mat) <- sub("C", "", colnames(mat))

mat <- mat[, as.character(c(1,2,3,5,4,6,8,7,13,10,11,9))]
plot.confusion.mat(mat = mat, plot.out = paste0(wd, "/confusion_mat.pdf"), 
                   column_names_rot = 0, column_names_centered = T)


clusters <- as.character(c(0:7, 9:11, 14, 16, 20))

work.dir <- "sync_all/sth/atac/merge_v1/"
getGroupBW.fanc(ArchRProj = aoi, groupBy = "scc", groups = clusters, split.by = "Sample",
                out.dir = paste0(work.dir, "/bw/scc/"), tileSize = 50)
