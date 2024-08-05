source("objects.R")
work.dir <- "sync_all/sth/atac/merge_v1/"
plot.dir <- paste0(work.dir, "/Plots/")
system(paste0("mkdir -p ", work.dir))

samples <- SAMPLES

arrow.files <- paste0("sync_all/sth/atac/per_sample/", samples, "/", samples, ".arrow")

aoi <- ArchRProject(ArrowFiles = arrow.files, outputDirectory = work.dir,
                    copyArrows = F)

soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")
cells <- intersect(aoi$cellNames, colnames(soi))
soi <- soi[, cells]
aoi <- aoi[cells, ]

saveRDS(aoi, paste0(work.dir, "/ao.Rds"))
saveRDS(soi, "sync_all/sth/rna/int_rpca20/soi.Rds")

aoi <- archr.add.seurat(ao = aoi, so = soi, meta = "seurat_clusters", as.factor = T)
aoi <- archr.cluster.pipe(ao = aoi, minTSS = 10, bc.metrics.file.list = NULL, dimsToUse = 1:15,
                          work.dir = work.dir, force = T, do.harmony = F, sample.order = samples,
                          so = soi)
