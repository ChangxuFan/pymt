source("objects.R")
work.dir <- "sync_all/sth/rna/int_rpca20"
plot.dir <- paste0(work.dir, "/plots/")

samples <- INFO.filtered$sample

soi <- readRDS(paste0(work.dir, "/soi.Rds"))
sol <- readRDS(paste0(work.dir, "/sol_v1.Rds"))

soi <- cluster.pipe(
  soi = soi, assay = "integrated",
  pc.dim = 1:30, cluster.resolution = 0.7, cluster.seed = 0,
  umap.seed = 10,
  work.dir = work.dir, plot.dir = plot.dir,
  project.name = "rpca", metas.include = "sample",
  save.rds = T, do.sct = F,
  plot.common.markers = T, split.by = "sample", split.order = samples
)

int.corr.plot(
  soi = soi, sol = sol, cluster.ident = "seurat_clusters", order = F,
  plot.dir = paste0(work.dir, "/int_corr_plot/")
)

int.corr(q.vec = sol, s = soi, meta = "seurat_clusters", 
         out.file = paste0(work.dir, "/int_corr_plot/int_corr.tsv"))
