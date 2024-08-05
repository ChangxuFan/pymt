source("objects.R") 
library(future)
work.dir <- "sync_all/sth/rna/int_rpca20"
plot.dir <- paste0(work.dir, "/plots/")
dir.create(plot.dir, showWarnings = F, recursive = T)

samples <- SAMPLES
sol <- utilsFanc::safelapply(paste0("fast_check/sth/per_sample/", samples, "/soi.Rds"), readRDS, threads = 4)
names(sol) <- samples

sol <- utilsFanc::safelapply(samples, function(sample) {
    so <- sol[[sample]]
    so$sample <- sample
    so <- RenameCells(so, new.names = colnames(so) %>% sub(".+#", "", .))
    so <- RenameCells(so, new.names = get.cell.names.seurat(so, style = "ArchR"))
}, threads = 8)
names(sol) <- samples
lapply(sol, function(so) head(colnames(so)))

saveRDS(sol, paste0(work.dir, "/sol_v1.Rds"))

features <- SelectIntegrationFeatures(object.list = sol, nfeatures = 3000)
sol <- PrepSCTIntegration(object.list = sol, anchor.features = features)
sol <- utilsFanc::safelapply(X = sol, FUN = RunPCA, features = features, assay = "SCT", threads = 8)
options(future.globals.maxSize = 40 * 1000 * 1024^2)
plan("sequential")

immune.anchors <- FindIntegrationAnchors(
    object.list = sol, normalization.method = "SCT",
    anchor.features = features, dims = 1:30,
    reduction = "rpca", k.anchor = 20,
    reference = c(1, 7))

soi.rpca20 <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:30)

saveRDS(soi.rpca20, paste0(work.dir, "/soi_v1.Rds"), compress = F)

soi.rpca20 <- cluster.pipe(
  soi = soi.rpca20, assay = "integrated",
  pc.dim = 1:30, cluster.resolution = 0.8,
  work.dir = work.dir, plot.dir = plot.dir,
  project.name = "rpca", metas.include = "sample",
  save.rds = T, do.sct = F,
  plot.common.markers = T, split.by = "sample", split.order = samples
) 

int.corr.plot(
  soi = soi.rpca20, sol = sol, cluster.ident = "seurat_clusters", order = F,
  plot.dir = paste0(work.dir, "/int_corr_plot/")
)

int.corr(q.vec = sol, s = soi.rpca20, meta = "seurat_clusters", 
         out.file = paste0(work.dir, "/int_corr_plot/int_corr.tsv"))

