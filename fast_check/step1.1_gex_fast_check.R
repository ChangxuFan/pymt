source("objects.R")
work.dir <- "fast_check/sth/"
plot.dir <- "fast_check/plots/"

paste0("mkdir -p ", plot.dir) %>% system()

sample.info <- data.frame(
  sample = SAMPLES,
  dir = paste0(EXON.ONLY.DIR, "/", SAMPLES, "/outs/filtered_feature_bc_matrix"),
  type = stringr::str_extract(SAMPLES, "tumor|ctrl"),
  tissue = stringr::str_extract(SAMPLES, "BM|SP")
)

write.table(sample.info, "fast_check/sample_info.tsv",
            col.names = T, row.names = F, quote = F, sep = "\t")

sol <- sc.qc.construct.so(sample.info = "fast_check/sample_info.tsv",
                          project.name = "fast_check",
                          mt.pattern = "mt-")

sol <- sc.qc.metadata.find.elbow.core(sol = sol, meta.term = "nFeature_RNA",
                                      override = T)

t <- sc.rank.qc(sol = sol, feature = "nFeature_RNA", ymax = 10000,
                plot.dir = plot.dir, project.name = "ori")

elbow.list <- list(
  BM_ctrl_ba1rep1 = 0,
  BM_tumor_ba1rep1 = 0,
  BM_ctrl_ba3rep1 = 0,
  BM_tumor_ba3rep1 = 0,
  SP_tumor_ba4rep1 = 750,
  SP_tumor_ba4rep2 = 750,
  SP_tumor_ba5rep1 = 800,
  SP_tumor_ba5rep2 = 800
)

sol <- sc.qc.metadata.find.elbow.core(sol = sol, meta.term = "nFeature_RNA",
                                      elbow.list = elbow.list,
                                      override = T)

t <- sc.rank.qc(sol = sol, feature = "nFeature_RNA", ymax = 10000,
                plot.dir = plot.dir, project.name = "lower_end")

sol <- sc.qc.elbow.filter(sol = sol,metas = c("nFeature_RNA"),
                          project.name = "fast_check", take.lower = F)

t <- sc.rank.qc(sol = sol, feature = "nFeature_RNA", ymax = 10000,
                plot.dir = plot.dir, project.name = "lower_end_post")

elbow.list <- list(
  BM_ctrl_ba1rep1 = 3800,
  BM_tumor_ba1rep1 = 4200,
  BM_ctrl_ba3rep1 = 3000,
  BM_tumor_ba3rep1 = 3000,
  SP_tumor_ba4rep1 = 3500,
  SP_tumor_ba4rep2 = 3500,
  SP_tumor_ba5rep1 = 4200,
  SP_tumor_ba5rep2 = 3800
)
sol <- sc.qc.metadata.find.elbow.core(sol = sol, meta.term = "nFeature_RNA",
                                      elbow.list = elbow.list,
                                      override = T)

t <- sc.rank.qc(sol = sol, feature = "nFeature_RNA", ymax = 10000,
                plot.dir = plot.dir, project.name = "higher_end")

sol <- sc.qc.metadata.find.elbow.core(sol = sol, meta.term = "percent.mt",
                                      override = T)
t <- sc.rank.qc(sol = sol, feature = "percent.mt", ymax = 100, plot.dir = plot.dir)

elbow.list <- list(
  BM_ctrl_ba1rep1 = 25,
  BM_tumor_ba1rep1 = 25,
  BM_ctrl_ba3rep1 = 25,
  BM_tumor_ba3rep1 = 38,
  SP_tumor_ba4rep1 = 25,
  SP_tumor_ba4rep2 = 25
)

sol <- sc.qc.metadata.find.elbow.core(sol = sol, meta.term = "percent.mt",
                                      elbow.list = elbow.list,
                                      override = T)

sol <- sc.qc.elbow.filter(sol = sol,metas = c("percent.mt", "nFeature_RNA"),
                          project.name = "fast_check")

t <- sc.rank.qc(sol = sol, feature = "nFeature_RNA", ymax = 10000, plot.dir = plot.dir, project.name = "post_filter")
t <- sc.rank.qc(sol = sol, feature = "percent.mt", ymax = 100, plot.dir = plot.dir, project.name = "post_filter")

saveRDS(sol, paste0(work.dir, "/sol_v1.Rds"), compress = F)

sol <- readRDS( paste0(work.dir, "/sol_v1.Rds"))

sol <- utilsFanc::safelapply(sol, function(so) {
  root.name <- so@meta.data$sample[1]
  utilsFanc::t.stat(root.name)
  so <- cluster.pipe(soi = so, assay = "SCT", pc.dim = 1:30, cluster.resolution = 0.7, 
                     work.dir = paste0(work.dir, "/per_sample/", root.name),
                     plot.dir = paste0(work.dir, "/per_sample/", root.name, "/plots/"), 
                     project.name = root.name, metas.include = "sample", 
                     sample.tsv = "fast_check/sample_info.tsv",
                     save.rds = T, do.sct = T, sct.n.vars = 3000, hm.dims = 11:30, 
                     plot.common.markers = T, bc.metrics.file.list = BC.METRICS.FILE.LIST.exonOnly)
  return(so)
}, threads = 4)

