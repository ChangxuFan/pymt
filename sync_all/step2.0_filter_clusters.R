source("objects.R")

samples <- SAMPLES[1:4]
sol <- paste0("fast_check/sth/per_sample/", samples, "/soi.Rds") %>% 
  lapply(readRDS)
names(sol) <- samples

rm.list <- list(
  c(5, 7, 17, 11, 10, 6, 13),
  c(13, 7, 10, 5, 8, 9, 19, 20, 14, 12),
  c(0, 10,5, 15, 13),
  c(4, 10, 12, 13, 14, 15)
) %>% lapply(as.character) %>% `names<-`(samples)

sol <- lapply(samples, function(sample) {
  so <- sol[[sample]]
  clusters.rm <- rm.list[[sample]]
  so <- so[, ! so$seurat_clusters %in% clusters.rm]
  return(so)
})
names(sol) <- samples

pl <- lapply(sol, function(so) {
  plot.panel.list("Ly6a", so, assay = "RNA", order = F) %>% return()
})
trash <- wrap.plots.fanc(pl, plot.out = "sync_all/step2.0_plots/validate_filter.png")

lapply(samples, function(sample) {
  saveRDS(sol[[sample]], paste0("fast_check/sth/per_sample/", samples, "/soi.Rds"))
})
