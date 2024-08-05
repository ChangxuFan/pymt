source("objects.R")
work.dir <- "fast_check/sth/atac/"
samples <- SAMPLES.BM
aol <- lapply(samples, function(x) {
  readRDS(paste0(work.dir, "/", x, "/ao_IterativeLSI.Rds")) %>% 
    return()
})
names(aol) <- samples

semi.clean.rmClus.list <- list(
  paste0("C", c(3, 4, 12)),
  paste0("C", c(8, 12, 13, 4, 1,2, 3)),
  paste0("C", c(8, 9, 13)),
  paste0("C", c(10))
)

names(semi.clean.rmClus.list) <- samples

semi.clean.rmSeuratClus.list <- list(
  "9", 
  "1",
  "4",
  "5"
)
names(semi.clean.rmSeuratClus.list) <- samples

round.2.pct <- list(
  c(0.02, 0.03, 0.04),
  c(0.02, 0.03, 0.04),
  c(0.02, 0.03, 0.04),
  c(0.02, 0.03, 0.04)
)

names(round.2.pct) <- samples

stats.df <- utilsFanc::safelapply(SAMPLES.good, function(sample) {
  ao <- aol[[sample]]
  semi.dir <- paste0(work.dir, "/", sample, "/semiClean/")
  dir.create(semi.dir, showWarnings = F, recursive = T)
  
  stats.df <- data.frame(sample = sample, pre.remove = ao$cellNames %>% length())
  
  ao <- ao[! ao$Clusters %in% semi.clean.rmClus.list[[sample]],]
  ao <- ao[! ao$seurat_clusters %in% semi.clean.rmSeuratClus.list[[sample]], ]
  
  saveRDS(ao, paste0(semi.dir, "/ao.Rds"))
  stats.df$post.remove <- ao$cellNames %>% length()
  write.table(stats.df, paste0(semi.dir, "/stats.tsv"), quote = F, col.names = T,
              row.names = F, sep = "\t")
  archr.titrate.doublets(ao = ao, plot.dir = semi.dir, root.name = "ao_semiClean",
                         label.vars = c("Clusters", "seurat_clusters"),
                         remove.frac = round.2.pct[[sample]])
  return(stats.df)
}, threads = 4) %>% Reduce(rbind, .)

write.table(stats.df, paste0(work.dir, "/semiClean_stats.tsv"), quote = F, col.names = T,
            row.names = F, sep = "\t")

round.2.pct <- c(0.04, 0.04, 0.04, 0.04)
names(round.2.pct) <- samples

utilsFanc::safelapply(samples, function(sample) {
  cleaned.v1.dir <- paste0(work.dir, "/", sample, "/cleaned_v1/")
  dir.create(cleaned.v1.dir, showWarnings = F, recursive = T)
  ao.source <- paste0(work.dir,"/", sample, "/semiClean/ao_semiClean_",round.2.pct[sample], "_ao.Rds")
  cmd <- paste0("cp ", ao.source, " ", cleaned.v1.dir, "/ao.Rds")
  
  print(cmd); system(cmd)
}, threads = 4)

