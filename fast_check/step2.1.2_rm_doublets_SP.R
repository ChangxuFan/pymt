source("objects.R")
work.dir <- "fast_check/sth/atac/"
samples <- SAMPLES.SP

aol <- lapply(samples, function(x) {
  readRDS(paste0(work.dir, "/", x, "/ao_IterativeLSI.Rds")) %>% 
    return()
})
names(aol) <- samples

lapply(samples, function(x) {
  archr.titrate.doublets(ao = aol[[x]],
                         plot.dir = paste0(work.dir, "/", x, "/doublets/"),
                         root.name = "ao", 
                         label.vars = c("Clusters", "seurat_clusters"),
                         remove.frac = c(0.04, 0.08, 0.10))
})

semi.clean.rmClus.list <- list(
  c(1, 2, 5, 8, 10, 13),
  c(5, 8, 9, 14, 15),
  c(1, 4, 8, 11, 12, 13),
  c(4, 7, 8, 9, 12)
) %>% lapply(function(x) paste0("C", x))
names(semi.clean.rmClus.list) <- samples

semi.clean.rmSeuratClus.list <- list(
  c(9, 10, 13:18),
  c(9, 11, 12:14),
  c(9, 10, 12, 14, 15, 16),
  c(8, 10, 11, 12, 14, 15) 
) %>% lapply(as.character)
names(semi.clean.rmSeuratClus.list) <- samples


round.2.pct <- list(
  c(0.02, 0.03, 0.04),
  c(0.02, 0.03, 0.04),
  c(0.02, 0.03, 0.04),
  c(0.02, 0.03, 0.04)
)

names(round.2.pct) <- samples

stats.df <- utilsFanc::safelapply(samples, function(sample) {
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
}, threads = 1) %>% Reduce(rbind, .)

write.table(stats.df, paste0(work.dir, "/semiClean_stats.tsv"), quote = F, col.names = T,
            row.names = F, sep = "\t")

round.2.pct <- c(0.03, 0.04, 0.02, 0.04)
names(round.2.pct) <- samples

utilsFanc::safelapply(samples, function(sample) {
  cleaned.v1.dir <- paste0(work.dir, "/", sample, "/cleaned_v1/")
  dir.create(cleaned.v1.dir, showWarnings = F, recursive = T)
  ao.source <- paste0(work.dir,"/", sample, "/semiClean/ao_semiClean_",round.2.pct[sample], "_ao.Rds")
  cmd <- paste0("cp ", ao.source, " ", cleaned.v1.dir, "/ao.Rds")
  
  print(cmd); system(cmd)
}, threads = 1)
