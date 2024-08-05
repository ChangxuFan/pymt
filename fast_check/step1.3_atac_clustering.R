source("objects.R")
addArchRThreads(threads = 4) 
work.dir <- "fast_check/sth/atac/"

work.dir.rna <- "fast_check/sth/per_sample/"
system("mkdir -p " %>% paste0(work.dir))

samples <- SAMPLES

utilsFanc::safelapply(samples, function(sample) {
  work.dir <- paste0(work.dir, "/", sample)
  work.dir.rna <- paste0(work.dir.rna, "/", sample)
  ao <- readRDS(paste0(work.dir, "/ao.Rds")) %>% recoverArchRProject()
  so <- readRDS(paste0(work.dir.rna, "/soi.Rds"))
  ao <- ao[ao$cellNames %in% get.cell.names.seurat(so = so, style = "ArchR"),]
  ao <- archr.cluster.pipe(ao = ao, so = so, minTSS = 10, work.dir = work.dir, do.harmony = F,
                           force = T, dimsToUse = 1:15)
  return()

}, threads = 4)

utilsFanc::safelapply(samples, function(sample) {
  work.dir <- paste0(work.dir, "/", sample)
  work.dir.rna <- paste0(work.dir.rna, "/", sample)
  ao <- readRDS(paste0(work.dir, "/ao_IterativeLSI.Rds")) %>% recoverArchRProject()
  so <- readRDS(paste0(work.dir.rna, "/soi.Rds"))
  archr.crosscheck.pipe(ao = ao, so = so, study.dir = paste0(work.dir, "/study/"),
                        samples = sample,
                        bc.metrics.file.list = BC.METRICS.FILE.LIST.exonOnly)
  return()
  
}, threads = 4)
