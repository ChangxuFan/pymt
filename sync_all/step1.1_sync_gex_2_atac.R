source("objects.R")
rna.from.master <- "fast_check/sth/rna/per_sample/"
atac.from.master <- "fast_check/sth/atac/"
atac.type <- "cleaned_v1"

samples <- SAMPLES
rna.dir <- "sync_all/sth/rna"
atac.dir <- "sync_all/sth/atac"
threads <- length(samples)
system(paste0("mkdir -p ", rna.dir, " ", atac.dir))

utilsFanc::safelapply(samples, function(sample) {
  rna.dir <- paste0(rna.dir, "/", sample)
  atac.dir <- paste0(atac.dir, "/", sample)
  system(paste0("mkdir -p ", rna.dir, " ", atac.dir))
  so <- readRDS(paste0(rna.from.master, "/", sample, "/soi.Rds"))
  ao <- readRDS(paste0(atac.from.master, "/", sample,"/", atac.type, "/ao.Rds")) %>% recoverArchRProject()
  synco <- seurat.sync.archr(so = so, ao = ao)
  so <- synco$so
  ao <- synco$ao
  so <- seurat.add.archr.meta(so = so, ao = ao)
  saveRDS(so, paste0(rna.dir, "/soi.Rds"), compress = F)
  ao <- save.archr.fanc(ao = ao, new.work.dir = atac.dir, copy.arrow = F, over.write = T)
  return()
}, threads = threads)