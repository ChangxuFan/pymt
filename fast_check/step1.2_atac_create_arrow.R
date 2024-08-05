source("objects.R")
addArchRThreads(threads = 4) 
work.dir <- "fast_check/sth/atac/"

system("mkdir -p " %>% paste0(work.dir))

frag.files <- FRAG.FILES.exonOnly

aol <- utilsFanc::safelapply(names(frag.files), function(sample) {
  frag.file <- frag.files[sample]
  sub.work.dir <- paste0(work.dir, "/", sample)
  ao <- ao.gen(frag.files = frag.file, minTSS = 10, minFrags = 1000,
               work.dir = sub.work.dir)
  saveRDS(ao, paste0(sub.work.dir, "/ao.Rds"))
  return(ao)
}, threads = 4)