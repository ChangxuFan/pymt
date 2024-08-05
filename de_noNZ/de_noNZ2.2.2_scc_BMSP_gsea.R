source("objects.R")
wd <- "de_noNZ/sth/de_noNZ2.2_scc_BMSP/"

de <- readRDS(paste0(wd, "/bulk.list.Rds"))

rnk <- de.2.rnk(pbl = de, microarray.mode = F, pbl.slot = "res",
                rank.by = "log_p", logp.max = 8, to.upper = T,
                pbl.topn = 12000, pbl.samples = SAMPLES, comp = "SP:BM", 
                out.dir = paste0(wd, "/rnk/"))

gsea.fanc(rnk.vec = rnk, gmt.vec = GMT.FILES.HALL.GOBP,
          out.dir = paste0(wd, "/gsea/"), thread.rnk = 1, thread.gmt = 12,
          n.plot = 80, run = T, force = F, parse = F, parse.only = F, zip = F)

