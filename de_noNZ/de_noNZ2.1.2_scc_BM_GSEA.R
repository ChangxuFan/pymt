source("objects.R")
wd <- "de_noNZ/sth/de_noNZ2.1_scc_BM/"

de <- readRDS(paste0(wd, "/bulk.list.Rds"))

rnk <- de.2.rnk(pbl = de, microarray.mode = F, pbl.slot = "res",
                rank.by = "log_p", logp.max = 8, to.upper = T,
                pbl.topn = 12000, pbl.samples = SAMPLES.BM, comp = "tumor:ctrl", 
                out.dir = paste0(wd, "/rnk/"))

gsea.fanc(rnk.vec = rnk, 
          gmt.vec = c(GMT.FILES.HALL.GOBP, 
                      GMT.FILES[c("c7.all.v7.2", "c8.all.v7.2")]),
          out.dir = paste0(wd, "/gsea/"), thread.rnk = 1, thread.gmt = 12,
          n.plot = 80, run = T, force = F, parse = F, parse.only = F, zip = F)

