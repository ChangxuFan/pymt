source("objects.R")
wd <- paste0("publication/sth/pub2.2_DC_DEG_DAR/")

de <- readRDS("sth/de_noNZ/de_noNZ2.1_scc_BM/bulk.list.Rds")
de <- deseq2.summary(de)

deseq2.xyplot(pbl = de["scc_6"], publication = T,
              comp.vec = "tumor:ctrl", comp.meta = "type",
              transformation = function(x) log2(x + 1),
              plot.dir = paste0(wd, "/"), 
              root.name = "xy_DEG_scc6", 
              device = "pdf", 
              plot.each.comp = T)

da <- readRDS("sth/da_noNZ/da_noNZ2.1_scc_BM//ao.bulk.list.Rds")
da <- deseq2.summary(da, log2fc.cutoff = 0.5, padj.cutoff = 0.1)

deseq2.xyplot(pbl = da["scc_6"], publication = T,
              comp.vec = "tumor:ctrl", comp.meta = "type",
              transformation = function(x) log2(x + 1),
              plot.dir = paste0(wd, "/"),
              root.name = "xy_DAR_scc6_log2fc0.5padj0.1",
              device = "pdf",
              plot.each.comp = T)
