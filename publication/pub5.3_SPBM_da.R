source("objects.R")
wd <- "publication/sth/pub5.3_SPBM_da/"
dir.create(wd)

da <- readRDS("da_noNZ/sth/da_noNZ2.2_scc_BMSP/ao.bulk.list.Rds")
da.dir <- "da_noNZ/sth/da_noNZ2.2_scc_BMSP"

clusters <- paste0("scc_", c("1", "9"))
names(clusters) <- clusters

da <- da[clusters]
da <- deseq2.summary(da, log2fc.cutoff = 0.5, padj.cutoff = 0.1)

lapply(clusters, function(cluster) {
  deseq2.xyplot(pbl = da[cluster], 
                comp.vec = "SP:BM", comp.meta = "tissue",
                transformation = function(x) log2(x + 1),
                plot.dir = paste0(wd, "/plot_xy/"), 
                root.name = paste0("xy_DAR", "_", cluster),
                label.list = markers, text.size = 6 * 0.36,
                device = "pdf", publication = T,
                plot.each.comp = F)
})

homer.plot.table(da.homer.dir = paste0(da.dir, "/homer/"), use.regex = F, col.num = 3,
                 plot.logo = T,
                 motif.map = paste0(da.dir, "/motif_cluster_map/scc1.tsv"), 
                 regex.include = paste0("scc_", c(1), "/res.exp_padj_top_1000"),
                 regex.from = "scc_", 
                 out.dir = wd)
