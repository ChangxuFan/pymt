source("objects.R")

wd <- "publication/sth/pub5.2_SPBM_de/"
dir.create(wd)

clusters <- paste0("scc_", c("1", "9"))
names(clusters) <- clusters

de.dir <- "de_noNZ/sth/de_noNZ2.2_scc_BMSP/"
de <- readRDS(paste0(de.dir, "/bulk.list.Rds"))

de <- de[clusters]
de <- deseq2.summary(de)

de$scc_1$summary$up.genes
de$scc_1$summary$down.genes

markers <- c("Selp", "Slc35d3", "Hes1")

lapply(clusters, function(cluster) {
  deseq2.xyplot(pbl = de[cluster], comp.vec = "SP:BM", comp.meta = "tissue",
                transformation = function(x) log2(x + 1),
                plot.dir = paste0(wd, "/plot_xy/"), 
                root.name = paste0("xy_DEG", "_", cluster),
                add.label = T, repel.direction = 'x',
                hjust = 0.5, vjust = 0.5, use.repel = T, italic.label = T,
                label.list = markers, text.size = 6 * 0.36,
                device = "pdf", publication = T,
                plot.each.comp = F)
})

# UMAPs of gene expression
lapply(c(T), function(bOrder) {
  lapply(c("Rbpj", "Rbpjl", "Gzmb", "Slc35d3", "Hes1", "Col18a", "Scin"), function(marker) {
    pl <- plot.panel.list(soi, hide.legend = F,
                          panel.list = marker, 
                          order = bOrder, 
                          split.by = "tissue",
                          assay = "RNA", invisible = T, 
                          plot.out = paste0(wd, "/panel_DEG_", bOrder, "_", marker, ".pdf"),
                          raster = T, 
                          label = F, label.size = 1.8,
                          publication = T, 
                          italic.title = T,
                          n.col = 2,
                          pt.size = 3, auto.adjust.raster.pt.size = F,
                          sub.width = 1, sub.height = 1)
  })
})
