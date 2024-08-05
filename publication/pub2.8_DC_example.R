source("objects.R")
wd <- "publication/sth/pub2.8_DC_example/"
dir.create(wd)

de <- readRDS("de_noNZ/sth/de_noNZ2.1_scc_BM/bulk.list.Rds")
de <- deseq2.summary(de)
de$scc_6$summary$up.genes
"Socs3" %in% de$scc_6$summary$up.genes # TRUE

soi$sample.pub <- SAMPLE.MAP[soi$sample]

lapply(c("sample.pub", "type"), function(split.by) {
  lapply(c(T, F), function(bOrder) {
    plot.panel.list(soi, panel.list = c("Socs3"),
                    order = bOrder, assay = "RNA", 
                    sample = SAMPLES.BM,
                    invisible = T, 
                    plot.out = paste0(wd, "/Socs3/split_by_", split.by , "_", bOrder , ".pdf"),
                    split.by = split.by, use.split.as.title = T,
                    raster = T, label = F, publication = T,
                    auto.adjust.raster.pt.size = F, pt.size = 5,
                    ymax = NULL, italic.title = F,
                    sub.width = 1.2, sub.height = 1.2)
  })
})

soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")

lapply(c("type"), function(split.by) {
  lapply(c(T, F), function(bOrder) {
    plot.panel.list(soi, panel.list = c("Tgfbi"),
                    order = bOrder, assay = "RNA", 
                    sample = SAMPLES.BM,
                    invisible = T, 
                    plot.out = paste0(wd, "/Tgfbi/Tgfbi_split_by_", split.by , "_", bOrder , ".pdf"),
                    split.by = split.by, use.split.as.title = T,
                    raster = T, label = F, publication = T,
                    auto.adjust.raster.pt.size = F, pt.size = 1,
                    ymax = NULL, italic.title = F,
                    sub.width = 1.2, sub.height = 1.2, n.col = 1)
  })
})
