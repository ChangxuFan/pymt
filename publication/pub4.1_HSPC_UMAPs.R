source("objects.R")
# we'll come back to this.
wd <- "publication/sth/pub4.1_HSPC_UMAPs/"
dir.create(wd)
soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")

clusters <- c("1", "9")

lapply(clusters, function(cluster) {
  pl <- plot.panel.list(soi, sample = SAMPLES.BM, hide.legend = T,
                        panel.list = "scc", 
                        binarize.panel = T, binarize.items = cluster,
                        order = F, 
                        assay = "RNA", invisible = T, 
                        plot.out = paste0(wd, "/highlight_c", cluster, "_BM.pdf"),
                        raster = T, label = T, label.size = 1.8,
                        publication = T, 
                        italic.title = F,
                        sub.width = 1.5, sub.height = 1.5)
  
})

