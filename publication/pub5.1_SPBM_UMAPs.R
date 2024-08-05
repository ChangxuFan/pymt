source("objects.R")

wd <- "publication/sth/pub5.1_SPBM_UMAPs/"
dir.create(wd)
soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")

clusters <- c("1", "9")
soi$scc[is.na(soi$scc)] <- "22"

lapply(clusters, function(cluster) {
  pl <- plot.panel.list(soi, hide.legend = T, 
                        panel.list = "scc", 
                        binarize.panel = T, binarize.items = cluster,
                        order = F, 
                        assay = "RNA", invisible = T, 
                        split.by = "tissue", use.split.as.title = T,
                        plot.out = paste0(wd, "/highlight_c", cluster, "_SPvsBM.pdf"),
                        raster = T, label = T, label.size = 1.8,
                        publication = T, 
                        italic.title = F,
                        sub.width = 1.5, sub.height = 1.5)
  
})
