source("objects.R")
wd <- "publication/sth/pub1.1_umaps/"
dir.create(wd)

soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")
#>>>>>>>> Plot cluster UMAP
p <- umap.fanc(soi, group.by = "seurat_clusters", 
          show.legends = F, axis.type = "nothing", 
          pt.size = 0.2, label.size = 2.5,
          label.groups = T,
          plot.out = paste0(wd, "/RNA_UMAP.pdf"))

color.map <- data.frame(cell_type = as.character(sort(unique(soi$seurat_clusters))),
                        color = utilsFanc::gg_color_hue(n = length(unique(soi$seurat_clusters))))
color.map$cat <- NA

color.map$cat[color.map$cell_type %in% c("1", "9", "20")] <- "HSPC"
color.map$cat[color.map$cell_type %in% c("0", "2", "11", "17")] <- "Ery"
color.map$cat[color.map$cell_type %in% c("7", "14")] <- "Mega"
color.map$cat[color.map$cell_type %in% c("4", "18")] <- "Neu"
color.map$cat[color.map$cell_type %in% c("5")] <- "Myeloid"
color.map$cat[color.map$cell_type %in% c("3")] <- "Mono/Mac/DC"
color.map$cat[color.map$cell_type %in% c("10")] <- "Mono/Mac"
color.map$cat[color.map$cell_type %in% c("6", "19")] <- "DC"
color.map$cat[color.map$cell_type %in% c("16")] <- "Lympho"
color.map$cat[color.map$cell_type %in% c("8", "15")] <- "Cycling"
color.map$cat[is.na(color.map$cat)] <- "Other"

cell.types <- c("HSPC", "Ery", "Mega", "Myeloid", "Neu", "Mono/Mac/DC", "Mono/Mac", "DC", "Lympho", "Cycling", "Other")

utilsFanc::check.intersect(color.map$cat, "x", cell.types, "y")

color.map <- color.map[utilsFanc::sort.by(x = color.map$cat, y = cell.types, return.order = T),]
write.table(color.map, paste0(wd, "/colorMap_wCat.tsv"), col.names = T, row.names = F,
            quote = F, sep = "\t")

color.map$cat <- factor(color.map$cat, levels = cell.types)

color.map <- color.map %>% split(f = color.map$cat) %>% 
  lapply(function(df) {
    cat <- df$cat[1] %>% as.character()
    df <- rbind(data.frame(cell_type = cat, color = ""), df[, 1:2])
  }) %>% do.call(rbind, .)

rownames(color.map) <- NULL


legend.plot(df = color.map, pt.size = 1.5, text.size = 6,
            out.file = paste0(wd, "/RNA_UMAP_legend.pdf"),
            width = 1, height = 4)

#>>>>>>>> Split cluster UMAP by sample:
samples <- INFO.filtered$sample
sample.map <- c("BM ctrl rep1", "BM tumor rep1", "BM ctrl rep2", "BM tumor rep2",
                "SP tumor rep1", "SP tumor rep2", "SP tumor rep3", "SP tumor rep4")
names(sample.map) <- samples

soi$sample.display <- sample.map[soi$sample]
soi$tissue.type <- paste0(soi$tissue, ' ', soi$type)
table(soi$tissue.type)
types <- soi$tissue.type %>% unique()

sol <- lapply(samples, function(sample) return(soi[, soi$sample == sample]))
names(sol) <- samples


pl <- lapply(samples, function(sample) {
  p <- umap.fanc(sol[[sample]], group.by = "seurat_clusters", 
                 show.legends = F, axis.type = "nothing", 
                 pt.size = 0.4, label.size = 2.5,
                 label.groups = T, 
                 plot.title = sample.map[sample],
                 plot.out = paste0(wd, "/RNA_UMAP_by_sample/", sample, ".pdf"))

})


pl <- lapply(types, function(type) {
  p <- umap.fanc(soi[, soi$tissue.type == type], group.by = "seurat_clusters", 
                 show.legends = F, axis.type = "nothing", 
                 pt.size = 0.1, label.size = 2.5,
                 label.groups = T, 
                 plot.title = type,
                 plot.out = paste0(wd, "/RNA_UMAP_by_sample/", type, ".pdf"))
  
})

#>>>>>>>> Plot marker genes
markers <- c("Hlf", "Cebpe", "Irf8", "Lyz2", "Epor","Pf4", "Dntt")
lapply(c(T, F), function(bOrder) {
  plot.panel.list(panel.list = markers, obj = soi, order = bOrder, assay = "RNA", 
                  plot.out = paste0(wd, "/gex_markers_", bOrder, ".pdf"),
                  raster = T, 
                  label = F, # label.size = 1.8,
                  publication = T, invisible = T, n.col = 2, 
                  sub.width = 1, sub.height = 1,
                  pt.size = 0.42, auto.adjust.raster.pt.size = F,
                  return.list = F, threads = 1)
  
})

#>>>>>>>> bubble plot:
m.bub <- c(
  "Hlf", "Procr", "Ly6a", "Myct1",
  "Klf1", "Gata1", "Epor", "Pf4", "Itga2b", 
  "Cebpe", "Gfi1", "Elane", "Gstm1", 
  "Irf8", "Csf1r", "Ly86", "Ccr2",
  "Flt3", "Dntt", "Il7r", 
  "Mki67", "Top2a", "Cdk1"
)
soi$cell_type <- soi$seurat_clusters
color.map.cat <- read.table(paste0(wd, "colorMap_wCat.tsv"), header = T, comment.char = "")

#>>>>>>>> We first plot by cluster
lapply(c("h", "v"), function(x) {
  text.size = 10
  if (x == "h") {
    soi$cell_type <- factor(soi$cell_type, levels = rev(color.map.cat$cell_type))
    width <- 8
    height <- 8
  } else {
    soi$cell_type <- factor(soi$cell_type, levels = color.map.cat$cell_type)
    width <- 8
    height <- 8
  }
  
  if (x == "v") {
    m.bub <- rev(m.bub)
  }
  
  p <- DotPlot(soi, assay = "RNA", features = m.bub, cluster.idents = F, group.by = "cell_type", 
               dot.scale = 5, dot.min = 0.01) +
    theme() + xlab(NULL) + ylab(NULL)
  
  p <- p %>% utilsFanc::theme.fc.1(
    rotate.x.45 = T, italic.x = F, rm.x.ticks = F, text.size = text.size)
  
  if (x == "v") {
    p <- p + theme(axis.text.y = element_text(size = text.size, family = "Arial", face = "italic",
                                              color = "black"))
    p <- p + coord_flip()
  }
  
  p <- p + theme(legend.key.size = unit(0.2, "in"))
  ggsave(paste0(wd, "/bubblePlot_", x,".pdf"), p,  
         device = cairo_pdf, width = width, height = height, 
         dpi = 300)
  
})

####### ATAC-UMAP
aoi <- readRDS("sync_all/sth/atac/merge_v1/ao.Rds")
aoi$Clusters <- sub("C", "", aoi$Clusters)
aoi$Clusters <- factor(aoi$Clusters, levels = gtools::mixedsort(unique(aoi$Clusters)))

p <- umap.fanc(aoi, group.by = "Clusters", 
               show.legends = F, axis.type = "nothing", 
               pt.size = 0.2, label.size = 2.5,
               label.groups = T,
               plot.out = paste0(wd, "/ATAC_UMAP.pdf"))

soi <- seurat.add.archr.embed(so = soi, ao = aoi, ao.embedding = "UMAP")

markers <- c("Hlf", "Cebpe", "Irf8", "Lyz2", "Epor","Pf4", "Dntt")
lapply(c(T, F), function(bOrder) {
  plot.panel.list(panel.list = markers, obj = soi, order = bOrder, assay = "RNA", 
                  reduction = "ArchRUMAP",
                  plot.out = paste0(wd, "/gex_markers_", bOrder, "_on_ATAC_UMAP.pdf"),
                  raster = T, 
                  label = F, # label.size = 1.8,
                  pt.size = 0.42, auto.adjust.raster.pt.size = F,
                  publication = T, invisible = T, n.col = 2, 
                  sub.width = 1, sub.height = 1,
                  return.list = F, threads = 1)
  
})

