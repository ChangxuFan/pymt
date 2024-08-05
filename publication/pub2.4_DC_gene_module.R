rm(list = ls())
source("objects.R")
wd <- "publication/sth/pub2.4_DC_gene_module/"
dir.create(wd)

samples <- SAMPLES.BM

de <- readRDS("de_noNZ/sth/de_noNZ2.1_scc_BM/bulk.list.Rds")
de <- deseq2.summary(de)

de.denardo <- readRDS("~/hmtp/scAR/reference/denardo2018/sth/denardo1.1_geo2R/bulk.list.Rds")

de.denardo <- deseq2.summary(de.denardo, padj.cutoff = 1, log2fc.cutoff = 0)
de.print.summary(de.denardo)
# simply dividing genes into up and down
# [1] "cluster: CDP"
# [1] "n.up:"
# [1] 8206
# [1] "n.down:"
# [1] 8014
# [1] "n.de:"
# [1] 16220

de.murphy <- readRDS("~/hmtp/scAR/reference/kim_murphy2023/sth/km23_2.1_deFromRaw/bulk.list.Rds")
de.murphy <- deseq2.summary(de.murphy, padj.cutoff = 1, log2fc.cutoff = 0)
de.print.summary(de.murphy)
# [1] "cluster: CDP"
# [1] "n.up:"
# [1] 8376
# [1] "n.down:"
# [1] 8289
# [1] "n.de:"
# [1] 16665

#>>>>>>>> Filter to only contain DEGs with the same trend in both reference datasets:

up.genes <- de$scc_6$summary$up.genes %>%  #start from 77
  .[. %in% de.denardo$CDP$summary$up.genes] %>% 
  .[. %in% de.murphy$CDP$summary$up.genes]
length(up.genes) # 51 genes!

saveRDS(up.genes, paste0(wd, "/scc6_up_genes_post_filter.Rds"))


#>>>>>>>> Plot a heatmap to cluster the genes based on their behaviour in scc6 vs scc1
genes.shared <- intersect(de$scc_1$bulkNorm$gene, de$scc_6$bulkNorm$gene)
raw.mat <- cbind(de$scc_1$bulk.mat[genes.shared,], de$scc_6$bulk.mat[genes.shared,])
colnames(raw.mat) <- paste0(c(rep("scc1_",4), rep("scc_6_",4)), 
                            colnames(raw.mat))

col.data <- data.frame(sample = colnames(raw.mat)) %>% 
  dplyr::mutate(type = stringr::str_extract(sample, "tumor|ctrl"))
rownames(col.data) <- col.data$sample

dds <- DESeqDataSetFromMatrix(raw.mat, colData = col.data, design = ~type)
dds <- estimateSizeFactors(dds)

norm.mat <- DESeq2::counts(dds, normalized = T)

sum(!up.genes %in% rownames(norm.mat)) # only 4
not.in.scc1 <- up.genes[!up.genes %in% rownames(norm.mat)]

deg.mat <- norm.mat[intersect(up.genes, rownames(norm.mat)),]
deg.mat.scaled <- deg.mat %>% t() %>% scale() %>% t()

plot.mat.auto.clustering(mat = deg.mat.scaled, cluster_rows = T, cluster_columns = F, 
                         height = 1400, k.m = 3,
                         show_row_names = T, plot.out = paste0(wd, "/deg_hm_clustering/test.png"))

# The heatmap is actually very helpful: some of the genes are clearly higher in scc_1, while others
# are clearly higher in scc_6. 

#>>>>>>>> divide up.genes into those higher in scc_6 or those higher in scc_1.
##########
# plot in figure
##########
deg.mean <- data.frame(scc1 = rowMeans(deg.mat[, 1:4]), scc6 = rowMeans(deg.mat[, 5:8]))

module <- list(
  scc1hi = rownames(deg.mean)[deg.mean$scc1 > deg.mean$scc6],
  scc6hi = c(not.in.scc1, rownames(deg.mean)[deg.mean$scc1 < deg.mean$scc6])
)

soi <- AddModuleScore(soi, features = module, name = "DEGsplit", assay = "RNA")

lapply(c("combineReplicate", "eachSample"), function(split.type) {
  if (split.type == "combineReplicate") {
    split.by <- "type.pub"
    split.order <- c("Ctrl", "Tumor")
  } else {
    split.by <- "sample"
    split.order <- samples
  }
  
  plot.panel.list(soi, panel.list = c("DEGsplit1", "DEGsplit2"),
                  order = T, assay = "RNA", invisible = T, 
                  plot.out = paste0(wd, "/seurat_module_sc6up_split_", split.type, ".pdf"),
                  sample = samples,
                  split.by = split.by, split.order = split.order, use.split.as.title = T,
                  raster = T, label = F, publication = T, ymax = NULL, ymin = 0.1, italic.title = F,
                  sub.width = 1.2, sub.height = 1.2)
  
})

length(module$scc6hi) # 37
saveRDS(module, paste0(wd, "/scc1_vs_6.Rds"))
#<<<<<<<< 