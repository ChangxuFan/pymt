source("objects.R")
soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")
wd <- paste0("publication/sth/pub3.1_c20_HSPC_module")

#>>>>>>>> We first identify HSPC markers
markers <- readRDS("sync_all/sth/rna/int_rpca20/markers.Rds")
markers <- markers %>% filter(cluster == 1, pct.2 < 0.15, pct.1 > 0.25) %>% pull(gene)
# [1] "Vldlr"  "Hlf"    "Gcnt2"  "Ly6a"   "Selp"   "Rbpms"  "Dipk1b" "Calcrl" "Myct1"  "Gpr171"

samples.sp <- SAMPLES.SP
samples <- samples.sp
samples.display <- uname(SAMPLE.MAP[samples])
soi <- soi[, soi$sample %in% samples]

lapply(c("eachSample"), function(split.type) {
  if (split.type == "combineReplicate") {
    # split.by <- "tissue"
    # split.order <- c("Ctrl", "Tumor")
  } else {
    split.by <- "sample"
    split.order <- samples
  }
  
  plot.panel.list(soi, panel.list = markers,
                  order = T, assay = "RNA", invisible = T, 
                  plot.out = paste0(wd, "/marker_genes/marker_panel_", split.type, ".pdf"),
                  sample = samples,
                  split.by = split.by, split.order = split.order, use.split.as.title = F,
                  raster = T, label = F, publication = T, ymax = NULL, ymin = 0.1, italic.title = F,
                  sub.width = 1.2, sub.height = 1.2)
  
})


#>>>>>>>> See if we can cluster c20 together with c1 using marker genes from various clusters:
# first get a gene by cluster matrix:
soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")
soi$sample.display <- SAMPLE.MAP[soi$sample]

groupings <- paste0(soi$sample.display, "..", soi$seurat_clusters)
names(groupings) <- colnames(soi)

rawmat <- aggregate.fanc(mat = soi@assays$RNA@counts, margin = 2, groupings = groupings,
                         binarize = F, take.mean = F, sort = T)

saveRDS(rawmat, "sync_all/sth/rna/int_rpca20/gene_by_sampleAndCluster/rawmat.Rds")

soi$groupings <- groupings
coldata <- soi@meta.data[, c("groupings", "sample.display", "seurat_clusters", "type")] %>% 
  unique() %>% dplyr::rename(sample = groupings) %>% arrange(sample)

rownames(coldata) <- coldata$sample
saveRDS(coldata, "sync_all/sth/rna/int_rpca20/gene_by_sampleAndCluster/coldata.Rds")

dds <- DESeqDataSetFromMatrix(rawmat, colData = coldata, design = ~type)
dds <- estimateSizeFactors(dds)

norm.mat <- DESeq2::counts(dds, normalized = T)
saveRDS(norm.mat, "sync_all/sth/rna/int_rpca20/gene_by_sampleAndCluster/norm.mat.Rds")

# next we find genes that we trust as lineage markers:
markers <- readRDS("sync_all/sth/rna/int_rpca20/markers.Rds")

top10 <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_logFC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() %>% as.data.frame()

any(duplicated(top10$gene)) # TRUE

top10 <- top10 %>% dplyr::filter(!duplicated(gene))

# plot heatmap:
lapply(samples.display, function(sample) {
  mat <- norm.mat[, grepl(sample, colnames(norm.mat))]
  colnames(mat) <- sub(".+\\.\\.", "", colnames(mat))
  # remove genes that are zero in all clusters.
  mat <- mat[rowSums(mat) > 0, ]
  mat.scaled <- mat %>% t() %>% scale() %>% t()
  genes <- intersect(top10$gene, rownames(mat.scaled))
  
  plot.out <- paste0(wd, "/marker_gene_hm_clustering/top10_autoCluster_wKmean_", sample,".png")
  
  plot.mat.auto.clustering(mat = mat.scaled[genes,], 
                           cluster_rows = T, cluster_columns = T, 
                           height = 1000, 
                           show_row_names = F, show_column_dend = T, 
                           plot.out = plot.out)
})
# It works well for all samples

# now we plot a mean of all samples and present that one in the manuscript:
sp.norm.mat <- norm.mat[, grepl("SP", colnames(norm.mat))]
ave.groupings <- colnames(sp.norm.mat) %>% sub(".+\\.\\.", "", .)
names(ave.groupings) <- colnames(sp.norm.mat)

ave.mat <- aggregate.fanc(mat = norm.mat, margin = 2, groupings = ave.groupings, binarize = F, 
                          take.mean = T, sort = T)
ave.mat.scaled <- ave.mat %>% t() %>% scale() %>% t()

plot.mat.rank.row(mat = ave.mat.scaled[top10$gene,],
                  plot.out = paste0(wd, "/marker_gene_hm_clustering/all_SPsamples_top10_rowRanking.pdf"),
                  width = 2.5, height = 3, title = "z-score",
                  hm.colors = c("darkblue", "white" , "deeppink"), hm.values = c(-3, 0,  3))
 