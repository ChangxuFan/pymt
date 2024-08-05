source("objects.R")
wd <- "publication/sth/pub2.6_DC_ATAC_module/"
dir.create(wd)

da <- readRDS("da_noNZ/sth/da_noNZ2.1_scc_BM//ao.bulk.list.Rds")
da <- deseq2.summary(da, padj.cutoff = 0.1, log2fc.cutoff = 0.5)
scc_6.upDAR <- da$scc_6$summary$up.genes

peakmat <- readRDS("sync_all/sth/atac/merge_v1/peakmat_full.Rds")
aoi <- readRDS("sync_all/sth/atac/merge_v1/ao.Rds")

samples <- SAMPLES.BM

cell.list <- archr.get.cells.grid(ao = aoi, cluster.ident = "scc", clusters = c("1", "6"),
                                  group.ident = "sample", groups = samples, melt = T)
peaks <- c(da$scc_1$bulkNorm$gene, da$scc_6$bulkNorm$gene) %>% unique()

gr <- vectorization.core(cell.list = cell.list, ao = aoi, mat = assay(peakmat), mat.name = "PeakMatrix",
                         deseq2.norm = F)
atac.mat <- mcols(gr)[, grepl("\\.\\.", colnames(mcols(gr)))] %>% as.matrix()
rownames(atac.mat) <- utilsFanc::gr.get.loci(gr)

atac.mat <- atac.mat[peaks, ]
atac.col.data <- data.frame(sample = colnames(atac.mat)) %>% 
  mutate(type = stringr::str_extract(sample, "tumor|ctrl"))

atac.dds <- DESeqDataSetFromMatrix(countData = atac.mat, colData = atac.col.data, design = ~type)
atac.dds <- estimateSizeFactors(atac.dds)
atac.norm <- counts(atac.dds, normalize = T)

atac.scaled <- atac.norm %>% t() %>% scale() %>% t()

plot.mat.auto.clustering(mat = atac.scaled[scc_6.upDAR, ], cluster_rows = T, cluster_columns = F, 
                         height = 2000, k.m = 3,
                         show_row_names = T, plot.out = paste0(wd, "/dar_hm_clustering/dar_hm.png"))
saveRDS(atac.dds, paste0(wd, "/scc1_vs_6_dds.Rds"))
dar.mat <- atac.norm[scc_6.upDAR,]
head(dar.mat)
dar.mean <- data.frame(scc1 = rowMeans(dar.mat[, 1:4]), scc6 = rowMeans(dar.mat[, 5:8]))

module <- list(
  scc1hi = rownames(dar.mean)[dar.mean$scc1 > dar.mean$scc6],
  scc6hi = rownames(dar.mean)[dar.mean$scc1 < dar.mean$scc6]
)

plot.mat.auto.clustering(mat = atac.scaled[unlist(module), ], cluster_rows = F, cluster_columns = F, 
                         height = 2000, 
                         show_row_names = F,
                         plot.out = paste0(wd, "/dar_hm_clustering/dar_modules_hm.png"))
saveRDS(module, paste0(wd, "/scc1_vs_6_modules.Rds"))


regions <- lapply(module, utilsFanc::loci.2.gr)
aoi <- addPeakAnnotations(aoi, regions = regions,
                          name = "DCupDARsplit", force = T)

aoi <- addDeviationsMatrix(
  ArchRProj = aoi, 
  peakAnnotation = "DCupDARsplit",
  force = TRUE
)

chromvar.mat <- ArchR::getMatrixFromProject(aoi, useMatrix = "DCupDARsplitMatrix")
assays(chromvar.mat) # deviations z

assayNames(chromvar.mat) <- c("d", "z")
soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")
soi <- seurat.add.archr.matrix(so = soi, mat.se = chromvar.mat, assay.root.name = "DCupDAR", overwrite = F)

identical(soi@assays$DCupDAR_d@data[1,], assays(chromvar.mat)[['d']][1, colnames(soi)]) # TRUE

pl <- lapply(c("d", "z"), function(type) {
  pl <- plot.panel.list(soi[, soi$sample %in% samples], 
                        panel.list = c("scc1hi", "scc6hi"), order = F, 
                        assay = paste0("DCupDAR_", type), invisible = T, 
                        plot.out = paste0(wd, "/scc1_vs_scc6_modules/scc_6_upDAR_chromVar_", type, "_all_samples.pdf"),
                        split.by = "sample", split.order = samples, use.split.as.title = T,
                        raster = T, label = F, publication = T, ymax = 2, ymin = 0.1,
                        sub.width = 1.5, sub.height = 1.5)
})

soi$type.pub <- soi$type %>% stringr::str_to_title()

pl <- lapply(c("d", "z"), function(type) {
  pl <- plot.panel.list(soi[, soi$sample %in% samples], 
                        panel.list = c("scc1hi", "scc6hi"), order = T, 
                        assay = paste0("DCupDAR_", type), invisible = T, 
                        plot.out = paste0(wd, "/scc1_vs_scc6_modules/scc_6_upDAR_chromVar_", type, "_combineReplicate.pdf"),
                        split.by = "type.pub", use.split.as.title = T, split.order = c("Ctrl", "Tumor"),
                        raster = T, label = F, publication = T, ymax = 2, ymin = 0.1,
                        italic.title = F,
                        sub.width = 1.2, sub.height = 1.2)
})

#####
