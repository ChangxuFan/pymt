source("objects.R")
de.dir <- "de_noNZ/sth/de_noNZ2.1_scc_BM/"
dir.create(de.dir)
samples <- INFO.filtered$sample[1:4]
# [1] "BM_ctrl_ba1rep1"  "BM_tumor_ba1rep1" "BM_ctrl_ba3rep1"  "BM_tumor_ba3rep1"

soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")
clusters <- as.character(c(0:7, 9:11, 14, 16))

soi@meta.data[, c("scc", "Clusters", "seurat_clusters")] %>% unique() %>% 
  filter(seurat_clusters %in% "6")

de <- bulk.list(so = soi, assay = "RNA", slot = "counts", 
                group.by = "sample", groups = samples, 
                sample.order = samples, coldata.columns = c("type", "batch"), 
                threads = 1, design.formula = ~ batch + type, contrast = c("type", "tumor", "ctrl"), 
                cluster.ident = "scc", clusters = clusters, 
                filter.nz = F, filter.size = c(1, 0.999), sequential.filter = T, single.sample = F, 
                do.plot = F, save.rds = T,  stop.on.error = T,
                work.dir = de.dir, plot.dir = paste0(de.dir, "/plots/"))

de <- deseq2.summary(de, save = T, gene.out.dir = paste0(de.dir, "/summary/"))

deseq2.xyplot(pbl = de, comp.vec = "tumor:ctrl", comp.meta = "type",
              transformation = function(x) log2(x + 1),
              plot.dir = paste0(de.dir, "/plot_xy/"), 
              root.name = "log2fc1",
              device = "png", theme.text.size = 12)

#>>>>>>>>
options(java.parameters = "- Xmx1024m")
de.write.xls.all.genes(de = de, out.file = "de_noNZ/sth/de_noNZ2.1_scc_BM/xls/all_genes.xlsx")
