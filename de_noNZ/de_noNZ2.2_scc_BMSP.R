source("objects.R")
de.dir <- "de_noNZ/sth/de_noNZ2.2_scc_BMSP/"
dir.create(de.dir)
samples <- SAMPLES

soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")

#####
# key: remove mt- genes
soi <- soi[!grepl("^mt\\-", rownames(soi)), ]

clusters <- as.character(c(0:7, 9:11, 14, 16))

de <- bulk.list(so = soi, assay = "RNA", slot = "counts",
                group.by = "sample", groups = samples, 
                sample.order = samples, coldata.columns = c("type", "batch", "tissue"), 
                threads = 1, design.formula = ~tissue, contrast = c("tissue", "SP", "BM"), 
                cluster.ident = "scc", clusters = clusters, 
                filter.nz = F, filter.size = c(1, 0.999), 
                sequential.filter = T, single.sample = F, 
                do.plot = F, save.rds = T,  stop.on.error = T,
                work.dir = de.dir, plot.dir = paste0(de.dir, "/plots/"))

de <- deseq2.summary(de, save = T, gene.out.dir = paste0(de.dir, "/summary/"))

deseq2.xyplot(pbl = de, comp.vec = "SP:BM", comp.meta = "tissue",
              transformation = function(x) log2(x + 1),
              plot.dir = paste0(de.dir, "/plot_xy/"), 
              root.name = "log2fc1",
              device = "png", theme.text.size = 12)


deseq2.plot.panel(pbl = de['scc_1'], so = soi, max.plots = 100, summary.slot = "summary", 
                  order = F, assay = "RNA",
                  cluster.ident = "scc", page.limit = 16,
                  n.col = 4, plot.type = "up.genes", sample.order = samples,
                  panel.label = T, threads = 1, 
                  plot.dir = paste0(de.dir, "/plot_panel"))


