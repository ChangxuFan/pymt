source("objects.R")

samples <- SAMPLES
peakmat <- readRDS("sync_all/sth/atac/merge_v1/peakmat_full.Rds")
peakmat <- assay(peakmat)

da.dir <- "da_noNZ/sth/da_noNZ2.2_scc_BMSP"
dir.create(da.dir, showWarnings = F, recursive = T)
aoi <- readRDS("sync_all/sth/atac/merge_v1/ao.Rds")

clusters <- as.character(c(0:7, 9:11, 14, 16))

da <- ao.2.bulk.list(ao = aoi, peakmat = peakmat, pseudo = F,
                     quantile.norm = T, ## key line
                     # filter.samples = filter.samples,
                     filter.nz = F, filter.size = c(6, 0.999), sequential.filter = T,
                     independentFiltering = F,
                     cell.prep.dir = paste0(da.dir, "/cell_prep/"),
                     coldata.columns = c("batch", "sample", "type", "tissue"),
                     cluster.ident = "scc", clusters = clusters,
                     group.ident = "sample", groups = samples,
                     threads = 1, design.formula = ~ tissue,
                     contrast = c("tissue", "SP", "BM"),
                     work.dir = da.dir, do.plot = F)

da <- deseq2.summary(da, save = T, gene.out.dir = paste0(da.dir, "/summary/"),
                     log2fc.cutoff = 1)

deseq2.xyplot(pbl = da, comp.vec = "SP:BM", comp.meta = "tissue",
              transformation = function(x) log2(x + 1),
              plot.dir = paste0(da.dir, "/plot_xy/"),
              root.name = "log2fc1",
              device = "png", threads = 1)

da <- deseq2.summary(da, save = T, gene.out.dir = paste0(da.dir, "/summary_log2fc0.5_padj0.1/"),
                     log2fc.cutoff = 0.5, padj.cutoff = 0.1)

deseq2.xyplot(pbl = da, comp.vec = "SP:BM", comp.meta = "tissue",
              transformation = function(x) log2(x + 1),
              plot.dir = paste0(da.dir, "/plot_xy/"),
              root.name = "log2fc0.5_padj0.1",
              device = "png", threads = 1)

samples.bg <- SAMPLES.BM
scatter.df <- data.frame(x = SAMPLES.BM, y = SAMPLES.SP)

a2bl.homer.pipe(a2bl = da[paste0("scc_", c(0, 1, 2, 9, 11, 14))],
                work.dir = paste0(da.dir, "/homer/"),
                slot = "res.exp", rank.key = "padj", trends = c("down"), 
                desc = F, top.n.vec = c(1000), n.bg.each = 25, no.replace = T,
                samples.for.bg = samples.bg, 
                scatter.xy.df = scatter.df, 
                threads.a2bl = 3, threads.titrate = 1, threads.homer = 4, 
                genome = "mm10", size = "given", denovo = F, homer.run = T)

