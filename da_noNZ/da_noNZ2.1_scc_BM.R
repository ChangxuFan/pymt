source("objects.R")

samples <- SAMPLES.BM
peakmat <- readRDS("sync_all/sth/atac/merge_v1/peakmat_full.Rds")
peakmat <- assay(peakmat)

da.dir <- "da_noNZ/sth/da_noNZ2.1_scc_BM"
dir.create(da.dir, showWarnings = F, recursive = T)
aoi <- readRDS("sync_all/sth/atac/merge_v1/ao.Rds")

clusters <- as.character(c(0:7, 9:11, 14, 16))

da <- ao.2.bulk.list(ao = aoi, peakmat = peakmat, pseudo = F,
                     quantile.norm = F, 
                     filter.nz = F, filter.size = c(6, 0.999), sequential.filter = T,
                     independentFiltering = F,
                     cell.prep.dir = paste0(da.dir, "/cell_prep/"),
                     coldata.columns = c("batch", "sample", "type"),
                     cluster.ident = "scc", clusters = clusters,
                     group.ident = "sample", groups = samples,
                     threads = 1, design.formula = ~ batch + type,
                     contrast = c("type", "tumor", "ctrl"),
                     work.dir = da.dir, do.plot = F)

da <- deseq2.summary(da, save = T, gene.out.dir = paste0(da.dir, "/summary/"),
                     log2fc.cutoff = 1)

deseq2.xyplot(pbl = da, comp.vec = "tumor:ctrl", comp.meta = "type",
              transformation = function(x) log2(x + 1),
              plot.dir = paste0(da.dir, "/plot_xy/"),
              root.name = "log2fc1",
              device = "png", threads = 1)


da <- deseq2.summary(da, log2fc.cutoff = 0.5)

deseq2.xyplot(pbl = da, comp.vec = "tumor:ctrl", comp.meta = "type",
              transformation = function(x) log2(x + 1),
              plot.dir = paste0(da.dir, "/plot_xy/"),
              root.name = "log2fc0.5",
              device = "png", threads = 1)

da <- deseq2.summary(da, log2fc.cutoff = 0.5, padj.cutoff = 0.1)

deseq2.xyplot(pbl = da, comp.vec = "tumor:ctrl", comp.meta = "type",
              transformation = function(x) log2(x + 1),
              plot.dir = paste0(da.dir, "/plot_xy/"),
              root.name = "log2fc0.5padj0.1",
              device = "png", threads = 1)


samples.bg <- samples[c(1,3)]
scatter.df <- data.frame(x = samples[c(2,4)], y = samples[c(1,3)])

a2bl.homer.pipe(a2bl = da[paste0("scc_", c(3,6,10))], work.dir = paste0(da.dir, "/homer/"),
                slot = "res.exp", rank.key = "padj", trends = c("up"), 
                desc = F, top.n.vec = c(1000), n.bg.each = 25, no.replace = T,
                samples.for.bg = samples.bg, 
                scatter.xy.df = scatter.df, 
                threads.a2bl = 3, threads.titrate = 1, threads.homer = 4, 
                genome = "mm10", size = "given", denovo = F, homer.run = T)

da <- readRDS(paste0(da.dir, "/ao.bulk.list.Rds"))
da <- deseq2.summary(da, log2fc.cutoff = 0.5, padj.cutoff = 0.1, 
                     gene.out.dir = paste0(da.dir, "/summary_log2fc0.5padj0.1/"))
