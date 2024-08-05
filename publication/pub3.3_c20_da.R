rm(list = ls())
source("objects.R")

wd <- paste0("publication/sth/pub3.3_c20_da")
dir.create(wd)
da.dir <- wd

aoi <- readRDS("publication/sth/pub3.2_c20_de/aoi_c20.Rds")

peakmat <- readRDS("sync_all/sth/atac/merge_v1/peakmat_full.Rds")
peakmat <- assay(peakmat)

da <- ao.2.bulk.list(ao = aoi, peakmat = peakmat, pseudo = F,
                     quantile.norm = F, ## key line
                     # filter.samples = filter.samples,
                     filter.nz = F, filter.size = c(6, 0.999), sequential.filter = T,
                     independentFiltering = F,
                     cell.prep.dir = paste0(da.dir, "/cell_prep/"),
                     coldata.columns = c("sample", "sc"),
                     cluster.ident = "pd", # clusters = clusters,
                     group.ident = "sample20", # groups = samples,
                     threads = 1, design.formula = ~ sample + sc,
                     contrast = c("sc", "sc20", "sc1"),
                     work.dir = da.dir, do.plot = F)

da <- deseq2.summary(da, save = T, gene.out.dir = paste0(da.dir, "/summary/"))

samples <- da$pd_sc20$coldata$sample20
samples.bg <- samples[1:4]
# [1] "1..SP_tumor_rep1" "1..SP_tumor_rep2" "1..SP_tumor_rep3" "1..SP_tumor_rep4"

scatter.df <- data.frame(x = samples[1:4], y = samples[5:8])

a2bl.homer.pipe(a2bl = da, work.dir = paste0(da.dir, "/homer/"),
                slot = "res.exp", rank.key = "padj", trends = c("up"), 
                desc = F, top.n.vec = c(1000), n.bg.each = 25, no.replace = T,
                samples.for.bg = samples.bg, 
                scatter.xy.df = scatter.df, 
                threads.a2bl = 1, threads.titrate = 1, threads.homer = 6, 
                genome = "mm10", size = "given", denovo = F, homer.run = T)

# plot a scatter plot:
da <- deseq2.summary(pbl = da)
deseq2.xyplot(pbl = da, comp.vec = "sc20:sc1", comp.meta = "sc",
              transformation = function(x) log2(x + 1),
              plot.dir = paste0(da.dir, "/plot_xy/"), 
              root.name = "xy_DAR",
              device = "pdf", publication = T,
              plot.each.comp = F)

# plot motif enrichment:
homer.plot.table(da.homer.dir = paste0(da.dir, "/homer/"), use.regex = F, col.num = 3,
                 plot.logo = T,
                 motif.map = paste0(da.dir, "/motif_cluster_map/sc20.tsv"), 
                 regex.include = paste0("sc", c(20), "/res.exp_padj_top_1000"),
                 regex.from = "pd_sc", 
                 out.dir = paste0(wd, "/homer/"), logo.width = 300)
