source("objects.R")
da <- readRDS("da_noNZ/sth/da_noNZ2.1_scc_BM/ao.bulk.list.Rds")
#>>>>>>>> we can calculate a correlation value of the log2 transformed 
da <- deseq2.summary(da)
da.dir <- "da_noNZ/sth/da_noNZ2.1_scc_BM/"
wd <- "publication/sth/pub4.3_HSPC_why_no_da"
dir.create(wd)

deseq2.xyplot(pbl = da, comp.vec = "tumor:ctrl", comp.meta = "type",
              transformation = function(x) log2(x + 1),
              plot.dir = paste0(da.dir, "/plot_xy/"), 
              root.name = "log2fc1_wCorr",
              device = "png", theme.text.size = 12, add.corr = T,
              plot.each.comp = F)

#>>>>>>>> Examine the log2fc of peaks around DEGs.
de <- readRDS("de_noNZ/sth/de_noNZ2.1_scc_BM/bulk.list.Rds")

gtf <- rtracklayer::import("~/genomes/mm10/gencode/gencode.default.gtf")

clusters <- paste0("scc_", c("1", "6"))
names(clusters) <- clusters

de <- de[clusters]
da <- da[clusters]

de <- deseq2.summary(de)

de <- de.add.promoter(de = de, gtf.gr = gtf, slot = "summary")

# pndg: peaks near DEG.
ext <- c(10^4, 10^5)
names(ext) <- utilsFanc::so.formatter(ext)

top.n <- de$scc_1$summary$n.up # 16

pndg <- lapply(ext, function(ext) {
  lapply(clusters, function(cluster) {
    genes <- de[[cluster]]$summary$up.genes[1:top.n]
    regions <- de[[cluster]]$summary$promoters$up + ext
    regions <- regions[regions$gene %in% genes]
    
    peakset <- da[[cluster]]$bulkNorm$gene %>% utilsFanc::loci.2.gr()
    o <- findOverlaps(peakset, regions)
    pndg <- peakset[sort(unique(queryHits(o)))] 
    utilsFanc::write.zip.fanc(
      pndg, out.file = paste0(wd, "/pndg/ProExt", utilsFanc::so.formatter(ext), "_", cluster, ".bed"),
      bed.shift = T)
    pndg <- pndg %>% utilsFanc::gr.get.loci()
    return(pndg)
  })
})

saveRDS(pndg, paste0(wd, "/pndg/pndg.Rds"))


res.shrink <- lapply(ext, function(ext) {
  ext.name <- utilsFanc::so.formatter(ext)
  lapply(clusters, function(cluster) {
    peaks <- pndg[[ext.name]][[cluster]]
    res <- da[[cluster]]$res.shrink.ashr.exp
    utilsFanc::check.intersect(peaks, "peaks", res$gene, "res.shrink.ashr.exp$gene")
    res <- res %>% filter(gene %in% peaks)
    rank.plot(df = res, vars = "log2FoldChange",
              outfile = paste0(wd, "/pndg_log2fc_shrink/pndg_log2fc_shrink_rankplot_", ext.name, "_", cluster, ".pdf"), 
              y.limit = c(-3, 3), publication = T, width = 1, height = 1)
    
    res <- res %>% arrange(desc(log2FoldChange))
    return(res)
  })
})

#>>>>>>>> investigate the single DAR in scc 9: are there any DEGs around it?
da <- deseq2.summary(da, log2fc.cutoff = 0.5, padj.cutoff = 0.1)

clusters <- paste0("scc_", c(1, 9))
names(clusters) <- clusters

lapply(clusters, function(cluster) {
  deseq2.xyplot(pbl = da[cluster], comp.vec = "tumor:ctrl", comp.meta = "type",
                transformation = function(x) log2(x + 1),
                plot.dir = paste0(wd, "/plot_xy/"), 
                root.name = paste0("xy_DAR", "_", cluster),
                device = "pdf", publication = T,
                plot.each.comp = F)
})

dar <- da$scc_9$summary$up.genes

dar.100k <- utilsFanc::loci.2.gr(dar) + 1000000

nearby.genes <- subsetByOverlaps(gtf, dar.100k, ignore.strand = T)$gene_name %>% unique()

de$scc_9$res.exp %>% filter(gene %in% nearby.genes, log2FoldChange > 0) %>% View()
# padj are all basically 1.