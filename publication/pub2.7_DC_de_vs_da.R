rm(list = ls())
source("objects.R")
wd <- "publication/sth/pub2.7_de_intersect_da"
dir.create(wd)

da <- readRDS("da_noNZ/sth/da_noNZ2.1_scc_BM//ao.bulk.list.Rds")
de <- readRDS("de_noNZ/sth/de_noNZ2.1_scc_BM//bulk.list.Rds")
gtf <- rtracklayer::import("~/genomes/mm10/gencode/gencode.default.gtf")

clusters <- paste0("scc_", c(6))
de <- de[clusters]
da <- da[clusters]

de <- deseq2.summary(de, gene.out.dir = paste0(wd, "/de_summary/"))

gene.sets <- list(
  upDEG = de$scc_6$summary$up.genes,
  core.genes = readRDS("publication/sth/pub2.4_DC_gene_module/scc6_up_genes_post_filter.Rds"),
  myelo.genes = readRDS("publication/sth/pub2.4_DC_gene_module/scc1_vs_6.Rds")$scc6hi
)

de$scc_6$summary <- NULL
samples <- INFO.filtered$sample[1:4]
samples.for.bg <- samples[c(1,3)]
scatter.xy.df <- data.frame(x = samples[c(1,3)], y = samples[c(2,4)])
n.fold.bg <- 25

des <- lapply(names(gene.sets), function(gs.name) {
  wd <- paste0(wd, "/", gs.name)
  gs <- gene.sets[[gs.name]]
  print(gs)
  de$scc_6[[gs.name]] <- list(up.genes = gs, n.up = length(gs))
  de <- de.add.promoter(de = de, gtf.gr = gtf, slot = gs.name)

  de <- de.add.bg(de = de, slot = gs.name, n.bg.each = n.fold.bg,
                  samples.for.bg = samples.for.bg, scatter.xy.df = scatter.xy.df, split.into.folds = T,
                  work.dir = paste0(wd, "/bg/"),
                  no.replace = T)

  de <- de.bg.expand(de = de, n.folds = n.fold.bg, slot = gs.name)
  for (i in 1:n.fold.bg) {
    de <- de.add.promoter(de, gtf.gr = gtf, slot = paste0("fold", i))
  }
  return(de)
})
names(des) <- names(gene.sets)

lapply(names(gene.sets), function(gs.name) {
  de <- des[[gs.name]]
  wd <- paste0(wd, "/", gs.name)
  lapply(c(10^4, 10^5), function(dist) {
    de.da.intersect.titrate(de = de, da = da, slot = gs.name, mode = "promoter", genome.size = 2.7 * 10^9,
                            padj.cutoffs = c(0.1, 0.2), log2FC.cutoffs = c(0.25, 0.5),
                            promoter.ext.1side = dist, directions = c("up"),
                            use.background = T, n.folds.bg = n.fold.bg,
                            out.dir = paste0(wd, "/promoter_distal_", utilsFanc::so.formatter(dist), "/"),
                            print.details = T)
    return()
  })
})

tsvs <- Sys.glob(paste0(wd, "/*/*/promoter_distal_*.tsv"))

lapply(tsvs, function(tsv) {
  de.da.intersect.titrate.plot.venn(in.tsv = tsv)
  return()
})


