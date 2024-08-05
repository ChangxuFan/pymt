rm(list = ls())
source("objects.R")
wd <- "publication/sth/pub3.4_de_da_intersect"
dir.create(wd)
da <- readRDS("publication/sth/pub3.3_c20_da/ao.bulk.list.Rds")
de <- readRDS("publication/sth/pub3.2_c20_de/bulk.list.Rds")
gtf <- rtracklayer::import("~/genomes/mm10/gencode/gencode.default.gtf")

de <- deseq2.summary(de, gene.out.dir = paste0(wd, "/de_summary/"))

gene.sets <- list(
  upDEG = de$pd_sc20$summary$up.genes
)

de$pd_sc20$summary <- NULL
ss <- de$pd_sc20$coldata %>% rownames()

samples.for.bg <- ss[1:4]
scatter.xy.df <- data.frame(x = ss[1:4], y = ss[5:8])
n.fold.bg <- 25

des <- lapply(names(gene.sets), function(gs.name) {
  wd <- paste0(wd, "/", gs.name)
  gs <- gene.sets[[gs.name]]
  # print(gs)
  de[[1]][[gs.name]] <- list(up.genes = gs, n.up = length(gs))
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
                            padj.cutoffs = c(0.05, 0.2), log2FC.cutoffs = c(0.25, 1),
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


#>>>>>>>> Plot some examples
# Candidates for examples:
#* Ahnak
#* Lgals3
#* Lmna
#* Emp1: promoter!!

markers <- c("Ahnak", "Emp1", "S100a4", "S100a6")
soi.small$sample.display <- SAMPLE.MAP[soi.small$sample]

# matching # of cells between sc20 and sc1:
cell.list <- get.cell.list(soi.small, n.cells.each = 250, 
                           group.by = "sample.display", split.by = "sc")

plot.panel.list(obj = soi.small, panel.list = markers, 
                cells = unlist(cell.list),
                assay = "RNA", # sample = samples.sp, 
                ident = "sample.display", 
                split.by = "sc.display", add.median = F,
                violin = T, n.col = 3, 
                sub.width = 1, sub.height = 1,
                publication = T, hide.legend = T, italic.title = T, pt.size = 0.01,
                plot.out = paste0(wd, "/violin/violin_browser_genes_match_cellNumber.pdf"),
                violin.adjust = 1)

soi.small$tmp <- paste0(soi.small$sc, "..", soi.small$sample.display)
tmp.levels <- soi.small$tmp %>% unique() %>% sort()
tmp.levels <- c(tmp.levels[5:8], tmp.levels[1:4])
soi.small$tmp <- factor(soi.small$tmp, tmp.levels)

color.map <- c(rep("#117780", 4), rep("#F35F5B", 4))
names(color.map) <- tmp.levels

lapply(markers, function(marker) {
  plot.panel.list(obj = soi.small, panel.list = marker, 
                  cells = unlist(cell.list),
                  assay = "RNA", # sample = samples.sp, 
                  ident = "tmp", 
                  add.median = F,
                  violin = T, n.col = 3, 
                  sub.width = 0.7, sub.height = 0.5, 
                  violin.remove.x = T, violin.color.map = color.map,
                  publication = T, hide.legend = T, italic.title = T, pt.size = 0.01,
                  plot.out = paste0(wd, "/violin/violin_browser_genes_match_cellNumber_", marker,".pdf"),
                  violin.adjust = 1)
  
})


