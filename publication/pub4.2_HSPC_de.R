source("objects.R")
wd <- "publication/sth/pub4.2_HSPC_de/"
dir.create(wd)
de <- readRDS("de_noNZ/sth/de_noNZ2.1_scc_BM/bulk.list.Rds")
de <- deseq2.summary(de)

clusters <- paste0("scc_", c(1, 9))
names(clusters) <- clusters

# Figure out what genes we want to label:
upDEGs <- de$scc_1$summary$up.genes
markers <- c("Vldlr", "Gda", "Csf2rb", "Il18rap")

lapply(clusters, function(cluster) {
  deseq2.xyplot(pbl = de[cluster], comp.vec = "tumor:ctrl", comp.meta = "type",
                transformation = function(x) log2(x + 1),
                plot.dir = paste0(wd, "/plot_xy/"), 
                root.name = paste0("xy_DEG", "_", cluster),
                add.label = T, repel.direction = 'x',
                hjust = 0.5, vjust = 0.5, use.repel = T, italic.label = T,
                label.list = markers, text.size = 6 * 0.36,
                device = "pdf", publication = T,
                plot.each.comp = F)
})

deseq2.xyplot(pbl = de['scc_9'], comp.vec = "tumor:ctrl", comp.meta = "type",
              transformation = function(x) log2(x + 1),
              plot.dir = paste0(wd, "/plot_xy/"), 
              root.name = paste0("xy_DEG", "_scc_9_w_scc1DEG"),
              highlight.genes = upDEGs,
              highlight.color = utilsFanc::gg_color_hue(2)[2],
              add.label = T, repel.direction = 'x',
              hjust = 0.5, vjust = 0.5, use.repel = T, italic.label = T,
              label.list = markers, text.size = 6 * 0.36,
              device = "pdf", publication = T,
              plot.each.comp = F)

# also plot for de.sz:
de.sz <- readRDS("~/hmtp/scAR/reference/schuettpelz2014/sth/sz1.1_geo2R/bulk.list.Rds")
de.sz <- microarray.quantile.norm(de.sz)
de.sz$szGCSF36h$res <- microarray.read.geo2r("~/hmtp/scAR/reference/schuettpelz2014/sth/sz1.1_geo2R/GSE55095_36h_GCSF_vs_PBS.top.table.tsv")
de.sz$szGCSF7d$res <- microarray.read.geo2r("~/hmtp/scAR/reference/schuettpelz2014/sth/sz1.1_geo2R/GSE55095_7d_GCSF_vs_PBS.top.table.tsv")
de.sz <- deseq2.summary(de.sz, log2fc.cutoff = 1, use.topn.p = 25)
de.sz <- deseq2.summary(de.sz, log2fc.cutoff = 1, use.topn.p = 50, summary.slot = "summary50")

saveRDS(de.sz, "~/hmtp/scAR/reference/schuettpelz2014/sth/sz1.1_geo2R/bulk.list.Rds")
de.sz <- de.sz['szGCSF7d']
deseq2.xyplot(pbl = de.sz, comp.vec = "GCSF:PBS", comp.meta = "type",
              pbl.slot = "bulkNorm",
              transformation = NULL, 
              plot.dir = paste0(wd, "/plot_xy/"), 
              root.name = "xy_DEG_on_sz",
              highlight.genes = de$scc_1$summary$up.genes, highlight.color = utilsFanc::gg_color_hue(2)[2],
              add.label = T,
              hjust = 0.5, vjust = 0.5, use.repel = T, italic.label = T,
              label.list = markers, text.size = 6 * 0.36,
              device = "pdf", publication = T,
              repel.direction = "x", highlight.ptsize = 0.08,
              plot.each.comp = F)

des <- list(fanc = de['scc_1'], sz = de.sz)
des$fanc$scc_1$fgsea.meta <- list(
  gene.set.name.display.up = "Gene Set: upDEGs in HSPC[1]\n(Tumor > Ctrl)",
  gene.set.name.display.down = "Gene Set: downDEGs in HSPC[1]\n(Tumor < Ctrl)",
  rank.name.display = "Genes in HSPC[1]:",
  up.name.display = "Tumor > Ctrl",
  down.name.display = "Tumor < Ctrl"
)

des$sz$szGCSF7d$fgsea.meta <- list(
  gene.set.name.display.up = paste0("Gene Set: Schuettpelz 2014\nG-CSF > PBS"),
  gene.set.name.display.down = paste0("Gene Set: Schuettpelz 2014\nG-CSF < PBS"),
  rank.name.display = paste0("Genes in Schuettpelz 2014:"),
  up.name.display = paste0("G-CSF > PBS"),
  down.name.display = paste0("G-CSF < PBS")
)

de.reciprocal.gsea(des = des, 
                   slot = "summary", out.dir = paste0(wd, "/gsea_against_reference/"), 
                   run = T, force = T, n.par.de = 1,
                   fgsea.plot = T, fgsea.width = 1.6)


