source("objects.R")
wd <- "publication/sth/pub4.2.2_HSPC_de_ref"
dir.create(wd)

de <- readRDS("de_noNZ/sth/de_noNZ2.1_scc_BM/bulk.list.Rds") %>% deseq2.summary()

des <- list(fanc = de['scc_1'], mr = readRDS("~/hmtp/scAR/reference/MR2019/sth/MR1.1_de/bulk.list.Rds"))
des$mr <- deseq2.summary(des$mr)

des$fanc$scc_1$fgsea.meta <- list(
  gene.set.name.display.up = "Gene Set: upDEGs in HSPC[1]\n(Tumor > Ctrl)",
  gene.set.name.display.down = "Gene Set: downDEGs in HSPC[1]\n(Tumor < Ctrl)",
  rank.name.display = "Genes in HSPC[1]:",
  up.name.display = "Tumor > Ctrl",
  down.name.display = "Tumor < Ctrl"
)

des$mr$HSC$fgsea.meta <- list(
  gene.set.name.display.up = paste0("Gene Set: MR 2019\nMy-HSC > Ly-HSC"),
  gene.set.name.display.down = paste0("Gene Set: MR 2019\nMy-HSC < Ly-HSC"),
  rank.name.display = paste0("Genes in MR 2019:"),
  up.name.display = paste0("My-HSC > Ly-HSC"),
  down.name.display = paste0("My-HSC < Ly-HSC")
)

de.reciprocal.gsea(des = des, 
                   slot = "summary", out.dir = paste0(wd, "/gsea_against_reference/"), 
                   run = T, force = T, n.par.de = 1,
                   fgsea.plot = T, fgsea.width = 1.6)

soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")

module <- des$mr$HSC$summary[c("up.genes", "down.genes")]

soi <- Seurat::AddModuleScore(object = soi, features = module, name = "MR", 
                              assay = "RNA")

samples <- SAMPLES.BM
lapply(c("eachSample"), function(split.type) {
  if (split.type == "combineReplicate") {
    split.by <- "type.pub"
    split.order <- c("Ctrl", "Tumor")
  } else {
    split.by <- "sample"
    split.order <- samples
  }
  
  plot.panel.list(soi, panel.list = c("MR1", "MR2"),
                  order = T, assay = "RNA", invisible = T, 
                  plot.out = paste0(wd, "/seurat_module_MR_DEG_", split.type, ".pdf"),
                  sample = samples,
                  split.by = split.by, split.order = split.order, use.split.as.title = F,
                  raster = T, label = F, publication = T, ymax = NULL, ymin = 0, italic.title = F,
                  sub.width = 1.2, sub.height = 1.2)
  
})
