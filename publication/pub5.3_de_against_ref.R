source("objects.R")
wd <- "publication/sth/pub5.4_de_against_ref"
dir.create(wd)

de <- readRDS("de_noNZ/sth/de_noNZ2.2_scc_BMSP/bulk.list.Rds")
de <- deseq2.summary(de)

des <- list(fanc = de['scc_1'], mr = readRDS("~/hmtp/scAR/reference/MR2019/sth/MR1.1_de/bulk.list.Rds"))
des$mr <- deseq2.summary(des$mr)

des$fanc$scc_1$fgsea.meta <- list(
  gene.set.name.display.up = "Gene Set: upDEGs in HSPC[1]\n(SP > BM)",
  gene.set.name.display.down = "Gene Set: downDEGs in HSPC[1]\n(SP < BM)",
  rank.name.display = "Genes in HSPC[1]:",
  up.name.display = "SP > BM",
  down.name.display = "SP < BM"
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
