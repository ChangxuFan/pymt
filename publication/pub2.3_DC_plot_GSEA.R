source("objects.R")

wd <- "publication/sth/pub2.3_DC_GSEA/"
dir.create(wd, showWarnings = F, recursive = T)

de <- readRDS("de_noNZ/sth/de_noNZ2.1_scc_BM/bulk.list.Rds")
de <- deseq2.summary(de)

de.denardo <- readRDS("reference/denardo2018/sth/denardo1.1_geo2R/bulk.list.Rds")

de.murphy <- readRDS("reference/kim_murphy2023/sth/km23_2.1_deFromRaw/bulk.list.Rds")
de.murphy <- deseq2.summary(de.murphy)

de$scc_6$fgsea.meta <- list(gene.set.name.display = "DEGs in DC[6]",
                          rank.name.display = "Genes in DC[6]: ",
                          up.name.display = "Tumor > Ctrl",
                          down.name.display = "Tumor < Ctrl")
de.murphy$CDP$fgsea.meta <- list(gene.set.name.display = "DEGs in Kim 2023",
                                 rank.name.display = "Genes in CDP (Kim 2023):",
                                 up.name.display = "IL-6 > Ctrl",
                                 down.name.display = "IL-6 < Ctrl")

de.denardo$CDP$fgsea.meta <- list(gene.set.name.display = "DEGs in Meyer 2018",
                                  rank.name.display = "Genes in CDP (Meyer 2018):",
                                  up.name.display = "Tumor > Ctrl",
                                  down.name.display = "Tumor < Ctrl")


de.reciprocal.gsea(des = list(fanc = de['scc_6'], murphy = de.murphy, denardo = de.denardo), 
                   slot = "summary", out.dir = paste0(wd, "/all/"), 
                   run = T,force = T, n.par.de = 1,
                   fgsea.plot = T)

