source("objects.R")
wd <- "publication/sth/pub1.2_DEG_bar_plots/"
dir.create(wd)

de <- readRDS("de_noNZ/sth/de_noNZ2.1_scc_BM//bulk.list.Rds")
# saveRDS(de, "de_noNZ/sth/de_noNZ1.1_tumor_vs_ctrl_sc/bulk.list.Rds")
aliases <- c("Ery[0]", "HSPC[1]", "Ery[2]", "Mono/Mac/DC[3]", 
             "Neu[4]", "Myelo[5]", "DC[6]", "Mega[7]", "HSPC[9]", 
             "Mono/Mac[10]", "Ery[11]", "Mega[14]", "Lympho[16]")

identical(stringr::str_extract(aliases, "\\d+"), stringr::str_extract(names(de), "\\d+"))

names(de) <- aliases
de <- de.sync.name(de)

de <- deseq2.summary(de)

de.deg.plot.bar(de = de, out.file = paste0(wd, "/deg_bar.pdf"), 
                up.name = "# upDEGs (Tumor > Ctrl)",
                down.name = "# downDEGs (Tumor < Ctrl)",
                width = 2.5, height = 2.5)

###
de <- readRDS("de_noNZ/sth/de_noNZ2.2_scc_BMSP///bulk.list.Rds")
aliases <- c("Ery[0]", "HSPC[1]", "Ery[2]", "Mono/Mac/DC[3]", 
             "Neu[4]", "Myelo[5]", "DC[6]", "Mega[7]", "HSPC[9]", 
             "Mono/Mac[10]", "Ery[11]", "Mega[14]", "Lympho[16]")

identical(stringr::str_extract(aliases, "\\d+"), stringr::str_extract(names(de), "\\d+"))

names(de) <- aliases
de <- de.sync.name(de)

de <- deseq2.summary(de)

# notably, we are not plotting Mono/Mac/DCs.
de.deg.plot.bar(de = de[!grepl("DC|Mono", names(de))], out.file = paste0(wd, "/deg_bar_SPBM.pdf"), 
                up.name = "# upDEGs (Spleen > BM)",
                down.name = "# downDEGs (Spleen < BM)",
                width = 2.5, height = 2.2)


######

clusters <- stringr::str_extract(aliases, "\\d+")
da <- readRDS("da_noNZ/sth/da_noNZ2.1_scc_BM//ao.bulk.list.Rds")

da <- da[paste0("scc_", c(0:7, 9:11, 14, 16))]
identical(stringr::str_extract(aliases, "\\d+"), stringr::str_extract(names(da), "\\d+"))


names(da) <- aliases
da <- de.sync.name(da)

da <- deseq2.summary(da)

de.deg.plot.bar(de = da, out.file = paste0(wd, "/dar_bar.pdf"), 
                up.name = "# upDARs (Tumor > Ctrl)",
                down.name = "# downDARs (Tumor < Ctrl)",
                width = 2.5, height = 2.5, palette.fc = "R4.fc2")

####
clusters <- stringr::str_extract(aliases, "\\d+")
da <- readRDS("da_noNZ/sth/da_noNZ2.2_scc_BMSP///ao.bulk.list.Rds")

da <- da[paste0("scc_", c(0:7, 9:11, 14, 16))]
identical(stringr::str_extract(aliases, "\\d+"), stringr::str_extract(names(da), "\\d+"))


names(da) <- aliases
da <- de.sync.name(da)

da <- deseq2.summary(da)

de.deg.plot.bar(de = da[!grepl("DC|Mono", names(da))], out.file = paste0(wd, "/dar_bar_SPBM.pdf"), 
                up.name = "# upDARs (Spleen > BM)",
                down.name = "# downDARs (Spleen < BM)",
                width = 2.5, height = 2.2, palette.fc = "R4.fc2")
