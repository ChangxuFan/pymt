source("objects.R")
wd <- "denardo2018/sth/denardo1.1_geo2R/"
dir.create(wd, showWarnings = F, recursive = T)

sample.info <- read.table("denardo2018/sample_info.tsv", header = T)
marray.mat <- GEO.microarray.2.de(GSE = "GSE99467", sample.map = sample.info, collapse.method = "mean")

marray.mat <- marray.mat[!grepl("Gm.+$", rownames(marray.mat)),]
nrow(marray.mat) # 16205
sample.info
#           GSM      sample celltype sample2
# 1  GSM2644432 MDP_Normal2      MDP Normal2
# 2  GSM2644433 MDP_Normal3      MDP Normal3
# 3  GSM2644434 MDP_Normal4      MDP Normal4
# 4  GSM2644435 MDP_Normal5      MDP Normal5
# 5  GSM2644436   MDP_PyMT2      MDP   PyMT2
# 6  GSM2644437   MDP_PyMT3      MDP   PyMT3
# 7  GSM2644438   MDP_PyMT4      MDP   PyMT4
# 8  GSM2644439   MDP_PyMT5      MDP   PyMT5
# 9  GSM2644440 CDP_Normal2      CDP Normal2
# 10 GSM2644441 CDP_Normal3      CDP Normal3
# 11 GSM2644442 CDP_Normal4      CDP Normal4
# 12 GSM2644443 CDP_Normal5      CDP Normal5
# 13 GSM2644444   CDP_PyMT2      CDP   PyMT2
# 14 GSM2644445   CDP_PyMT3      CDP   PyMT3
# 15 GSM2644446   CDP_PyMT4      CDP   PyMT4
# 16 GSM2644447   CDP_PyMT5      CDP   PyMT5

saveRDS(marray.mat, paste0(wd, "/marray.mat.Rds"))

sample.info$type <- stringr::str_extract(sample.info$sample, "Normal|PyMT")
de <- diffbind.2.raw.a2bl(sample.info = sample.info, skip.deseq2 = T,
                          mat = marray.mat, use.genewise.disperisons = T,
                          cluster.ident = "celltype",sample.col = "sample", rename.sample.to = "sample2",
                          filter.nz = F, filter.size = c(0, 1), sequential.filter = T,
                          quantile.norm = F,
                          independentFiltering = F, 
                          single.sample = F,
                          coldata.columns = c("type"), 
                          design.formula = ~ type, 
                          contrast = c("type", "PyMT", "Normal"), 
                          threads = 1,
                          work.dir = wd)

de <- microarray.quantile.norm(de)
de <- de["CDP"]

de$CDP$res <- microarray.read.geo2r("reference/denardo2018/sth/denardo1.1_geo2R/GSE99467.top.table.CDP_PyMT_vs_Normal.tsv")

de <- deseq2.summary(de, padj.cutoff = 0.1)
saveRDS(de, paste0(wd, "/bulk.list.Rds"))

