source("objects.R")
wd <- "schuettpelz2014/sth/sz1.1_geo2R/"
dir.create(wd, showWarnings = F, recursive = T)

GSE="GSE55095"
sample.map <- read.table("schuettpelz2014/sample.map.tsv", header = T, sep = "\t")
marray.mat <- GEO.microarray.2.de(GSE = GSE, sample.map = sample.map, collapse.method = "mean")
marray.mat <- marray.mat[!grepl("Gm.+$", rownames(marray.mat)),]
nrow(marray.mat) # 19493

sample.info <- sample.map %>% 
  mutate(timepoint = stringr::str_extract(sample, "36h|7d"),
         type = stringr::str_extract(sample, "GCSF|PBS"),
         sample2 = sub("_36h_|_7d_", "_", sample))
sample.info <- sample.info %>% filter(!grepl("_4", sample))

#           GSM     sample timepoint type sample2
# 1  GSM1329461  PBS_36h_1       36h  PBS   PBS_1
# 2  GSM1329462  PBS_36h_2       36h  PBS   PBS_2
# 3  GSM1329463  PBS_36h_3       36h  PBS   PBS_3
# 4  GSM1329465 GCSF_36h_1       36h GCSF  GCSF_1
# 5  GSM1329466 GCSF_36h_2       36h GCSF  GCSF_2
# 6  GSM1329467 GCSF_36h_3       36h GCSF  GCSF_3
# 7  GSM1329469   PBS_7d_1        7d  PBS   PBS_1
# 8  GSM1329470   PBS_7d_2        7d  PBS   PBS_2
# 9  GSM1329471   PBS_7d_3        7d  PBS   PBS_3
# 10 GSM1329472  GCSF_7d_1        7d GCSF  GCSF_1
# 11 GSM1329473  GCSF_7d_2        7d GCSF  GCSF_2
# 12 GSM1329474  GCSF_7d_3        7d GCSF  GCSF_3

write.table(sample.info, "schuettpelz2014/sample_info.tsv", col.names = T, row.names = F,
            sep = "\t", quote = F)
saveRDS(marray.mat, paste0(wd, "/marray.mat.full.Rds"))
marray.mat <- marray.mat[, !grepl("_4", colnames(marray.mat))]
saveRDS(marray.mat, paste0(wd, "/marray.mat.Rds"))
