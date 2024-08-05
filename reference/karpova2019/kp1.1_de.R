source("objects.R")
wd <- "karpova2019/sth/kp1.1_geo2R/"
dir.create(wd)

sample.info <- read.table("karpova2019/sample_info.tsv", header = T)
marray.mat <- GEO.microarray.2.de(GSE = "GSE123505", sample.map = sample.info, collapse.method = "mean")

marray.mat <- marray.mat[!grepl("Gm.+$", rownames(marray.mat)),]
nrow(marray.mat) # 15673

saveRDS(marray.mat, paste0(wd, "/marray.mat.Rds"))
sample.info$celltype <- "LSK"

de <- diffbind.2.raw.a2bl(sample.info = sample.info, skip.deseq2 = T,
                          mat = marray.mat, use.genewise.disperisons = T,
                          cluster.ident = "celltype",sample.col = "sample", 
                          filter.nz = F, filter.size = c(0, 1), sequential.filter = T,
                          quantile.norm = F,
                          independentFiltering = F, 
                          single.sample = F,
                          coldata.columns = c("treatment"), 
                          design.formula = ~ treatment, 
                          contrast = c("treatment", "GCSF", "ctrl"), 
                          threads = 1,
                          work.dir = wd)

treatments <- sample.info$treatment %>% unique()
# [1] "ctrl"    "Grob"    "GCSF"    "AMD3100"
comps <- paste0(treatments[2:4], ":", treatments[1])

for (comp in comps) {
  treatment <- sub(":.+$", "", comp)
  print(paste0("Processing: ", treatment))
  if(treatment == "GCSF") {
    padj.cutoff <- 0.01
  } else {
    padj.cutoff <- 0.1
  }
  stats.file <- paste0(wd, "/GSE123505.top.table.", sub(":", "_vs_", comp), ".tsv")
  if (!file.exists(stats.file)) {
    stop(paste0(stats.file, " doesn't exist"))
  }
  res.slot <- paste0("res_", treatment) 
  summ.slot <- paste0("summary_", treatment)
  de$LSK[[res.slot]] <- microarray.read.geo2r(stats.file)
  de <- deseq2.summary(de, padj.cutoff = padj.cutoff, summary.slot = summ.slot, res.slot = res.slot)
  
  deseq2.xyplot(pbl = de, comp.vec = comp, comp.meta = "treatment", 
                pbl.slot = "bulkNorm", summ.slot = summ.slot,
                highlight.ptsize = 0.2,
                plot.dir = paste0(wd, "/plot_xy/"), root.name = paste0("GEO2R_DEG_", comp),
                device = "png", theme.text.size = 15, 
                plot.each.comp = F)
}

saveRDS(de, paste0(wd, "/bulk.list.Rds"))
