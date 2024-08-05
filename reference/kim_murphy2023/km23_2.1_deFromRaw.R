source("objects.R")
wd <- work.dir <- de.dir <- "kim_murphy2023/sth/km23_2.1_deFromRaw/"
dir.create(wd, showWarnings = F, recursive = T)

samples.ori <- c("cdp_1", "cdp_2", "cdp_flt3l_1", "cdp_flt3l_2", "cdp_flt3l_3",
                 "cdp_flt3l_il_6_1", "cdp_flt3l_il_6_2", "cdp_flt3l_il_6_3")
samples <- samples.ori %>% sub("cdp", "CDP", .) %>% sub("flt3l", "Flt3l", .) %>% sub("il_6", "Il6", .)
raw.df <- lapply(samples.ori, function(sample) {
  dir <- "~/hmtp/scAR/reference/kim_murphy2023/sth/GEO/"
  file <- Sys.glob(paste0(dir, "/*", sample, "*"))
  # these files were downloaded from GSE215751
  if (length(file) != 1) stop("length(file) != 1")
  df <- read.table(file, header = T)
  df <- df[, c("Feature", "Count")]
  colnames(df) <- c("Feature_ID", sample)
  return(df)
}) %>% Reduce(full_join, .)

any(is.na(raw.df)) # False

gene.map <- read.table("kim_murphy2023/sth/GEO/GSE215751_Merged_differential_expression_results.tsv", header = T,
                       quote = "", sep = "\t")

gene.map <- gene.map[, c("Feature_ID", "external_gene_name")]
colnames(gene.map) <- c("Feature_ID", "gene")
nrow(gene.map) # 16750. 

raw.df <- left_join(gene.map, raw.df)
any(is.na(raw.df)) # False
raw.df$Feature_ID <- NULL
sum(duplicated(raw.df$gene)) # 59

raw.df <- raw.df[!duplicated(raw.df$gene),]

rownames(raw.df) <- raw.df$gene
raw.df$gene <- NULL
filtered.mat <- raw.df %>% as.matrix()


identical(colnames(filtered.mat), samples.ori) # TRUE
colnames(filtered.mat) <- samples

sample.info <- data.frame(
  sample = samples,
  celltype = "CDP",
  type = stringr::str_extract(samples, "Flt3l(_Il6)*")
)
sample.info$type[is.na(sample.info$type)] <- "untreated"


de <- diffbind.2.raw.a2bl(sample.info = sample.info,
                          mat = filtered.mat,
                          cluster.ident = "celltype",sample.col = "sample",
                          filter.nz = F, filter.size = c(0, 0.999), sequential.filter = T,
                          independentFiltering = T,
                          coldata.columns = c("type", "sample"), 
                          design.formula = ~ type, 
                          contrast = c("type", "Flt3l_Il6", "Flt3l"), 
                          threads = 1,
                          work.dir = work.dir)

de <- deseq2.summary(de, save = T, gene.out.dir = paste0(wd, "/summary/"))

deseq2.xyplot(pbl = de, comp.vec = "Flt3l_Il6:Flt3l", comp.meta = "type",
              transformation = function(x) return(log2(x + 1)),
              add.label = F, 
              plot.dir = paste0(de.dir, "/plot_xy/"), root.name = "log2",
              device = "png", 
              plot.each.comp = F)

saveRDS(de, paste0(wd, "/bulk.list.Rds"))
