source("objects.R")
de.dir <- "MR2019/sth/MR1.1_de/"
dir.create(de.dir)

expr <- read.table("MR2019/sth/GEO/GSE112769_Ly_My_HSCs.genesfpkm.txt.gz", sep = "\t", header = T)

sum(duplicated(expr$Gene.name)) # only 2!

# we pretend that fpkm is raw. Since it's mostly a linear scaling anyway...
expr <- expr[!duplicated(expr$Gene.name),]
# 14k genes, seems reasonable. 
rawmat <- expr[, c(3:ncol(expr))] %>% as.matrix()
rownames(rawmat) <- expr$Gene.name

colSums(rawmat) # about 200k each sample. We need to scale this up a little bit. 
# the actual sequencing depth was 30 million.

de.my <- readRDS("~/hmtp/scAR/bulkVerify/sth/de/de1.1_BM/bulk.list.Rds")
de.my$KSL$bulk.mat %>% colSums() # about 20million

# Let's keep it conservative and scale it up to 10 million.
scale.factor <- 10 * 10^6 / 200000
# 50

rawmat <- floor(rawmat * 50)

cutoff <- 5
n.pass <- 3
bMat <- rawmat > cutoff
print(sum(rowSums(bMat) >= n.pass)) #11949
filtered.mat <- rawmat[rowSums(bMat) >= n.pass,]
filtered.mat <- filtered.mat[!grepl("^Gm", rownames(filtered.mat)), ]
nrow(filtered.mat) #11936

samples <- colnames(filtered.mat)

sample.info <- data.frame(
  sample = samples,
  bias = stringr::str_extract(samples, "Ly|My"),
  age = stringr::str_extract(samples, "You|Old"),
  rep = paste0("rep", stringr::str_extract(samples, "\\d")),
  celltype = "HSC"
)

# > sample.info
#    sample bias age  rep
# 1  LyYou1   Ly You rep1
# 2  LyYou2   Ly You rep2
# 3  LyYou3   Ly You rep3
# 4  LyOld1   Ly Old rep1
# 5  LyOld2   Ly Old rep2
# 6  LyOld3   Ly Old rep3
# 7  MyYou1   My You rep1
# 8  MyYou2   My You rep2
# 9  MyYou3   My You rep3
# 10 MyOld1   My Old rep1
# 11 MyOld2   My Old rep2
# 12 MyOld3   My Old rep3

# we get rid of LyYou1, because it is very weird on a PCA plot
samples <- samples[samples != "LyYou1"]
sample.info <- sample.info %>% filter(sample %in% samples)

write.table(sample.info, paste0(de.dir, "/sample_info.tsv"), sep = "\t", 
            row.names = F, col.names = T, quote = F)

de <- diffbind.2.raw.a2bl(sample.info = sample.info,
                          mat = filtered.mat,
                          cluster.ident = "celltype",sample.col = "sample", 
                          filter.nz = F, filter.size = c(0, 0.999), sequential.filter = T,
                          independentFiltering = T,
                          coldata.columns = c("bias", "age", "rep"), 
                          design.formula = ~ age + bias, 
                          contrast = c("bias", "My", "Ly"), 
                          threads = 1,
                          work.dir = de.dir)

# We didn't see a strong dependece on replicate on the PCA plot. Not regressing on it.

de <- deseq2.summary(de, gene.out.dir = paste0(de.dir, "/summary"))

deseq2.xyplot(pbl = de, comp.vec = "My:Ly", comp.meta = "bias",
              transformation = function(x) return(log2(x + 1)),
              plot.dir = paste0(de.dir, "/plot_xy/"), root.name = "log2",
              device = "png", theme.text.size = 15, 
              plot.each.comp = F)

deseq2.hm(pbl = de, use.bubble = F, 
          plot.dir = paste0(de.dir, "/plot_hm/"), dense.hm = T)

