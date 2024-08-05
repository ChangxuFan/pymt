rm(list = ls())
source("objects.R")
soi <- readRDS("sync_all/sth/rna/int_rpca20/soi.Rds")
aoi <- readRDS("sync_all/sth/atac/merge_v1/ao.Rds")
wd <- paste0("publication/sth/pub3.2_c20_de")
dir.create(wd)
de.dir <- wd

soi$sample <- SAMPLE.MAP[soi$sample]
soi <- soi[, soi$tissue == "SP" & soi$seurat_clusters %in% c("1", "20")]
soi$sc <- paste0("sc", as.character(soi$seurat_clusters))
soi$sample20 <- paste0(soi$sc, "..", soi$sample)
soi$pd <- "sc20"

aoi <- aoi[aoi$cellNames %in% colnames(soi),]
identical(sort(colnames(soi)), sort(aoi$cellNames)) # TRUE

aoi <- archr.add.seurat.m(ao = aoi, so = soi, metas = c("sc", "sample", "sample20", "pd"), as.factor = F)
saveRDS(aoi, paste0(wd, "/aoi_c20.Rds"))
saveRDS(soi, paste0(wd, "/soi_c20.Rds"))

de <- bulk.list(so = soi, assay = "RNA", slot = "counts", 
                group.by = "sample20", # groups = samples, 
                # sample.order = samples, 
                coldata.columns = c("sample", "sc"), 
                threads = 1, design.formula = ~ sample + sc, 
                contrast = c("sc", "sc20", "sc1"), 
                cluster.ident = "pd", #  clusters = clusters, 
                filter.nz = F, filter.size = c(1, 0.999), sequential.filter = T, single.sample = F, 
                do.plot = F, save.rds = T,  stop.on.error = T,
                work.dir = de.dir, plot.dir = paste0(de.dir, "/plots/"))

de <- deseq2.summary(de, save = T, gene.out.dir = paste0(de.dir, "/summary/"))

deseq2.xyplot(pbl = de, comp.vec = "sc20:sc1", comp.meta = "sc",
              transformation = function(x) log2(x + 1),
              plot.dir = paste0(de.dir, "/plot_xy/"), 
              root.name = "xy_DEG",
              add.label = T,
              hjust = 0.5, vjust = 0.5, use.repel = T, italic.label = T,
              label.list = c("S100a4", "S100a6", "Ahnak"), text.size = 6 * 0.36,
              device = "pdf", publication = T,
              plot.each.comp = T)

# violin plot:
markers <- c("Ahnak", "S100a4", "S100a6")
samples.sp <- SAMPLE.MAP[SAMPLES.SP]
soi$sc.display <- paste0("HSPC[", soi$seurat_clusters, "]")
plot.panel.list(obj = soi, panel.list = markers, assay = "RNA", sample = samples.sp, 
                ident = "sc.display", add.median = F,
                violin = T, n.col = 3, 
                sub.width = 1, sub.height = 1,
                publication = T, hide.legend = T, italic.title = T, pt.size = 0.01,
                plot.out = paste0(de.dir, "/violin/markers_violin.pdf"))

