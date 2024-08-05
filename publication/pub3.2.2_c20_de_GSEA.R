source("objects.R")
wd <- "publication/sth/pub3.2_c20_de/"
de <- readRDS(paste0(wd, "/bulk.list.Rds"))
ss <- rownames(de$pd_sc20$coldata)
# [1] "sc1..SP_tumor_rep1"  "sc1..SP_tumor_rep2"  "sc1..SP_tumor_rep3"  "sc1..SP_tumor_rep4"  "sc20..SP_tumor_rep1"
# [6] "sc20..SP_tumor_rep2" "sc20..SP_tumor_rep3" "sc20..SP_tumor_rep4"

rnk <- de.2.rnk(pbl = de, microarray.mode = F, pbl.slot = "res",
                rank.by = "log_p", logp.max = 8, to.upper = T,
                pbl.topn = 12000, pbl.samples = ss, comp = "sc20:sc1", 
                out.dir = paste0(wd, "/rnk/"))

# gsea.fanc(rnk.vec = rnk, gmt.vec = GMT.FILES.HALL.GOBP, 
#           out.dir = paste0(wd, "/gsea/"), thread.rnk = 1, thread.gmt = 12,
#           n.plot = 80, run = T, force = F, parse = F, parse.only = F, zip = F)

gsea.fanc(rnk.vec = rnk, gmt.vec = GMT.FILES["c2.cgp.v7.2"], 
          out.dir = paste0(wd, "/gsea/"), thread.rnk = 1, thread.gmt = 12,
          n.plot = 80, run = T, force = F, parse = F, parse.only = F, zip = F)

# We first plot cell cycle and cell migration:
### Line Mark: went into figures
jobs <- read.table(paste0(wd,"/gsea/plot_curve/plot_jobs.tsv"), header = T, sep = "\t")
lapply(1:nrow(jobs), function(i) {
  print(paste0("Processing job ", i))
  job <- jobs[i, ]
  params <- as.list(job)
  params$plot.out <- paste0(wd, "/gsea/plot_curve/", basename(job$rnk), "..", job$gmt.gene.set.name, ".pdf")
  params$gene.set.name.display <- paste0("Gene Set:\n", params$gene.set.name.display)
  do.call(gsea.plot.curve, params)
})


# now we use fgsea to plot a table of GSEA results.
# We focus on GOBP only. 
rnk <- paste0(wd, "/rnk/cluster_pd_sc20.rnk")
fg <- fgsea.fanc(rnk = rnk, gmt = GMT.FILES.HALL.GOBP["c5.go.bp.v7.2"])

top.pathways <- lapply(c("up", "down"), function(direction) {
  fgsea.top.genesets(fgsea.res = fg$res, n = 20, direction = direction)
})

names(top.pathways) <- c("up", "down")

lapply(c("up", "down"), function(direction) {
  label.style <- list(size = 6)
  p <- fgsea::plotGseaTable(pathways = fg$genesets[top.pathways[[direction]]],
                            stats = fg$ranks, 
                            fgseaRes = fg$res, 
                            gseaParam=1, 
                            colwidths = c(7, 3, 0.8, 1.2, 1.2),
                            pathwayLabelStyle = list(size = 5), headerLabelStyle = label.style, 
                            valueStyle = label.style, axisLabelStyle = label.style)
  ggsave(paste0(wd, "/gsea/plot_table/GOBP_top20_", direction,".pdf"), p, 
         device = cairo_pdf, height = 3, width = 5)

})
