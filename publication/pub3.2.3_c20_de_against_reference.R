source("objects.R")
wd <- paste0("publication/sth/pub3.2_c20_de//")

des <- list()
des$fanc <- readRDS("publication/sth/pub3.2_c20_de/bulk.list.Rds")
des$fanc <- deseq2.summary(des$fanc)

de.karpova <- readRDS("~/hmtp/scAR/reference/karpova2019/sth/kp1.1_geo2R/bulk.list.Rds")
treatments <- de.karpova$LSK$coldata$treatment %>% unique()
treatments <- treatments[-1]
# [1] "Grob"    "GCSF"    "AMD3100"

for (i in treatments) {
  des[[i]] <- de.karpova
  des[[i]][[1]]$summary <- de.karpova[[1]][[paste0("summary_", i)]]
  des[[i]][[1]]$res <- de.karpova[[1]][[paste0("res_", i)]]
  
  if (i == "Grob") {
    tr <- "Gro-β"
  } else if (i == "GCSF") {
    tr <- "G-CSF"
  } else {
    tr <- i
  }
  
  des[[i]][[1]]$fgsea.meta <- list(
    gene.set.name.display.up = paste0("Gene Set: Karpova 2019\n", tr, " > ctrl"),
    gene.set.name.display.down = paste0("Gene Set: Karpova 2019\n", tr, " < ctrl"),
    rank.name.display = paste0("Genes in Karpova 2019:"),
    up.name.display = paste0(tr, " > ctrl"),
    down.name.display = paste0(tr, " < ctrl")
  )
}

des$fanc$pd_sc20$fgsea.meta <- list(
  gene.set.name.display.up = "Gene Set: upDEGs\nHSPC[20] > HSPC[1]",
  gene.set.name.display.down = "Gene Set: downDEGs\nHSPC[20] < HSPC[1]",
  rank.name.display = "Genes in HSPC clusters:",
  up.name.display = "20 > 1",
  down.name.display = "20 < 1"
)

de.reciprocal.gsea(des = des, 
                   slot = "summary", out.dir = paste0(wd, "/gsea_against_reference/"), 
                   run = F, force = T, n.par.de = 1,
                   fgsea.plot = T, fgsea.width = 1.6)

#>>>>>>>> We next generate a xy plot for GCSF:
markers <- c("S100a4", "S100a6", "Ahnak")
top10 <- des$fanc$pd_sc20$summary$up.genes[1:21] # one gene doesn't exist in Karpova et al. add one more

deseq2.xyplot(pbl = des$GCSF, comp.vec = "GCSF:ctrl", comp.meta = "treatment",
              pbl.slot = "bulkNorm",
              transformation = NULL, 
              plot.dir = paste0(wd, "/plot_xy/"), 
              root.name = "xy_DEG_on_GCSF",
              highlight.genes = top10, highlight.color = utilsFanc::gg_color_hue(2)[2],
              add.label = T,
              hjust = 0.5, vjust = 0.5, use.repel = T, italic.label = T,
              label.list = markers, text.size = 6 * 0.36,
              device = "pdf", publication = T,
              repel.direction = "x", highlight.ptsize = 0.08,
              plot.each.comp = T)

