source("objects.R")
wd <- "publication/sth/pub2.5_DC_homer/"
dir.create(wd)

da.dir <- "da_noNZ/sth/da_noNZ2.1_scc_BM//"
deseq2.homer.gather(paste0(da.dir, "/homer/"), 
                    regex.include = paste0("scc_", c(6), "/"),
                    regex.from = "scc_")

homer.plot.table(da.homer.dir = paste0(da.dir, "/homer/"), use.regex = F, col.num = 3,
               plot.logo = T,
               motif.map = paste0(da.dir, "/motif_cluster_map/scc6.tsv"), 
               regex.include = paste0("scc_", c(6), "/res.exp_padj_top_1000"),
               regex.from = "scc_", 
               out.dir = wd)

# see if there is a significant overlap between peaks with Stat3 and CEBP motifs

da <- readRDS(paste0(da.dir, "/ao.bulk.list.Rds"))
da <- da['scc_6']
da <- deseq2.summary(da, log2fc.cutoff = 0.5, padj.cutoff = 0.1)
dar <- da$scc_6$summary$up.genes %>% utilsFanc::loci.2.gr()

motifs <- list(
  Stat3 = rtracklayer::import("~/motifs/homer/sth/each_motif/Stat3(Stat).bed"),
  CEBP = rtracklayer::import("~/motifs/homer/sth/each_motif/CEBP(bZIP).bed")
)

motif.df <- lapply(motifs, function(motifs) {
  o <- findOverlaps(dar, motifs)
  out <- rep(F, length(dar))
  out[(1:length(dar) %in% queryHits(o))] <- T
  return(out)
}) %>% as.data.frame()

View(motif.df)
stats <- list(n = nrow(motif.df))
stats$n.Stat3 <- sum(motif.df$Stat3)
stats$n.CEBP <- sum(motif.df$CEBP)
stats$n.both <- sum(motif.df$CEBP & motif.df$Stat3)
stats$pct.Stat3 <- stats$n.Stat3/stats$n
stats$pct.CEBP <- stats$n.CEBP/stats$n
stats$pct.both <- stats$n.both/stats$n
stats$pct.both.expected <- stats$pct.Stat3 * stats$pct.CEBP
stats$enrich <- stats$pct.both / stats$pct.both.expected
stats <- stats %>% as.data.frame() %>% t()
stats <- round(stats, digits = 2)
write.table(stats, paste0(wd, "/CEBP_Stat3_co-occurence.tsv"), col.names = F,
            row.names = T, sep = "\t", quote = F)
# Conclusion: no significant overlap at all.