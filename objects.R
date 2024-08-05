source("~/R_packages/common/base.R")
source("~/R_packages/common/sc.R")
suppressMessages(addArchRGenome("mm10"))
suppressMessages(library(DESeq2))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(GenomicInteractions))
suppressMessages({
  R.files <- Sys.glob("R_*/*.R")
  for (R.file in R.files) {
    source(R.file)
  }
})

SAMPLES <- c("BM_ctrl_ba1rep1", "BM_tumor_ba1rep1", "BM_ctrl_ba3rep1", "BM_tumor_ba3rep1", 
             "SP_tumor_ba4rep1", "SP_tumor_ba4rep2", "SP_tumor_ba5rep1", "SP_tumor_ba5rep2")
SAMPLES.SP <- c("SP_tumor_ba4rep1", "SP_tumor_ba4rep2", "SP_tumor_ba5rep1", "SP_tumor_ba5rep2")
SAMPLES.BM <- c("BM_ctrl_ba1rep1", "BM_tumor_ba1rep1", "BM_ctrl_ba3rep1", "BM_tumor_ba3rep1")
SAMPLE.MAP <- c("BM_ctrl_rep1", "BM_tumor_rep1", "BM_ctrl_rep2", "BM_tumor_rep2", 
                "SP_tumor_rep1", "SP_tumor_rep2", "SP_tumor_rep3", "SP_tumor_rep4")
names(SAMPLE.MAP) <- SAMPLES

EXON.ONLY.DIR <- "~/hmtp/scAR/spbm2/sth/count_exonOnly/"

BC.METRICS.FILE.VEC.exonOnly <- paste0(EXON.ONLY.DIR, "/",  
                                       SAMPLES,
                                       "/outs/per_barcode_metrics.csv")

names(BC.METRICS.FILE.VEC.exonOnly) <- SAMPLES

BC.METRICS.FILE.LIST.exonOnly <- as.list(BC.METRICS.FILE.VEC.exonOnly)

FRAG.FILES.exonOnly <- paste0(EXON.ONLY.DIR, "/", SAMPLES, "/outs/atac_fragments.tsv.gz")

names(FRAG.FILES.exonOnly) <- SAMPLES

GMT.FILES <- Sys.glob("/bar/cfan/genomes/MSigDB/v7.2/msigdb_v7.2_GMTs/*symbols.gmt")
names(GMT.FILES) <- basename(GMT.FILES) %>% sub(".symbols.gmt", "", .)
GMT.FILES.core <- GMT.FILES[grepl(paste0(c("^c2", "^c5", "^h.all"), collapse = "|"),
                                  names(GMT.FILES))]
GMT.FILES.HALL.GOBP <- GMT.FILES[grepl(paste0(c("^c5.go.bp", "^h.all"), collapse = "|"),
                                       names(GMT.FILES))]

GSEA.DEFAULT <- paste0(" -collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 ",
                       " -scoring_scheme weighted -create_svgs false -include_only_symbols true",
                       " -make_sets true -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false")


de.fast.alias <- function(de, prefix = "seurat_clusters_") {
  if (any(!grepl(prefix, names(de)))) {
    stop()
  }
  aliases <- c("Ery[0]", "HSPC[1]", "Ery[2]", "Mono/Mac/DC[3]", 
               "Neu[4]", "Myelo[5]", "DC[6]", "Mega[7]", "HSPC[9]", 
               "Mono/Mac[10]", "Ery[11]", "Mega[14]", "Lympho[16]")
  clusters <- paste0(prefix, stringr::str_extract(aliases, "\\d+"))
  
  utilsFanc::check.intersect(names(de), "names(de)", clusters, "clusters")
  
  utilsFanc::check.dups(clusters, x.name = "clusters")
  
  names(aliases) <- clusters
  names(de) <- aliases[names(de)]
  de <- de.sync.name(de)
  return(de)
}
