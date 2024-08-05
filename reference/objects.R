source("~/R_packages/common/base.R")
source("~/R_packages/common/sc.R")
source("~/R_packages/common/bioc.R")
suppressMessages(library(DESeq2))
suppressMessages({
  R.files <- Sys.glob("R_*/*.R")
  for (R.file in R.files) {
    source(R.file)
  }
})

GMT.FILES <- Sys.glob("/bar/cfan/genomes/MSigDB/v7.2/msigdb_v7.2_GMTs/*symbols.gmt")
names(GMT.FILES) <- basename(GMT.FILES) %>% sub(".symbols.gmt", "", .)
GMT.FILES.core <- GMT.FILES[grepl(paste0(c("^c2", "^c5", "^h.all"), collapse = "|"),
                                  names(GMT.FILES))]
GMT.FILES.HALL.GOBP <- GMT.FILES[grepl(paste0(c("^c5.go.bp", "^h.all"), collapse = "|"),
                                       names(GMT.FILES))]

GSEA.DEFAULT <- paste0(" -collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 ",
                       " -scoring_scheme weighted -create_svgs false -include_only_symbols true",
                       " -make_sets true -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false")

