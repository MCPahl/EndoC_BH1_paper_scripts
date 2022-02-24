library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(MotifDb)
library(TFBSTools)
library(tidyverse)

#load("v2g/T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_simplified.RData")
load("v2g/T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_simplified_22XY.RData")
proxylist = t2d_v2g_simplified$proxy
proxylist = unlist(strsplit(proxylist, "\\|"))
proxylist = unique(proxylist[grepl("\\.",proxylist)==FALSE])
#### extract variants from dbSNP144
variants <- snps.from.rsid(rsid = proxylist,
                           dbSNP = SNPlocs.Hsapiens.dbSNP151.GRCh38,
                           search.genome = BSgenome.Hsapiens.UCSC.hg38)

#### subset total MotifDb
jaspar2020=subset (
        MotifDb, dataSource %in% c("jaspar2020")
        )
jaspar2020=subset(
        jaspar2020, organism %in% c("Mmusculus","Hsapiens")
)

#### run motifbreakR
motifbreakr.results <- motifbreakR(snpList = variants, pwmList = jaspar2020,
                                   filterp = TRUE,
                                   threshold = 1e-4,
                                   method = "ic",
                                   bkg = c(A=0.25, C=0.25, G=0.25, T=0.25))

motifbreakr.results = data.frame(motifbreakr.results )
names(motifbreakr.results)[6] = "proxy"

save(motifbreakr.results, file ="data/motifbreakR_results_22XY.RData")



