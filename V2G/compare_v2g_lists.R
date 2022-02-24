library(tidyverse)
library(viridis)
load("v2g/T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_simplified_22XY.RData")
load("v2g/T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_promoterSNPs_22XY.Rdata")

islet_v2g <- read.delim("data/gene_lists/miguelescalada_islet_pchic_t2d_list.txt")
endoc_screen <- read.delim("data/gene_lists/Thomsen_2016_endocbh1_screen_hits.txt")

endoc_bh1_promoter = promoter_ocr_v2g[promoter_ocr_v2g$celltype=="EndoC_BH1",]
endoc_genes = unique(c(endoc_bh1_promoter$gene_name, t2d_v2g_simplified[t2d_v2g_simplified$celltype=="EndoC_BH1",]$implicated_gene_name))


high_confidence_list  = unique(unlist(strsplit(islet_v2g$Active.gene.s..connected.with.associated.variant.through.pcHi.C.high.confidence.interaction.s..CHiCAGO.5., split=",")))
endoc_genes[endoc_genes %in% high_confidence_list]

