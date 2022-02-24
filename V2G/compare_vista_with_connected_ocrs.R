library(GenomicRanges)
library(tidyverse)
library(rtracklayer)

vista_enhancers <- import("data/references/vista/VISTA.hg38.bed")
vista.tissue = read.delim("data/references/vista/VISTA_enhancer_positive.tissueExpression.txt",header=F)
names(vista.tissue) = c("name", "tissue", "rep")

load("data/intersections/capture_c_atac_seq_celltypes_22XY.RData")
intersection_by_celltype = do.call("rbind",intersection_by_celltype)

ocr_dat = intersection_by_celltype %>% select(contains("ocr")) %>% distinct()

ocr.gr = GRanges(seqnames= ocr_dat$ocr_chr, ranges = IRanges(ocr_dat$ocr_start+1, ocr_dat$ocr_end), id = ocr_dat$ocr_id)

index = as.data.frame(findOverlaps(ocr.gr, vista_enhancers))

data.frame(ocr.gr[index[,1],], vista_enhancers$name[index[,2]])

ocr_match = data.frame(ocr_dat[index[,1],], name = vista_enhancers$name[index[,2]])

ocr_match = left_join(ocr_match, vista.tissue)