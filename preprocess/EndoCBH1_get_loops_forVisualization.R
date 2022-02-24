library(tidyverse)
library(GenomicRanges)

data_dir = "/mnt/isilon/sfgi/pahlm/projects/T2D_endocBH1/data/pr_captureC/"
files = list.files(data_dir)
files = files[grepl("ibed", files)]

gene = "ABCC8"

dat = lapply(files, function(file){
	celltype = gsub("_.frag.ibed", "", file)
	loop_ocr = read.delim(
		paste0(data_dir, "/", file)
	)
	loop_ocr = loop_ocr[grepl(gene, loop_ocr$bait_name),]
	loop_ocr %>% 
	select(bait_chr, bait_start, bait_end, otherEnd_chr,otherEnd_start, otherEnd_end, score) %>%
	unique() %>% 
	mutate(celltype = celltype)
})

dat = do.call("rbind", dat)

dat = split(dat, dat$celltype)

lapply(dat, function(dat_celltype){
	new.name = paste(unique(dat_celltype$celltype), "_loops.arcs", sep="")
	dat_celltype = dat_celltype[order(dat_celltype$bait_chr, dat_celltype$bait_start, dat_celltype$bait_end, dat_celltype$otherEnd_chr, dat_celltype$otherEnd_start, dat_celltype$otherEnd_end),]
	write.table(dat_celltype[,c(1:7)], file = new.name, sep = "\t", quote=F, row.names=F, col.names=F)
	})