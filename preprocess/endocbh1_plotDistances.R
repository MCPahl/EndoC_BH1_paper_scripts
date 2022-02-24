library(GenomicFeatures)
library(tidyverse)
setwd("/mnt/isilon/sfgi/pahlm/projects/T2D_endocBH1")
load("/mnt/isilon/sfgi/pahlm/projects/T2D_endocBH1/data/intersections/capture_c_atac_seq_celltypes_22XY.RData")


ocr_interactions = lapply(intersection_by_celltype, function(intsect_cell){
	intsect_cell$dist = ifelse(intsect_cell$bait_chr == intsect_cell$otherEnd_chr,
		floor(abs((intsect_cell$bait_start+intsect_cell$bait_end)/2-(intsect_cell$otherEnd_start+intsect_cell$otherEnd_end)/2)),
		NA
	)
	intsect_cell
	})
ocr_interactions = do.call("rbind", ocr_interactions)

pdf("plots/promoter_ocr_distances.pdf", useDingbats=FALSE)
ggplot(ocr_interactions, aes(color=celltype, x= log10(dist) ))+
stat_ecdf(position="identity")+
theme_minimal()+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
#scale_x_continuous(limits = c(-3, 2))
dev.off()





gwas.matrix = list_to_matrix()
gwas.matrix = make_comb_mat(gwas.matrix)

pdf("plots/Upset_T2D_GWAS_genes.pdf", width=12, height=7,)
UpSet(gwas.matrix,
	#set_order = c("AaM","Bipolar"BMI", "HEIGHT", "MDD", "SLEEP"),
	pt_size = unit(5, "mm"), lwd = 1,
	comb_col = c("red", "blue", "black", "darkgreen", "orange", "purple")[comb_degree(gwas.matrix)],
	comb_order = order(comb_size(gwas.matrix)),
	#set_on_rows
	)
dev.off()
