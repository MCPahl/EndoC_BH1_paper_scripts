library(GenomicFeatures)
library(tidyverse)

setwd("/mnt/isilon/sfgi/pahlm/projects/T2D_endocBH1")

txdb <- makeTxDbFromGFF("/mnt/isilon/sfgi/pahlm/annotationFiles/genecodeV30/gencode.v30.annotation.gtf", format="gtf")

make_location_anno = function(x){
	promoters = promoters(x,upstream=1500, downstream=500)
	promoters$class = "promoter"
	exon = unlist(exonsBy(x))
	exon$class = "exon"
	fiveUTR = unlist(fiveUTRsByTranscript(x))
	fiveUTR$class = "fiveUTR"
	intron = unlist(intronsByTranscript(x))
	intron$class = "intron"
	threeUTR = unlist(threeUTRsByTranscript(x))
	threeUTR$class= "threeUTR"
	list(promoters, exon, fiveUTR, intron, threeUTR)
}

genomic_anno = make_location_anno(txdb)



annotate_feature = function(df, annoList){
	df.gr = GRanges(df)
	df.feature.anno = lapply(annoList, function(feature){
		index = unique(as.data.frame(findOverlaps(df.gr, feature),)[,1])
		ifelse(1:nrow(df) %in% index, unique(feature$class), "." )
	})
	df.feature.anno = do.call("cbind", df.feature.anno)
	df.feature.anno = data.frame(df.feature.anno)
	df$feature = ifelse(df.feature.anno$X1 == "promoter", "promoter", 
		ifelse(df.feature.anno$X2 == "fiveUTR", "fiveUTR", 
			ifelse(df.feature.anno$X3 == "threeUTR", "threeUTR",
				ifelse(df.feature.anno$X4 == "exon", "exon",
					ifelse(df.feature.anno$X4 == "intron", "intron", "intergenic")))))
	
	df
}

load("v2g/T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_22XY.RData")
load("gwas_data/proxy_snpset/T2D/T2D_proxies_ldlink_hg38.RData")

t2d_proxy = t2d_v2g %>% separate_rows(proxy, sep="\\|") %>% select(proxy, celltype) %>% left_join(T2D_proxies_ldlink_hg38) %>% select(proxy, celltype, chr, start, end) %>%  unique()
annotated_snps = annotate_feature(t2d_proxy, genomic_anno)
annotated_snps_tally = annotated_snps %>% group_by(celltype, feature) %>% unique() %>% tally() %>% ungroup()
annotated_snps_totals = annotated_snps_tally %>% group_by(celltype) %>% summarize(sum = sum(n))
annotated_snps_tally = left_join(annotated_snps_tally, annotated_snps_totals) %>% mutate(ratio = n/sum)
annotated_snps_tally$celltype = factor(annotated_snps_tally$celltype, levels= c("EndoC_BH1", "HEPG2", "hMSC_Adipocytes", "hMSC_osteoblasts", "SGBS_Diff", "SGBS_Undiff"))



pdf("plots/implicated_proxies_genomic_annotation__.pdf")
ggplot(annotated_snps_tally, aes(x= celltype, y=ratio, fill= feature))+
geom_bar(position="stack", stat="identity")+
 theme_minimal()+
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top")+
  coord_flip()+
  scale_x_discrete(limits = rev(levels(annotated_snps_tally$celltype)))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()




