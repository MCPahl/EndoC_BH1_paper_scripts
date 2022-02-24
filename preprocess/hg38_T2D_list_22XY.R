library("GenomicRanges")
library("tidyverse")

load("capture_c_atac_seq_celltypes_22XY.RData") #intersection_by_celltype
load("T2D_proxies_ldlink_hg38.RData") #T2D_proxies_ldlink_hg38

gene_features = read.delim("../../../annotationFiles/genecodeV30/gencode.v30.gene_type.txt")

hla = read.delim("/mnt/isilon/sfgi/pahlm/annotationFiles/TOPMED_hg38_references/blacklist/hg38/HLA_region.bed",header=F)
names(hla) = c("chr", "start", "end")
hla$start = hla$start+1

snps.gr = GRanges(seqnames=T2D_proxies_ldlink_hg38$chr, ranges=IRanges(T2D_proxies_ldlink_hg38$start, T2D_proxies_ldlink_hg38$end))
hla.gr = GRanges(hla)

hla_index = as.data.frame(findOverlaps(snps.gr, hla.gr))[,1]

T2D_proxies_ldlink_hg38 = T2D_proxies_ldlink_hg38[-hla_index, ]



snps.gr = GRanges(seqnames=T2D_proxies_ldlink_hg38$chr, ranges=IRanges(T2D_proxies_ldlink_hg38$start, T2D_proxies_ldlink_hg38$end))


#intersect
out = lapply(intersection_by_celltype, function(intersection_set){
	intersection.gr = GRanges(seqnames=intersection_set$ocr_chr, ranges=IRanges(intersection_set$ocr_start, intersection_set$ocr_end))
	index = as.data.frame(findOverlaps(intersection.gr, snps.gr))
	data.frame(intersection_set[index[,1],], T2D_proxies_ldlink_hg38[index[,2],])
})

out = do.call("rbind", out)
out = out[out$proxy!=".",]
out = out %>% 
	separate_rows(bait_name, sep="\\|") %>% 
	mutate(bait_gene_name=gsub("\\+.*","",bait_name), bait_transcript_id=gsub(".*\\+","",bait_name)) %>%
	mutate(otherEnd_gene_name=gsub("\\+.*","",otherEnd_name), otherEnd_transcript_id=gsub(".*\\+","",otherEnd_name)) %>% 
	rename(gene_name = bait_gene_name, transcript_id = bait_transcript_id) %>%
	left_join(gene_features) %>% 
	rename(bait_gene_name = gene_name, bait_transcript_id = transcript_id, bait_gene_id = gene_id, bait_biotype = biotype, bait_gene_type = gene_type) %>% 
	rename(gene_name = otherEnd_gene_name, transcript_id = otherEnd_transcript_id) %>%
	left_join(gene_features) %>% 
	rename(otherEnd_gene_name = gene_name, otherEnd_transcript_id = transcript_id, otherEnd_gene_id = gene_id, otherEnd_biotype = biotype, otherEnd_gene_type = gene_type) %>% 
	select(-c(bait_name,bait_transcript_id, otherEnd_name, otherEnd_transcript_id, MAF, Dprime, R2, coordinate_hg19, Correlated_Alleles, RegulomeDB, chr, start, end, otherEnd_gene_type, bait_gene_type)) %>%
	mutate(proxy_located = ifelse(bait_chr == ocr_chr & abs((ocr_start+ocr_end)/2 - (bait_start+out$bait_end)/2) < abs((ocr_start+ocr_end)/2 - (otherEnd_start+otherEnd_end)/2), "bait", "otherEnd"))%>%
	unique() %>%
	group_by(across(c(-proxy))) %>%
	summarize(proxy=paste(unique(proxy), collapse="|")) %>% 
	ungroup() %>%
	mutate(loop_dist = ifelse(bait_chr != otherEnd_chr, Inf, floor(abs((otherEnd_start + otherEnd_end)/2 - (bait_start + bait_end)/2))))%>%
	filter( (is.na(bait_gene_name) & proxy_located=="otherEnd")==F) %>% 
	filter( (is.na(otherEnd_gene_name) & proxy_located=="bait")==F) %>% 
	relocate(sentinel, proxy, ancestry, bait_chr, bait_start, bait_end, bait_gene_name, bait_gene_id, otherEnd_chr, otherEnd_start, otherEnd_end, otherEnd_gene_name, otherEnd_gene_id, loop_dist,proxy_located, N_reads, score, resolution, celltype, ocr_chr, ocr_start, ocr_end, ocr_id)

out = out[(is.na(out$bait_gene_name) & is.na(out$otherEnd_gene_name)) == F, ]

t2d_v2g = out
save(t2d_v2g, file = "T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_22XY.RData")
write.table(t2d_v2g, file = "T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_22XY.txt", sep="\t", quote=F, row.names=F)


out_simplified = out %>%
 ungroup() %>%
 mutate(implicated_gene_name = ifelse(proxy_located=="otherEnd", out$bait_gene_name, out$otherEnd_gene_name)) %>% 
 mutate(implicated_gene_id = ifelse(proxy_located=="otherEnd",  out$bait_gene_id, out$otherEnd_gene_id)) %>% 
 select(sentinel, proxy, ancestry, celltype, implicated_gene_name, implicated_gene_id, resolution) %>%
 unique() 
 out_simplified = out_simplified[complete.cases(out_simplified),]

t2d_v2g_simplified = out_simplified
save(t2d_v2g_simplified, file = "T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_simplified_22XY.RData")
write.table(t2d_v2g_simplified, file = "T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_simplified_22XY.txt", sep="\t", quote=F, row.names=F)




promoters = read.delim("../../../annotationFiles/genecodeV30/gencode.v30.promoter.bed", header=F)
names(promoters) = c("chr", "start", "end", "transcript", "score", "strand")
promoters = promoters %>% mutate(gene_name = gsub("\\+.*", "", transcript), transcript_id = gsub(".*\\+", "", transcript)) %>% select(-transcript)
promoters = left_join(promoters, gene_features)

promoters.gr = GRanges(promoters)

index = as.data.frame(findOverlaps(promoters.gr, snps.gr))

promoter_snps = data.frame(promoters[index[,1],], T2D_proxies_ldlink_hg38[index[,2],])

promoter_snps = promoter_snps %>% select(sentinel, proxy, chr.1, end.1, ancestry, gene_name, gene_id, biotype) %>% distinct() %>% rename(proxy_chr = chr.1, proxy_pos = end.1)

write.table(promoter_snps , file = "T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_promoterSNPs_22XY.txt", sep="\t", quote=F, row.names=F)
