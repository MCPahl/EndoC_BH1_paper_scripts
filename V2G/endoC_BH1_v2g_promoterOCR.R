library("GenomicRanges")
library("tidyverse")

load("gwas_data/proxy_snpset/T2D/T2D_proxies_ldlink_hg38.RData") #T2D_proxies_ldlink_hg38

T2D_proxies_ldlink_hg38 = T2D_proxies_ldlink_hg38 %>% filter(proxy != ".")


gene_features = read.delim("../../annotationFiles/genecodeV30/gencode.v30.gene_type.txt")
gene_promoters = read.delim("../../annotationFiles/genecodeV30/gencode.v30.promoter.bed",header=FALSE)
names(gene_promoters) = c("chr", "start","end", "anno", "score", "strand")
gene_promoters$gene_name = gsub("\\+.*", "", gene_promoters$anno)
gene_promoters$transcript_id = gsub(".*\\+", "", gene_promoters$anno)
gene_promoters = left_join(gene_promoters, gene_features)
gene_promoters$anno = NULL
gene_promoters$start = gene_promoters$start+1

hla = read.delim("/mnt/isilon/sfgi/pahlm/annotationFiles/TOPMED_hg38_references/blacklist/hg38/HLA_region.bed",header=F)
names(hla) = c("chr", "start", "end")
hla$start = hla$start+1

snps.gr = GRanges(seqnames=T2D_proxies_ldlink_hg38$chr, ranges=IRanges(T2D_proxies_ldlink_hg38$start, T2D_proxies_ldlink_hg38$end))
hla.gr = GRanges(hla)

hla_index = as.data.frame(findOverlaps(snps.gr, hla.gr))[,1]

T2D_proxies_ldlink_hg38 = T2D_proxies_ldlink_hg38[-hla_index, ]

atac.files = list.files("data/atac/atac_peaks/")

atac_peaks = lapply(atac.files, function(atac.file){
	dat = read.delim(paste0("data/atac/atac_peaks/", atac.file))
	dat = dat[,c(1:4,11)]
	names(dat) = c("chr", "start", "end", "id", "celltype")
	dat$start = dat$start+1
	dat 
})
atac_peaks = do.call("rbind", atac_peaks)

index = as.data.frame(findOverlaps(GenomicRanges::GRanges(gene_promoters), GenomicRanges::GRanges(atac_peaks)))

promoter_ocr = unique(data.frame(gene_promoters[index[,1],], atac_peaks[index[,2],]))

promoter_ocr = promoter_ocr %>% select(gene_name, transcript_id, gene_id, biotype, gene_type, chr.1, start.1, end.1, id, celltype) %>%
	rename(chr = chr.1, start = start.1, end = end.1) %>%
	unique()


index = as.data.frame(findOverlaps(GenomicRanges::GRanges(promoter_ocr), snps.gr))

promoter_ocr_v2g = data.frame(promoter_ocr[index[,1],], T2D_proxies_ldlink_hg38[index[,2],])
promoter_ocr_v2g = promoter_ocr_v2g %>% 
select(gene_name, gene_id, biotype, gene_type,chr, start, end, id, celltype,sentinel, proxy, ancestry) %>%
unique()

save(promoter_ocr_v2g, file = "v2g/T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_promoterSNPs_22XY.Rdata")
write.table(promoter_ocr_v2g, file = "v2g/T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_promoterSNPs_22XY.txt", quote=F, row.names=F, sep="\t")