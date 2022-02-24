library(data.table)
library(tidyverse)

sentinels = read.delim("gwas_data/sentinel_snpset/sentinels_hg19.bed", header=F)
sentinels = sentinels %>% 
select(V1, V3, V4) %>% 
rename(chrom = V1, pos = V3, rsid = V4)

m_lead = read.csv("/mnt/isilon/sfgi/pahlm/projects/T2D_v2g/leadSNPs/leadSNPs_T2D_McCarthy.csv")
v_lead = read.csv("/mnt/isilon/sfgi/pahlm/projects/T2D_v2g/leadSNPs/leadSNPs_T2D_Voight.csv")

mlead = m_lead %>% 
	select(rsid, chr, position, risk_allele, other_allele,  study, Ancestry) %>%
	mutate(ancestry = paste(study, Ancestry, sep="_")) %>%
	select(-c(study, Ancestry)) %>%
	unique() %>%
	mutate(chr = paste("chr", chr, sep=""), position= as.integer(gsub(",", "", position)))

vlead = v_lead %>% 
	select(rsid, chr, position, risk.allele, other_allele, Study, Ancestry) %>%
	mutate(ancestry = paste(Study, Ancestry, sep="_")) %>%
	select(-c(Study, Ancestry)) %>%
	unique() %>%
	mutate(chr = paste("chr", chr, sep="")) %>%
	rename(risk_allele = risk.allele)

lead = rbind(mlead, vlead)

gcv19= read.delim("/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.gene_length.txt",header=F)
gcv19= gcv19 %>% select(V1, V2) %>% rename(gene_id = V1, gene_name = V2)

#GTEX eQTL
dir = "/mnt/isilon/sfgi/pahlm/datasets_downloaded/GTeX/eqtl_v7/GTEx_Analysis_v7_eQTL_sig_associations/"

files = list.files(dir)
eqtl_out = lapply(files, function(file){
	eqtl = fread(paste0(dir, "/", file))
	eqtl.m = eqtl %>% 
		separate(variant_id, sep="_", into = c("chr", "position", "allele1", "allele2", "genome")) %>%
		mutate(chr = paste0("chr", chr)) %>%
		select(-genome)
	eqtl.m$position = as.integer(eqtl.m$position)
	sentinels_eqtl_sig_risk = inner_join(lead, eqtl.m)
	sentinels_eqtl_sig_risk = left_join(sentinels_eqtl_sig_risk,gcv19)
})

for(i in seq_along(eqtl_out)){
	eqtl_out[[i]]$tissue = gsub(".v7.signif_variant_gene_pairs.txt.gz", "", files[[i]])
}

eqtl_out = do.call("rbind", eqtl_out)

write.csv(eqtl_out, file= "tables/gtex_eqtl_sentinel_variants_alltissues.txt", quote=F, row.names=F)

eqtl_out.m = eqtl_out %>% mutate(gene_id = gsub("\\..*", "", gene_id))

load("v2g/T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_22XY.RData")


t2d_v2g_conv = t2d_v2g %>% select(sentinel, proxy, ancestry, bait_gene_name, bait_gene_id, celltype) %>%
	unique() %>%
	mutate(bait_gene_id = gsub("\\..*", "", bait_gene_id)) %>% 
	rename(rsid = sentinel, ancestry = ancestry, gene_id = bait_gene_id, gene_name = bait_gene_name)

t2d_pcc_eqtl_intersect = inner_join(eqtl_out.m, t2d_v2g_conv ) %>% 
	rename(sentinel = rsid) %>%
	select(sentinel, proxy, celltype, ancestry, gene_name, gene_id, tissue) %>% 
	filter(
		(celltype == "hMSC_Adipocyte" | grepl("SGBS", celltype)) & grepl("Adipose", tissue) |
		(celltype == "HEPG2" & tissue == "Liver") | 
		(celltype == "EndoC_BH1" & tissue == "Pancreas")  
	)%>%
	unique() 
write.csv(t2d_pcc_eqtl_intersect, file= "tables/pcc_gtex_eqtl_sentinel_variants_cell_tissues.txt", quote=F, row.names=F)


load("v2g/T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_promoterSNPs_22XY.Rdata")
promoter_ocr_v2g_conv = promoter_ocr_v2g %>% select(sentinel, proxy, ancestry, gene_name, gene_id, celltype) %>%
	group_by(sentinel, ancestry, gene_name, gene_id, celltype) %>% 
	summarize(proxy = sort(unique(paste(proxy, collapse="|")))) %>%
	unique() %>% 
	mutate(gene_id = gsub("\\..*", "", gene_id))   %>% 
	rename(rsid = sentinel)

t2d_promoter_eqtl_intersect = inner_join(eqtl_out.m, promoter_ocr_v2g_conv ) %>% 
	select(rsid, proxy, celltype, ancestry, gene_name, gene_id, tissue) %>% 
	filter(
		(celltype == "hMSC_Adipocyte" | grepl("SGBS", celltype)) & grepl("Adipose", tissue) |
		(celltype == "HEPG2" & tissue == "Liver") | 
		(celltype == "EndoC_BH1" & tissue == "Pancreas")  
	)%>%
	rename(sentinel = rsid) %>%
	unique() 

write.csv(t2d_promoter_eqtl_intersect, file= "tables/promoter_ocr_gtex_eqtl_sentinel_variants_cell_tissues.txt", quote=F, row.names=F)


