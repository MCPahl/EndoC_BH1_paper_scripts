#Islet eqtl

library(data.table)
library(tidyverse)

#Get eqtl
islet_eqtl <- fread("data/eqtl/islet_meta.results.annotated.qvalue.postfreq.tbl.gz")

load("v2g/T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_22XY.RData")

sentinels = read.delim("gwas_data/sentinel_snpset/sentinels_hg19.bed", header=F)
sentinels = sentinels %>% 
select(V1, V3, V4) %>% 
rename(chrom = V1, pos = V3, rsid = V4)

sentinels_eqtl = inner_join(sentinels, islet_eqtl)

sentinels_eqtl_sig = sentinels_eqtl %>%
 filter(q_storey < 0.05) %>%
 select(rsid, Allele1, Allele2, Zscore, Direction, q_storey, gene_name, gene, gene_type) %>% 
 rename(gene_id = gene , eqtl_Allele1 = Allele1, eqtl_Allele2 = Allele2, eqtl_Zscore = Zscore, eqtl_Direction = Direction, eqtl_q_storey = q_storey) %>%
 unique()



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

sentinels_eqtl_sig_risk = inner_join(lead, sentinels_eqtl_sig)


sentinels_eqtl_sig_risk$increased_expression_snp = toupper(ifelse(sentinels_eqtl_sig_risk$eqtl_Zscore >0, sentinels_eqtl_sig_risk$eqtl_Allele1, sentinels_eqtl_sig_risk$eqtl_Allele2))
sentinels_eqtl_sig_risk$decreased_expression_snp = toupper(ifelse(sentinels_eqtl_sig_risk$eqtl_Zscore >0, sentinels_eqtl_sig_risk$eqtl_Allele2, sentinels_eqtl_sig_risk$eqtl_Allele1))

sentinels_eqtl_sig_risk$eqtl_class = ifelse(sentinels_eqtl_sig_risk$increased_expression_snp == sentinels_eqtl_sig_risk$risk_allele &
	   sentinels_eqtl_sig_risk$decreased_expression_snp == sentinels_eqtl_sig_risk$other_allele,
		"risk increased expression", 
		ifelse(sentinels_eqtl_sig_risk$decreased_expression_snp == sentinels_eqtl_sig_risk$risk_allele &
	   sentinels_eqtl_sig_risk$increased_expression_snp == sentinels_eqtl_sig_risk$other_allele,
	   "risk decreased expression",
	   "remove"
		)
)

sentinels_eqtl_sig_risk

t2d_v2g_conv = t2d_v2g %>% select(sentinel, proxy, ancestry, bait_gene_name, bait_gene_id, celltype) %>%
	filter(celltype == "EndoC_BH1") %>% 
	unique() %>%
	mutate(bait_gene_id = gsub("\\..*", "", bait_gene_id)) %>% 
	rename(rsid = sentinel, ancestry = ancestry, gene_id = bait_gene_id, gene_name = bait_gene_name)

t2d_pcc_eqtl_intersect = inner_join(sentinels_eqtl_sig_risk, t2d_v2g_conv ) %>% 
	select(rsid, proxy, ancestry, celltype, risk_allele, other_allele, eqtl_Zscore, eqtl_Direction, eqtl_class, eqtl_q_storey, gene_name, gene_id, gene_type ) %>%
	group_by(rsid, proxy, celltype, risk_allele, other_allele, eqtl_Zscore, eqtl_Direction, eqtl_class, eqtl_q_storey, gene_name, gene_id, gene_type ) %>%
	summarize(ancestry = paste(sort(unique(ancestry)), collapse="|")) %>% 
	unique() %>% 
	rename(sentinel = rsid)


write.csv(sentinels_eqtl_sig_risk, file= "tables/islet_eqtl_sentinel_variants.txt", quote=F, row.names=F)

write.csv(t2d_pcc_eqtl_intersect , file= "tables/pcc_islet_eqtl_sentinel_variants.txt", quote=F, row.names=F)

