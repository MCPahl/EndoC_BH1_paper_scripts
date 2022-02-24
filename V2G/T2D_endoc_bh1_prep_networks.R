library(tidyverse)
load("v2g/T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_promoterSNPs_22XY.Rdata")
load("v2g/T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_simplified_22XY.RData")

promoter_v2g_simplified = promoter_ocr_v2g %>% 
	select(sentinel, proxy, gene_name, gene_id, celltype) %>%
	mutate(class = "promoterOCR") %>% 
	unique()


PIR_v2g_simplifed = t2d_v2g_simplified %>% 
	separate_rows(proxy) %>%
	select(sentinel, proxy, implicated_gene_name, implicated_gene_id,celltype) %>%
	rename(gene_name = implicated_gene_name, gene_id = implicated_gene_id) %>%
	mutate(class = "PIROCR") %>% 
	unique()


v2g = bind_rows(promoter_v2g_simplified, PIR_v2g_simplifed) %>% 
	group_by(sentinel, proxy, gene_name, gene_id, celltype) %>%
	summarize(class = paste(unique(class), collapse="|"))

v2g$celltype = v2g$celltype %>% 
	str_replace("EndoC_BH1", "EndoCBH1") %>%
	str_replace("hMSC_Adipocytes", "Adipose") %>%
	str_replace("HEPG2", "Hepg2") %>%
	str_replace("hMSC_Osteoblasts", "Osteoblasts") %>%
	str_replace("hMSC_osteoblasts", "Osteoblasts") %>%
	str_replace("SGBS_diff", "SGBS_Diff") %>%
	str_replace("SGBS_undiff", "SGBS_Undiff") 


load("data/rnaseq/tpm_data.Rdata")
tpm_dat$celltype = gsub("_rep.*", "",tpm_dat$sample)

tpm_dat = tpm_dat %>% group_by(gene_name, gene_id, celltype) %>% summarize(tpm = mean(tpm))

v2g = left_join(v2g, tpm_dat)
v2g = v2g[complete.cases(v2g),]

v2g.l = split(v2g, v2g$celltype)

lapply(v2g.l, function(v){
	write.table(v, file = paste0("tables/anno_cytoscape/", unique(v$celltype), "_variant2gene_mapping.tsv"), 
		quote=F, row.names=F, sep="\t")
	})

t2d_mah = read.csv("/mnt/isilon/sfgi/pahlm/projects/T2D_v2g/leadSNPs/leadSNPs_T2D_McCarthy.csv")
t2d_mah  = t2d_mah %>% mutate(sentinel = rsid, study = "Mahajan2018") %>% 
select(sentinel, chr, P, study, Ancestry)

t2d_vuj = read.csv("/mnt/isilon/sfgi/pahlm/projects/T2D_v2g/leadSNPs/leadSNPs_T2D_Voight.csv")
t2d_vuj  = t2d_vuj %>% mutate(sentinel = rsid, study = "Vujkovic2020") %>% 
select(sentinel, chr, P, study, Ancestry)

sent_list = rbind(t2d_vuj, t2d_mah)

v2g = left_join(v2g, sent_list)