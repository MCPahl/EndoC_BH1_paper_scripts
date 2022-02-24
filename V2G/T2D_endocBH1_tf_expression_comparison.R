library(tidyverse)
library(viridis)
load("v2g/T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_simplified_22XY.RData")
load("data/motifbreakR_results.RData")
load("data/rnaseq/tpm_data.Rdata")

tpm_mean = tpm_dat %>% mutate(celltype = gsub("_rep.*", "", sample)) %>%
	 group_by(gene_id, gene_name, celltype) %>%
	 summarize(tpm_mean = mean(tpm)) %>% 
	 rename(gene_name = gene_name, gene_id= gene_id)


 tpm_mean$celltype = ifelse(tpm_mean$celltype == "EndoCBH1", "EndoC_BH1", 
 	ifelse(tpm_mean$celltype =="Adipose", "hMSC_Adipocytes", 
 		ifelse(tpm_mean$celltype == "Hepg2", "HEPG2",
 			ifelse(tpm_mean$celltype == "Osteoblasts", "hMSC_osteoblasts", 
 				ifelse(tpm_mean$celltype == "SGBS_Diff", "SGBS_Diff", "SGBS_Undiff")))))


anno = read.delim("/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/jasper2020/TF_motif.anno.txt",stringsAsFactors=FALSE)

anno = anno %>% rename(geneSymbol=TF , providerName=motif_id) %>% select(geneSymbol, gene_name, TF_family, providerName) %>% unique()

motif_targets = motifbreakr.results %>% left_join(anno) %>% left_join(t2d_v2g_simplified) %>% left_join(tpm_mean)

motif_targets_dat = motif_targets %>% select(geneSymbol, proxy, alleleDiff, celltype, implicated_gene_name, implicated_gene_id, tpm_mean) %>% unique()


motif_targets_dat = motif_targets_dat[complete.cases(motif_targets_dat), ]


pdf("plots/tf_motif_proxy_disruption_connected_gene_expression.pdf", useDingbats=FALSE)
ggplot(motif_targets_dat, aes(y=log2(tpm_mean+1), x=alleleDiff))+
geom_hex()+
scale_fill_continuous(type = "viridis") +
theme_minimal()+
facet_wrap(vars(celltype))
dev.off()


tf_targets_dat = motif_targets_dat %>% group_by(geneSymbol, celltype) %>% summarize(alleleDiff_abs = mean(alleleDiff), tpm_mean = mean(tpm_mean), proxy_count = length(unique(proxy)), gene_count = length(unique(implicated_gene_id))) 

pdf("plots/tf_motif_proxy_disruption_connected_gene_expression_tf_mean.pdf", useDingbats=FALSE, width=10, height=10)
ggplot(tf_targets_dat , aes(y=gene_count, x=alleleDiff_abs))+
geom_point(aes(color = log2(tpm_mean+1), size= proxy_count))+
scale_color_continuous(type = "viridis") +
  geom_text(aes(label=ifelse(gene_count > quantile(gene_count, probs=0.95),as.character(geneSymbol),'')),hjust=0,vjust=0)+
theme_minimal()+
facet_wrap(vars(celltype))+
theme(aspect.ratio=1)
dev.off()

write.csv(tf_targets_dat, file = "tables/tf_expression_disruption.csv", quote=F, row.names=F)