library(tidyverse)

dir = "data/gene_lists/"
files = c("gwas_catalog_t2d.txt", "miguelescalada_islet_pchic_t2d_list.txt", "open_targets_t2d_list.txt", "Thomsen_2016_endocbh1_screen_hits.txt")


gwas_cat = read.delim(paste(dir, files[1], sep="/"))
gwas_cat = unique(unlist(strsplit(gsub(" ", "", gwas_cat$Mapped.gene), ",")))
gwas_cat = gwas_cat[(gwas_cat == "-")==F]
gwas_cat = gwas_cat[grepl("x", gwas_cat)==F]


miguelescalada_islet_pchic = read.delim(paste(dir, files[2], sep="/"))

miguelescalada_islet_pchic = c(miguelescalada_islet_pchic$Active.gene.s..connected.with.associated.variant.through.pcHi.C.high.confidence.interaction.s..CHiCAGO.5., 
miguelescalada_islet_pchic$Active.gene.s..connected.with.associated.variant.through.pcHi.C.moderate.confidence.interaction.s..CHiCAGO.2.5.,
miguelescalada_islet_pchic$Active.gene.s..assigned.to.variant.containing.enhancer.through.imputation,
miguelescalada_islet_pchic$Active.gene.s..connected.with.variant.through.hub,
miguelescalada_islet_pchic$Active.gene.within.the.same.bait.fragment.as.promoter.variant..or.is.in.close.linear.proximity...10.kb..to.variant.containing.enhancer)

miguelescalada_islet_pchic = unique(unlist(strsplit(miguelescalada_islet_pchic, split=" ")))
miguelescalada_islet_pchic = unique(unlist(strsplit(miguelescalada_islet_pchic, split=",")))

miguelescalada_islet_pchic = miguelescalada_islet_pchic[(miguelescalada_islet_pchic != "None")]

open_targets_t2d = read.delim(paste(dir, files[3], sep="/"))
open_targets_t2d = unique(unlist(strsplit(gsub(" ", "", open_targets_t2d$L2G), split = ",")))


Thomsen_2016_endocbh1_screen_hits = read.delim(paste(dir, files[4], sep="/"))
Thomsen_2016_endocbh1_screen_hits = unique(Thomsen_2016_endocbh1_screen_hits$Gene)



Thomsen_2016_endocbh1_screen_hits = data.frame(source = "Thomsen_2016_endocbh1_screen_hits", gene_name = Thomsen_2016_endocbh1_screen_hits)

miguelescalada_islet_pchic = data.frame(source = "miguelescalada_islet_pchic", gene_name = miguelescalada_islet_pchic)

open_targets_t2d = data.frame(source = 'open_targets_t2d', gene_name = open_targets_t2d)

gwas_cat = data.frame(source = "gwas_cat", gene_name = gwas_cat)

gene_lists = rbind(Thomsen_2016_endocbh1_screen_hits, rbind(miguelescalada_islet_pchic, rbind(gwas_cat, open_targets_t2d) ))


load("v2g/T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_promoterSNPs_22XY.Rdata")
load("v2g/T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_simplified_22XY.RData")


endoc_bh1_genes = unique(c(t2d_v2g_simplified$implicated_gene_name[t2d_v2g_simplified$celltype=="EndoC_BH1"],
	promoter_ocr_v2g[promoter_ocr_v2g$celltype == "EndoC_BH1",]$gene_name))


gene_annot_overlap = gene_lists %>% group_by(source) %>% summarize(count = length(gene_name), overlap = length(intersect(gene_name,endoc_bh1_genes))) %>% mutate(overlap_percent = overlap/count * 100)

pdf("plots/compare_genelists.pdf", useDingbats=FALSE, height = 2, width = 3)
ggplot(gene_annot_overlap, aes(x = source, y = overlap_percent))+
	geom_bar(stat="identity", width=0.5)+
	theme_minimal()+
	ylim(0,100)+
	theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()

gene_list_table = gene_lists[gene_lists$gene_name %in% intersect(gene_lists$gene_name,endoc_bh1_genes),]

write.csv(gene_list_table, file="tables/endoc_bh1_overlapping_geneset_table.csv", quote=F, row.names=F)

klaus_pancreas_v2g = read.csv("data/gene_lists/chun_T2D_Mahajan_NatGenet2018b.v2g_mapping.csv")
klaus_pancreas_v2g = unique(data.frame(gene_name= klaus_pancreas_v2g$gene_name, source = paste("Su2021_", klaus_pancreas_v2g$cell, sep="")))

klaus_pancreas_v2g_gene_annot_overlap = klaus_pancreas_v2g  %>% group_by(source) %>% summarize(count = length(gene_name), overlap = length(intersect(gene_name,endoc_bh1_genes))) %>% mutate(overlap_percent = overlap/count * 100)

pdf("plots/klaus_compare_genelists.pdf", useDingbats=FALSE, height = 2, width = 3)
ggplot(klaus_pancreas_v2g_gene_annot_overlap , aes(x = source, y = overlap_percent))+
	geom_bar(stat="identity", width=0.5)+
	theme_minimal()+
	ylim(0,100)+
	theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()


panc_cells_hic_gene_list_table = klaus_pancreas_v2g[klaus_pancreas_v2g$gene_name %in% intersect(klaus_pancreas_v2g$gene_name,endoc_bh1_genes),]
write.csv(panc_cells_hic_gene_list_table, file="tables/endoc_bh1_overlapping_klauspanccells_table.csv", quote=F, row.names=F)
