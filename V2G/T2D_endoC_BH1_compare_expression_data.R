library(tidyverse)

load("data/rnaseq/tpm_data.Rdata")
load("data/intersections/capture_c_atac_seq_celltypes_22XY.RData")

gene = read.delim("/mnt/isilon/sfgi/pahlm/annotationFiles/genecodeV30/gencode.v30.gene_type.txt")
gene = gene %>% select(gene_name, gene_id, transcript_id)

intersection_by_celltype = do.call('rbind', intersection_by_celltype)

intersection_spread = intersection_by_celltype %>%
	separate_rows(bait_name, sep = "\\|") %>%
	separate_rows(otherEnd_name, sep = "\\|") %>% 
	mutate(bait_gene_name = gsub("\\+.*", "", bait_name),
			bait_transcript_id = gsub(".*\\+", "", bait_name), 
			bait_gene_name = gsub("\\+.*", "", bait_name),
			otherEnd_transcript_id = gsub(".*\\+", "", otherEnd_name), 
			otherEnd_gene_name = gsub("\\+.*", "", otherEnd_name),
			) %>% select(-bait_name, otherEnd_name)

gene = gene %>% rename(bait_gene_name = gene_name, bait_transcript_id = transcript_id, bait_gene_id = gene_id)

intersection_spread = left_join(intersection_spread, gene)

gene = gene %>% rename(otherEnd_gene_name = bait_gene_name, otherEnd_transcript_id = bait_transcript_id, otherEnd_gene_id = bait_gene_id)

intersection_spread = left_join(intersection_spread, gene)

intersection_spread = intersection_spread %>% select(-c(bait_transcript_id, otherEnd_transcript_id)) %>% unique()

ocr_connected_genes = intersection_spread %>%
	select(bait_gene_id, celltype) %>%
	unique() %>%
	rename(gene_id = bait_gene_id) %>%
	mutate(class = "contacted")


 tpm_mean = tpm_dat %>% mutate(celltype = gsub("_rep.*", "", sample)) %>% group_by(gene_id, gene_name, celltype) %>% summarize(tpm_mean = mean(tpm))


 tpm_mean$celltype = ifelse(tpm_mean$celltype == "EndoCBH1", "EndoC_BH1", 
 	ifelse(tpm_mean$celltype =="Adipose", "hMSC_Adipocytes", 
 		ifelse(tpm_mean$celltype == "HEPG2", "Hepg2",
 			ifelse(tpm_mean$celltype == "Osteoblasts", "hMSC_osteoblasts", 
 				ifelse(tpm_mean$celltype == "SGBS_Diff", "SGBS_Diff", "SGBS_Undiff")))))

tpm_mean = left_join(tpm_mean, ocr_connected_genes)

tpm_mean$class[is.na(tpm_mean$class)] = "not_contacted"

pdf("plots/contacted_genes_expression.pdf", useDingbats=FALSE)
ggplot(tpm_mean, aes(x = celltype, y = log2(tpm_mean+1), fill = class))+
	geom_boxplot(outlier.shape = NA)+
	#geom_violin()+
  theme_minimal()+
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()

load("/mnt/isilon/sfgi/pahlm/annotationFiles/gtex/v8/gtex_median_tissue_expression_tpm.RData")
gtex_relevant = gtex.tissue.med %>%
 gather(key = tissue, value = tpm, Adipose...Subcutaneous, Adipose...Visceral..Omentum., Adrenal.Gland, Artery...Aorta, Artery...Coronary, Artery...Tibial, Bladder, Brain...Amygdala, Brain...Anterior.cingulate.cortex..BA24., Brain...Caudate..basal.ganglia., Brain...Cerebellar.Hemisphere, Brain...Cerebellum, Brain...Cortex, Brain...Frontal.Cortex..BA9., Brain...Hippocampus, Brain...Hypothalamus, Brain...Nucleus.accumbens..basal.ganglia., Brain...Putamen..basal.ganglia., Brain...Spinal.cord..cervical.c.1., Brain...Substantia.nigra, Breast...Mammary.Tissue, Cells...Cultured.fibroblasts, Cells...EBV.transformed.lymphocytes, Cervix...Ectocervix, Cervix...Endocervix, Colon...Sigmoid, Colon...Transverse, Esophagus...Gastroesophageal.Junction, Esophagus...Mucosa, Esophagus...Muscularis, Fallopian.Tube, Heart...Atrial.Appendage, Heart...Left.Ventricle, Kidney...Cortex, Kidney...Medulla, Liver, Lung, Minor.Salivary.Gland, Muscle...Skeletal, Nerve...Tibial, Ovary, Pancreas, Pituitary, Prostate, Skin...Not.Sun.Exposed..Suprapubic., Skin...Sun.Exposed..Lower.leg., Small.Intestine...Terminal.Ileum, Spleen, Stomach, Testis, Thyroid, Uterus, Vagina, Whole.Blood) %>%
  filter(tissue %in% c("Adipose...Subcutaneous", "Adipose...Visceral..Omentum", "Pancreas", "Liver"))

gtex_relevant

key_df = data.frame(tissue = c("Adipose...Subcutaneous", "Adipose...Visceral..Omentum.", "Adipose...Subcutaneous", "Adipose...Visceral..Omentum.","Adipose...Subcutaneous", "Adipose...Visceral..Omentum.", "Liver", "Pancreas"), celltype = c("hMSC_Adipocytes","hMSC_Adipocytes", "SGBS_Diff","SGBS_Diff", "SGBS_Undiff", "SGBS_Undiff", "HEPG2", "EndoC_BH1"))



gtex_connected = gtex_relevant %>% left_join(key_df) %>% left_join(ocr_connected_genes)

gtex_connected = gtex_connected[gtex_connected$tissue %in% key_df$tissue,]

gtex_connected$class[is.na(gtex_connected$class)] = "not_connected"


pdf("plots/gtex_connected_genes_expression.pdf", useDingbats=FALSE)
ggplot(gtex_connected, aes(x = tissue, y = log2(tpm+1), fill = class))+
	geom_boxplot(outlier.shape = NA)+
  theme_minimal()+
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()