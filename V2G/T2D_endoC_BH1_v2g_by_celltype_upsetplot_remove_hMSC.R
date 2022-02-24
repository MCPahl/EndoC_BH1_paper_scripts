library(GenomicFeatures)
library(tidyverse)
library(ComplexHeatmap)
setwd("/mnt/isilon/sfgi/pahlm/projects/T2D_endocBH1")

load("v2g/T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_simplified_22XY.RData")


t2d_v2g_connections = t2d_v2g_simplified %>% separate_rows(proxy, sep="\\|") %>% select(implicated_gene_name, celltype) %>% unique

t2d_v2g_connections = t2d_v2g_connections[t2d_v2g_connections$celltype %in% c("hMSC_Adipocytes","hMSC_osteoblasts")==F, ]

t2d_v2g_connections = split(t2d_v2g_connections, t2d_v2g_connections$celltype)

t2d_v2g_connections_upset = lapply(t2d_v2g_connections, function(x){
	x$implicated_gene_name
	})

names(t2d_v2g_connections_upset) = names(t2d_v2g_connections)

gwas.matrix = list_to_matrix(t2d_v2g_connections_upset)
gwas.matrix = make_comb_mat(gwas.matrix)

pdf("plots/Upset_T2D_GWAS_genes_no_hMSC.pdf", width=12, height=4,)
UpSet(gwas.matrix,
	#set_order = c("AaM","Bipolar"BMI", "HEIGHT", "MDD", "SLEEP"),
	pt_size = unit(5, "mm"), lwd = 1,
	comb_col = c("red", "blue", "black", "darkgreen", "orange", "purple")[comb_degree(gwas.matrix)],
	comb_order = order(comb_size(gwas.matrix)),
	#set_on_rows
	)
dev.off()


