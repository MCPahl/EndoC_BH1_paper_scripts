library("tidyverse")

load("data/motifbreakR_results.RData")
load("gwas_data/proxy_snpset/T2D/T2D_proxies_ldlink_hg38.RData")

#Add sentinel info back in
MCCarthy= read.csv("/mnt/isilon/sfgi/pahlm/projects/T2D_v2g/leadSNPs/leadSNPs_T2D_McCarthy.csv")
MCCarthy = MCCarthy[,c(1,7,8,9)]
names(MCCarthy)[3] = "Study"

Voight = read.csv("/mnt/isilon/sfgi/pahlm/projects/T2D_v2g/leadSNPs/leadSNPs_T2D_Voight.csv")
Voight = Voight[,c(2,8,9,10)]
sentinels = rbind(MCCarthy,Voight)
names(sentinels)= c("sentinel", "GWAS_sentinel_P", "study", "ancestry")

sentinels$ancestry = paste(sentinels$study,sentinels$ancestry,sep="_")
sentinels$study=NULL
sentinels$ancestry = gsub("Voight_TRANSE", "Voight_transEth", sentinels$ancestry)
EndoC_BH1_16358

tmp = sentinels[is.na(as.numeric(sentinels$GWAS_sentinel_P)),]$GWAS_sentinel_P
sentinels$GWAS_sentinel_P <- as.numeric(sentinels$GWAS_sentinel_P)
tmp = c(1.e-8, 1.8e-8, 3.7e-8, 5.3e-9, 8.2e-10, 2.3e-11, 5.7e-14, 1.1e-9, 8.5e-9, 2.7e-15, 1.1e-8, 5.4e-11, 5.5e-16, 1.3e-9,2.8e-8, 5.6e-60, 2.2e-9, 2.2e-9, 3.3e-12, 4.2e-9, 6.1e-13 )
sentinels$GWAS_sentinel_P[is.na(sentinels$GWAS_sentinel_P)] <- tmp


T2D_proxies_ldlink_hg38 = left_join(T2D_proxies_ldlink_hg38, sentinels)


motifbreakr.results = left_join(motifbreakr.results,

motifbreakr.results.simplified = motifbreakr.results %>% select(proxy, geneSymbol, providerName, alleleDiff, sentinel, ancestry, GWAS_sentinel_P)


 T2D_proxies_ldlink_hg38)


pdf("plots/TF_motifBreakR_vs_GWAS_Pvalue.pdf")
ggplot(motifbreakr.results.simplified %>% select(GWAS_sentinel_P, alleleDiff) %>% unique(), aes(x= -log10(GWAS_sentinel_P), y = abs(alleleDiff)))+
	geom_point()
dev.off()


load("v2g/T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight.RData")
t2d_v2g  = t2d_v2g %>% separate_rows(proxy, sep="\\|")

left_join(t2d_v2g, )