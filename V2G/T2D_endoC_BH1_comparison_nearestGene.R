#R get hg19 positions

library(tidyverse)

#Add nearest gene
Voight = read.csv("../T2D_v2g/leadSNPs/leadSNPs_T2D_Voight.csv")
voight.sent.positions = data.frame(sentinel = Voight$rsid, chr = paste("chr", Voight$chr,sep=""), start = Voight$position, end = Voight$position)
MCCarthy= read.csv("../T2D_v2g/leadSNPs/leadSNPs_T2D_McCarthy.csv")
mccarthy.sent.positions = data.frame(sentinel = MCCarthy$rsid, chr = paste("chr", MCCarthy$chr,sep=""), start = as.integer(gsub(",", "", MCCarthy$position)), end = as.integer(gsub(",", "", MCCarthy$position)))
sent_positions = unique(rbind(voight.sent.positions, mccarthy.sent.positions))
sent_positions$start= sent_positions$start-1
sent_positions = sent_positions[order(sent_positions$chr, sent_positions$start, sent_positions$end),]
sent_positions = sent_positions[,c(2,3,4,1)]

write.table(sent_positions, file="gwas_data/sentinel_snpset/sentinels_hg19.bed", quote=F,row.name=F, col.names=F, sep="\t")


#Liftover bash hg19->hg38
cd /mnt/isilon/sfgi/pahlm/projects/T2D_endocBH1/gwas_data/sentinel_snpset
liftOver \
<(sort-bed sentinels_hg19.bed) \
/mnt/isilon/sfgi/crossMapChainFiles/hg19ToHg38.over.chain.gz sentinels_hg38.bed sentinels.unmapped


#Map nearest_gene
library(tidyverse)
library(GenomicRanges)
sent_liftover = read.delim("sentinels_hg38.bed",header=F)
names(sent_liftover ) = c("chr", "start", "end", "sentinel")
sent_liftover$start = sent_liftover$start+1

load("/mnt/isilon/sfgi/pahlm/annotationFiles/genecodeV30//gene_tss.Rdata")


gene_tss$start = gene_tss$pos
gene_tss$end = gene_tss$pos


gene_tss = gene_tss[gene_tss$gene_type == "protein_coding" | gene_tss$gene_type == "noncoding_RNA_long",]

nearest = gene_tss[nearest(GRanges(sent_liftover), GRanges(gene_tss)),]


sent_liftover = data.frame(sent_liftover, nearest_gene_name = nearest$gene_name, nearest_gene_id = nearest$gene_id, nearest_gene_type = nearest$gene_type)

write.csv(sent_liftover, file= "sentinels_hg38_nearestgene.csv", quote=F, row.names=F)


library(tidyverse)
load("v2g/T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_simplified_22XY.RData"
nearest = read.csv("gwas_data/sentinel_snpset/sentinels_hg38_nearestgene.csv")

comp = left_join(t2d_v2g_simplified, nearest)

comp_sent_stats = comp %>% group_by(celltype, sentinel) %>% 
summarize(nearest_implicated = Reduce("|", implicated_gene_id %in% nearest_gene_id), implicated_gene_number= length(unique(implicated_gene_id))) %>%
ungroup() %>%
mutate(class = ifelse(implicated_gene_number == 1 & nearest_implicated == TRUE, "nearest_only", 
                    ifelse(implicated_gene_number > 1 & nearest_implicated == TRUE, "multiple_including_nearest", "not_inc_nearest"))) %>%
select(celltype, class) %>%
group_by(celltype, class) %>%
tally()

pdf("plots/sentinel_nearest_gene_comparison.pdf", useDingbats=FALSE, width =3, height=2)
ggplot(comp_sent_stats, aes(x= celltype, y = n, fill = class))+
geom_bar(position="stack", stat="identity") +
coord_flip()+
theme_minimal()+
scale_x_discrete(limits = rev(levels(comp_sent_stats$celltype)))+
theme(legend.position="top",
    panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()