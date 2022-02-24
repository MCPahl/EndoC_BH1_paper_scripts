library(data.table)
setwd("/mnt/isilon/sfgi/pahlm/gwas_summary_stats/common_metabolic_disease/PCOS")
gwas = fread("PCOS_summary_data_19092018.txt")
snps.anno = fread("/mnt/isilon/sfgi/commonSnps/hg19/snp150Common.bed")
names(gwas)[1] = "SNP"
gwas$SNP = gsub(":ID", "", gwas$SNP)
snps.anno$V1 = gsub("chr","", snps.anno$V1)
names(snps.anno) = c("chr", "start", "end", "rsid")

snps.anno$SNP = paste(snps.anno$chr, snps.anno$end,sep=":")

setkey(snps.anno,SNP)
setkey(gwas, SNP)

out = snps.anno[gwas]

 out$SNP = ifelse(is.na(out$rsid),out$SNP, out$rsid)
# out = out[,c(5:14)]
 out2 = out[complete.cases(out),]
write.table(out2, file="PCOS_Day2018_rsid_added_summary_stats.txt", quote=F, row.names=FALSE,sep="\t", col.names=TRUE)

