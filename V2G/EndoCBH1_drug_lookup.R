library("rDGIdb")

endoc_bh1_genes

result <- queryDGIdb(endoc_bh1_genes)

results_table = resultSummary(result)
results_table_detailed = detailedResults(result)

gene_results_table =  byGene(result)

write.table(results_table, file = "tables/v2g_all_drug_info_results_table.txt", sep="\t", quote=F, row.names=F)
write.table(results_table_detailed ,file = "tables/v2g_all_drug_info_results_table_details.txt", sep="\t", quote=F, row.names=F)
write.table(gene_results_table, file = "tables/v2g_all_drug_gene_results_table.txt", sep="\t", quote=F, row.names=F)
