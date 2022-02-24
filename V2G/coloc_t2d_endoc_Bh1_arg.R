library(data.table)
library(tidyverse)
library(coloc)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=1) {
  stop("One argument (index) must be supplied", call.=FALSE)
}

index = args[1]

#helper function
run_coloc = function(gwas_roi, eqtl_roi, gene){
	merged_ref = inner_join(gwas_roi, eqtl_roi)
	merged_ref = merged_ref[(merged_ref$frequency > 0.1 &  merged_ref$frequency < 0.9 ), ]
	if(nrow(merged_ref)> 0 ){
	coloc_out = coloc.abf(
		dataset1 = list(
			pvalues = merged_ref$pvalue, 
			beta = (merged_ref$effect_size), 
			varbeta=(merged_ref$standard_error)^2,
			N=merged_ref$sample_size, 
			type="cc",
			s = unique(merged_ref$n_cases/merged_ref$sample_size),
			MAF = as.numeric(merged_ref$frequency)
		), 
		dataset2= list(pvalues = merged_ref$eqtl_pval_nominal,
		  #beta=merged_ref$eqtl_slope,
		  varbeta= (merged_ref$eqtl_slope_se)^2,
		   N=merged_ref$eqtl_ma_samples,
		   type="quant",
		   MAF=merged_ref$eqtl_maf
		 )
	)
	coloc.out = data.frame(leadSNP = unique(gwas_roi$leadSNP), gene_id = gene, eqtl_tissue = gsub(".allpairs.txt.gz", "", gsub("/mnt/isilon/sfgi/pahlm/datasets_downloaded/GTeX/eqtl_v7/GTEx_Analysis_v7_eQTL_all_associations/", "", eqtl_path)), rbind(coloc_out$summary))
}
}

#Give lead SNPs (or snps to check for coloc) and summary stats
leadSNPs = fread("/mnt/isilon/sfgi/pahlm/projects/T2D_v2g/leadSNPs/leadSNPs_T2D_McCarthy.csv")
sumstats = fread("/mnt/isilon/sfgi/pahlm/projects/T2D_endocBH1/gwas_data/harmonize/Mahajan.NatGenet2018b.T2D.European.harmonized.txt")


#Get the path to the eqtl files
eqtl_dir = "/mnt/isilon/sfgi/pahlm/datasets_downloaded/GTeX/eqtl_v7/GTEx_Analysis_v7_eQTL_all_associations"
eqtl_files = list.files(eqtl_dir)
eqtl_paths = paste(eqtl_dir, eqtl_files, sep = "/")

#Identify and subset the summary stats by regions around the lead SNP (1MB window centered at each lead SNP)
leadSNPs$position = as.numeric(gsub(",", "", leadSNPs$position))


regions_of_interest = lapply(leadSNPs$rsid,  function(leadSNP_rsid){
	leadSNP_panel = sumstats[sumstats$variant_id == leadSNP_rsid,]
	leadSNP_panel_of_interest = sumstats[sumstats$chromosome == leadSNP_panel$chromosome & sumstats$position > leadSNP_panel$position-500000 & sumstats$position < leadSNP_panel$position+500000,]
	leadSNP_panel_of_interest$leadSNP = leadSNP_rsid
	leadSNP_panel_of_interest
})

regions_of_interest = do.call("rbind",regions_of_interest)
regions_of_interest = split(regions_of_interest, regions_of_interest$leadSNP)


#Loop through each tissue, get snps within region to coloc
eqtl_path = eqtl_paths[index]
tissue = gsub(".allpairs.txt.gz", "", gsub("/mnt/isilon/sfgi/pahlm/datasets_downloaded/GTeX/eqtl_v7/GTEx_Analysis_v7_eQTL_all_associations/", "", eqtl_path)
	eqtl_data = fread(eqtl_path)
	#Loop through each region to get SNPs within region
	coloc_out_genes_signals= lapply(regions_of_interest, function(region_of_interest){
		eqtl_region_data = eqtl_data[eqtl_data$variant_id %in% region_of_interest$panel_variant_id,]
		eqtl_region_data = eqtl_region_data[complete.cases(eqtl_region_data),]
		if(length(eqtl_region_data$gene_id)>0){
			gene_test = unique(eqtl_region_data$gene_id)
			#subset and rename eqtl genes, loop through each gene
			coloc_out_genes = lapply(gene_test, function(gene){
			eqtl_region_gene_data =	 eqtl_region_data %>% 
			filter(gene_id == gene) %>%
			select(gene_id,
				 variant_id,
				 tss_distance,
				 ma_samples,
				 maf,
				 pval_nominal,
				 slope,
				 slope_se) %>% 
			rename(panel_variant_id = variant_id,
				  eqtl_tss_distance = tss_distance,
				  eqtl_ma_samples = ma_samples, 
				  eqtl_maf = maf, 
				  eqtl_pval_nominal = pval_nominal, 
				  eqtl_slope = slope, 
				  eqtl_slope_se = slope_se)
			print(paste(unique(region_of_interest$leadSNP), gene, sep="_"))
			coloc_out = run_coloc(region_of_interest, eqtl_region_gene_data, gene)
			})
		coloc_out_genes = do.call("rbind", coloc_out_genes )
		}
})
coloc_out_genes_signals = do.call("rbind",tissue, coloc_out_genes_signals)

save(coloc_out_genes_signals, file = paste("data/coloc", "_coloc_results.RData", sep=""))




