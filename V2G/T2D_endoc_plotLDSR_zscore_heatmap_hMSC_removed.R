#promoter OCR and 
#
setwd("/mnt/isilon/sfgi/pahlm/projects/T2D_endocBH1/ldsc")
library(ggplot2)
files = list.files("out")
files= files[grepl(".results", files)]

pldsr = lapply(files,function(file){
	x = read.delim(paste0("out/",file))
	x$Category = gsub("_regElements.results","", file)
	x = x[1,]
	x
})
pldsr = do.call("rbind", pldsr)

pldsr$Category = gsub("_rsid_added_summary_stats.txt","",gsub(".results", "", pldsr$Category))

pldsr$celltype <- ifelse(grepl("EndoC_BH1", pldsr$Category), "EndoC_BH1",
    ifelse(grepl("HEPG2",pldsr$Category), "HEPG2",
    ifelse(grepl("hMSC_Adipocyte",pldsr$Category), "hMSC_Adipocyte",
    ifelse(grepl("hMSC_BMP2", pldsr$Category), "hMSC_BMP2",
    ifelse(grepl("iPSC_Hepatocytes", pldsr$Category), "iPSC_Hepatocytes",
    ifelse(grepl("SGBS_Diff", pldsr$Category), "SGBS_Diff", "SGBS_Undiff"
))))))

pldsr = pldsr[pldsr$celltype %in% c("hMSC_BMP2", "hMSC_Adipocyte", "iPSC_Hepatocytes")==F,]

pldsr$fdr = p.adjust(pldsr$Enrichment_p, method="BH")
pldsr$significant = ifelse(pldsr$fdr <0.05, "Sig", ifelse(pldsr$Enrichment_p < 0.05,"nominal", "NS"))


pldsr$trait = gsub("SGBS_Diff_", "",
    gsub("SGBS_Undiff_","",
    gsub("EndoC_BH1_", "", 
    gsub("HEPG2_", "", 
    pldsr$Category))))

pldsr$trait = factor(pldsr$trait, levels= (unique(pldsr$trait)))
pldsr$celltype = factor(pldsr$celltype, levels = (unique(pldsr$celltype)))

pldsr$z = -qnorm(pldsr$Enrichment_p)


write.csv(pldsr, file = "../tables/pldsr_output_MSC_dervived_cells_removed.csv",quote=F, row.names=F )

pdf("plots/ldsc_zscore_heatmap_MSC_dervived_cells_removed.pdf", useDingbats=FALSE)
ggplot(pldsr, aes(x=celltype, y=trait, fill=z))+
 geom_tile(color="black")+
 scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000")+
 #geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error*1.96, ymax=Enrichment+Enrichment_std_error*1.96), width=0)+
 #scale_colour_manual(values = c("orange","black","red"))+
 #geom_hline(yintercept = 1, linetype="dotted", color="blue", size=1)+
 #coord_flip()+
 theme_minimal()+
 theme(axis.title.y = element_blank(),
  axis.ticks.y = element_blank(),
  legend.position = "top")+
  scale_y_discrete(limits = rev(levels(pldsr$trait)))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
 dev.off()




pdf("/plots/", useDingbats=FALSE)
ggplot(pldsr, aes(x=Category, y=Enrichment, col=significant, fill=significant))+
 geom_point(size=5)+
 geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error*1.96, ymax=Enrichment+Enrichment_std_error*1.96), width=0)+
 scale_colour_manual(values = c("orange","black","red"))+
 geom_hline(yintercept = 1, linetype="dotted", color="blue", size=1)+
 coord_flip()+
 theme_minimal()+
 theme(axis.title.y = element_blank(),
  axis.ticks.y = element_blank(),
  legend.position = "top")+
  scale_x_discrete(limits = rev(levels(pldsr$Category)))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
 dev.off()




pdf("../plots/pldsr_enrichment_results.pdf", useDingbats=FALSE)
ggplot(pldsr, aes(x=Category, y=Enrichment, col=significant, fill=significant))+
 geom_point(size=5)+
 geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error*1.96, ymax=Enrichment+Enrichment_std_error*1.96), width=0)+
 scale_colour_manual(values = c("orange","black","red"))+
 geom_hline(yintercept = 1, linetype="dotted", color="blue", size=1)+
 coord_flip()+
 theme_minimal()+
 theme(axis.title.y = element_blank(),
  axis.ticks.y = element_blank(),
  legend.position = "top")+
  scale_x_discrete(limits = rev(levels(pldsr$Category)))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
 dev.off()

##################################################################################################

#promoter OCR and covid

library(ggplot2)
files = list.files()
files= files[grepl("allOCRs.results", files)]

pldsr = lapply(files,function(file){
	x = read.delim(file)
	x$Category = gsub("_allOCRs.results","", file)
	x = x[1,]
	x
})
pldsr = do.call("rbind", pldsr)

pldsr$Category =factor(pldsr$Category, levels=c("hESC", "Monocyte", "NaiveB", "GCB", "NaiveT", "TFH"))
pldsr$fdr = p.adjust(pldsr$Enrichment_p, method="BH")
pldsr$significant = ifelse(pldsr$fdr <0.05, "Sig", ifelse(pldsr$Enrichment_p < 0.05,"nominal", "NS"))
pdf("../../plots/pldsr_enrichment_results_allOCRs.pdf", useDingbats=FALSE)
ggplot(pldsr, aes(x=Category, y=Enrichment, col=significant, fill=significant))+
 geom_point(size=5)+
 geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error*1.96, ymax=Enrichment+Enrichment_std_error*1.96), width=0)+
 scale_colour_manual(values = c("orange","black","red"))+
 geom_hline(yintercept = 1, linetype="dotted", color="blue", size=1)+
 coord_flip()+
 theme_minimal()+
 theme(axis.title.y = element_blank(),
  axis.ticks.y = element_blank(),
  legend.position = "top")+
  scale_x_discrete(limits = rev(levels(pldsr$Category)))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
 dev.off()