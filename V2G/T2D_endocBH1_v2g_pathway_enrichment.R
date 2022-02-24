library(tidyverse)

load("/mnt/isilon/sfgi/pahlm/annotationFiles/msigdb_v7.0_GMTs/c2.cp.reactome.v7.0.symbols.Rdata")

c2.cp.reactome = lapply(c2.cp.reactome, function(x){ x[grepl("HIST",x)==F] })

HG.test = function(grp,anno){
  universe = unique(unlist(anno)) #Set of genes that have been annotated, number of balls 
  gset = grp
  gset = gset[gset %in% universe]
  tmp = lapply(anno,function(j){
      p = length(j[j %in% gset]) #white balls drawn from an urn
      m = length(j) #white balls in the urn
      n = length(universe) - m #Number of black balls in the urn
      k = length(gset) #Number of balls drawn from the urn
      expected = ceiling(length(j)/length(universe)*length(gset)) #Number of white balls expected to be drawn, rounding up (replace ceiling with floor to round down, or remove to leave the decimal)
      hg.pval = phyper(p,m,n,k,lower.tail=FALSE) #Calculate p value
      genes = paste(j[j %in% gset], collapse=":") #Combine the list of genes
      data.frame(p.val = hg.pval,genes= genes, observed=p, expected=expected) #Put the data together in a data.frame
    })
   out = do.call("rbind",tmp) #Merge the list of dataframes into a single dataframe
   out$Pathway = names(tmp) #Convert row.names to a column
   row.names(out) = NULL #ls() the row.names
   out$p.adjust = p.adjust(out$p.val,method="fdr") #Calculate p.value after FDR correction
   out #Print results
}

load("v2g/T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_simplified_22XY.RData")
load("v2g/T2D_v2g_mapping_TOPMED_Hg38_McCarthy_Voight_promoterSNPs_22XY.Rdata")

promoter_ocr_v2g$celltype[promoter_ocr_v2g$celltype=="hMSC_Osteoblasts"] = "hMSC_osteoblasts"

enriched_pathways = lapply(unique(promoter_ocr_v2g$celltype), function(celltype){
genes = unique(c(promoter_ocr_v2g[promoter_ocr_v2g$celltype==celltype,]$gene_name, t2d_v2g_simplified[t2d_v2g_simplified$celltype==celltype,]$implicated_gene_name))
geneset.hg <- HG.test(genes, c2.cp.reactome)
  geneset.hg = geneset.hg[order(geneset.hg$p.val),]
  group = geneset.hg
  #group = geneset.hg[geneset.hg$p.adjust<0.05 & geneset.hg$observed>1,]
  group$Pathway = factor(group$Pathway,levels=group$Pathway)
  group$celltype = celltype
  group
})
enriched_pathways = do.call("rbind", enriched_pathways)

enriched_pathways_plot = enriched_pathways %>% 
  mutate(logPValue = -log10(p.val)) %>%
  group_by(celltype) %>%
  mutate(rank = rank(p.val, ties.method = "first")) %>% 
  filter(rank <= 10) %>% 
  mutate(Pathway = fct_reorder(paste(celltype, Pathway,sep=":"), rank))



#enriched_pathways = enriched_pathways[enriched_pathways$p.adjust < 0.05,]

pdf("plots/reactome_enrichment.pdf", useDingbats=FALSE)
ggplot(enriched_pathways_plot, aes(x= Pathway))+
  geom_bar(aes(y = observed/expected), stat="identity", fill = "lightblue")+
  geom_point(aes(y = logPValue), color = "orange")+
  theme_minimal()+
  theme(aspect.ratio=1)+
    theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
    facet_wrap(vars(celltype), scales= "free", nrow= 2)
dev.off()