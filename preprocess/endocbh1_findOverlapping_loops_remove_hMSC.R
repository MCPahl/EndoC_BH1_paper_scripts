library(GenomicFeatures)
library(tidyverse)
library(data.table)
setwd("/mnt/isilon/sfgi/pahlm/projects/T2D_endocBH1")

dir = "data/pr_captureC/"
files = list.files(dir)
files = files[grepl(".ibed", files)]
loops = lapply(files, function(file){
	path=paste0(dir, "/", file)
	dat <- fread(path)
	dat$class = gsub(".ibed", "", file)
	dat
})

loops = do.call("rbind", loops)

loops_simple = loops %>% 
	mutate(bait = paste0(bait_chr, ":", bait_start, "-", bait_end),
		otherEnd = paste0(otherEnd_chr, ":", otherEnd_start, "-", otherEnd_end),
		loop_coords = paste0(bait_chr, ":", bait_start, "-", bait_end,"~",otherEnd_chr, ":", otherEnd_start, "-", otherEnd_end),
		celltype = gsub("_.frag", "", class)) %>%
	select(bait, otherEnd, loop_coords, celltype) %>%
	filter(celltype %in% c("hMSC_Adipocytes", "hMSC_osteoblasts") == F) %>% 
	unique()


celltypes = unique(loops_simple$celltype)

celltype_loop_jaccard <- lapply(celltypes, function(celltype_i){
	set_i <- lapply(celltypes, function(celltype_j){
		loop_set_i = loops_simple[loops_simple$celltype == celltype_i,]$loop_coords
		loop_set_j = loops_simple[loops_simple$celltype == celltype_j,]$loop_coords
		loop_jaccard = length(intersect(loop_set_i, loop_set_j))/length(union(loop_set_i, loop_set_j))
		data.frame(celltype_i, celltype_j, loop_jaccard )
	})
	do.call("rbind", set_i)
})

celltype_loop_jaccard = do.call("rbind", celltype_loop_jaccard)

celltype_loop_jaccard$comp = apply(celltype_loop_jaccard ,1, function(x){
	paste(sort(t(c(x[1], x[2]))), collapse=":")
	})

dups = base::duplicated(celltype_loop_jaccard$comp)

celltype_loop_jaccard = celltype_loop_jaccard[dups==F, ]

pdf("plots/loop_by_celltype_comparison_heatmap_hMSC_removed.pdf")
ggplot(celltype_loop_jaccard, aes(x = celltype_i, y = factor(celltype_j, levels=rev(unique(celltype_j))), fill = loop_jaccard))+
geom_tile(color="black")+
scale_fill_gradient(low = "white", high = "red") +
geom_text(aes(label = round(loop_jaccard, digits=2)), color = "black", size = 4) +
theme_minimal()+
guides(fill = guide_colourbar(barwidth = 0.5,
                                barheight = 20))+
theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
theme(aspect.ratio=1)
dev.off()




gwas.matrix = list_to_matrix()
gwas.matrix = make_comb_mat(gwas.matrix)

pdf("plots/Upset_T2D_GWAS_genes.pdf", width=12, height=7,)
UpSet(gwas.matrix,
	#set_order = c("AaM","Bipolar"BMI", "HEIGHT", "MDD", "SLEEP"),
	pt_size = unit(5, "mm"), lwd = 1,
	comb_col = c("red", "blue", "black", "darkgreen", "orange", "purple")[comb_degree(gwas.matrix)],
	comb_order = order(comb_size(gwas.matrix)),
	#set_on_rows
	)
dev.off()





library(ComplexHeatmap)



loop_simple_sets =loops_simple %>% select(loop_coords,celltype)
loop_simple_sets = split(loop_simple_sets, loop_simple_sets$celltype)
loop_simple_set = lapply(loop_simple_sets, function(loop_simple_set){
	loop_simple_set$loop_coords
})


loop_matrix = list_to_matrix(loop_simple_set )
loop_matrix = make_comb_mat(loop_matrix)

pdf("plots/Upset_loop_celltype.pdf", width=12, height=7,)
UpSet(loop_matrix,
	#set_order = c("AaM","Bipolar"BMI", "HEIGHT", "MDD", "SLEEP"),
	pt_size = unit(5, "mm"), lwd = 1,
	#comb_col = c("red", "blue", "black", "darkgreen", "orange", "purple")[comb_degree(gwas.matrix)],
	comb_order = order(comb_size(loop_matrix)),
	#set_on_rows
	)
dev.off()


