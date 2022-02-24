#Intersect ATAC-seq with CaptureC (build reference for V2G)
#Use GenomicInteractions Package
#Selected interactions where otherEnd overlaps with ATAC-seq reproducible peak or bait overlaps with ATAC-seq peak and other end is annotated (open bait 2 baits)

library(GenomicInteractions)
library(data.table)

setwd("/mnt/isilon/sfgi/pahlm/projects/T2D_endocBH1/data")

#helper functions
chicago_ibed2genomicInteraction = function(ibed){
	ancor1=GRanges(seqnames=ibed$bait_chr, ranges=IRanges(ibed$bait_start,ibed$bait_end))
	ancor2=GRanges(seqnames=ibed$otherEnd_chr, ranges=IRanges(ibed$otherEnd_start,ibed$otherEnd_end))
	ibed.gi = GenomicInteractions::GenomicInteractions(ancor1, ancor2)
	ibed.gi$bait_name = ibed$bait_name
	ibed.gi$otherEnd_name = ibed$otherEnd_name
	ibed.gi$counts = ibed$N_reads
	ibed.gi$score = ibed$score
	ibed.gi$resolution = ibed$resolution
	ibed.gi$celltype = ibed$celltype
	ibed.gi
}

chiapet_micc_bedpe2genomicInteraction = function(chia){
	names(chia)[1:6] = c("chrA", "startA", "endA", "chrB", "startB", "endB")
	ancor1=GRanges(seqnames=chia$chrA, ranges=IRanges(chia$startA+1, chia$endA))
	ancor2=GRanges(seqnames=chia$chrB, ranges=IRanges(chia$startB+1, chia$endB))
	chia.gi = GenomicInteractions::GenomicInteractions(ancor1, ancor2)
	chia.gi$peakA = chia$peakAr
	chia.gi$peakB = chia$peakB
	chia.gi$cA = chia$cA
	chia.gi$cB = chia$cB
	chia.gi$cAB = chia$cAB
	chia.gi$prob = chia$"-log10(1-PostProb)"
	chia.gi$fdr = chia$fdr
	chia.gi
}

#load capture c and atac data

pr_capturec_files = list.files("pr_captureC")
pr_capturec_files = pr_capturec_files[grepl("ibed", pr_capturec_files)] 

atac_files = list.files("atac/atac_peaks/")

celltypes = unique(gsub("_.frag.ibed", "", pr_capturec_files))

celltypes = celltypes[celltypes!="iHEP"]

#load Capture C and ATAC-seq data
capturec_dat = lapply(celltypes, function(celltype){
	pr_capturec_file_frags = pr_capturec_files[grepl(celltype, pr_capturec_files)]
	pr_capturec_1frag = fread(paste0("pr_captureC/", pr_capturec_file_frags[grepl("1frag", pr_capturec_file_frags)]))
	pr_capturec_1frag$resolution = "1frag"
	pr_capturec_4frag = fread(paste0("pr_captureC/", pr_capturec_file_frags[grepl("4frag", pr_capturec_file_frags)]))
	pr_capturec_4frag$resolution = "4frag"
	pr_capturec = rbind(pr_capturec_1frag, pr_capturec_4frag)
	pr_capturec$celltype = celltype
	pr_capturec
})
capturec_dat = do.call("rbind", capturec_dat)

atacseq_dat = lapply(celltypes, function(celltype){
	atac_file = atac_files[grepl(celltype, atac_files)]
	atac_dat = fread(paste0("atac/atac_peaks/", atac_file))
	names(atac_dat) = c("chr", "start", "end", "id", "score", "strand", "enrichment", "logP", "logQ", "center", "celltype", "rep1", "rep2", "rep3", "sum","ratio")
	atac_dat$celltype = celltype
	atac_dat = atac_dat[atac_dat$ratio>0.5,]
	atac_dat$celltype
	atac_dat$start = atac_dat$start+1
	atac_dat$center = atac_dat$center+1
	atac_dat[,c(1:4,11)]
})
atacseq_dat = do.call("rbind", atacseq_dat)

#intersect by side
intersection_by_celltype = lapply(celltypes, function(celltype){
	ct = celltype
	atac.cell = atacseq_dat[celltype == ct,]
	capt.cell = capturec_dat[capturec_dat$celltype == ct,]
	capt.gi = chicago_ibed2genomicInteraction(capt.cell)
	atac.gr = GRanges(atac.cell)
	index = as.data.frame(findOverlaps(capt.gi, atac.gr, use.region="second",))
	oeOpen = data.frame(capt.cell[index[,1],], atac.cell[index[,2],])

	index2 = as.data.frame(findOverlaps(capt.gi, atac.gr, use.region="first",))
	baitOpen = data.frame(capt.cell[index[,1],], atac.cell[index[,2],])
	baitOpen = baitOpen[baitOpen$otherEnd_name==F,]

	out = rbind(oeOpen, baitOpen)
	names(out) = c("bait_chr", "bait_start","bait_end", "bait_name", "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name", "N_reads", "score", "resolution", "celltype", "ocr_chr", "ocr_start", "ocr_end", "ocr_id", "ocr_celltype")
	out$ocr_celltype=NULL
	out
})

#prep for saving
dir.create("intersections")
dir.create("intersections/tables")

#save rdata
save(intersection_by_celltype, file = "intersections/capture_c_atac_seq_celltypes.RData")

#write output
lapply(intersection_by_celltype, function(intersection_dat){
	intersection_dat$start = intersection_dat$bait_start-1
	intersection_dat$end = intersection_dat$otherEnd_start-1
	intersection_dat$ocr_start = intersection_dat$otherEnd_start-1
	write.table(intersection_dat, file = paste0("intersections/", unique(intersection_dat$celltype), "_captureC_atac_seq_intersection.txt"), quote=FALSE, row.names=F, sep="\t")
})



