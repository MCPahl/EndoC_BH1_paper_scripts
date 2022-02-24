dir="/mnt/isilon/sfgi/pahlm/annotationFiles/TOPMED_hg38_references/chicago_baitmap"
cd $dir


ln -s /mnt/isilon/sfgi/pahlm/annotationFiles/TOPMED_hg38_references/Digest_hg38_TOPMED_DpnII_Arima_None_15-12-15_11-08-2021.txt TOPMED_hg38_Arima.txt
ln -s /mnt/isilon/sfgi/pahlm/annotationFiles/TOPMED_hg38_references/Digest_hg38_TOPMED_DpnII_None_15-08-12_11-08-2021.txt TOPMED_hg38_DpnII.txt
chgrp sfgi-admin TOPMED_hg38_Arima.txt
chgrp sfgi-admin TOPMED_hg38_DpnII.txt
chmod 775 TOPMED_hg38_DpnII.txt
chmod 775 TOPMED_hg38_Arima.txt

##### CREATE rmap
sed 1,2d /mnt/isilon/sfgi/programs/hicup_v0.7.4/digested_genomes/TOPMED_hg38_DpnII.txt| cut -f 1-3 | awk 'BEGIN{OFS="\t"; FS="\t"}{print $0,NR}' >  chicago_baitmap_22XY/DpnII_1frag_TOPMED_hg38.rmap
cat  chicago_baitmap_22XY/DpnII_1frag_TOPMED_hg38.rmap|grep -v chrM| grep -v chrUn | grep -v ERCC | grep -v chrEBV | grep -v random > chicago_baitmap_22XY/DpnII_1frag_TOPMED_hg38.rmap.tmp
mv chicago_baitmap_22XY/DpnII_1frag_TOPMED_hg38.rmap.tmp chicago_baitmap_22XY/DpnII_1frag_TOPMED_hg38.rmap
cd chicago_baitmap_22XY

##### CREATE baitmap
## liftOver probes.bed
liftOver \
<(sort-bed /mnt/isilon/sfgi/suc1/analyses/wells/captureC/atypical_promoterome/ProbesFromAle/parse/probes.bed) \
/mnt/isilon/sfgi/crossMapChainFiles/hg19ToHg38.over.chain.gz probes_hg38.bed probes.unmapped

wc -l probes.unmapped
# 34 probes.unmapped (17 probes not mapped)

wc -l probes_hg38.bed
# 127456 probes_hg38.bed

ln -s ../probes_hg38.filter.bed .

bedtools intersect \
-a probes_hg38.filter.bed \
-b <(awk 'BEGIN{OFS="\t"; FS="\t"}{$2=$2-1; print $1,$2,$3,$4}' DpnII_1frag_TOPMED_hg38.rmap) \
-wao > DpnII_1frag_TOPMED_hg38.vs.probe.txt

cut -f 1-4 DpnII_1frag_TOPMED_hg38.vs.probe.txt | sort -u | wc -l
# 127448 probes mapped

#Run in command line
summarize_probe2bait.R

wc -l DpnII_1frag_hg38_baits.bed 
# 36687 baits

## annotate baits
bedtools intersect -a DpnII_1frag_hg38_baits.bed -b /mnt/isilon/sfgi/pahlm/annotationFiles/genecodeV30/gencode.v30.promoter.bed -wa -wb > DpnII_1frag_TOPMED_hg38_baits.vs.promoter.txt

###
closest-features --closest --dist --delim "\t" \
DpnII_1frag_hg38_baits.bed \
<(sed "/#/d" /mnt/isilon/sfgi/pahlm/annotationFiles/genecodeV30/gencode.v30.TSS.bed) \
> bait2closestTSS.txt

dir="/mnt/isilon/sfgi/pahlm/annotationFiles/TOPMED_hg38_references/chicago_baitmap"
bait_bed_file="$dir/chicago_baitmap_22XY/DpnII_1frag_hg38_baits.bed"
output="$dir/chicago_baitmap_22XY/DpnII_1frag_TOPMED_hg38_baits.baitmap"
promoter2bait_file="$dir/chicago_baitmap_22XY/DpnII_1frag_TOPMED_hg38_baits.vs.promoter.txt"
bait2closestTSS_file="/mnt/isilon/sfgi/pahlm/annotationFiles/TOPMED_hg38_references/chicago_baitmap/chicago_baitmap_22XY/bait2closestTSS.txt"
gene_tss_rdata="/mnt/isilon/sfgi/pahlm/annotationFiles/genecodeV30/gene_tss.Rdata"

Rscript ~/scripts/R/annotate_bait2gene.R $bait_bed_file $output $promoter2bait_file $bait2closestTSS_file $gene_tss_rdata

##### CREATE chicago reference
mkdir -p $dir/chicago_baitmap_22XY/chicago/1frag
cd $dir/chicago_baitmap_22XY/chicago/1frag
ln -s /mnt/isilon/sfgi/pahlm/annotationFiles/TOPMED_hg38_references/chicago_baitmap/chicago_baitmap_22XY/DpnII_1frag_TOPMED_hg38_baits.baitmap DpnII_1frag_TOPMED_hg38.baitmap
ln -s /mnt/isilon/sfgi/pahlm/annotationFiles/TOPMED_hg38_references/chicago_baitmap/chicago_baitmap_22XY/DpnII_1frag_TOPMED_hg38.rmap DpnII_1frag_TOPMED_hg38.rmap

echo "source ~/.bashrc
conda activate py27
python /mnt/isilon/sfgi/programs/CHiCAGO/chicagoTools/makeDesignFiles.py --baitmapfile=DpnII_1frag_TOPMED_hg38.baitmap --rmapfile=DpnII_1frag_TOPMED_hg38.rmap --binsize=2500 --outfilePrefix=chicago_1frag
" > makeDesignFiles.sh
sed -i '1i#!/bin/bash' makeDesignFiles.sh

sbatch --mem=20G -t 12:00:00 makeDesignFiles.sh

echo -e "binsize\t2500" > settingsFile.txt

############################ 4frag############################
mkdir -p $dir/chicago_baitmap_22XY/chicago/4frag
cd $dir/chicago_baitmap_22XY/chicago/4frag
##### CREATE rmap
perl ~/scripts/concatenate_neighbor_frag.pl 4 $dir/chicago_baitmap_22XY/DpnII_1frag_TOPMED_hg38.rmap  > DpnII_4frag_TOPMED_hg38.rmap

##### CREATE baitmap
## generate bait.bed
bedtools intersect \
-a <(awk 'BEGIN{OFS="\t"}{$2=$2-1; print $0}' DpnII_4frag_TOPMED_hg38.rmap) \
-b $dir/chicago_baitmap_22XY/DpnII_1frag_hg38_baits.bed \
-wa -u > $dir/chicago_baitmap_22XY/chicago/4frag/DpnII_4frag_TOPMED_hg38.bed

wc -l DpnII_4frag_TOPMED_hg38.bed
# 34339 DpnII_4frag_hg38_baits.bed

## annotate baits
bedtools intersect -a DpnII_4frag_TOPMED_hg38.bed -b /mnt/isilon/sfgi/pahlm/annotationFiles/genecodeV30/gencode.v30.promoter.bed -wa -wb > DpnII_4frag_TOPMED_hg38_baits.vs.promoter.txt


closest-features --closest --dist --delim "\t" \
DpnII_4frag_TOPMED_hg38.bed \
<(sed "/#/d" /mnt/isilon/sfgi/pahlm/annotationFiles/genecodeV30/gencode.v30.TSS.bed) \
> bait2closestTSS_4frag.txt

dir="/mnt/isilon/sfgi/pahlm/annotationFiles/TOPMED_hg38_references/chicago_baitmap/"
bait_bed_file="$dir/DpnII_4frag_TOPMED_hg38.bed"
output=$dir"chicago_baitmap_22XY/DpnII_4frag_TOPMED_hg38_baits.baitmap"
promoter2bait_file="$dir/chicago_baitmap_22XY/DpnII_4frag_TOPMED_hg38_baits.vs.promoter.txt"
bait2closestTSS_file="$dir/chicago_baitmap_22XY/chicago/4frag/bait2closestTSS_4frag.txt"
gene_tss_rdata="/mnt/isilon/sfgi/pahlm/annotationFiles/genecodeV30/gene_tss.Rdata"

Rscript ~/scripts/R/annotate_bait2gene.R $bait_bed_file $output $promoter2bait_file $bait2closestTSS_file $gene_tss_rdata

##### CREATE chicago reference
mkdir -p $dir/chicago
cd $dir/chicago
ln -s $dir/DpnII_4frag_TOPMED_hg38_baits.baitmap DpnII_4frag_TOPMED_hg38_baits.baitmap 
ln -s $dir/DpnII_4frag_hg38.rmap DpnII_4frag_hg38.rmap

echo "source ~/.bashrc
conda activate py27
python /mnt/isilon/sfgi/programs/CHiCAGO/chicagoTools/makeDesignFiles.py --baitmapfile=$dir/chicago_baitmap_22XY/chicago/4frag//DpnII_4frag_TOPMED_hg38_baits.baitmap  --rmapfile=$dir/chicago_baitmap_22XY/chicago/4frag/DpnII_4frag_TOPMED_hg38.rmap --binsize=10000 --outfilePrefix=chicago_4frag --removeAdjacent=False
" > makeDesignFiles_4frag.sh
sed -i '1i#!/bin/bash' makeDesignFiles_4frag.sh

sbatch --mem=20G -t 12:00:00 makeDesignFiles_4frag.sh

echo -e "binsize\t10000" > settingsFile.txt
echo -e "removeAdjacent\tFalse" >> settingsFile.txt

######################### ln -s to chicago_map
cd /mnt/isilon/sfgi/suc1/chicago_map
ln -s /mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV35/AleDpnII_1frag/DpnII_1frag_hg38.rmap hg38_DpnII.1frag.rmap
ln -s /mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV35/AleDpnII_1frag/DpnII_1frag_hg38_baits.baitmap hg38_DpnII_baits.1frag.baitmap

ln -s /mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV35/AleDpnII_4frag/DpnII_4frag_hg38.rmap hg38_DpnII.4frag.rmap
ln -s /mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV35/AleDpnII_4frag/DpnII_4frag_hg38_baits.baitmap hg38_DpnII_baits.4frag.baitmap


############################ 8frag############################
mkdir -p $dir/8frag
##### CREATE rmap
perl ~/scripts/concatenate_neighbor_frag.pl 8 $dir/DpnII_1frag_TOPMED_hg38.rmap  > DpnII_8frag_TOPMED_hg38.rmap

##### CREATE baitmap
## generate bait.bed
bedtools intersect \
-a <(awk 'BEGIN{OFS="\t"}{$2=$2-1; print $0}' DpnII_8frag_TOPMED_hg38.rmap) \
-b $dir/DpnII_1frag_hg38_baits.bed \
-wa -u > $dir/DpnII_frag_TOPMED_hg38.bed

wc -l DpnII_4frag_TOPMED_hg38.bed
#  32909 DpnII_8frag_TOPMED_hg38.bed

## annotate baits
bedtools intersect -a DpnII_8frag_TOPMED_hg38.bed -b /mnt/isilon/sfgi/pahlm/annotationFiles/genecodeV30/gencode.v30.promoter.bed -wa -wb > DpnII_8frag_TOPMED_hg38_baits.vs.promoter.txt


closest-features --closest --dist --delim "\t" \
DpnII_8frag_TOPMED_hg38.bed \
<(sed "/#/d" /mnt/isilon/sfgi/pahlm/annotationFiles/genecodeV30/gencode.v30.TSS.bed) \
> bait2closestTSS_8frag.txt

bait_bed_file="$dir/DpnII_8frag_TOPMED_hg38.bed"
output="$dir/DpnII_8frag_TOPMED_hg38_baits.baitmap"
promoter2bait_file="$dir/DpnII_8frag_TOPMED_hg38_baits.vs.promoter.txt"
bait2closestTSS_file="$dir/bait2closestTSS_8frag.txt"
gene_tss_rdata="/mnt/isilon/sfgi/pahlm/annotationFiles/genecodeV30/gene_tss.Rdata"

Rscript ~/scripts/R/annotate_bait2gene.R $bait_bed_file $output $promoter2bait_file $bait2closestTSS_file $gene_tss_rdata

##### CREATE chicago reference
mkdir -p $dir/chicago
cd $dir/chicago
ln -s $dir/DpnII_8frag_TOPMED_hg38_baits.baitmap DpnII_8frag_TOPMED_hg38_baits.baitmap 
ln -s $dir/DpnII_8frag_TOPMED_hg38.rmap DpnII_8frag_TOPMED_hg38.rmap

echo "module load python27
python /mnt/isilon/sfgi/programs/CHiCAGO/chicagoTools/makeDesignFiles.py --baitmapfile=$dir/chicago/DpnII_8frag_TOPMED_hg38_baits.baitmap  --rmapfile=$dir/chicago/DpnII_8frag_TOPMED_hg38.rmap --binsize=10000 --outfilePrefix=chicago_8frag --removeAdjacent=False
" > makeDesignFiles_8frag.sh

sbatch --mem=20G -t 6:00:00 makeDesignFiles_8frag.sh

echo -e "binsize\t10000" > settingsFile.txt
echo -e "removeAdjacent\tFalse" >> settingsFile.txt

######################### ln -s to chicago_map
cd /mnt/isilon/sfgi/suc1/chicago_map
ln -s /mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV35/AleDpnII_1frag/DpnII_1frag_hg38.rmap hg38_DpnII.1frag.rmap
ln -s /mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV35/AleDpnII_1frag/DpnII_1frag_hg38_baits.baitmap hg38_DpnII_baits.1frag.baitmap

ln -s /mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV35/AleDpnII_4frag/DpnII_4frag_hg38.rmap hg38_DpnII.4frag.rmap
ln -s /mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV35/AleDpnII_4frag/DpnII_4frag_hg38_baits.baitmap hg38_DpnII_baits.4frag.baitmap

