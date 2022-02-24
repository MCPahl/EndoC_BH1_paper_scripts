dir="/mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome_TOPMED_hg38/chicago" # change (chicago run parent directory)
chicago_run="hMSC_Adipocytes_2" # change (your chicago run)
atacSeq_file="/mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome_TOPMED_hg38/chicago/EndoC_BH1/dummy_atac.bed" # change (optional - better use your corresponding atacseq bed file)
hicup_dir="/mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome_TOPMED_hg38" # change (hicup parent dir)
hicup_names=("hMSC_Adipocyte_7078" "hMSC_Adipocyte_8002" "hMSC_Adipocyte_8004") # change (hicup replicate dir name)

mkdir $dir/$chicago_run
mkdir $dir/$chicago_run/1frag
mkdir $dir/$chicago_run/4frag

######### FEATURES PREP
# MAKE SIMBOYLIC OF ATAC-SEQ CONSERVATIVE PEAK TO CURRENT DIR (optional)
cd $dir/$chicago_run
echo -e "DUMMY\t$atacSeq_file" > feature.list # if use your own ATACseq, can change DUMMY to ATAC

#### 1FRAG CHICAGO
cd $dir/$chicago_run/1frag/
# MAKE SIMBOYLIC LINK OF HICUP CHICAGO INPUT FILE TO CURRENT DIR
for hicup_name in ${hicup_names[@]}; do
        file=$(ls $hicup_dir/$hicup_name/hicup/*_1frag/*_1frag.chinput)
        cp $file .
done


Rscript ~/SFGI_scripts_mp/R/remove_non1to22XY_chr_fragments_from_chinput.R

# RUN 1FRAG CHICAGO
files=$(ls $dir/$chicago_run/1frag/*_filtered.chinput | tr "\n" "," | sed "s/,$//")

sbatch --mem-per-cpu 80G -c 4 \
~/captureC/scripts/bash/runChicago.sh \
/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_1frag_design2/settingsFile.txt \
/mnt/isilon/sfgi/pahlm/annotationFiles/TOPMED_hg38_references/chicago_baitmap/chicago_baitmap_22XY/chicago/1frag \
$dir/$chicago_run/feature.list \
$files chicagoRes



#### 4FRAG CHICAGO
cd $dir/$chicago_run/4frag
# MAKE SIMBOYLIC LINK OF HICUP CHICAGO INPUT FILE TO CURRENT DIR
for hicup_name in ${hicup_names[@]}; do
        file=$(ls $hicup_dir/$hicup_name/hicup/*_4frag/*_4frag.chinput)
        cp $file .
done

# RUN 4FRAG CHICAGO
Rscript ~/SFGI_scripts_mp/R/remove_non1to22XY_chr_fragments_from_chinput.R

# RUN 4FRAG CHICAGO
files=$(ls $dir/$chicago_run/4frag/*_filtered.chinput | tr "\n" "," | sed "s/,$//")

sbatch --mem-per-cpu 60G -c 4 \
~/captureC/scripts/bash/runChicago.sh \
/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_4frag_design2/settingsFile.txt \
/mnt/isilon/sfgi/pahlm/annotationFiles/TOPMED_hg38_references/chicago_baitmap/chicago_baitmap_22XY/chicago/4frag \
$dir/$chicago_run/feature.list \
$files chicagoRes


cd $dir/$chicago_run/1frag
sbatch --mem 32G ~/captureC/scripts/bash/exportIbed.sh


cd $dir/$chicago_run/4frag
sbatch --mem 32G ~/captureC/scripts/bash/exportIbed.sh

dir="/mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome_TOPMED_hg38/chicago" # change (chicago run parent directory)
chicago_run="SGBS_Undiff_2" # change (your chicago run)
atacSeq_file="/mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome_TOPMED_hg38/chicago/EndoC_BH1/dummy_atac.bed" # change (optional - better use your corresponding atacseq bed file)
hicup_dir="/mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome_TOPMED_hg38" # change (hicup parent dir)
hicup_names=("SGBS_Undiff_rep1" "SGBS_Undiff_rep2" "SGBS_Undiff_rep3") # change (hicup replicate dir name)

