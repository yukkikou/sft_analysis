#!/bin/bash
cosmic="/media/db/genomes/Hsapiens/hg38_gatk4/Cosmic.v90.vcf.gz"
pon="/media/db/genomes/Hsapiens/hg38_gatk4/panel_of_normal.vcf.gz"
dbsnp="/media/db/genomes/Hsapiens/hg38_gatk4/dbSNP156.hg38.vcf.gz"
germline="/media/db/genomes/Hsapiens/hg38_gatk4/gnomad_WGS_v4.vcf.gz"
germline2="/media/db/genomes/Hsapiens/hg38_gatk4/gnomad_exomes_v4.vcf.gz"
sentieon="/home/panda/software/sentieon/bin/sentieon"
hg38="/media/db/genomes/Hsapiens/hg38_gatk4/bwa/Homo_sapiens_assembly38.fasta"
samtools="/home/panda/software/samtools-1.9/samtools"
tmp_path="/media/ssdtmp"

project=$1
tumor_bam=$2
#normal_bam=$3
NCPU=$3
if [[ "$#" -ne 3 ]]; then
	echo "$(tput bold; tput setaf 2)Usage: sentieon.TN project_name tumor_bam threads
	This script is use for somatic mutation filtering of tumor bam without normal pairing control
	only use hg38 as reference
	$(tput sgr0)"
	exit 1
fi

tumor_recal="$(dirname $tumor_bam)/recal_data.table"
#normal_recal="$(dirname $normal_bam)/recal_data.table"
mkdir -p $(pwd)/$project
rootdir=$(pwd)/$project
# disk reading is not the limiting step for TN algorithym
#:<<comment
#mkdir -p $tmp_path/sentieon_tmp
#rsync -avPhW $tumor_bam $normal_bam $tumor_recal $normal_recal $tumor_bam".bai" $normal_bam".bai" $tmp_path/sentieon_tmp
#tumor_bam="$tmp_path/sentieon_tmp/$(basename $tumor_bam)"
#normal_bam="$tmp_path/sentieon_tmp/$(basename $normal_bam)"
#tumor_recal="$/sentieon_tmp/recal_data.table"
#romal_recal="$tmp_path/sentieon_tmp/recal_data.table"
#comment
tumor_name=$($samtools view -H $tumor_bam | \
grep "@RG" | \
grep -v "@PG" | \
awk -F "\t" '{print $6}' | \
awk -F ":" '{print $2}')
#normal_name=$(samtools view -H $normal_bam | \
#grep "@RG" | \
#grep -v "@PG" | \
#awk -F "\t" '{print $6}' | \
#awk -F ":" '{print $2}')
  
tn_cmd="time $sentieon driver -t $NCPU \
-r $hg38 \
-i $tumor_bam \
-q $tumor_recal \
--algo TNhaplotyper2 \
--tumor_sample $tumor_name \
--germline_vcf $germline \
--germline_vcf $germline2 \
--pon $pon \
$rootdir/$project.tmp.vcf.gz \
--algo OrientationBias --tumor_sample $tumor_name \
$rootdir/ORIENTATION_DATA \
--algo ContaminationModel --tumor_sample $tumor_name \
--vcf $germline \
--vcf $germline2 \
--tumor_segments $rootdir/CONTAMINATION_DATA.segments \
$rootdir/CONTAMINATION_DATA"

tn_filter="$sentieon driver -t $NCPU \
-r $hg38 \
--algo TNfilter \
--tumor_sample $tumor_name \
-v $rootdir/$project.tmp.vcf.gz \
--contamination $rootdir/CONTAMINATION_DATA \
--tumor_segments $rootdir/CONTAMINATION_DATA.segments \
--orientation_priors $rootdir/ORIENTATION_DATA \
$rootdir/$project.TN.vcf.gz"

echo $tn_cmd
echo $tn_filter

$tn_cmd

$tn_filter
