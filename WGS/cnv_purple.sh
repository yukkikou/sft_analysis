#!/bin/bash
set -x

cobalt="/home/panda/software/cobalt-1.8.jar"
hg38="/media/db/genomes/Hsapiens/hg38_gatk4/bwa/Homo_sapiens_assembly38.fasta"
cobalt_GC="/media/db/genomes/cobalt/GC_profile.hg38.1000bp.cnp"
amber="/home/panda/software/amber-3.2.jar"
amber_pon="/media/db/genomes/amber3/GermlineHetPon.hg38.vcf.gz"
blacklist_hg38="/home/panda/software/Blacklist/lists/hg38-blacklist.v2.bed"
purple="/home/panda/software/purple-2.43.jar"

if [[ "$#" -ne 3 ]]; then
    echo "$(tput bold; tput setaf 2)Usage: purple-copy-number project_name tumor_bam gender
    $(tput sgr0)"
    exit 1
fi

project=$1
tumor_bam=$2
gender=$3
NCPU=10

mkdir -p $project
output_dir=$(pwd)/$project
gender=${gender^^}

if [ $gender == "FEMALE" ]; then
	ref_bam="/media/data/PON/pon.female/pon.female.deduped.bam"
	echo "Use female PON bam $ref_bam"
elif [ $gender == "MALE" ]; then
	ref_bam="/media/data/PON/pon.male/pon.male.deduped.bam"
	echo "use male PON bam $ref_bam"
fi

#get sample ID from bam
ref_id=$(samtools view -H $ref_bam | \
grep "@RG" | \
grep -v "@PG" | \
awk -F "\t" '{print $6}' | \
awk -F ":" '{print $2}')

tumor_id=$(samtools view -H $tumor_bam | \
grep "@RG" | \
grep -v "@PG" | \
awk -F "\t" '{print $6}' | \
awk -F ":" '{print $2}')
:<<comment1
comment1
# run COBALT
time java -jar -Xmx30G $cobalt \
-reference $ref_id \
-reference_bam $ref_bam \
-tumor $tumor_id \
-tumor_bam $tumor_bam \
-output_dir $output_dir \
-threads $NCPU \
-gc_profile $cobalt_GC \
1>$output_dir/$project".cobalt.stdout" 2>&1

# run AMBER
time java -Xmx64G -jar $amber \
-tumor_only \
-tumor $tumor_id \
-tumor_bam $tumor_bam \
-output_dir $output_dir \
-threads $NCPU \
-loci $amber_pon \
1>$output_dir/$project".amber.stdout" 2>&1

# run gridss
#time gridss -t $NCPU \
#-b $blacklist_hg38 \
#-r $hg38 \
#-w $output_dir \
#-o $output_dir/$project".vcf" \
#-a $output_dir/$project".assembly.bam" \
#--jvmheap 40g \
#$tumor_bam 1>$output_dir/$project".gridss.stdout" 2>&1

#grep "PASS" $output_dir/$project".vcf" | \
#cat <(grep "^#" $output_dir/$project".vcf") - \
#> $output_dir/$project".filtered.vcf"

#bgzip $output_dir/$project".vcf"
#tabix -p vcf $output_dir/$project".vcf.gz"
#bgzip $output_dir/$project".filtered.vcf"
#tabix -p vcf $output_dir/$project".filtered.vcf.gz"

#deleting temp dir
#rm $output_dir/$project".deduped.bam.gridss.working" -rf
#rm $output_dir/$project".assembly.bam.gridss.working" -rf
#rm $output_dir/$project".vcf.gridss.working" -rf

# run purple
time java -jar $purple \
-threads $NCPU \
-tumor_only \
-tumor $tumor_id \
-output_dir $output_dir \
-amber $output_dir \
-cobalt $output_dir \
-gc_profile $cobalt_GC \
-ref_genome $hg38 \
-reference $ref_id \
-circos /home/panda/software/circos/bin/circos 
:<<comment2
comment2



