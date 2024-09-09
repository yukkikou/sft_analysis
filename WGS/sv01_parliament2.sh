#!/bin/bash

sample=$1
bam_file=$sample".deduped.bam"
bai_file=$sample".deduped.bam.bai"
working_dir=$(pwd)
source_dir="/media/data/SFT_WGS_BAM/"
fasta_file="Homo_sapiens_assembly38.fasta"
fasta_idx="Homo_sapiens_assembly38.fasta.fai"

#prepare input dir
mkdir -p $working_dir/$sample/input
rsync -avPhW $source_dir/$sample/$bam_file $working_dir/$sample/input/
rsync -avPhW $source_dir/$sample/$bai_file $working_dir/$sample/input/
rsync -avPhW $fasta_file $working_dir/$sample/input
rsync -avPhW $fasta_idx $working_dir/$sample/input
docker_cmd="docker run \
    -v $working_dir/$sample/input/:/home/dnanexus/in \
    -v $working_dir/$sample/parliament2:/home/dnanexus/out dnanexus/parliament2:latest \
    --bam $bam_file \
    --bai $bai_file \
    --fai $fasta_idx \
    -r $fasta_file \
    --prefix $sample \
    --manta --cnvnator --lumpy \
    --delly_deletion --delly_insertion \
    --delly_inversion --delly_duplication \
    --genotype  --breakdancer  \
    --svviz --svviz_only_validated_candidates"
$docker_cmd
wait
echo "Parliament task end"
## get ride of temporal files
docker system prune -f

# check annot env 
#if [[ "$CONDA_DEFAULT_ENV" == "annotsv" ]]; then
#  echo "Conda environment "AnnotSV" is already active."
#else
#  conda activate annotsv
#fi

#out_put_vcf="$working_dir/$sample/parliament2/$sample".combined.genotyped.vcf""
#AnnotSV -SVinputFile $out_put_vcf -outputDir "$working_dir/$sample/parliaments/AnnotSV"

## filter SV
#filter_input="$working_dir/$sample/parliaments/AnnotSV/$sample".combined.genotyped.annotated.tsv""
#filter_output="$working_dir/$sample/parliaments/AnnotSV/$sample".combined.filtered.tsv""
#head -n 1 $filter_input > $filter_output
#rg \"SUPP=3\" $filter_input >> $filter_output
#gzip "$working_dir/$sample/parliaments/AnnotSV/*.tsv"

## remove input file
rm "$working_dir/$sample/input/" -rf
#rm "$working_dir/$sample/input/*.fasta*" -f




