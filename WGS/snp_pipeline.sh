#!/usr/bin/env bash
#SBATCH -p single-ib,single-em,single-fc,compute,gpu
#SBATCH --nodes 1
#SBATCH --cpus-per-task=10
#SBATCH -J snv
#SBATCH --mem 30G

# environment
workdir=project/SFT
outdir=project/SFT/WGS
genome=reference/Homo_sapiens/GATKvcf/GRCh38_hg38/Homo_sapiens_assembly38.fasta
samlist=$outdir/snv_sample.list

thread=$SLURM_CPUS_PER_TASK

# array infomation
i=$SLURM_ARRAY_TASK_ID
sample_array=(`cat $samlist`)
file=${sample_array[$i]}
prefix=$(basename $file ".TN.vcf.gz")
mvcf=$outdir/1_somaticMerge/Somatic_SFT_merge_norm.vcf.gz
mfcvf=$outdir/1_somaticMerge/Somatic_SFT_merge_norm_filter.vcf.gz

# singularity
sifDir=~/singularity
gatkSif="singularity exec $sifDir/gatk/gatk-4.sif gatk "
bcftoolsSif="singularity exec $sifDir/bcftools/bcftools-1.19.sif bcftools"
aCon="perl /home/hxy/software/annovar/convert2annovar.pl"
aTab="perl /home/hxy/software/annovar/table_annovar.pl"
humandb="/home/hxy/software/annovar/humandb"

# main
echo "*** slurm node is $SLURMD_NODENAME ***"
echo "*** jod id is $SLURM_JOB_ID ***"

# vcf left align
mkdir -p $workdir/$prefix && cd $workdir/$prefix

echo "*** bcftools start: `date` ***" && $bcftoolsSif norm \
    --atomize -d none -f $genome \
    -Oz -o ${prefix}.TN.norm.vcf.gz \
    --threads $thread \
    --write-index $file

# data transfor
mkdir -p $outdir/tmp_leftAlign/$prefix
cp $workdir/$prefix/* $outdir/tmp_leftAlign/$prefix

# merge
mkdir -p $workdir/somaticMerge && cd $workdir/somaticMerge

echo "*** bcftools start: `date` ***" && $bcftoolsSif merge \
    --force-samples \
    -l $samlist \
    -Oz -o Somatic_SFT_merge_norm.vcf.gz \
    --threads $thread \
    --write-index

# data transfor
mkdir -p $outdir/1_somaticMerge
cp $workdir/somaticMerge/* $outdir/1_somaticMerge/

# vcf merge
mkdir -p $workdir/filter && cd $workdir/filter

echo "*** bcftools filter 1 start: `date` ***" && $bcftoolsSif view \
    -i 'FILTER="PASS" & AF>0.05' \
    -Oz -o $mfvcf \
    --threads $thread \
    --write-index \
    $mvcf

echo "*** bcftools stats start: `date` ***" && $bcftoolsSif stats \
    --threads $thread \
    -s - \
    -I \
    $mfvcf > Somatic_SFT_merge_norm_filter.stats

# data transfor
mkdir -p $outdir/1_annovar
cp $workdir/annovar/* $outdir/1_annovar


# annotation
mkdir -p $workdir/annovar && cd $workdir/annovar

echo "*** Annovar start: `date` ***" && $aTab --outfile Somatic_SFT_merge_norm_filter_anno \
    --buildver hg38 \
    --checkfile \
    --remove \
    --nastring . \
    --thread $thread \
    -protocol refGene,knownGene,ensGene,EUR.sites.2015_08,AFR.sites.2015_08,ALL.sites.2015_08,AMR.sites.2015_08,EAS.sites.2015_08,dbscsnv11,intervar_20180118,exac03,exac03nonpsych,exac03nontcga,gene4denovo201907,gnomad40_exome,gnomad40_genome,avsnp150,clinvar_20221231,wgRna \
    -operation g,g,g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r \
    $outdir/1_annovar/Somatic_SFT_merge_norm_filter.avinput $humandb

# work space release
rm -r $workdir/$prefix

echo "*** pipeline end: `date` ***"
