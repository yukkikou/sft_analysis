#!/usr/bin/env bash
#SBATCH -p mix
#SBATCH --nodes 1
#SBATCH --cpus-per-task=50
#SBATCH -J cjrnaseq
#SBATCH --mem 130G
#SBATCH --array 0-96

# environment
workdir=project/SFT
outdir=project/SFT/3_RNAseq
rnadir=raw/SFT/RNA-seq
genome=reference/homo_sapiens/GRCh38/genome/GRCh38.primary_assembly.genome.fa
index=reference/homo_sapiens/GRCh38/index/hisat2/GRCh38.CHR.fa
samlist=$outdie/rnaseq_sample_info3.tsv
thread=$SLURM_CPUS_PER_TASK

# array infomation
i=$SLURM_ARRAY_TASK_ID
sample_array=(`cat $samlist`)
prefix=${sample_array[$i]}
raw1=$rnadir/${prefix}_1.fq.gz
#raw1=$rnadir/${prefix}_R1.fq.gz
raw2=$rnadir/${prefix}_2.fq.gz
#raw2=$rnadir/${prefix}_R2.fq.gz
hlog=${prefix}.hisat2.summary

# singularity
sifDir=/media/share/data4/container
trimSif="singularity exec $sifDir/trimGalore/trimGalore.sif"
hisat2Sif="singularity exec $sifDir/hisat2/hisat2.sif"
samtoolsSif="singularity exec $sifDir/samtools/samtools.sif"

# main
echo "*** slurm node is $SLURMD_NODENAME ***"
echo "*** jod id is $SLURM_JOB_ID ***"
echo "*** sample is ${prefix} *** "

# trim and mapping
mkdir -p $workdir/$prefix && cd $workdir/$prefix

clean1=${prefix}_1_val_1.fq.gz
clean2=${prefix}_2_val_2.fq.gz
echo "*** trimGalore start: `date` ***" && $trimSif trim_galore \
   -q $thread --paired \
   --length 15 \
   --stringency 3 \
   -o ./ $raw1 $raw2

echo "*** hisat2 start: `date` ***" && $hisat2Sif hisat2 \
    -p $thread \
    --sensitive \
    --rna-strandness RF \
    -x $index \
    --new-summary \
    --summary-file $hlog \
    -1 $clean1 -2 $clean2 \
    -S ${prefix}.hisat2.sam
  
echo "*** samtools start: `date` ***"
$samtoolsSif samtools view -@ $thread \
    -b -o ${prefix}.hisat2.bam \
    ${prefix}.hisat2.sam

$samtoolsSif samtools sort \
    -@ $thread \
    -o ${prefix}.hisat2.sorted.bam \
    ${prefix}.hisat2.bam

# data transfor
rm ${prefix}.hisat2.sam
rm ${prefix}.hisat2.bam

mkdir -p $outdir/1_bam/$prefix
cp $workdir/$prefix/* $outdir/1_bam/$prefix/

# work space release
rm -r $workdir/$prefix

echo "*** pipeline end: `date` ***"
