#!/usr/bin/env bash
#SBATCH -p mix
#SBATCH --nodes 1
#SBATCH --cpus-per-task=20
#SBATCH -J stringtie
#SBATCH --mem 50G

# environment
workdir=project/SFT
outdir=project/SFT/3_RNAseq
samlist=$outdir/rnaseq_sample_info.tsv
gtf=reference/homo_sapiens/GRCh38/annotation/gencode.v45.basic.annotation.gtf

# array infomation
i=$SLURM_ARRAY_TASK_ID
sample_array=(`cat $samlist`)
prefix=${sample_array[$i]}
bam=$outdir/1_bam/${prefix}/${prefix}.hisat2.sorted.bam

thread=$SLURM_CPUS_PER_TASK

# singularity
sifDir=/media/share/data4/container
stSif="singularity exec $sifDir/stringtie/stringtie-2.2.1.sif stringtie"
samtoolsSif="singularity exec $sifDir/samtools/samtools.sif"

# main
echo "*** slurm node is $SLURMD_NODENAME ***"
echo "*** jod id is $SLURM_JOB_ID ***"
echo "*** sample is ${prefix} *** "

# trim and mapping
mkdir -p $workdir/stringtie/${prefix} && cd $workdir/stringtie/${prefix}

echo "*** stringtie start: `date` ***" && $stSif \
    -e --rf -G $gtf \
    -p $thread \
    -A ${prefix}_gene_abundance.tsv \
    -o ${prefix}.stringtie2.gtf \
    $bam

# data transfor
mkdir -p $outdir/3_stringtie
cp -r $workdir/stringtie/${prefix} $outdir/3_stringtie/

# work space release
rm -r $workdir/stringtie/${prefix}

echo "*** pipeline end: `date` ***"
