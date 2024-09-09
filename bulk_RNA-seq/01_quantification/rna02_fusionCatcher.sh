#!/usr/bin/env bash
#SBATCH -p mix
#SBATCH --nodes 1
#SBATCH --cpus-per-task=40
#SBATCH -J fusion
#SBATCH --mem 120G

# environment
workdir=project/SFT
outdir=project/SFT/3_RNAseq
rnadir=raw/SFT/RNA-seq
samlist=$outdir/rnaseq_sample_info2.tsv
index=reference/homo_sapiens/GRCh38/index/fusioncather/human_v102
thread=$SLURM_CPUS_PER_TASK

# array infomation
i=$SLURM_ARRAY_TASK_ID
sample_array=(`cat $samlist`)
prefix=${sample_array[$i]}
raw1=$rnadir/${prefix}_1.fq.gz
raw2=$rnadir/${prefix}_2.fq.gz

# singularity
sifDir=/media/share/data4/container
fcSif="singularity exec $sifDir/fusioncatcher/fusioncatcher-1.33.sif fusioncatcher"

# main
echo "*** slurm node is $SLURMD_NODENAME ***"
echo "*** jod id is $SLURM_JOB_ID ***"

# trim and mapping
mkdir -p $workdir/fusion/${prefix} && cd $workdir/fusion/${prefix}

echo "*** fusioncatcher start: `date` ***" && $fcSif \
    -p $thread \
    --skip-fastqtk \
    --Xmx=100g \
    -d $index \
    -i $raw1,$raw2 \
    -o ./

# data transfor
mkdir -p $outdir/2_fusion
cp -r $workdir/fusion/$prefix $outdir/2_fusion/

# work space release
rm -r $workdir/fusion/$prefix

echo "*** pipeline end: `date` ***"
