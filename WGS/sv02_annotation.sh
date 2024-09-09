#!/usr/bin/env bash
#SBATCH -p compute,gpu,single-ib,single-em,single-fc
#SBATCH --nodes 1
#SBATCH --cpus-per-task=8
#SBATCH -J AnnotSV
#SBATCH -o /home/hxy/log/%x-%A_%a.out
#SBATCH -e /home/hxy/log/%x-%A_%a.err
#SBATCH --mem 20G
#SBATCH --array 92-93

# environment
workdir=/media/work/hxy/SFT
outdir=/media/share/data3/hxy/CJ_SFT/SV
samlist=/media/share/data3/hxy/project/SFT/0_data/1_CJ/sv_vcf_list.tsv

thread=$SLURM_CPUS_PER_TASK

# array infomation
i=$SLURM_ARRAY_TASK_ID
sample_array=(`cat $samlist`)
file=${sample_array[$i]}
prefix=$(basename $file ".combined.genotyped.vcf")

# conda
if [[ "$CONDA_DEFAULT_ENV" == "AnnotSV" ]]; then
  echo "Conda environment "AnnotSV" is already active."
else
  mamba activate AnnotSV && echo "Activate AnnotSV "
fi

annotsv="/home/hxy/software/AnnotSV-master/bin/AnnotSV"

# main
echo "*** slurm node is $SLURMD_NODENAME ***"
echo "*** jod id is $SLURM_JOB_ID ***"

# annotation
mkdir -p $workdir/$prefix && cd $workdir/$prefix

echo "*** AnnotSV start: `date` ***" && time $annotsv \
    -SVinputFile $file \
    -outputDir AnnotSV

# test
echo "********"
cat AnnotSV/${prefix}.combined.genotyped.annotated.tsv | grep "NAB2;STAT6"
# compress
gzip AnnotSV/*.tsv

# data transfor
cp -r $workdir/$prefix/AnnotSV $outdir/$prefix/parliament2/

# work space release
rm -r $workdir/$prefix

echo "*** pipeline end: `date` ***"
