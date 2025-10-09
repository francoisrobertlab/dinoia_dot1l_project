#!/bin/bash
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=90G
#SBATCH --output=filterbam-%A_%a.out

# exit when any command fails
set -e

index="${1:-1}"
threads="${SLURM_CPUS_PER_TASK:-1}"

if [[ -n "$CC_CLUSTER" ]]
then
  module purge
  module load StdEnv/2020
  module load samtools/1.17
  module load robtools/core/2.0

  index=${SLURM_ARRAY_TASK_ID:-0}
  index=$((index+1))
fi

sample=$(awk -v sample_index="$index" \
    '$0 !~ /[ \t]*#/ {ln++} ln == sample_index {print $1}' samples.txt)
sample="${sample%%[[:cntrl:]]}"

bam="${sample}_m_musculus.sorted.marked.bam"
filtered="${sample}_m_musculus.sorted.marked.filtered.bam"
samtools view -b -F 2048 -F 256 -f 2 --threads "$threads" -o "$filtered" "$bam"
samtools index -@ "$threads" "$filtered"
