#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=90G
#SBATCH --output=trimmomatic-%A_%a.out

# exit when any command fails
set -e

index="${1:-1}"
threads=${SLURM_CPUS_PER_TASK:-1}

if [[ -n "$CC_CLUSTER" ]]
then
  module load StdEnv/2018.3
  module load trimmomatic/0.36

  index=${SLURM_ARRAY_TASK_ID:-0}
  index=$((index+1))
fi

sample=$(awk -v sample_index="$index" \
    '$0 !~ /[ \t]*#/ {ln++} ln == sample_index {print $1}' samples.txt)
sample="${sample%%[[:cntrl:]]}"

echo "Processing sample $sample"

java -jar trimmomatic.jar PE \
    "${sample}_R1.fastq.gz" "${sample}_R2.fastq.gz" \
    "${sample}-paired_R1.fastq.gz" "${sample}-unpaired_R1.fastq.gz" \
    "${sample}-paired_R2.fastq.gz" "${sample}-unpaired_R2.fastq.gz" \
    ILLUMINACLIP:adapters.fa:2:30:10: TRAILING:3 MINLEN:25 \
    -threads "$threads"

fastqc --threads "$threads" "$sample"*.fastq.gz
