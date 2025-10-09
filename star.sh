#!/bin/bash
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=90G
#SBATCH --output=star-%A_%a.out

# exit when any command fails
set -e

index="${1:-1}"
threads=${SLURM_CPUS_PER_TASK:-1}
tmpdir=${SLURM_TMPDIR:-${PWD}}

if [[ -n "$CC_CLUSTER" ]]
then
  module purge
  module load StdEnv/2020
  module load star/2.7.9a
  module load samtools/1.17
  module load picard/2.26.3
  module load robtools/core/2.0

  index=${SLURM_ARRAY_TASK_ID:-0}
  index=$((index+1))
fi

genomes=(m_musculus s_cerevisiae)
sample=$(awk -v sample_index="$index" \
    '$0 !~ /[ \t]*#/ {ln++} ln == sample_index {print $1}' samples.txt)
sample="${sample%%[[:cntrl:]]}"

for genome in "${genomes[@]}"
do
  fastq_1="${sample}-paired_R1.fastq.gz"
  fastq_2="${sample}-paired_R2.fastq.gz"
  star_index="${genome}/star_index"
  
  rm -rf "${tmpdir}/${sample}"

  STAR \
    --runThreadN "$threads" \
    --runMode alignReads \
    --genomeDir "${star_index}" \
    --readFilesIn "$fastq_1" "$fastq_2" \
    --readFilesCommand zcat \
    --quantMode TranscriptomeSAM GeneCounts \
    --twopassMode Basic \
    --outSAMunmapped None \
    --outSAMattrRGline "ID:${sample}" "PU:${sample}" "SM:${sample}" LB:unknown PL:illumina \
    --outSAMtype BAM Unsorted \
    --outTmpDir "${tmpdir}/${sample}" \
    --outFileNamePrefix "${sample}_${genome}."

  samtools sort \
    --threads "$threads" \
    -o "${sample}_${genome}.sorted.bam" \
    "${sample}_${genome}.Aligned.out.bam"

  java -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
    INPUT="${sample}_${genome}.sorted.bam" \
    OUTPUT="${sample}_${genome}.sorted.marked.bam" \
    METRICS_FILE="${sample}_${genome}.sorted.marked.metrics" \
    REMOVE_DUPLICATES=false \
    ASSUME_SORTED=true \
    MAX_RECORDS_IN_RAM=2000000 \
    VALIDATION_STRINGENCY=LENIENT \
    TMP_DIR="${tmpdir}/${sample}"

  samtools index "${sample}_${genome}.sorted.marked.bam"

  rm "${sample}_${genome}.sorted.bam"
done
