#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=90G
#SBATCH --output=star-index-%A.out

# exit when any command fails
set -e

if [[ -n "$CC_CLUSTER" ]]
then
  # load required modules
  module purge
  module load StdEnv/2020 star/2.7.9a
fi

threads=${SLURM_CPUS_PER_TASK:-1}

echo "Generating index for Mus Musculus"

STAR --runMode genomeGenerate \
  --genomeDir m_musculus/star_index \
  --genomeFastaFiles m_musculus/GRCm38.primary_assembly.customCH12genome.fa \
  --genomeSAindexNbases 12 \
  --sjdbGTFfile m_musculus/gencode.vM23.primary_assembly.customCH12.annotation.gff3 \
  --runThreadN "$threads"

echo "Generating index for S. Cerevisiae"

STAR --runMode genomeGenerate \
  --genomeDir s_cerevisiae/star_index \
  --genomeFastaFiles s_cerevisiae/s_cerevisiae.R64-1-1.fa \
  --genomeSAindexNbases 10 \
  --sjdbGTFfile s_cerevisiae/s_cerevisiae.R64-1-1.gtf \
  --runThreadN "$threads"

