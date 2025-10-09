#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=sample-genomecoverage-%A.out

# exit when any command fails
set -e

threads="${SLURM_CPUS_PER_TASK:-1}"
tmpdir="${SLURM_TMPDIR:-/tmp}"

if [[ -n "$CC_CLUSTER" ]]
then
  module purge
  module load StdEnv/2020
  module load samtools/1.17
  module load deeptools/3.5.1
fi

bin_size=10

genome_coverage() {
  local sample=$1
  local scaleFactor=$2
  local suffix=$3
  local bam="${sample}_m_musculus.sorted.marked.filtered.bam"

  echo "Calling bam coverage for sample ${sample}${suffix} with scale factor ${scaleFactor}"
  bamCoverage --scaleFactor "$scaleFactor" --binSize "$bin_size" -p "$threads" \
      -b "$bam" -o "${sample}_m_musculus${suffix}-cov.bw"

  echo "Calling bam coverage for forward reads on sample ${sample}${suffix} with scale factor ${scaleFactor}"
  local bam_forward="${tmpdir}/${sample}-for.bam"
  local bam_forward1="${tmpdir}/${sample}-for1.bam"
  local bam_forward2="${tmpdir}/${sample}-for2.bam"
  samtools view -b -f 128 -F 16 --threads "$threads" "$bam" > "$bam_forward1"
  samtools view -b -f 80 --threads "$threads" "$bam" > "$bam_forward2"
  samtools merge --threads "$threads" -f "$bam_forward" "$bam_forward1" "$bam_forward2"
  samtools index "$bam_forward"
  bamCoverage --scaleFactor "$scaleFactor" --binSize "$bin_size" -p "$threads" \
      -b "$bam_forward" -o "${sample}_m_musculus${suffix}-cov-for.bw"

  echo "Calling bam coverage for reverse reads on sample ${sample}${suffix} with scale factor -${scaleFactor}"
  local bam_reverse="${tmpdir}/${sample}-rev.bam"
  local bam_reverse1="${tmpdir}/${sample}-rev1.bam"
  local bam_reverse2="${tmpdir}/${sample}-rev2.bam"
  samtools view -b -f 144 --threads "$threads" "$bam" > "$bam_reverse1"
  samtools view -b -f 64 -F 16 --threads "$threads" "$bam" > "$bam_reverse2"
  samtools merge --threads "$threads" -f "$bam_reverse" "$bam_reverse1" "$bam_reverse2"
  samtools index "$bam_reverse"
  bamCoverage --scaleFactor "-${scaleFactor}" --binSize "$bin_size" -p "$threads" \
      -b "$bam_reverse" -o "${sample}_m_musculus${suffix}-cov-rev.bw"
}

genome_coverage WTG-1 0.902055939
genome_coverage WTG-2 1.124646102
genome_coverage WTG-3 1.014961548
genome_coverage Dot1l-KOB-1 1.128120282
genome_coverage Dot1l-KOB-2 0.930284654
genome_coverage Dot1l-KOB-3 0.896388353

genome_coverage WTG-1 0.091777052 -spike
genome_coverage WTG-2 0.102371434 -spike
genome_coverage WTG-3 0.083561247 -spike
genome_coverage Dot1l-KOB-1 0.134205958 -spike
genome_coverage Dot1l-KOB-2 0.099105788 -spike
genome_coverage Dot1l-KOB-3 0.105260809 -spike
