# Dot1l project code

This repo contains code for the TT-seq pipeline used in *Subramani et al.* (2025) article: **DOT1L activity limits transcription elongation velocity and favors RNAPII pausing to facilitate mutagenesis by Activation Induced Deaminase**.

### Requirements

* A system running Linux or similar
* Python version 3.12.4 or later
* Trimmomatic version 0.36
* bedtools version 2.31.0
* FastQC version 0.12.1
* Bowtie 2 version 2.5.4
* SAMtools version 1.20
* Picard version 3.1.0
* deepTools version 3.5.1
* MultiQC version 1.25.2
* R version 4.4 with package DESeq2 version 1.44.0

## Running the pipeline

### Download genomes

```shell
wget https://github.com/francoisrobertlab/dinoia_dot1l_project/releases/download/1.0/GRCm38.primary_assembly.customCH12genome.fa.gz
wget https://github.com/francoisrobertlab/dinoia_dot1l_project/releases/download/1.0/gencode.vM23.primary_assembly.customCH12.annotation.gff3.gz
wget https://ftp.ensembl.org/pub/release-112/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-112/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.112.gtf.gz
gunzip *.gz
mkdir m_musculus
cp GRCm38.primary_assembly.customCH12genome.fa m_musculus
cp gencode.vM23.primary_assembly.customCH12.annotation.gff3 m_musculus
mkdir s_cerevisiae
cp Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa s_cerevisiae
cp Saccharomyces_cerevisiae.R64-1-1.112.gtf s_cerevisiae
```

### Download FASTQ files

> [!WARNING]  
> Add SRA identifiers to samples.txt to download FASTQ.

### Trimming reads and alignment

```shell
seq 1 6 | xargs -L 1 -I % bash trimmomatic.sh %
bash star-index.sh
seq 1 6 | xargs -L 1 -I % bash star.sh %
seq 1 6 | xargs -L 1 -I % bash filterbam.sh %
```

### Merge STAR counts

```shell
python merge_star_counts.py -o m_musculus_counts.txt --name "(.*)_m_musculus.ReadsPerGene.out.tab" *_m_musculus.ReadsPerGene.out.tab
sed -i.bak 's/-/_/g' m_musculus_counts.txt
```

### Obtain counts for mouse ribosomal genes

```shell
python extract_go_genes.py --go mouse_ribosomal_genes.txt --genes gencode.vM23.primary_assembly.customCH12.annotation.gff3 m_musculus_counts.txt > m_musculus_counts_ribosomes.txt
```

### Scale factors

DESeq2 scale factors are computed using the estimated size factors from DESeq2 using this formula: `1 / estimated size factor`

To obtain the estimated size factors from DESeq2, run

```shell
Rscript deseq2-estimateSizeFactor.R
```

### Genome coverage tracks

```shell
bash sample-genomecoverage.sh
```
