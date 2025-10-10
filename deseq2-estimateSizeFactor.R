library("DESeq2")

count_matrix = read.csv("m_musculus_counts.txt", sep="\t", row.names="Gene")

coldata <- read.csv("samples-condition.txt", sep="\t", row.names="Sample")
coldata$Condition <- factor(coldata$Condition)

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ Condition)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
