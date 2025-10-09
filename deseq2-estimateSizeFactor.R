library("DESeq2")

count_matrix = read.csv("count_matrix.txt", sep="\t", row.names="Gene_ID")

coldata <- read.csv("samples-condition.txt", sep="\t", row.names="Sample")
coldata$Condition <- factor(coldata$Condition)

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ Condition)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
