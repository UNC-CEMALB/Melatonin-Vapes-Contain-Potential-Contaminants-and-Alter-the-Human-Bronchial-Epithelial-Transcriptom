library(shinyBS)
library(shiny)
library(pheatmap)
library(magrittr)
library(DESeq2)
library(ggplot2)
library("stringr")
library(tidyverse)
library(janitor)
library(EnhancedVolcano)
library(RUVSeq)

#Reading in data
setwd("C:\\Users\\kschi\\OneDrive\\Desktop\\MT papers")
countData <- read.csv('raw_counts_allsamples.csv', header = TRUE, sep = ",")
head(countData)
metaData <- read.csv('metadata_allsamples.csv', header = TRUE, sep = ",")
Dups <- duplicated(countData[,1])
summary(Dups)
countData <- countData %>% column_to_rownames("Geneid")

metaData$samples <- gsub("-", ".", metaData$samples, fixed=TRUE)


#QA/QC
#Background filter
nsamp <- ncol(countData)

total_median <- median(as.matrix(countData))


countData <- countData %>% rownames_to_column("Gene")

genes_above_background <- countData %>% # Start from the 'countdata' dataframe
  pivot_longer(cols=!Gene, names_to = "sampleID", values_to="expression") %>% # Melt the data so that we have three columns: gene, exposure condition, and expression counts
  mutate(above_median=ifelse(expression>total_median,1,0)) %>% # Add a column that indicates whether the expression of a gene for the corresponding exposure condition is above (1) or not above (0) the median of all count data
  group_by(Gene) %>% # Group the dataframe by the gene
  summarize(total_above_median=sum(above_median)) %>% # For each gene, count the number of exposure conditions where the expression was greater than the median of all count data
  filter(total_above_median>=.2*nsamp) %>% # Filter for genes that have expression above the median in at least 20% of the samples
  dplyr::select(Gene) # Select just the genes that pass the filter

# Then filter the original 'countdata' dataframe for only the genes above background. 
countData <- left_join(genes_above_background, countData, by="Gene")

#Sample filter
countData_T <- countData %>% 
  pivot_longer(cols=!Gene, names_to="sampleID",values_to="expression") %>% 
  pivot_wider(names_from=Gene, values_from=expression)

countData_T$rowsum <- rowSums(countData_T[2:ncol(countData_T)])

countData_T <- countData_T %>% filter(rowsum!=0)

countData_T <- countData_T %>% dplyr::select(!rowsum) 


countData <- countData_T %>%
  pivot_longer(cols=!sampleID, names_to = "Geneid",values_to="expression") %>% 
  pivot_wider(names_from = sampleID, values_from = "expression")

countData <- countData %>% column_to_rownames("Geneid")
metaData <- metaData %>% column_to_rownames("samples")

#RUVseq

groups <- as.factor(metaData$batch)
differences <- makeGroups(groups)
exprSet <- newSeqExpressionSet(as.matrix(countData),phenoData = 
                                 data.frame(groups,row.names=colnames(countData)))


ruv_set <- RUVs(exprSet, rownames(countData), k=1, differences)










#Differential analysis
dds <- DESeqDataSetFromMatrix(countData=counts(ruv_set), 
                              colData=pData(ruv_set), 
                              design=~0+groups+W_1)

dds <- DESeq(dds)
res <- results(dds)
head(results(dds, tidy=TRUE))
summary(res)
res <- res[order(res$padj),]
head(res)
write.csv(as.data.frame(res[order(res$padj),] ),"C:\\Users\\kschi\\OneDrive\\Desktop\\MT papers\\ctrls_degs.csv")


trt_groups_noIL <- list('vapenoIL', 'MTnoIL')
ctrlnoIL <- "noMTnoIL"

trt_groups_IL <- list('vapeIL', 'MTIL')
ctrlIL <- "noMTIL"


# Run experiment


sigdegs <- data.frame(matrix(ncol = 0, nrow = 10))

# Loop through and extract and export results for all contrasts (treatments vs. control)
for (trt in trt_groups_noIL){ # Iterate for each of the treatments listed in 'trt_groups'
  cat(trt) # Print which treatment group we are on in the loop 
  res <- results(dds, pAdjustMethod = "BH", contrast = c("groups",trt,ctrlnoIL)) # Extract the results of the DESeq2 analysis specifically for the comparison of the treatment group for the current iteration of the loop with the control group
  summary(res) # Print out a high-level summary of the results
  ordered <- as.data.frame(res[order(res$padj),]) # Make a dataframe of the results and order them by adjusted p-value from lowest to highest
  top10 <- head(ordered, n=10) # Make dataframe of the first ten rows of the ordered results
  cat("\nThe 10 most significantly differentially expressed genes by adjusted p-value:\n\n")
  print(top10) # View the first ten rows of the ordered results
  pfilt.05 <- nrow(ordered %>% filter(padj<0.05)) # Get the number of genes that are significantly differentially expressed where padj < 0.05
  cat("\nThe number of genes showing significant differential expression where padj < 0.05 is ", pfilt.05)
  pfilt.10 <- nrow(ordered %>% filter(padj<0.1)) # Get the number of genes that are significantly differentially expressed where padj < 0.10
  cat("\nThe number of genes showing significant differential expression where padj < 0.10 is ", pfilt.10,"\n\n")
  write.csv(ordered, paste0("Module2_4_Output_StatisticalResults_",trt ,".csv")) # Export the full dataframe of ordered results as a csv
  sigdegs[trt]<- row.names(top10)
  assign(paste0("Comparison", trt), ordered)
}
for (trt in trt_groups_IL){ # Iterate for each of the treatments listed in 'trt_groups'
  cat(trt) # Print which treatment group we are on in the loop 
  res <- results(dds, pAdjustMethod = "BH", contrast = c("groups",trt,ctrlIL)) # Extract the results of the DESeq2 analysis specifically for the comparison of the treatment group for the current iteration of the loop with the control group
  summary(res) # Print out a high-level summary of the results
  ordered <- as.data.frame(res[order(res$padj),]) # Make a dataframe of the results and order them by adjusted p-value from lowest to highest
  top10 <- head(ordered, n=10) # Make dataframe of the first ten rows of the ordered results
  cat("\nThe 10 most significantly differentially expressed genes by adjusted p-value:\n\n")
  print(top10) # View the first ten rows of the ordered results
  pfilt.05 <- nrow(ordered %>% filter(padj<0.05)) # Get the number of genes that are significantly differentially expressed where padj < 0.05
  cat("\nThe number of genes showing significant differential expression where padj < 0.05 is ", pfilt.05)
  pfilt.10 <- nrow(ordered %>% filter(padj<0.1)) # Get the number of genes that are significantly differentially expressed where padj < 0.10
  cat("\nThe number of genes showing significant differential expression where padj < 0.10 is ", pfilt.10,"\n\n")
  write.csv(ordered, paste0("Module2_4_Output_StatisticalResults_",trt ,".csv")) # Export the full dataframe of ordered results as a csv
  sigdegs[trt]<- row.names(top10)
  assign(paste0("Comparison", trt), ordered)
}


#Getting normalized counst with hgnc symbols for GSEA

norm_counts <- as.data.frame(counts(dds, normalized=TRUE))
nc <-as.data.frame(normalized_counts)
nc$ensembl <- row.names(nc)
signc <- nc[rownames(top100),] #normalized counts of top 100 degs
rownames(norm_counts)>NULL


norm_counts <- norm_counts %>% rownames_to_column("ensembl")

library( "biomaRt" )
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = norm_counts$ensembl,
                  mart = ensembl )
idx <- match( norm_counts$ensembl, genemap$ensembl_gene_id )
norm_counts$entrez <- genemap$entrezgene[ idx ]
norm_counts$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

norm_counts <- subset(norm_counts, norm_counts$hgnc_symbol != "")
norm_counts = norm_counts[order(norm_counts[,'Gene'],-norm_counts[,'KDS.HBEC.HC1']),]
norm_counts = norm_counts[!duplicated(norm_counts$Gene),]
norm_counts[,1] <- norm_counts[,21]
norm_counts <- subset(norm_counts, select= -c(19,21))
colnames(norm_counts)[1] <- 'Gene'
rownames(norm_counts) <- NULL

norm_counts <- norm_counts %>% column_to_rownames("Gene")

write.csv(as.data.frame(norm_counts),"C:\\Users\\kschi\\OneDrive\\Desktop\\MT papers\\normalized_counts_all.csv")


#Volcano Plots
#MTnoIL
ComparisonMTnoIL <- ComparisonMTnoIL %>% rownames_to_column("ensembl")
library( "biomaRt" )
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = ComparisonMTnoIL$ensembl,
                  mart = ensembl )
idx <- match( ComparisonMTnoIL$ensembl, genemap$ensembl_gene_id )
ComparisonMTnoIL$entrez <- genemap$entrezgene[ idx ]
ComparisonMTnoIL$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

ComparisonMTnoIL <- subset(ComparisonMTnoIL, ComparisonMTnoIL$hgnc_symbol != "")
ComparisonMTnoIL = ComparisonMTnoIL[!duplicated(ComparisonMTnoIL$hgnc_symbol),]
rownames(ComparisonMTnoIL) <- NULL
ComparisonMTnoIL <- ComparisonMTnoIL %>% column_to_rownames("hgnc_symbol")

EnhancedVolcano(ComparisonMTnoIL, 
                lab = rownames(ComparisonMTnoIL), # Label information from dataset (can be a column name)
                x = 'log2FoldChange', # Column name in dataset with l2fc information
                y = 'padj', # Column name in dataset with adjusted p-value information
                ylab = "-Log(FDR-adjusted p)", # Y-axis label
                pCutoff= 0.1, # Set p-value cutoff
                ylim=c(0,5), # Limit y-axis for better plot visuals
                xlim=c(-2,2), # Limit x-axis (similar to in MA plot y-axis)
                title="Volcano Plot", # Set title
                subtitle = "MTnoIL", # Set subtitle
                legendPosition = 'bottom') # Put legend on bottom






heatmap_samp <- norm_counts[1:50,]


#plotting and saving heatmap
pdf("CBDq24_HM.pdf", width = 6, height = 7)

dd = subset(metaData, select = -c(samples))
rownames(dd) <- metaData$samples
dd <- as.data.frame(dd)

heatmap <- pheatmap(
  norm_counts,
  cluster_rows = TRUE, # Cluster the rows of the heatmap (genes in this case)
  cluster_cols = TRUE, # Cluster the columns of the heatmap (samples),
  show_rownames = FALSE, # There are too many genes to clearly show the labels
  annotation_col = metaData,
  main = "IL-13 Treatment",
  colorRampPalette(c(
    "red",
    "white",
    "blue"
  ))(25
  ),
  scale = "row", # Scale values in the direction of genes (rows)
  fontsize = 20
)
png("topgenes_heatmap_allsamples.png", width = 1280, height = 1440, bg = NA)
heatmap
dev.off()









#a <- filter(hmtable,hmtable$Gene %in% aframe[degsorted$Symbol])

##Plotting GSEA pathways

library(enrichplot)
library(PPInfer)
library(DOSE)
setwd("C:\\Users\\kschi\\gsea_home\\output\\nov08\\hall_no13_sol.Gsea.1668003850706")
#WikiPathways top 10

pos <- read.table(file = 'gsea_report_for_0_1668003850706.tsv', sep = '\t', header = TRUE)
pos <- data.frame(pos)
neg <- read.table(file = 'gsea_report_for_4_1668003850706.tsv', sep = '\t', header = TRUE)
neg <- data.frame(neg)
gsea <- rbind.data.frame(pos,neg)
gsea <- gsea[c(1,6,8)]
gsea <- gsea[sort(abs(gsea$NES),decreasing=T,index.return=T)[[2]],]
gsea <- gsea[1:10,]
gsea <- gsea[order(gsea$NES, decreasing=TRUE),]
rownames(gsea) <- NULL
gsea$NAME <- gsub('_', ' ', gsea$NAME)
gsea$NAME <- gsub('HALLMARK','',gsea$NAME)
gsea$NAME <- stringr::str_wrap(gsea$NAME, 35)

ggplot(gsea, aes(NES, y= reorder(NAME, +NES), fill = FDR.q.val)) +               # ggplot2 barplot with single color
  geom_col() +
  scale_fill_continuous(low = "red2",  high = "mediumblue", limits=c(0,0.37)) +
  theme(plot.title = element_text(hjust = 0.5)) + labs(x="Normalized Enrichment Score", y=' ') +
  scale_x_continuous(limits = c(-3.5,3.5)) +
  geom_text(aes(label = NAME,
                hjust=ifelse(NES>=0,1, 0), x = ifelse(NES < 1, .1, -.1)),color="black", size = 4.5, fontface = 'bold') +
  theme(axis.ticks.y = element_blank(),
         axis.ticks = element_blank(),
         legend.position = 'none',
          axis.text.y = element_blank())+
  geom_vline(xintercept = 0, linewidth = 2)+
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(face = 'bold', size = 14, color = 'black'),axis.line.x.bottom = element_line(linewidth=2),axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.x.bottom = element_text(size = 20), panel.background = element_rect(fill='transparent'), axis.ticks.y = element_blank(), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), legend.background = element_rect(fill = "transparent"),axis.ticks.x.bottom = element_line(linewidth=2,color='black'))
ggsave('legend.png', width = 8, height = 8)   
dev.off()
