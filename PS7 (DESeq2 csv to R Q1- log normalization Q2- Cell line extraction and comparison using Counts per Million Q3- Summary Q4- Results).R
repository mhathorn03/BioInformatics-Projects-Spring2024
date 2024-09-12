## BCH 339N
## RNA-Seq Analysis - 21pt

## We learned about NGS and RNA-Seq in class.
## In this assignment, we will analyze RNA-Seq data using DESeq2. 
## If you want to learn more about this package, check out the vignette
## https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
## We are going to carry out exploratory and differential expression analyses. 

# install.packages("BiocManager")
# BiocManager::install("DESeq2")
library('DESeq2')

## In this homework, we will use the fread function from the data.table package. 
## fread allows us to extract .csv.gz files directly from the web. 
## If you are interested in learning more, data.table is a very popular R package.
install.packages('data.table')
library('data.table')

## We are going to analyze the following dataset
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65912
## Note that there is a link to .csv.gz file with RNA Expression on this page.
## https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65912/suppl/GSE65912_RNASeq.csv.gz

rnaseq_counts = fread ( 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65912/suppl/GSE65912_RNASeq.csv.gz', 
                        data.table = F)


## We will assign gene names to the row.names attribute below. 
## The column headers correspond to samples and look like: "GM12878_Rep1_Counts"
## GM12878 corresponds to the cell line identifier. 
## Rep1 corresponds to the replicate number
row.names(rnaseq_counts) = rnaseq_counts$V1
rnaseq_counts = as.matrix (rnaseq_counts[,-1] ) 

## Q1 - 6pt 
## Data Exploration: 

## A typical exploratory analysis is to inspect a scatterplot of read counts
## Plot two scatterplots using log10 of the gene read counts. 
## For the first plot compare GM12878_Rep1_Counts to GM12878_Rep2_Counts
## For the second plot compare GM12878_Rep1_Counts to GM12891_Rep1_Counts
## Describe your observations regarding these two plots in a few sentences. 

counts_subset <- data.frame(
  GM12878_Rep1 = rnaseq_counts[, "GM12878_Rep1_Counts"],
  GM12878_Rep2 = rnaseq_counts[, "GM12878_Rep2_Counts"],
  GM12891_Rep1 = rnaseq_counts[, "GM12891_Rep1_Counts"]
)

log_counts_subset <- log10(counts_subset)


plot(log_counts_subset$GM12878_Rep1, log_counts_subset$GM12878_Rep2,
     xlab = "GM12878_Rep1 Counts (log10)",
     ylab = "GM12878_Rep2 Counts (log10)",
     main = "Scatterplot: GM12878_Rep1 vs GM12878_Rep2 (log10)")

plot(log_counts_subset$GM12878_Rep1, log_counts_subset$GM12891_Rep1,
     xlab = "GM12878_Rep1 Counts (log10)",
     ylab = "GM12891_Rep1 Counts (log10)",
     main = "Scatterplot: GM12878_Rep1 vs GM12891_Rep1 (log10)")

# Both plots show a positive correlation, but the first plot (GM12878_Rep1 vs GM12878_Rep2)
# has a stronger correlation and less variability (or deviation from the diagonal line), than 
# plot 2 (GM12878_Rep1 vs GM12891_Rep1). Both plots have more spread towards the origin
# suggesting that there are many genes with low expression levels between the samples. Overall, 
# the second plot shows a weaker correlation because of the wider spread and larger dispersion of
# points.

## Q2 - 3pt
## We will  extract the cell line names from the column names. 
## We will then define a factor with two levels of length equal to the number of columns of rnaseq_counts
##  the reference level of this factor will match the cell lines:
## c("GM19238", "GM19240", "GM12892", "GM12878" )
## The alternative level of the factor will match: c("GM19239", "GM12891" ) 

cell_line_ids = sapply (strsplit(colnames(rnaseq_counts), "_"), "[[" , 1 ) 
Factor_level_1 = c("GM19238", "GM19240", "GM12892", "GM12878" )
Factor_level_2 =  c("GM19239", "GM12891" )
factor_of_interest = cell_line_ids %in% Factor_level_1
factor_of_interest = as.factor(factor_of_interest)
factor_of_interest = relevel(factor_of_interest, ref = "TRUE")


## Use DESeqDataSetFromMatrix to prepare the data for differential expression analysis.
dds <- DESeqDataSetFromMatrix(countData = rnaseq_counts,
                               colData = DataFrame(Variable = factor_of_interest ),
                               design = ~ Variable)


counts_per_million <- function (count_matrix) { 
  return (t(apply(count_matrix, 1, function(x){1000000*x/colSums(rnaseq_counts)})))
  
}

## Next, we will filter out low expressed genes from our DESeq object
keep <- rowSums(counts_per_million(counts(dds)) > 1) > 3
dds <- dds[keep,]
## Describe in your own words what the previous two lines of code achieve and how.  

# The first line of code extracts the raw counts from dds (using counts), normalizes and scales the counts by a 
# factor of one million using the counts_per_million function, then it filters the genes for expression
# levels with a CPM greater than 1. The rowSums() function then sums the number samples that have a CPM
# greater than 1. If the sum is greater than 3, it's stored in the keep variable whereas the ones that aren't are 
# filtered out. The dds is then subsetted to the keep values in the second line of code. 



## Q3 - 6pt
dds <- DESeq(dds)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef= "Variable_FALSE_vs_TRUE")
resLFC <- resLFC[order(resLFC$pvalue),]
summary(resLFC, alpha = 0.05)
## Explain in your own words the output of the summary().

# Summary() gives a summary of the DESeq2 analysis. Out of 11563 genes (the number of genes with a nonzero read count),
# 3.6% of genes had a higher expression when the condition was false, and 4.1% of genes had a decrease of expression
# when the varaible was false, with an adjusted p-value <0.05. 10 outliers were identified, and
# no genes had a low count, defined as having a mean count less than 5. 

## Generate an MA-plot using resLFC
?DESeq2::plotMA

plotMA(resLFC)

## Q4 -6pt
## Next we will create an interactive html to explore our results
BiocManager::install("ReportingTools")
library("ReportingTools")
des2Report <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2',
                         title = 'RNA-seq analysis of differential expression using DESeq2',
                         reportDirectory = "./Class Hw")
## This might take a while. 
publish(dds,des2Report, pvalueCutoff=0.05,
        annotation.db="org.Hs.eg", factor = factor_of_interest,
        reportDir="./Class Hw")
finish(des2Report)
## If everything went well, you should be able to inspect an html file containing your results.

## Even if the above code fails, you can see the most differentially expressed genes with the following
resLFC

## You will notice that the top five differentially expressed genes by adjusted p-value should be: 
## KDM5C, PSMA8, KDM6A, ZFX, SMC1A
## Search what each of these genes are. 
## Do you notice any shared feature(s)? If so, can you speculate what our factor_of_interest was?

# All these genes are involved in chromatin structure or modifying histone proteins, all regulate gene expression
# by modulating chromatin structure. The factor of interest in differential expression analysis is probably related
# to processes affecting chromatin structure and gene regulation. 
