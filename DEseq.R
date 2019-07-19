##  AUTHOR: SHWETHA HARA SRIDHAR
##  TECHDEV TEAM, DEPARTMENT OF GENETICS AND GENOMICS
##  ICAHN SCHOOL OF MEDICINE AT MOUNT SINAI, NY
##  VERSION: 19.7.18
##  METHOD:
##        THIS RSCRIPT, USES THE FINAL COUNTS COMBINED OUTPUT CSV FILE TO PERFORM DESEQ2 ANALYSIS
##        IT FILTERS FOR A P-VALUE CUT OFF GIVEN BY THE USER
##        ORDERS THE DEGENES BASED ON THE LOG2FOLDCHANGE VALUES CALCULATED BY THE DESEQ METHOD
##        IT ALSO WRITE A PRE-RANKED .RNK FILE WITH THE DIFFERENTIALLY EXPRESSED GENES WITH THEIR 
##        RESPECTIVE LOG2FOLDCHANGE VALUES, WHICH CAN BE FURTHER USED FOR GSEA ANALYSIS
##        FINALLY, IT PLOTS A HEATMAP IDENTIFYING THE UPREGULATION AND DOWNREGULATION OF THE GENESET


### USING THE COMBINED COUNTS TO PERFORM DESEQ AND IDENTIFY DIFFERENTIALLY EXPRESSED GENES ###

library(DESeq2)
D_express <- function(in_matrix, condition, p_val=".", count_cutoff = 0, outDEseq, outrank){
  countdata <- in_matrix
  
  (coldata <- data.frame(row.names=colnames(countdata), 
                         condition))
  
  dds <- DESeqDataSetFromMatrix(countData=countdata, 
                                colData=coldata, 
                                design=~condition)
  keep <- rowSums(counts(dds)) >= count_cutoff
  dds <- dds[keep,]
  dds <- DESeq(dds)
  res<-results(dds)
  if (isTRUE(p_val != ".")){
    res_filt <- res[ which(res$padj < p_val), ]  
    # filtering based on p-val cut off of 1e-10 and saving it into a variable
    resOrdered <- res_filt[order(-res_filt$log2FoldChange),]             # ordering based on fold change values
    print("Filtered DEseq:")
    write.csv(as.data.frame(resOrdered),file=outDEseq)                                          # uncomment to save the DEseq results as a csv file
    f <- as.data.frame(resOrdered)
    f$lfcSE<- f$stat<- f$pvalue<- f$padj<- f$baseMean <- NULL
    write.table(as.data.frame(f),file=outrank, sep = "\t", col.names = FALSE, quote = FALSE)                                          # uncomment to save the DEseq results as a csv file
  }
  else{
    print("DEseq RESULTS:")
    f <- as.data.frame(res)
    f$lfcSE<- f$stat<- f$pvalue<- f$padj<- f$baseMean <- NULL
    write.table(as.data.frame(f), file=outrank, sep = "\t", col.names = FALSE, quote = FALSE)                                          # uncomment to save the DEseq results as a csv file
    write.csv(as.data.frame(res), file=outDEseq)                                          # uncomment to save the DEseq results as a csv file
    
  }
  
  
  dds
}

### USING THE DESEQ RESULTS TO PLOT A HEATMAP ###

library("pheatmap")
library("genefilter")

heat_mapping <- function(dds, cutoff=0, p_val= ".", output_heatmap){
  
  if (isTRUE(p_val != ".")){
    dds <- dds[which(results(dds)$padj < p_val)]
    dds <- dds[order(results(dds)$padj),]
  }else{
    dds <- dds[order(results(dds)$padj),]
  }
  if (cutoff>0){
    dds <- dds[1:cutoff]
    dds <- dds[order(-results(dds)$log2FoldChange),] 
  }
  else{
    dds <- dds[order(-results(dds)$log2FoldChange),] 
  }
  select <- order(rowMeans(counts(dds,normalized=TRUE)), 
                  decreasing=TRUE)
  nt <- normTransform(dds)                                             # defaults to log2(x+1)
  log2.norm.counts <- assay(nt)[select,]
  df <- as.data.frame(colData(dds)[,c("condition", "sizeFactor")]) 
  df$sizeFactor <- NULL                                                # removing sizeFactor column out of colData(dds)
  log2.norm.counts 
  
  pheatmap(log2.norm.counts, 
           cluster_rows=FALSE, 
           show_rownames=TRUE,
           cluster_cols=FALSE, 
           annotation_col=df, 
           cellwidth = 10, cellheight = 4,
           filename = output_heatmap, 
           fontsize = 3,
           border_color = "black")
}


################### CANDIDA #######################
file_format = "*.csv"                                                  # path to final combined counts .csv
run_path = "/Candida/Candida_featureCounts/"

setwd(run_path)
skipit = 1
cutoff <- 50
files <- list.files(run_path, file_format)
length(files)
data_counts <- read.csv(files[1], 
                        sep = ",", 
                        header = TRUE)
data_matrix<-as.matrix(data_counts[,-1])
rownames(data_matrix)<-data_counts[,1]
datain <- data_matrix
colnames(datain)
colnames(datain) <- factor(c("18_1",
                             "17_1",
                             "18_2",
                             "16_1"))
datain <- datain[ , c("16_1", "17_1","18_1", "18_2")]         #Ordering the columns based on conditions
condition <- factor(c("Sensitive","Resistant","Resistant","Resistant")) 

p_val_cutoff <- 0.05 
Dexpress_res <- D_express(datain, condition, 
                          p_val = p_val_cutoff,
                          outDEseq = "DEgenes.csv", 
                          outrank = "Preranked_DEgenes.rnk")
heat_mapping(Dexpress_res, cutoff = cutoff, p_val = p_val_cutoff,
             output_heatmap = "16_vs_17_18_1.2.pdf")  

heat_mapping(Dexpress_res, cutoff = cutoff, 
             output_heatmap = "withoutpvalcutoff_16_vs_17_18_1.2.pdf")  

