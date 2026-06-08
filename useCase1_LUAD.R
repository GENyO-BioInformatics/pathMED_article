# Load necessary libraries
library(pathMED)
library(biomaRt)
library(ggplot2)
library(WriteXLS)
library(readxl)
suggest <- c("ada","AUCell", "Biobase", "BiocGenerics", "BiocStyle", "doMC", "fgsea", "gam", "GSEABase", "import", "kernlab", "klaR", "knitr", "mboost", "MLeval",
             "randomForest", "ranger", "rmarkdown", "RUnit", "SummarizedExperiment", "utils", "xgboost")
for(i in suggest){library(i,character.only = T)}

# Save in this variable the number of CPUs to be used (default = 10)
ncores <- 10

# Download data from https://kb.linkedomics.org/download#LUAD and save into 
# LinkedOmicsKB_LUAD folder

# Create output folders
dir.create("SupplTables")
dir.create("figures")

# Read and prepare the metadata
metadata <- read.delim("LinkedOmicsKB_LUAD/LUAD_meta.txt")
metadata <- metadata[-1,]
rownames(metadata) <- metadata[,1]

# Read and process omics data
exprTumor <- read.delim("LinkedOmicsKB_LUAD/LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt", row.names = 1)
protTumor <- read.delim("LinkedOmicsKB_LUAD/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt", row.names = 1)
exprNormal <- read.delim("LinkedOmicsKB_LUAD/LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt", row.names = 1)
protNormal <- read.delim("LinkedOmicsKB_LUAD/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt", row.names = 1)


## Discard genes in pseudoautosomal regions
exprTumor <- exprTumor[grep("PAR", rownames(exprTumor), invert = T),]
exprNormal <- exprNormal[grep("PAR", rownames(exprNormal), invert = T),]

# Remove ENSEMBL version numbers
rownames(exprTumor) <- gsub("\\..*", "", rownames(exprTumor))
rownames(protTumor) <- gsub("\\..*", "", rownames(protTumor))
rownames(exprNormal) <- gsub("\\..*", "", rownames(exprNormal))
rownames(protNormal) <- gsub("\\..*", "", rownames(protNormal))

# Translate ENSEMBL IDs to Gene Symbols

## Transcriptomics tumor
connect_ensembl()
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                values=rownames(exprTumor),mart = mart)
G_list <- G_list[!duplicated(G_list[,1]),]
rownames(G_list) <- G_list[,1]

exprTumor <- exprTumor[rownames(exprTumor) %in% G_list[,1],]
exprTumor$GeneSymbol <- G_list[rownames(exprTumor), 2]
exprTumor_symbol <- aggregate(exprTumor, list(exprTumor$GeneSymbol), median)
exprTumor_symbol <- na.omit(exprTumor_symbol)
rownames(exprTumor_symbol) <- exprTumor_symbol[,1]
exprTumor_symbol <- exprTumor_symbol[,-1]
exprTumor_symbol <- exprTumor_symbol[,-ncol(exprTumor_symbol)]

## Transcriptomics normal
exprNormal <- exprNormal[rownames(exprNormal) %in% G_list[,1],]
exprNormal$GeneSymbol <- G_list[rownames(exprNormal), 2]
exprNormal_symbol <- aggregate(exprNormal, list(exprNormal$GeneSymbol), median)
exprNormal_symbol <- na.omit(exprNormal_symbol)
rownames(exprNormal_symbol) <- exprNormal_symbol[,1]
exprNormal_symbol <- exprNormal_symbol[,-1]
exprNormal_symbol <- exprNormal_symbol[,-ncol(exprNormal_symbol)]


## Proteomics tumor
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                values=rownames(protTumor),mart = mart)
G_list <- G_list[!duplicated(G_list[,1]),]
rownames(G_list) <- G_list[,1]

protTumor <- protTumor[rownames(protTumor) %in% G_list[,1],]
protTumor$GeneSymbol <- G_list[rownames(protTumor), 2]
protTumor_symbol <- aggregate(protTumor, list(protTumor$GeneSymbol), median)
protTumor_symbol <- na.omit(protTumor_symbol)
rownames(protTumor_symbol) <- protTumor_symbol[,1]
protTumor_symbol <- protTumor_symbol[,-1]
protTumor_symbol <- protTumor_symbol[,-ncol(protTumor_symbol)]


## Proteomics normal
protNormal <- protNormal[rownames(protNormal) %in% G_list[,1],]
protNormal$GeneSymbol <- G_list[rownames(protNormal), 2]
protNormal_symbol <- aggregate(protNormal, list(protNormal$GeneSymbol), median)
protNormal_symbol <- na.omit(protNormal_symbol)
rownames(protNormal_symbol) <- protNormal_symbol[,1]
protNormal_symbol <- protNormal_symbol[,-1]
protNormal_symbol <- protNormal_symbol[,-ncol(protNormal_symbol)]

resultsList <- list()

for (iteration in seq_len(10)) {
  cat("Iteration:",iteration,"\n")
  set.seed(iteration)
  
  #### Use different samples for training and validation
  samplesExprNormal <- sample(colnames(exprNormal_symbol), round(ncol(exprNormal_symbol) * 0.75), replace = F)
  samplesProtNormal <- colnames(protNormal_symbol)[!colnames(protNormal_symbol) %in% samplesExprNormal]
  
  samplesOnlyTumor <- setdiff(colnames(exprTumor_symbol), colnames(exprNormal_symbol))
  
  samplesExprTumor <- c(samplesExprNormal, sample(samplesOnlyTumor, 7, replace = F))
  samplesProtTumor <- colnames(protTumor_symbol)[!colnames(protTumor_symbol) %in% samplesExprTumor]
  
  exprNormal_symbol_it <- exprNormal_symbol[,samplesExprNormal]
  exprTumor_symbol_it <- exprTumor_symbol[,samplesExprTumor]
  
  protNormal_symbol_it <- protNormal_symbol[,samplesProtNormal]
  protTumor_symbol_it <- protTumor_symbol[,samplesProtTumor]
  
  
  ############## Merge data ######################
  colnames(exprTumor_symbol_it) <- paste0(colnames(exprTumor_symbol_it), "_Tumor")
  colnames(exprNormal_symbol_it) <- paste0(colnames(exprNormal_symbol_it), "_Normal")
  
  commonGenes <- intersect(rownames(exprTumor_symbol_it), rownames(exprNormal_symbol_it))[-1]
  exprMerged <- cbind(exprTumor_symbol_it[commonGenes,], exprNormal_symbol_it[commonGenes,])
  
  exprMeta <- data.frame(Group = c(rep("Tumor", ncol(exprTumor_symbol_it)), rep("Normal", ncol(exprNormal_symbol_it))),
                         Patient = c(colnames(exprTumor_symbol_it), colnames(exprNormal_symbol_it)),
                         row.names = colnames(exprMerged))
  
  colnames(protTumor_symbol_it) <- paste0(colnames(protTumor_symbol_it), "_Tumor")
  colnames(protNormal_symbol_it) <- paste0(colnames(protNormal_symbol_it), "_Normal")
  
  commonGenes <- intersect(rownames(protTumor_symbol_it), rownames(protNormal_symbol_it))
  protMerged <- cbind(protTumor_symbol_it[commonGenes,], protNormal_symbol_it[commonGenes,])
  
  protMeta <- data.frame(Group = c(rep("Tumor", ncol(protTumor_symbol_it)), rep("Normal", ncol(protNormal_symbol_it))),
                         Patient = c(colnames(protTumor_symbol_it), colnames(protNormal_symbol_it)),
                         row.names = colnames(protMerged))
  
  ########## CALCULATE SCORES ####################
  for (score in c("GSVA", "ssGSEA", "Z-score", "Plage", "singscore", "AUCell",
                  "MDT", "MLM", "ORA", "UDT", "ULM",
                  "FGSEA", "norm_FGSEA", "WMEAN", "norm_WMEAN",
                  "corr_WMEAN", "WSUM", "norm_WSUM", "corr_WSUM")) {
    cat("Calculating",score,"score \n")

    score_cores <- ncores
    if (.Platform$OS.type == "windows" && score == "AUCell") {
      score_cores <- 1
      message("Using 1 core for AUCell on Windows because doMC is not available")
    }
    
    exprScores <- getScores(exprMerged, "go_bp", method = score, cores = score_cores)
    protScores <- getScores(protMerged, "go_bp", method = score, cores = score_cores)
    
    # Center scores to 0 by omic type
    exprScores <- t(apply(exprScores, 1, function(x) x-mean(x)))
    protScores <- t(apply(protScores, 1, function(x) x-mean(x)))
    
    commonPathways <- intersect(rownames(na.omit(exprScores)), rownames(na.omit(protScores)))
    exprScores <- exprScores[commonPathways,]
    protScores <- protScores[commonPathways,]
    
    ########### Train with expression ##############
    set.seed(123)
    
    exprModel <- tryCatch(trainModel(exprScores, exprMeta, var2predict = "Group", 
                                     positiveClass = "Tumor", pairingColumn = "Patient",
                                     models = methodsML("svmLinear", "character"), repeatsCV = 1),
                          error = function(e) {NULL})
    
    if(!is.null(exprModel)) {
      prediction <- predictExternal(protScores, exprModel, realValues = protMeta$Group,
                                    positiveClass = "Tumor")
      
      rownames(prediction$stats) <- prediction$stats[,1]
      
      if (is.null(resultsList[[score]])) {
        resultsList[[score]] <- list()
      }
      
      if (is.null(resultsList[[score]][["svmLinear"]])) {
        resultsList[[score]][["svmLinear"]] <- list()
      }
      
      for (metric in prediction$stats[,1]){
        if (is.null(resultsList[[score]][["svmLinear"]][[metric]])) {
          resultsList[[score]][["svmLinear"]][[metric]] <- prediction$stats[metric,2]
        }
        else {
          resultsList[[score]][["svmLinear"]][[metric]] <- c(resultsList[[score]][["svmLinear"]][[metric]], prediction$stats[metric,2])
          
        }
      }
    }
    
    else {
      if (is.null(resultsList[[score]])) {
        resultsList[[score]] <- list()
      }
      
      if (is.null(resultsList[[score]][["svmLinear"]])) {
        resultsList[[score]][["svmLinear"]] <- list()
      }
      
      for (metric in c(
        "mcc", "balacc", "accuracy", "recall", "specificity",
        "npv", "precision", "fscore"
      )){
        if (is.null(resultsList[[score]][["svmLinear"]][[metric]])) {
          resultsList[[score]][["svmLinear"]][[metric]] <- NA
        }
        else {
          resultsList[[score]][["svmLinear"]][[metric]] <- c(resultsList[[score]][["svmLinear"]][[metric]], NA)
          
        }
      }
    }
  }
}


# Table 1
meansDF <- matrix(NA, nrow = length(resultsList[["GSVA"]][["svmLinear"]]), 
                  ncol = length(resultsList), dimnames = list(names(resultsList[["GSVA"]][["svmLinear"]]),
                                                              names(resultsList)))
sdsDF <- meansDF

for (method in colnames(meansDF)) {
  for (score in rownames(meansDF)) {
    meansDF[score, method] <- mean(resultsList[[method]][["svmLinear"]][[score]])
    sdsDF[score, method] <- sd(resultsList[[method]][["svmLinear"]][[score]])
  }
}

meansDF <- meansDF[,order(meansDF["mcc",], decreasing = T)]

meansDF2 <- lapply(seq_len(ncol(meansDF)), function(x) {paste0(as.character(round(unlist(meansDF[,x]), 2)), " (",
                                                               as.character(round(unlist(sdsDF[,x]), 2)), ")")})

meansDF2 <- t(do.call("cbind", meansDF2))

rownames(meansDF2) <- colnames(meansDF)
colnames(meansDF2) <- rownames(meansDF)

colnames(meansDF2) <- c("Accuracy", "Precision", "Recall", "Specificity", "BalAcc", "FScore", "MCC", "NPV")
WriteXLS(meansDF2, "Table 1.xlsx", row.names = T)



# Analysis with several gene sets databases
data(genesetsData)
geneSetsList <- genesetsData

resultsListGenesets <- list()

for (iteration in seq_len(10)) {
  cat(iteration)
  cat("\n")
  
  #### Use different samples for training and validation
  set.seed(iteration)
  
  samplesExprNormal <- sample(colnames(exprNormal_symbol), round(ncol(exprNormal_symbol) * 0.75), replace = F)
  samplesProtNormal <- colnames(protNormal_symbol)[!colnames(protNormal_symbol) %in% samplesExprNormal]
  
  samplesOnlyTumor <- setdiff(colnames(exprTumor_symbol), colnames(exprNormal_symbol))
  
  samplesExprTumor <- c(samplesExprNormal, sample(samplesOnlyTumor, 7, replace = F))
  samplesProtTumor <- colnames(protTumor_symbol)[!colnames(protTumor_symbol) %in% samplesExprTumor]
  
  exprNormal_symbol_it <- exprNormal_symbol[,samplesExprNormal]
  exprTumor_symbol_it <- exprTumor_symbol[,samplesExprTumor]
  
  protNormal_symbol_it <- protNormal_symbol[,samplesProtNormal]
  protTumor_symbol_it <- protTumor_symbol[,samplesProtTumor]
  
  
  ############## Merge data ######################
  colnames(exprTumor_symbol_it) <- paste0(colnames(exprTumor_symbol_it), "_Tumor")
  colnames(exprNormal_symbol_it) <- paste0(colnames(exprNormal_symbol_it), "_Normal")
  
  commonGenes <- intersect(rownames(exprTumor_symbol_it), rownames(exprNormal_symbol_it))[-1]
  exprMerged <- cbind(exprTumor_symbol_it[commonGenes,], exprNormal_symbol_it[commonGenes,])
  
  exprMeta <- data.frame(Group = c(rep("Tumor", ncol(exprTumor_symbol_it)), rep("Normal", ncol(exprNormal_symbol_it))),
                         Patient = c(colnames(exprTumor_symbol_it), colnames(exprNormal_symbol_it)),
                         row.names = colnames(exprMerged))
  
  colnames(protTumor_symbol_it) <- paste0(colnames(protTumor_symbol_it), "_Tumor")
  colnames(protNormal_symbol_it) <- paste0(colnames(protNormal_symbol_it), "_Normal")
  
  commonGenes <- intersect(rownames(protTumor_symbol_it), rownames(protNormal_symbol_it))
  protMerged <- cbind(protTumor_symbol_it[commonGenes,], protNormal_symbol_it[commonGenes,])
  
  protMeta <- data.frame(Group = c(rep("Tumor", ncol(protTumor_symbol_it)), rep("Normal", ncol(protNormal_symbol_it))),
                         Patient = c(colnames(protTumor_symbol_it), colnames(protNormal_symbol_it)),
                         row.names = colnames(protMerged))
  
  ########## CALCULATE SCORES ####################
  for (geneSet in names(geneSetsList)) {
    exprScores <- getScores(exprMerged, geneSetsList[[geneSet]], method = "Z-score")
    protScores <- getScores(protMerged, geneSetsList[[geneSet]], method = "Z-score")
    
    # Center scores to 0 by omic type
    exprScores <- t(apply(exprScores, 1, function(x) x-mean(x)))
    protScores <- t(apply(protScores, 1, function(x) x-mean(x)))
    
    commonPathways <- intersect(rownames(na.omit(exprScores)), rownames(na.omit(protScores)))
    exprScores <- exprScores[commonPathways,]
    protScores <- protScores[commonPathways,]
    
    ########### Train with expression ##############
    set.seed(123)
    exprModel <- tryCatch(trainModel(exprScores, exprMeta, var2predict = "Group", 
                                     positiveClass = "Tumor", pairingColumn = "Patient",
                                     models = methodsML("svmLinear", "character"), repeatsCV = 1),
                          error = function(e) {NULL})
    
    if(!is.null(exprModel)) {
      prediction <- predictExternal(protScores, exprModel, realValues = protMeta$Group,
                                    positiveClass = "Tumor")
      
      rownames(prediction$stats) <- prediction$stats[,1]
      
      if (is.null(resultsListGenesets[[geneSet]])) {
        resultsListGenesets[[geneSet]] <- list()
      }
      
      for (metric in prediction$stats[,1]){
        if (is.null(resultsListGenesets[[geneSet]][[metric]])) {
          resultsListGenesets[[geneSet]][[metric]] <- prediction$stats[metric,2]
        }
        else {
          resultsListGenesets[[geneSet]][[metric]] <- c(resultsListGenesets[[geneSet]][[metric]], prediction$stats[metric,2])
          
        }
      }
    }
    
    else {
      if (is.null(resultsListGenesets[[geneSet]])) {
        resultsListGenesets[[geneSet]] <- list()
      }
      
      if (is.null(resultsListGenesets[[geneSet]])) {
        resultsListGenesets[[geneSet]] <- list()
      }
      
      for (metric in c(
        "mcc", "balacc", "accuracy", "recall", "specificity",
        "npv", "precision", "fscore"
      )){
        if (is.null(resultsListGenesets[[geneSet]][[metric]])) {
          resultsListGenesets[[geneSet]][[metric]] <- NA
        }
        else {
          resultsListGenesets[[geneSet]][[metric]] <- c(resultsListGenesets[[geneSet]][[metric]], NA)
          
        }
      }
    }
  }
}


meansList <- list()
sdsList <- list()
for (metric in c(
  "mcc", "balacc", "accuracy", "recall", "specificity",
  "npv", "precision", "fscore"
)){
  meansList[[metric]] <- lapply(resultsListGenesets, function(x) {mean(x[[metric]])})
  sdsList[[metric]] <- lapply(resultsListGenesets, function(x) {sd(x[[metric]])})
}

meansDF <- do.call(cbind, meansList)
sdsDF <- do.call(cbind, sdsList)

mergedFormatted <- as.data.frame(
  mapply(function(meanCol, sdCol) {
    sprintf("%.2f (%.2f)", round(unlist(meanCol), 2), round(unlist(sdCol), 2))
  }, meansList, sdsList, SIMPLIFY = FALSE)
)
rownames(mergedFormatted) <- rownames(meansDF)
WriteXLS(mergedFormatted, ExcelFileName = "SupplTables/Supplementary Table 1.xlsx", row.names = T, col.names = T)



## External validation
# Download supplementary table 2 from Satpathy et al. (DOI: 10.1016/j.ccell.2025.07.011)
# https://www.cell.com/cms/10.1016/j.ccell.2025.07.011/attachment/0b7d3069-ccef-47d8-b983-cb7cd036e362/mmc3.xlsx


# Prepare external data
extData <- as.data.frame(read_excel("mmc3.xlsx", sheet=2))
extData <- extData[,c(8, 231:415, 629:802)]
extData <- extData[!duplicated(extData[,1]),]
rownames(extData) <- extData[,1]
extData <- extData[,-1]

extData[] <- lapply(extData, function(x) {
  if (is.character(x) || is.factor(x)) {
    as.numeric(as.character(x))
  } else {
    x
  }
})

colnames(exprTumor_symbol) <- paste0(colnames(exprTumor_symbol), "_Tumor")
colnames(exprNormal_symbol) <- paste0(colnames(exprNormal_symbol), "_Normal")

commonGenes <- intersect(rownames(exprTumor_symbol), rownames(exprNormal_symbol))[-1]
exprMerged <- cbind(exprTumor_symbol[commonGenes,], exprNormal_symbol[commonGenes,])

exprMeta <- data.frame(Group = c(rep("Tumor", ncol(exprTumor_symbol)), rep("Normal", ncol(exprNormal_symbol))),
                       Patient = c(colnames(exprTumor_symbol), colnames(exprNormal_symbol)),
                       row.names = colnames(exprMerged))



extMeta <- data.frame(Group = sub(".*\\.", "", colnames(extData)),
                      row.names = colnames(extData))
extMeta[extMeta[,1] == "T",1] <- "Tumor"
extMeta[extMeta[,1] == "N",1] <- "Normal"

extScores <- getScores(extData, "go_bp", "Z-score")


exprScores <- getScores(exprMerged, "go_bp", method = "Z-score")

# Center scores to 0 by omic type
exprScores <- t(apply(exprScores, 1, function(x) x-mean(x)))
extScores <- t(apply(extScores, 1, function(x) x-mean(x)))

commonPathways <- intersect(rownames(na.omit(exprScores)), rownames(na.omit(extScores)))
exprScores <- exprScores[commonPathways,]
extScores <- extScores[commonPathways,]

# Train with expression 
set.seed(123)

exprModel <- trainModel(exprScores, exprMeta, var2predict = "Group", 
                        positiveClass = "Tumor", pairingColumn = "Patient",
                        models = methodsML("svmLinear", "character"), repeatsCV = 1)

externalRealValues <- extMeta$Group
names(externalRealValues) <- rownames(extMeta)

# Predict with proteomics from external cohort

externalPrediction <- predictExternal(extScores, exprModel,
                                      realValues = externalRealValues,
                                      positiveClass = "Tumor"
                                      )

print(externalPrediction$stats)





