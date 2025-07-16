# Load necessary libraries
library(pathMED)
library(biomaRt)
library(ggplot2)
library(WriteXLS)

# Download data from https://kb.linkedomics.org/download#LUAD and save into 
# LinkedOmicsKB folder

# Create output folders
dir.create("SupplTables")
dir.create("figures")

# Read and prepare the metadata
metadata <- read.delim("LinkedOmicsKB_LUAD/LUAD_meta.txt")
metadata <- metadata[-1,]
rownames(metadata) <- metadata[,1]

pheno <- read.delim("LinkedOmicsKB_LUAD/LUAD_phenotype.txt")

# Read and process omics data
exprTumor <- read.delim("LinkedOmicsKB_LUAD/LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt", row.names = 1)
protTumor <- read.delim("LinkedOmicsKB_LUAD/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt", row.names = 1)

## Discard genes in pseudoautosomal regions
exprTumor <- exprTumor[grep("PAR", rownames(exprTumor), invert = T),]

# Remove ENSEMBL version numbers
rownames(exprTumor) <- gsub("\\..*", "", rownames(exprTumor))
rownames(protTumor) <- gsub("\\..*", "", rownames(protTumor))

# Translate ENSEMBL IDs to Gene Symbols

## Transcriptomics
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

## Proteomics
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
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


resultsList <- list()

for (iteration in seq_len(10)) {

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
    for (score in c("GSVA", "ssGSEA", "Z-score", "Plage", "singscore", "AUCell",
                    "MDT", "MLM", "ORA", "UDT", "ULM",
                    "FGSEA", "norm_FGSEA", "WMEAN", "norm_WMEAN",
                    "corr_WMEAN", "WSUM", "norm_WSUM", "corr_WSUM")) {

        exprScores <- getScores(exprMerged, "go_bp", method = score, cores = 10)
        protScores <- getScores(protMerged, "go_bp", method = score, cores = 10)
        
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


# Figure 1
dataBarplot <- matrix(NA, nrow = length(resultsList[["GSVA"]][["svmLinear"]]), 
                      ncol = length(resultsList), dimnames = list(names(resultsList[["GSVA"]][["svmLinear"]]),
                                                                  names(resultsList)))
for (method in colnames(dataBarplot)) {
    for (score in rownames(dataBarplot)) {
        dataBarplot[score, method] <- mean(resultsList[[method]][["svmLinear"]][[score]])
    }
}
dataBarplot <- apply(dataBarplot, 2, function(x) {x[is.na(x)] <- 0; return(x)})
dataBarplot <- data.frame(dataBarplot)

dataSuppl <- data.frame(t(dataBarplot))[]
colnames(dataSuppl) <- c("Accuracy", "Precision", "Recall", "Specificity", "BalAcc", "FScore", "MCC", "NPV")

WriteXLS(dataSuppl, "SupplTables/Supplementary Table 1.xlsx", row.names = T)

topMethods <- names(sort(dataBarplot["accuracy",], decreasing = T))[1:5]
dataBarplot <- dataBarplot[c("accuracy", "recall", "precision", "specificity", "fscore", "mcc"),topMethods]

dataBarplot <- stack(dataBarplot)
dataBarplot$metric <- factor(rep(c("Accuracy", "Recall", "Precision",  "Specificity", "F1", "MCC"), 5))

p <- ggplot(dataBarplot, aes(x=factor(metric, c("Accuracy", "Recall", "Precision",  "Specificity", "F1", "MCC")), y=values, fill=ind)) +
    geom_bar(position = "dodge", stat = "identity", color = "black") +
    xlab("Metric") +
    ylab("Value") +
    theme_classic() +
    scale_fill_jco() +
    ylim(0,1) +
    labs(fill="Method") +
    theme(legend.position = "top", legend.key.size = unit(0.02, "npc"),
          legend.text = element_text(size=8),
          legend.title = element_text(size=8))

ggsave("figures/Figure2.pdf", p, scale=1.2, width = 3.5, height = 2.4)



# Analysis with several gene sets databases
data(genesetsData)
geneSetsList <- genesetsData

resultsList <- list()

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
            
            if (is.null(resultsList[[geneSet]])) {
                resultsList[[geneSet]] <- list()
            }
            
            for (metric in prediction$stats[,1]){
                if (is.null(resultsList[[geneSet]][[metric]])) {
                    resultsList[[geneSet]][[metric]] <- prediction$stats[metric,2]
                }
                else {
                    resultsList[[geneSet]][[metric]] <- c(resultsList[[geneSet]][[metric]], prediction$stats[metric,2])
                    
                }
            }
        }
        
        else {
            if (is.null(resultsList[[geneSet]])) {
                resultsList[[geneSet]] <- list()
            }
            
            if (is.null(resultsList[[geneSet]])) {
                resultsList[[geneSet]] <- list()
            }
            
            for (metric in c(
                "mcc", "balacc", "accuracy", "recall", "specificity",
                "npv", "precision", "fscore"
            )){
                if (is.null(resultsList[[geneSet]][[metric]])) {
                    resultsList[[geneSet]][[metric]] <- NA
                }
                else {
                    resultsList[[geneSet]][[metric]] <- c(resultsList[[geneSet]][[metric]], NA)
                    
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
    meansList[[metric]] <- lapply(resultsList, function(x) {mean(x[[metric]])})
    sdsList[[metric]] <- lapply(resultsList, function(x) {sd(x[[metric]])})
}

meansDF <- do.call(cbind, meansList)
sdsDF <- do.call(cbind, sdsList)
meansDF2 <- lapply(seq_len(8), function(x) {paste0(as.character(round(unlist(meansDF[,x], 4))), " (",
                                                   as.character(round(unlist(sdsDF[,x], 4))), ")")})

mergedFormatted <- as.data.frame(
    mapply(function(meanCol, sdCol) {
        sprintf("%.2f (%.2f)", round(unlist(meanCol), 2), round(unlist(sdCol), 2))
    }, meansList, sdsList, SIMPLIFY = FALSE)
)
rownames(mergedFormatted) <- rownames(meansDF)
WriteXLS(mergedFormatted, ExcelFileName = "SupplTables/Supplementary Table 2.xlsx", row.names = T, col.names = T)

