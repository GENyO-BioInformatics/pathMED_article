# Clone the github https://github.com/cclab-brca/neoadjuvant-therapy-response-predictor
system("git clone https://github.com/cclab-brca/neoadjuvant-therapy-response-predictor")

# Process data (adapted from https://github.com/cclab-brca/neoadjuvant-therapy-response-predictor)

#load packages
library(data.table)
library(edgeR)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(readxl)
library(pathMED)
library(WriteXLS)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Create output folders
dir.create("figures")
dir.create("SupplTables")

# change the contents of variable baseDir to the root analysis folder
baseDir <- "neoadjuvant-therapy-response-predictor-master/"
# load directory structure downloaded from github site
source(paste0(baseDir,"/R/LoadDirectoryStructure.R"))
setwd("..")

# load clinical metadata (Supplementary Table 1)
metadata  <- data.frame(read_excel(paste0(suppTableDir, "/Supplementary Tables.xlsx"),sheet = 1))

# combine ER and HER2 status
metadata$ERHER2.status <- ifelse(metadata$ER.status=="POS","ER+ HER2-","ER- HER2-")
metadata$ERHER2.status <- ifelse(metadata$HER2.status=="POS","HER2+",metadata$ERHER2.status)
metadataFull <- metadata

# Perform analyses with cases that had an RCB assessment and received more than one cycle of therapy
# as detailed in Methods (Statistical testing)
metadata <- metadata[metadata$RCB.category!="NA",]
metadata <- metadata[metadata$Chemo.cycles>1 & metadata$aHER2.cycles>1,]

# load Gene Ensembl ID to Hugo ID dictionary
ensemblToHugo <- read.table(paste0(resourcesDir,"EnsemblID.to.Hugo.v87.tsv.gz"), header=T, stringsAsFactors = F,sep="\t")

# load RNA data (Supplementary Table 3)
rnadata  <- data.frame(read_excel(paste0(suppTableDir,"Supplementary Tables.xlsx"),sheet = 3))


p <- metadata
p$RCB.score <- as.numeric(p$RCB.score)

# load RNAseq raw counts (Methods)
transneo.counts <- data.frame(fread(paste0(dataDir,"transneo-diagnosis-RNAseq-rawcounts.tsv.gz"),header=T, sep="\t",stringsAsFactors = F),row.names = 1)
transneo.counts <- transneo.counts[,colnames(transneo.counts) %in% p$Donor.ID]

# update metadata - retain only samples that have RNAseq data
p <- p[p$Donor.ID %in% colnames(transneo.counts),]
rownames(p) <- p$Donor.ID

# Filtering and normalization
y <- DGEList(transneo.counts)
minCPM      <- 1
minNoToKeep <- 10
keep        <- rowSums(cpm(y)>minCPM)>=minNoToKeep
y <- y[keep , , keep.lib.sizes=FALSE]
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y, method = "TMM")

exprNormalized <- cpm(y, log=F)
logExpr <- log2(exprNormalized+1)

# Translate to Gene Symbols
entrez <- read.table(paste0(resourcesDir,"EnsemblID.to.Entrez.tsv.gz"), sep="\t",header = T,stringsAsFactors = F)
entrez <- entrez[!is.na(entrez$EntrezGene.ID),]
symbol    <- na.omit(mapIds(EnsDb.Hsapiens.v86, keys=row.names(logExpr), column="GENENAME", keytype="GENEID", multiVals="first"))
symbol <- symbol[!duplicated(symbol)]
logExpr <- logExpr[names(symbol),]
rownames(logExpr) <- symbol

# Gene sets to be used
geneSetsList <- c("go_bp", "go_mf", "reactome", "kegg")

resultsBRCA_binary <- list()
modelsBRCA_binary <- list()
scoresBRCA_binary <- list()


# Calculate scores and train/validate predictive models of the binary pCR/RD outcome
for (geneSet in geneSetsList) {

    # Calculate scores
    exprScores <- getScores(logExpr, geneSet, method = "GSVA", cores = 10)
    
    # Get the top 1000 by variance
    varianceScores <- na.omit(apply(exprScores, 1, function(x) {var(x, na.rm = T)}))
    nFeatures <- min(length(varianceScores), 1000)
    varFeatures <- sort(varianceScores, decreasing = T)[1:nFeatures]
    
    exprScores <- exprScores[names(varFeatures),]

    ########### Train with expression ##############
    set.seed(123)
    exprModel <- trainModel(exprScores, p, var2predict = "pCR.RD", 
                            positiveClass = "pCR", Koutter = 3, Kinner = 3,
                            models = methodsML("lda", 
                                "character", tuneLength = 100))
    
    # Save results
    if (is.null(resultsBRCA_binary[[geneSet]])) {
        resultsBRCA_binary[[geneSet]] <- list()
    }
    
    
    for (metric in rownames(exprModel$stats)){
        if (is.null(resultsBRCA_binary[[geneSet]][[metric]])) {
            resultsBRCA_binary[[geneSet]][[metric]] <- exprModel$stats[metric,1]
        }
        else {
            resultsBRCA_binary[[geneSet]][[metric]] <- c(resultsBRCA_binary[[geneSet]][[metric]], exprModel$stats[metric,1])
            
        }
    }
    
    if (is.null(modelsBRCA_binary[[geneSet]])) {
        modelsBRCA_binary[[geneSet]] <- list()
    }
    
    modelsBRCA_binary[[geneSet]] <- exprModel$model
    scoresBRCA_binary[[geneSet]] <- exprScores
}

# Estimate importance for each gene sets database
importanceList_binary <- list()

for (database in names(modelsBRCA_binary)) {
    importance <- varImp(modelsBRCA_binary[[database]], scale = F)[["importance"]]
    importance_df <- data.frame(Feature = rownames(importance), Importance = importance[,1])
    importance_df <- importance_df[order(-importance_df$Importance), ]
    if (database %in% c("go_bp", "go_mf")) {
        rownames(importance_df) <- gsub(".", ":", importance_df$Feature, fixed=T)
    }
    else if (database == "reactome") {
        rownames(importance_df) <- gsub(".", "-", importance_df$Feature, fixed=T)
    }
    else {
        rownames(importance_df) <- importance_df$Feature
    }
    importance_df$database <- database
    importance_annotated <- ann2term(importance_df)
    rownames(importance_annotated) <- importance_annotated$ID
    importance_df$term <- importance_annotated[rownames(importance_df), "term"]
    importanceList_binary[[database]] <- importance_df
}

mergedImportances_binary <- do.call(rbind, importanceList_binary)

corrected_terms <- c()
for (database in names(importanceList_binary)) {
    corrected_terms <- c(corrected_terms, rownames(importanceList_binary[[database]]))
}

rownames(mergedImportances_binary) <- corrected_terms

# Sort by importance
mergedImportances_binary <- mergedImportances_binary[order(-mergedImportances_binary$Importance), ]
mergedScores_binary <- do.call(rbind, scoresBRCA_binary)


# Heatmap of top terms
topFeatures <- unlist(lapply(importanceList_binary, function(x){return(rownames(x)[1:10])}))

dataHeatmap <- t(scale(t(mergedScores_binary[topFeatures,])))

columnAnnot <- "RCB.category"
columnAnnot2 <- "RCB.score"

# Order columns by the RCB score
colOrder <- rownames(p[order(p[[columnAnnot2]]),])

# Categorical annotation: RCB.category
annotCol <- factor(p[[columnAnnot]])
names(annotCol) <- rownames(p)
n <- length(levels(annotCol))
samples_colors <- setNames(brewer.pal(min(n, 8), "Set2")[1:n], levels(annotCol))

annotDB <- factor(mergedImportances_binary[topFeatures, "database"])
databases_colors <- setNames(brewer.pal(8, "Pastel1")[1:length(levels(annotDB))], levels(annotDB))

# Continuous annotation: RCB.score
annotCol2 <- as.numeric(p[[columnAnnot2]])
names(annotCol2) <- rownames(p)

RCBScore_range <- range(annotCol2, na.rm = TRUE)

RCBScore_col_fun <- colorRamp2(
    breaks = seq(RCBScore_range[2], RCBScore_range[1], length.out = 3),
    colors = plasma(3)
)

# Annotations (no legend for RCB.score here)
col_ha <- HeatmapAnnotation(
    RCB.category = annotCol,
    RCB.score = annotCol2,
    col = list(
        RCB.category = samples_colors,
        RCB.score = RCBScore_col_fun
    ),
    show_annotation_name = T,
    show_legend = F
)
row_ha <- rowAnnotation(
    Database = annotDB,
    col = list(
        Database = databases_colors
    ),
    show_annotation_name = F,
    show_legend = F
)

# Define a separate legend for the RCB.score
rcb_score_legend <- Legend(
    title = "RCB\nScore",
    at = RCBScore_range,
    col_fun = RCBScore_col_fun,
    direction = "vertical"
)

rcb_category_legend <- Legend(
    title = "RCB    \nCategory",
    at = levels(annotCol),
    legend_gp = gpar(fill = samples_colors)
)

database_legend <- Legend(
    title = "Database",
    at = levels(annotDB),
    legend_gp = gpar(fill = databases_colors)
)
# Create the heatmap

rng <- range(dataHeatmap, na.rm = TRUE)

col_fun <- colorRamp2(
    c(rng[1], 0, rng[2]),
    c("#4C72B0", "white", "#DD8452")
)

ht <- Heatmap(dataHeatmap,
              name = "GSVA\nZ-score",
              col = col_fun,
              column_order = colOrder,
              cluster_rows = TRUE,
              top_annotation = col_ha,
              right_annotation   = row_ha,
              show_column_names = FALSE,
              row_labels = mergedImportances_binary[rownames(dataHeatmap), "term"],
              show_row_dend = FALSE,
              row_names_gp = gpar(fontsize = 8)
)

# Draw heatmap and pack legends in one column
pdf("figures/Figure2a.pdf", 8, 5, pointsize = 12)
draw(
    ht,
    annotation_legend_list = list(rcb_category_legend, rcb_score_legend, database_legend),
    heatmap_legend_side = "left",
    annotation_legend_side = "left",
    merge_legends = TRUE,
    padding = unit(c(2, 2, 2, 15), "mm"),
    
)
dev.off()

# Write tables with the results
resultsMerged <- do.call(cbind, resultsBRCA_binary)
colnames(resultsMerged) <- names(resultsBRCA_binary)
resultsMerged <- as.data.frame(t(resultsMerged))
resultsMerged <- resultsMerged[,c(3, 7, 4, 5, 2, 8, 1, 6)]
colnames(resultsMerged) <- c("Accuracy", "Precision", "Recall", "Specificity", "BalAcc", "FScore", "MCC", "NPV")
WriteXLS(resultsMerged, "SupplTables/Supplementary Table 2.xlsx", row.names = T)
WriteXLS(mergedImportances_binary, "SupplTables/Supplementary Table 4.xlsx", row.names = F)



# Repeat the whole analysis predicting the continuous RCB score

resultsBRCA <- list()
modelsBRCA <- list()
scoresBRCA <- list()

for (geneSet in geneSetsList) {
    
    # Calculate scores
    exprScores <- getScores(logExpr, geneSet, method = "GSVA", cores = 10)
    
    # Get the top 1000 by variance
    varianceScores <- na.omit(apply(exprScores, 1, function(x) {var(x, na.rm = T)}))
    nFeatures <- min(length(varianceScores), 1000)
    varFeatures <- sort(varianceScores, decreasing = T)[1:nFeatures]
    
    exprScores <- exprScores[names(varFeatures),]
    
    ########### Train with expression ##############
    set.seed(123)
    exprModel <- trainModel(exprScores, p, var2predict = "RCB.score", 
                            Koutter = 3, Kinner = 3,
                            models = methodsML("rf", 
                                               "numeric", tuneLength = 100))
    
    # Save results
    if (is.null(resultsBRCA[[geneSet]])) {
        resultsBRCA[[geneSet]] <- list()
    }
    
    
    for (metric in rownames(exprModel$stats)){
        if (is.null(resultsBRCA[[geneSet]][[metric]])) {
            resultsBRCA[[geneSet]][[metric]] <- exprModel$stats[metric,1]
        }
        else {
            resultsBRCA[[geneSet]][[metric]] <- c(resultsBRCA[[geneSet]][[metric]], exprModel$stats[metric,1])
            
        }
    }
    
    if (is.null(modelsBRCA[[geneSet]])) {
        modelsBRCA[[geneSet]] <- list()
    }
    
    modelsBRCA[[geneSet]] <- exprModel$model
    scoresBRCA [[geneSet]] <- exprScores
}


importanceList <- list()

for (database in names(modelsBRCA)) {
    importance <- varImp(modelsBRCA[[database]], scale = F)[["importance"]]
    importance_df <- data.frame(Feature = rownames(importance), Importance = importance[,1])
    importance_df <- importance_df[order(-importance_df$Importance), ]
    if (database %in% c("go_bp", "go_cc", "go_mf", "hpo")) {
        rownames(importance_df) <- gsub(".", ":", importance_df$Feature, fixed=T)
    }
    else if (database %in% c("reactome", "lincs")) {
        rownames(importance_df) <- gsub(".", "-", importance_df$Feature, fixed=T)
    }
    else {
        rownames(importance_df) <- importance_df$Feature
    }
    importance_df$database <- database
    importance_annotated <- ann2term(importance_df)
    rownames(importance_annotated) <- importance_annotated$ID
    importance_df$term <- importance_annotated[rownames(importance_df), "term"]
    importanceList[[database]] <- importance_df
}

mergedImportances <- do.call(rbind, importanceList)
mergedScores <- do.call(rbind, scoresBRCA)

corrected_terms <- c()
for (database in names(importanceList)) {
    corrected_terms <- c(corrected_terms, rownames(importanceList[[database]]))
}

rownames(mergedImportances) <- corrected_terms

mergedImportances <- mergedImportances[order(-mergedImportances$Importance), ]
mergedScores <- do.call(rbind, scoresBRCA)



####### Heatmap of top terms
topFeatures <- unlist(lapply(importanceList, function(x){return(rownames(x)[1:10])}))

dataHeatmap <- t(scale(t(mergedScores[topFeatures,])))



columnAnnot <- "RCB.category"
columnAnnot2 <- "RCB.score"

# Order columns by the continuous score
colOrder <- rownames(p[order(p[[columnAnnot2]]),])

# Categorical annotation: RCB.category
annotCol <- factor(p[[columnAnnot]])
names(annotCol) <- rownames(p)
n <- length(levels(annotCol))
samples_colors <- setNames(brewer.pal(min(n, 8), "Set2")[1:n], levels(annotCol))

annotDB <- factor(mergedImportances[topFeatures, "database"])
databases_colors <- setNames(brewer.pal(8, "Pastel1")[1:length(levels(annotDB))], levels(annotDB))

# Continuous annotation: RCB.score
annotCol2 <- p[[columnAnnot2]]
names(annotCol2) <- rownames(p)

RCBScore_col_fun <- colorRamp2(
    breaks = seq(RCBScore_range[2], RCBScore_range[1], length.out = 3),
    colors = plasma(3)
)

# Annotations (no legend for RCB.score here)
col_ha <- HeatmapAnnotation(
    RCB.category = annotCol,
    RCB.score = annotCol2,
    col = list(
        RCB.category = samples_colors,
        RCB.score = RCBScore_col_fun
    ),
    show_annotation_name = T,
    show_legend = F
)
row_ha <- rowAnnotation(
    Database = annotDB,
    col = list(
        Database = databases_colors
    ),
    show_annotation_name = F,
    show_legend = F
)

# Define a separate legend for the RCB.score
rcb_score_legend <- Legend(
    title = "RCB\nScore",
    at = RCBScore_range,
    col_fun = RCBScore_col_fun,
    direction = "vertical"
)

rcb_category_legend <- Legend(
    title = "RCB    \nCategory",
    at = levels(annotCol),
    legend_gp = gpar(fill = samples_colors)
)

database_legend <- Legend(
    title = "Database",
    at = levels(annotDB),
    legend_gp = gpar(fill = databases_colors)
)
# Create the heatmap

rng <- range(dataHeatmap, na.rm = TRUE)

col_fun <- colorRamp2(
    c(rng[1], 0, rng[2]),
    c("#4C72B0", "white", "#DD8452")
)

ht <- Heatmap(dataHeatmap,
              name = "GSVA\nZ-score",
              col = col_fun,
              column_order = colOrder,
              cluster_rows = TRUE,
              top_annotation = col_ha,
              right_annotation   = row_ha,
              show_column_names = FALSE,
              row_labels = mergedImportances[rownames(dataHeatmap), "term"],
              show_row_dend = FALSE,
              row_names_gp = gpar(fontsize = 8)
)

# Draw heatmap and pack legends in one column
pdf("figures/Figure2b.pdf", 8, 5, pointsize = 12)
draw(
    ht,
    annotation_legend_list = list(rcb_category_legend, rcb_score_legend, database_legend),
    heatmap_legend_side = "left",
    annotation_legend_side = "left",
    merge_legends = TRUE,
    padding = unit(c(2, 2, 2, 10), "mm")
    )
dev.off()

resultsMerged <- do.call(cbind, resultsBRCA)
colnames(resultsMerged) <- names(resultsBRCA)
resultsMerged <- as.data.frame(t(resultsMerged))
resultsMerged <- resultsMerged[,-7]
colnames(resultsMerged)[1] <- c("R")
WriteXLS(resultsMerged, "SupplTables/Supplementary Table 3.xlsx", row.names = T)
WriteXLS(mergedImportances, "SupplTables/Supplementary Table 5.xlsx", row.names = F)
