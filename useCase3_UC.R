# Load dependencies
set.seed(12345)
library(GEOquery)
library(readxl)
library(tidyr)
library(parallel)
library(pathMED)
library(ggplot2)
library(pheatmap)

## Load and pre-process data ----

getGEOSuppFiles("GSE206171")
clin <- data.frame(read_xlsx("GSE206171/GSE206171_Array_data_GEO.xlsx", range="A2:EX11"))
clin <- clin[,-2]
rownames(clin) <- clin[,1]
clin <- clin[,-1]
clin <- as.data.frame(t(clin))

clin<-clin[colnames(exprData),c("Mayscore","Disease or control")]
colnames(clin)<-c("MayScore","Disease")

## Gene expression
exprData <- data.frame(read_xlsx("GSE206171/GSE206171_Array_data_GEO.xlsx", range="A12:EX17693"))
rownames(exprData) <- exprData[,1]
colnames(exprData)[1] <- "probeid"

exprData <- exprData[!(grepl(pattern = "NA", exprData[,1]) | is.na(exprData[,1])),]

genome <- exprData[,c("probeid","genesymbol")]
exprData <- exprData[,grepl("h0",colnames(exprData))]

## Annotation to gene symbol
genome <- genome %>% `colnames<-`(c("fromGenes","toGenes")) %>% 
  replace(.=="",NA) %>% drop_na()

exprData <- exprData[genome$fromGenes,]
exprData <- exprData[!duplicated(genome$toGenes),]
rownames(exprData) <- genome$toGenes[!duplicated(genome$toGenes)]


## Get Genesets database to use (Reactome with and without dissectDB) ----

dataList <- buildRefObject(data=list(exprData), metadata = list(clin),
                         groupVar = "Disease", controlGroup = "Control")

data(genesetsData)
REAC <- genesetsData[["reactome"]]

REAC.d <- dissectDB(refObject = dataList, geneSets = "reactome",
                minPathSize = 8, minSplitSize = 3, maxSplits = NULL,
                explainedVariance = 70, percSharedGenes = 90)

## Remove co-linear genesets
dbs <- list("reactome" = REAC, "reactome.f" = REAC.d)
dbs.f <- lapply(dbs,function(geneSets){
  net <- do.call("rbind",lapply(seq_len(length(geneSets)),function(x){
    res <- data.frame("source" = rep(names(geneSets)[[x]],length(geneSets[[x]])),
                      "target" = as.character(geneSets[[x]]),
                      "weight" = rep(1, length(geneSets[[x]])),
                      "mor" = rep(1, length(geneSets[[x]])))
    return(res)
  }))
  
  co.lin <- as.data.frame(decoupleR::check_corr(net,.source = "source",
                                                .target = "target",
                                                .mor = "mor"))
  
  net <- net[!net$source %in% as.character(co.lin[
    co.lin$correlation > 0.5, "source"]), ]
  
  res <- geneSets[names(geneSets) %in% net$source]
  return(res)
})

names(dbs.f) <- names(dbs)

rm(REAC,REAC.d,genesetsData)

## Get Single Sample scores ----

clin$Diagnosis <- ifelse(clin$Disease=="Control","Healthy","Case")

labels <- clin$Diagnosis
names(labels) <- rownames(clin)

algs <- c("M-Scores", "GSVA", "ssGSEA", "singscore", "Plage", "Z-score", "AUCell",
        "MDT", "MLM", "ORA", "UDT", "ULM", "FGSEA", "norm_FGSEA", "WMEAN", 
        "norm_WMEAN", "corr_WMEAN", "WSUM", "norm_WSUM", "corr_WSUM")

SCORES <- list()
SCORES.d <- list()
gc()

for(i in 1:length(algs)){
  alg <- algs[i]
  print(alg)
  
  if(alg=="MLM"){
    gnset <- dbs.f[[1]]
  } else{
    gnset <- dbs[[1]]
  }
  
  ncores <- ifelse(alg %in% c("M-Scores", "MDT", "UDT", "ULM"), detectCores()-2, 1)
  SCORES[[alg]] <- getScores(inputData = exprData, geneSets = gnset, method = alg,
                 labels = labels,cores = ncores)
  gc()
}

for(i in 1:length(algs)){
  alg <- algs[i]
  print(alg)
  
  if(alg=="MLM"){
    gnset <- dbs.f[[2]]
  } else{
    gnset <- dbs[[2]]
  }
  
  ncores <- ifelse(alg %in% c("M-Scores", "MDT", "UDT", "ULM"), detectCores()-2, 1)
  SCORES.d[[alg]] <- getScores(inputData = exprData, geneSets = gnset, method = alg,
                           labels = labels, cores = ncores)
  gc()
}

## Correlation with disease activity ----

ODB <- do.call("rbind", lapply(1:length(SCORES), function(i){
  tmp.scores <- SCORES[[i]]
  tmp.scores <- tmp.scores[,rownames(clin)]
  res <- do.call("rbind",lapply(1:nrow(tmp.scores), function(x){
    x <- as.numeric(tmp.scores[x,])
    
    return(c(cor.test(as.numeric(clin$MayScore), x, method = "pearson")$p.value,
             cor.test(as.numeric(clin$MayScore), x, method = "pearson")$estimate))
  }))
  colnames(res) <- c("pvalue","corr")
  
  res <- data.frame("geneset" = rownames(tmp.scores),
                  res,
                  "method" = names(SCORES)[i],
                  "approach" = "Original")
  return(res)
}))


DDB <- do.call("rbind",lapply(1:length(SCORES.d),function(i){
  tmp.scores <- SCORES.d[[i]]
  tmp.scores <- tmp.scores[,rownames(clin)]
  res <- do.call("rbind", lapply(1:nrow(tmp.scores), function(x){
    x <- as.numeric(tmp.scores[x,])
    
    return(c(cor.test(as.numeric(clin$MayScore), x, method = "pearson")$p.value,
             cor.test(as.numeric(clin$MayScore), x, method = "pearson")$estimate))
  }))
  colnames(res) <- c("pvalue","corr")
  
  res <- data.frame("geneset" = rownames(tmp.scores),
                  res,
                  "method"=names(SCORES.d)[i],
                  "approach" = "Dissected")
  return(res)
}))


## Filtering tables
M <- rbind(ODB, DDB)

methodsList <- unique(DDB$method)

df <- do.call("rbind",lapply(methodsList, function(meth){
  
  tmp <- M[!is.na(M$pvalue),]
  tmp <- tmp[tmp$pvalue<=0.05 & tmp$method==meth,]
  o.db <- tmp[tmp$approach=="Original",]
  d.db <- tmp[tmp$approach!="Original",]
  
  o.db <- o.db[order(o.db$pvalue,decreasing = F),]
  d.db <- d.db[order(d.db$pvalue,decreasing = F),]
  
  if(nrow(o.db) < 100){
      n <- nrow(o.db)
  } else{
      n <- 100
      }
  o.db <- o.db[1:n,]
  if(nrow(d.db) < 100){
      n<-nrow(d.db)
  } else{
          n<-100
          }
  d.db <- d.db[1:n,]
  
  res <- rbind(o.db,d.db)
  res$corr <- abs(res$corr)
  return(res)
}))

rownames(df) <- NULL



## Filtering approaches with mean lower than 0.6

means <- do.call("rbind", lapply(methodsList, function(meth){
  return(c(mean(as.numeric(df[df$method==meth & df$approach=="Original","corr"])),
    mean(as.numeric(df[df$method==meth & df$approach != "Original","corr"]))))
}))
rownames(means) <- methodsList
colnames(means) <- c("oDB","dDB")

selmethodsList <- rownames(means)[!apply(means,1,max) < 0.6]

df <- df[df$method %in% selmethodsList,]
df$method <- factor(df$method, levels = unique(df$method))
df$approach <- factor(df$approach, levels <- c("Original", "Dissected"))
p <- ggplot(df, aes(x = method, y = corr, fill = approach)) +
  geom_jitter(aes(color = approach), 
              position = position_jitterdodge(jitter.width = 1, dodge.width = 0.8), 
              alpha = 0.8, size = 0.75) +
  ylim(0.58,0.9)+
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.6) +
  labs(x = "",
    y = "Abs(top 100 correlations)",
    fill = "",
    color = "") +
  theme_light(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

ggsave("figures/Figure3a.pdf", p, scale=2, width = 4, height = 2.4)




DDB$genesetID <- sub("\\.spl.*", "", DDB$geneset)
ODB$genesetID <- sub("\\.spl.*", "", ODB$geneset)

df_merged <- merge(DDB[DDB$method=="Z-score",], ODB[ODB$method=="Z-score",], by = "genesetID", all.x = TRUE)
df_merged <- df_merged[,c(1:2, 4, 3, 9, 8)]

colnames(df_merged) <- c("genesetID", "subset", "Corr.subset", "PVal.subset", "Corr.original", "PVal.original") 

df_merged$diff <- df_merged$Corr.subset - df_merged$Corr.original
df_merged <- df_merged[order(abs(df_merged$diff),decreasing=T),]


pathwayInt <- "R-HSA-556833"
sets <- df_merged[df_merged$genesetID==pathwayInt,"subset"]
sets <- sets[!is.na(sets)]


tmp <- SCORES.d[["M-Scores"]]
tmp <- tmp[sets,]

clin$MayScore <- as.numeric(clin$MayScore)

p <- pheatmap(tmp[,rownames(clin[order(clin$MayScore,decreasing = F),])],
                   cluster_cols = F,
                   show_colnames = F,border_color = "black",
                   fontsize = 6,annotation_col = clin[,c("MayScore", "Disease")],
                   breaks = seq(-6, 6, length.out = 101),
              width = 4, height = 2.4,
                   main = "R-HSA-556833",
                    colorRampPalette(c("#4575b4", "white", "#d73027"))(100))

ggsave("figures/Figure3b.pdf", p$gtable, scale=2, width = 4, height = 2.4)


supplTable6 <- df_merged[df_merged$genesetID==pathwayInt, -1]
supplTable6$Genes.subset <- sapply(supplTable6$subset, function(x) {paste(REAC.d[[x]], collapse = ", ")})

WriteXLS(supplTable6, "SupplTables/Supplementary Table 6.xlsx", row.names = F)




