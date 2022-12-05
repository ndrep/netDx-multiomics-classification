cat("Load data")
cnv <- readRDS("~/Desktop/netDx/dataset/cnv.rds")
labels_pfi <- readRDS("~/Desktop/netDx/dataset/labels_pfi.rds")
mirna <- readRDS("~/Desktop/netDx/dataset/mirna.rds")
mrna <- readRDS("~/Desktop/netDx/dataset/mrna.rds")
proteins <- readRDS("~/Desktop/netDx/dataset/proteins.rds")

cat("Normalize mrna, mirna and proteins with min_max_norm")
source("~/Desktop/netDx/R/min_max_norm.R")
mirna <- min_max_norm(mirna)
mrna <- min_max_norm(mrna)
proteins <- min_max_norm(proteins)

cat("Create dataframe pheno")
pheno <- data.frame(ID = c(rownames(cnv)), STATUS = as.character(c(labels_pfi)))
rownames(pheno) <- rownames(cnv)

cat("Class balance")
fig <- plot_ly(x = c("0", "1"),                                    
        y = c(sum(pheno$STATUS == 0), sum(pheno$STATUS == 1)),
        color = c('0', '1'),
        type = "bar")

fig <- fig %>% layout(title = "Class Balance",
                      xaxis = list(title = "label"),
                      yaxis = list(title = "occurences"))

fig


cat("Removing the constants features")
cnv <- as.data.frame(cnv)
cnv <- removeConstantFeatures(cnv, 0, show.info = FALSE)

mirna <- as.data.frame(mirna)
mirna <- removeConstantFeatures(mirna, 0, show.info = FALSE)

mrna <- as.data.frame(mrna)
mrna <- removeConstantFeatures(mrna, 0, show.info = FALSE)

proteins <- as.data.frame(proteins)
proteins <- removeConstantFeatures(proteins, 0, show.info = FALSE)

cat("Tranform data in MultiAssayExperiment object")
brcaList <- list(t(mrna), t(mirna), t(proteins), t(cnv), pheno)
names(brcaList) <- c("mrna", "mirna", "proteins", "cnv", "pheno")
brca <- convertToMAE(brcaList)

cat("Grouping variables to define features")
groupList <- list();
pathFile <- sprintf("%s/extdata/pathways.gmt", path.package("netDx"))
pathwayList <- readPathways(pathFile,MIN_SIZE=5L,MAX_SIZE=300L)

gene_names=rownames(brcaList[["mrna"]])
netList <- foreach(k=1:length(pathwayList)) %do% {
  idx <- which(gene_names %in% pathwayList[[k]])
}

n_genes_in_pathway = sapply(netList,FUN=length)
pathwayList=pathwayList[n_genes_in_pathway > 1]
groupList[["mrna"]] <- pathwayList

for(i in 2:length(experiments(brca))){
  tmp <- list(rownames(experiments(brca)[[i]]));
  names(tmp) <- names(brca)[i]
  groupList[[names(brca)[[i]]]] <- tmp
}

cat("Define patient similarity for each network")
sims <- list(
  mrna="pearsonCorr",
  mirna="pearsonCorr",
  proteins="pearsonCorr",
  cnv="pearsonCorr"
)

cat("Holdout validation set")
set.seed(123)
dsets <- subsampleValidationData(brca,pctValidation=0.1)
brca <- dsets$trainMAE
holdout <- dsets$validationMAE

setwd("~/Desktop/netDx/")
outDir <- "/Users/andreapennati/Desktop/netDx/pred"
nco <- round(parallel::detectCores()*0.75)
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)

cat("Model training")
t0 <- Sys.time()
model <-buildPredictor(
  dataList=brca,          
  groupList=groupList,    
  sims=sims,
  outDir=outDir,          
  trainProp=0.8,          
  numSplits=10L,            
  featSelCutoff=1L,       
  featScoreMax=1L,    
  numCores=nco,            
  debugMode=FALSE,
  keepAllData=FALSE,     
  logging="default" ) 

t1 <- Sys.time()
print(t1-t0)

cat("Examine results")
results <- getResults(
  model,
  unique(colData(brca)$STATUS),
  featureSelCutoff=1L,
  featureSelPct=0.50
)

cat("Validate on independent samples")
predModel <- suppressMessages(
  predict(trainMAE=brca, testMAE=holdout, 
          groupList=groupList, 
          selectedFeatures=results$selectedFeatures, 
          sims=sims,
          outDir=outDir, verbose = FALSE)
)

cat("Plot results of validation")
perf <- getPerformance(predModel, 
                       unique(colData(brca)$STATUS))

plotPerf_multi(list(perf$prCurve), 
               plotType = "PR",
               plotTitle = sprintf(
                 "Validation with %i samples", 
                 nrow(colData(holdout))))

psn <- suppressMessages(getPSN(
  brca,
  groupList,
  sims=sims,
  selectedFeatures=results$selectedFeatures
))

tsne <- tSNEPlotter(
  psn$patientSimNetwork_unpruned, 
  colData(brca)
)
