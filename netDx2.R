source("./R/min_max_norm.R")
library("plotly")
library("mlr")
library("netDx")
library("curatedTCGAData")
library(Rtsne)

cat("Load data")
cnv <- readRDS("/Users/a.pennati/Desktop/netDx-multiomics-classification/dataset/cnv.rds")
mirna <- readRDS("/Users/a.pennati/Desktop/netDx-multiomics-classification/dataset/mirna.rds")
mrna <- readRDS("/Users/a.pennati/Desktop/netDx-multiomics-classification/dataset/mrna.rds")
proteins <- readRDS("/Users/a.pennati/Desktop/netDx-multiomics-classification/dataset/proteins.rds")

labels_pfi <- readRDS("/Users/a.pennati/Desktop/netDx-multiomics-classification/dataset/labels_pfi.rds")

################################################################################
cat("Create dataframe pheno")
pheno <- data.frame(ID = rownames(cnv), STATUS = as.character(labels_pfi))
rownames(pheno) <- rownames(cnv)

################################################################################
cat("Class balance")
fig <- plot_ly(x = c("0", "1"),                                    
               y = c(sum(pheno$STATUS == "0"), sum(pheno$STATUS == "1")),
               color = c('0', '1'),
               type = "bar")

fig <- fig %>% layout(title = "Class Balance",
                      xaxis = list(title = "label"),
                      yaxis = list(title = "occurences"))

fig

################################################################################
cat("Removing the constants features")
cnv <- as.data.frame(cnv)
cnv <- removeConstantFeatures(cnv, 0, show.info = FALSE)

mirna <- as.data.frame(mirna)
mirna <- removeConstantFeatures(mirna, 0, show.info = FALSE)

mrna <- as.data.frame(mrna)
mrna <- removeConstantFeatures(mrna, 0, show.info = FALSE)

proteins <- as.data.frame(proteins)
proteins <- removeConstantFeatures(proteins, 0, show.info = FALSE)

################################################################################
cat("Tranform data in MultiAssayExperiment object")
brcaList <- list(t(mrna), t(mirna), t(proteins), t(cnv), pheno)
names(brcaList) <- c("mrna", "mirna", "proteins", "cnv", "pheno")
brca <- convertToMAE(brcaList)

cat("Grouping variables to define features")
groupList <- list();

for(i in 1:length(experiments(brca))){
  tmp <- list(rownames(experiments(brca)[[i]]));
  names(tmp) <- names(brca)[i]
  groupList[[names(brca)[[i]]]] <- tmp
}

################################################################################
cat("Define patient similarity for each network")
sims <- list(
  mrna="sim.eucscale",
  mirna="sim.eucscale",
  proteins="sim.eucscale",
  cnv="sim.eucscale"
)

################################################################################
cat("Holdout validation set")
set.seed(123)
dsets <- subsampleValidationData(brca,pctValidation=0.1)
brca <- dsets$trainMAE
holdout <- dsets$validationMAE

setwd("/Users/a.pennati/Desktop/netDx-multiomics-classification/")
outDir <- "/Users/a.pennati/Desktop/netDx-multiomics-classification/pred"
nco <- round(parallel::detectCores())
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)

################################################################################
cat("Model training")
t0 <- Sys.time()
model <-buildPredictor(
  dataList=brca,          
  groupList=groupList,    
  sims=sims,
  outDir=outDir,          
  trainProp=0.8,          
  numSplits=2L,            
  featSelCutoff=1L,       
  featScoreMax=2L,
  numCores=nco,            
  debugMode=FALSE,
  keepAllData=FALSE,     
  logging="default" ) 

t1 <- Sys.time()
print(t1-t0)

################################################################################
cat("Confusion Matrix")
confMat <- confusionMatrix(model)

################################################################################
cat("Examine results")
results <- getResults(
  model,
  unique(colData(brca)$STATUS),
  featureSelCutoff=2L,
  featureSelPct=0.50
)

################################################################################
cat("Validate on independent samples")
outDir <- paste(tempdir(), randAlphanumString(), 
                sep = getFileSep())
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
dir.create(outDir)

predModel <- suppressMessages(
  predict(trainMAE=brca, testMAE=holdout, 
          groupList=groupList, 
          selectedFeatures=results$selectedFeatures, 
          sims=sims,
          outDir=outDir, verbose = FALSE)
)

################################################################################
cat("Plot results of validation")
perf <- getPerformance(predModel, 
                       unique(colData(brca)$STATUS))

plotPerf_multi(list(perf$rocCurve),
               plotTitle = sprintf(
                 "Validation: %i samples", 
                 nrow(colData(holdout))))

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
