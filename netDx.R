library("netDx")
library("curatedTCGAData")

######LOAD DATA#################################################################
cnv <- as.data.frame(readRDS("/Users/a.pennati/Desktop/netDx-multiomics-classification/dataset/cnv.rds"))
mirna <- as.data.frame(readRDS("/Users/a.pennati/Desktop/netDx-multiomics-classification/dataset/mirna.rds"))
mrna <- as.data.frame(readRDS("/Users/a.pennati/Desktop/netDx-multiomics-classification/dataset/mrna.rds"))
proteins <- as.data.frame(readRDS("/Users/a.pennati/Desktop/netDx-multiomics-classification/dataset/proteins.rds"))

labels_pfi <- readRDS("/Users/a.pennati/Desktop/netDx-multiomics-classification/dataset/labels_pfi.rds")


#################CROSS CORRELATION VISUALIZATION####################################
library(lares)

corr_cross(proteins, # name of dataset
           max_pvalue = 0.05, # display only significant correlations (at 5% level)
           top = 10 # display top 10 couples of variables (by correlation coefficient)
)

######NORMALIZE MRNA, MIRNA AND PROTEINS WITH MIN_MAX_NORM######################
source("./R/min_max_norm.R")
mirna <- as.data.frame(round(apply(mirna, 2, min_max_norm), 2))
mrna <- as.data.frame(round(apply(mrna, 2, min_max_norm), 2))
proteins <- as.data.frame(round(apply(proteins, 2, min_max_norm), 2))

################REMOVING FEATURES############################

library("mlr")
cnv <- removeConstantFeatures(cnv, 0, show.info = FALSE)
mirna <- removeConstantFeatures(mirna, 0, show.info = FALSE)
mrna <- removeConstantFeatures(mrna, 0, show.info = FALSE)
proteins <- removeConstantFeatures(proteins, 0, show.info = FALSE)

#########################DROP FEATURE CORRELATED################################
source("./R/drop_feature_correlated.R")
cnv <- drop_feature_correlated(cnv)
mirna <- drop_feature_correlated(mirna)
mrna <- drop_feature_correlated(mrna)
proteins <- drop_feature_correlated(proteins)


############CREATE DATAFRME PHENO##############################

pheno <- data.frame(ID = rownames(cnv), STATUS = as.character(labels_pfi))
rownames(pheno) <- rownames(cnv)

#########################CLASS BALANCE##########################################
library("plotly")
fig <- plot_ly(x = c("0", "1"),                                    
        y = c(sum(pheno$STATUS == "0"), sum(pheno$STATUS == "1")),
        color = c('0', '1'),
        type = "bar")

fig <- fig %>% layout(title = "Class Balance",
                      xaxis = list(title = "label"),
                      yaxis = list(title = "occurences"))

fig



####################TRANSFORM DATA IN MULTIASSAYEXPERIMENT####################################
brcaList <- list(t(mirna),t(cnv),t(proteins), pheno)
names(brcaList) <- c("mirna","cnv","proteins", "pheno")
brca <- convertToMAE(brcaList)

#pathFile <- fetchPathwayDefinitions("October",2020)
#pathwayList <- readPathways(pathFile)

####################GROUPING VARIABLE TO DEFINE FEATURE####################################
groupList <- list();

for(i in 1:length(experiments(brca))){
  tmp <- list(rownames(experiments(brca)[[i]]));
  names(tmp) <- names(brca)[i]
  groupList[[names(brca)[[i]]]] <- tmp
}

###############DEFINE PATIENT SIMILARITY FOR EACH NETWORK########################
sims <- list(
  mirna="sim.eucscale",
  cnv="sim.eucscale",
  proteins="sim.eucscale"
)



#####################HOLDOUT VALIDATION SET#############################
set.seed(123)
dsets <- subsampleValidationData(brca,pctValidation=0.1)
brca <- dsets$trainMAE
holdout <- dsets$validationMAE

setwd("/Users/a.pennati/Desktop/netDx-multiomics-classification/")
outDir <- "/Users/a.pennati/Desktop/netDx-multiomics-classification/pred"
nco <- round(parallel::detectCores())
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)

##################MODEL TRAINING#################################
t0 <- Sys.time()
model <-buildPredictor(
  dataList=brca,          
  groupList=groupList,    
  sims=sims,
  outDir=outDir,          
  trainProp=0.8,          
  numSplits=1L,            
  featSelCutoff=2L,       
  featScoreMax=10L,
  numCores=nco,            
  debugMode=FALSE,
  keepAllData=FALSE,     
  logging="default" ) 

t1 <- Sys.time()
print(t1-t0)

###########################CONFUSION MATRIX#################################
confMat <- confusionMatrix(model)

######################EXAMINE RESULTS##################################
results <- getResults(
  model,
  unique(colData(brca)$STATUS),
  featureSelCutoff=3L,
  featureSelPct=0.50
)

#################VALIDATE ON INDIPENDENT SAMPLES################################
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

######################PLOT RESULTS OF VALIDATION#################################
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

#####################DATA VISUALIZATION#########################################
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
