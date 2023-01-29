suppressWarnings(suppressMessages(require(netDx)))

## ----eval=TRUE----------------------------------------------------------------
suppressMessages(library(curatedTCGAData))

## ----eval=TRUE----------------------------------------------------------------
brca <- suppressMessages(
  curatedTCGAData("BRCA",
                  c("mRNAArray","miRNA*","Methylation_methyl27*"),
                  dry.run=FALSE,version="1.1.38"))



## ----eval=TRUE----------------------------------------------------------------
# setup brca data
prepareData <- function(dat, setBinary=FALSE) {
  ### clean up stage variable
  staget <- sub("[abcd]","",sub("t","",colData(dat)$pathology_T_stage))
  staget <- suppressWarnings(as.integer(staget))
  colData(dat)$STAGE <- staget
  
  ### remove NA PAM50 calls, remove normal samples
  tmp <- colData(dat)$PAM50.mRNA
  if (!setBinary){
    idx <- which(tmp %in% c("Normal-like","HER2-enriched"))
  } else {
    idx <- union(which(tmp %in% c("Normal-like","HER2-enriched","Luminal B")),
                 which(is.na(staget)))
  }
  idx <- union(idx, which(is.na(tmp)))
  pID <- colData(dat)$patientID
  tokeep <- setdiff(pID, pID[idx])
  dat <- dat[,tokeep,]
  pam50 <- colData(dat)$PAM50.mRNA
  
  ### where a patient has multiple instances of the same assay
  ### just keep the first instance encountered
  smp <- sampleMap(dat)
  expr <- assays(dat)
  for (k in 1:length(expr)) {
    samps <- smp[which(smp$assay==names(expr)[k]),]
    notdup <- samps[which(!duplicated(samps$primary)),"colname"]
    #message(sprintf("%s: %i notdup", names(expr)[k], length(notdup)))
    dat[[k]] <- suppressMessages(dat[[k]][,notdup])
  }
  
  ### create ID, STATUS columns, remove spaces/hyphens from patient labels
  pID <- colData(dat)$patientID
  colData(dat)$ID <- pID
  colData(dat)$STATUS <- pam50
  colData(dat)$STATUS <- gsub(" ",".",colData(dat)$STATUS)
  colData(dat)$STATUS <- gsub("-",".",colData(dat)$STATUS)
  
  if (setBinary){
    st <- colData(dat)$STATUS
    st[which(!st %in% "Luminal.A")] <- "other"
    colData(dat)$STATUS <- st
  }
  
  return(dat)
}
brca <- prepareData(brca,setBinary=TRUE)

## ----eval=TRUE----------------------------------------------------------------
pID <- colData(brca)$patientID
colData(brca)$ID <- pID

## ----eval=TRUE----------------------------------------------------------------
set.seed(123)
dsets <- subsampleValidationData(brca,pctValidation=0.1)
brca <- dsets$trainMAE
holdout <- dsets$validationMAE

## ----eval=TRUE----------------------------------------------------------------
groupList <- list()

# genes in mRNA data are grouped by pathways
pathFile <- sprintf("%s/extdata/pathway_ex3.gmt", path.package("netDx"))
pathList <- suppressMessages(readPathways(pathFile))
groupList[["BRCA_mRNAArray-20160128"]] <- pathList


groupList[["clinical"]] <- list(
  age="patient.age_at_initial_pathologic_diagnosis",
  stage="STAGE"
)

## ----eval=TRUE----------------------------------------------------------------
tmp <- list(rownames(experiments(brca)[[1]]));
names(tmp) <- names(brca)[1]
groupList[[names(brca)[[1]]]] <- tmp

tmp <- list(rownames(experiments(brca)[[2]]));
names(tmp) <- names(brca)[2]
groupList[[names(brca)[2]]] <- tmp

tmp <- list(rownames(experiments(brca)[[3]]));
names(tmp) <- names(brca)[3]
groupList[[names(brca)[3]]] <- tmp

## ----eval=TRUE----------------------------------------------------------------
sims <- list(
  a="pearsonCorr",
  b="normDiff",
  c="pearsonCorr",
  d="pearsonCorr"
)

# map layer names to sims
names(sims) <- names(groupList)

## ----eval=TRUE----------------------------------------------------------------
nco <- round(parallel::detectCores())
message(sprintf("Using %i of %i cores", nco, parallel::detectCores()))

outDir <- paste(tempdir(),"pred_output",sep=getFileSep()) # use absolute path
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
numSplits <- 2L

## ----eval=TRUE----------------------------------------------------------------
t0 <- Sys.time()
model <- 
  buildPredictor(
    dataList=brca,			## your data
    groupList=groupList,	## grouping strategy
    sims=sims,
    outDir=outDir, 			## output directory
    trainProp=0.8,			## pct of samples to use to train model in
    ## each split
    numSplits=2L,			 ## number of train/test splits
    featSelCutoff=1L,		## threshold for calling something
    ## feature-selected
    featScoreMax=2L,	## max score for feature selection
    numCores=nco,			  ## set higher for parallelizing
    debugMode=FALSE,
    keepAllData=FALSE,	    ## set to TRUE for debugging
    logging="default"     ## set to "default" for messages
  )
t1 <- Sys.time()
print(t1-t0)

## ----eval=TRUE----------------------------------------------------------------
results <- getResults(
  model,
  unique(colData(brca)$STATUS),
  featureSelCutoff=2L,
  featureSelPct=0.50
)

confMat <- confusionMatrix(model)

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

## ----eval=TRUE----------------------------------------------------------------
perf <- getPerformance(predModel, 
                       unique(colData(brca)$STATUS))
summary(perf)

## ----eval=TRUE----------------------------------------------------------------
plotPerf_multi(list(perf$rocCurve),
               plotTitle = sprintf(
                 "BRCA Validation: %i samples", 
                 nrow(colData(holdout))))
plotPerf_multi(list(perf$prCurve), 
               plotType = "PR",
               plotTitle = sprintf(
                 "BRCA Validation: %i samples", 
                 nrow(colData(holdout))))


## ---- class.source="codeblock",fig.width=8,fig.height=8, eval=TRUE------------
## this call doesn't work in Rstudio; for now we've commented this out and saved the PSN file. 
psn <- suppressMessages(getPSN(
  brca,
  groupList,
  sims=sims,
  selectedFeatures=results$selectedFeatures
))

## -----------------------------------------------------------------------------
library(Rtsne)
tsne <- tSNEPlotter(
  psn$patientSimNetwork_unpruned, 
  colData(brca)
)

## ----eval=TRUE----------------------------------------------------------------
sessionInfo()

