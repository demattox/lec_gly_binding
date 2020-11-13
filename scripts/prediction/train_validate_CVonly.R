#!/usr/bin/env Rscript

library(randomForest)
library(reshape)
library(ggplot2)
library(caret)
library(MLmetrics)
library(philentropy)
library(vcd)
library(PRROC)

OG_Dir = getwd() # Original directory from which script was called
dirInd = strsplit(OG_Dir, split = '/')[[1]]
dirInd = dirInd[length(dirInd)]

#randInd =  dirInd[length(dirInd)]


homeDir = '/dartfs-hpc/rc/home/y/f002tsy/cbklab/Mattox/lec_gly_binding/'
# homeDir = '/Users/dmattox/cbk/lec_gly_binding/'
setwd(homeDir)


#######################
# Functions
#######################

pullDiverseCases <- function(scaledData, startInds, thresh, prevSampled = NA){
  
  if (any(is.na(prevSampled))) {
    # if no inds being passed in that were sampled from the positive class
    if (length(startInds) == 1) {
      out = startInds
    } else if (length(startInds) == 2) {
      if (suppressMessages(distance(scaledFeats[startInds,])) >= thresh) {
        out = startInds # Two binding sites are greater than the threshold distance from each other, keep both
      } else{
        out = sample(startInds, size = 1) # Two binding sites are within the threshold distance from each other, pick one at random
      }
      
    } else {
      cnter = 1
      out = rep(0, length(startInds)) # Hold the set of diverse indices
    }
  } else{
    # if inds are passed in from the positive class
    
    out = c(prevSampled, rep(0, length(startInds)))
    cnter = length(prevSampled) + 1
    
    if (length(startInds) == 1) {
      # One new neg site to compare to previously smapled positive sites
      
      distMat = suppressMessages(distance(x = rbind(scaledData[out[out != 0], ], scaledData[startInds, ]))) # Find all pairwise distances b/w binding sites in cluster
      if (is.matrix(distMat)){ # if more than one previously sampled binding site
        distMat = distMat[1:sum(out != 0), -1 * (1:sum(out != 0))] # Get the top n rows of the distance matrix dropping the first n columns, where n = the number of indices already sampled
      }
      if (!any(distMat < thresh)) {
        out[cnter] = startInds
      }
    }
  }
  
  
  if (any(out == 0)){
    while( length(startInds) >= 2 ){
      
      out[cnter] = sample(startInds, size = 1) # Sample an index randomly
      cnter = cnter + 1 # increase the sample count
      
      startInds = startInds[! startInds %in% out] # Drop sampled indices from vector of remaining indices
      
      distMat = suppressMessages(distance(x = rbind(scaledData[out[out != 0],], scaledData[startInds,]))) # Find all pairwise distances b/w binding sites in cluster
      distMat = distMat[1:sum(out != 0),-1*(1:sum(out != 0))] # Get the top n rows of the distance matrix dropping the first n columns, where n = the number of indices already sampled
      
      if (is.matrix(distMat)){
        distMat = apply(X = distMat, MARGIN = 2, FUN = min) # For each of the remaining binding sites, take the minimum pairwise distance to any sampled binding site
      }
      
      dropInds = startInds[distMat < thresh ]
      startInds = startInds[! startInds %in% dropInds]
      
      if(length(startInds) == 1){
        out[cnter] = startInds
      }
      
    }
  }
  
  out = out[out != 0]
  return(out)
}

sampleDiverseSitesByLig <- function(clusterIDs, testClust, featureSet, ligandTag, distThresh, scaledFeatureSet = featureSet){
  # Sample diverse binding sites with ligand and without ligand [ligandTag] for each cluster in clusterIDs, except the cluster held out for LO(C)O validation indicated by testClust
  # Samples binding sites from the specified feature set, calculates Euclidean distance b/w binding sites from scaledFeatureSet
  # Binding sites sampled randomly if Euc. distance to any previously sampled binding sites is greater than distThresh (median pariwise distance between all binding sites)
  
  # Drop the excluded cluster
  uniClusts = unique(clusterIDs)
  uniClusts = uniClusts #[ ! uniClusts %in% testClust] 
  
  dat = as.data.frame(matrix(0, nrow = nrow(featureSet), ncol = ncol(featureSet)))
  colnames(dat) = colnames(predFeats)
  dat$bound = F
  dat$clus = 0
  
  j = 1 # index for writing to returned dataframe (dat)
  
  for (i in 1:length(uniClusts)) {
    inds =  (1:nrow(featureSet))[clusterIDs == uniClusts[i]] # all the row indices matching a specific cluster 
    negInds = inds[! ligandTag[inds]] # Indices of binding sites w/o ligand
    posInds = inds[ligandTag[inds]] # Indices of binding sites w/ ligand
    
    if (length(posInds) > 0){
      outInds = pullDiverseCases(scaledData = scaledFeatureSet, startInds = posInds, thresh = distThresh)
    } else {
      outInds = NA
    }
    if (length(negInds > 0)){
      outInds = pullDiverseCases(scaledData = scaledFeatureSet, startInds = negInds, thresh = distThresh, prevSampled = outInds)
    }
    
    dat[(j:(j+length(outInds) - 1)), (1:ncol(predFeats))] = predFeats[outInds, ] # set feature values for representative binding sites
    dat$bound[(j:(j+length(outInds) - 1))] = ligandTag[outInds] # Set bound variable
    dat$clus[(j:(j+length(outInds) - 1))] = uniClusts[i] # Set cluster ID
    
    j = j + length(outInds)
  }
  dat = dat[-1*(j:nrow(dat)), ]
  dat$bound = as.factor(dat$bound)
  return(dat)
}

f2 <- function (data, lev = NULL, model = NULL, beta = 2) {
  precision <- posPredValue(data$pred, data$obs, positive = "TRUE")
  recall  <- sensitivity(data$pred, data$obs, postive = "TRUE")
  f2_val <- ((1 + beta^2) * precision * recall) / (beta^2 * precision + recall)
  names(f2_val) <- c("F2")
  f2_val
}

# f3 <- function (data, lev = NULL, model = NULL, beta = 3) {
#   precision <- posPredValue(data$pred, data$obs, positive = "TRUE")
#   recall  <- sensitivity(data$pred, data$obs, postive = "TRUE")
#   f3_val <- ((1 + beta^2) * precision * recall) / (beta^2 * precision + recall)
#   names(f3_val) <- c("F3")
#   f3_val
# }

pCnt <- function(x){
  return(x/sum(x))
}

getKPRFb <- function(conMatDF){
  # conMatDF has rows of different confusion matrices, with the columns ordered as TP, TN, FP, FN
  # sums each column and finds performance metrics
  f2TestCon = apply(X = conMatDF, MARGIN = 2, FUN = sum)
  
  TP = f2TestCon[[1]]
  TN = f2TestCon[[2]]
  FP = f2TestCon[[3]]
  FN = f2TestCon[[4]]
  
  f2_validationRecall = TP / ( TP + FN)
  f2_validationPrec = TP / (TP + FP)
  f2_validationF2 = ((1+(2^2)) * f2_validationPrec * f2_validationRecall) / (2^2 * f2_validationPrec + f2_validationRecall)
  f3score = ((1+(3^2)) * f2_validationPrec * f2_validationRecall) / (3^2 * f2_validationPrec + f2_validationRecall)
  f4score = ((1+(4^2)) * f2_validationPrec * f2_validationRecall) / (4^2 * f2_validationPrec + f2_validationRecall)
  randAcc = ((TN+FP)*(TN+FN) + (FN+TP)*(FP+TP)) / sum(c(TP + TN + FP + FN))^2
  testAcc = (TP+TN)/(TP+TN+FP+FN)
  f2_validationKappa = (testAcc - randAcc) / (1 - randAcc)
  return(list(kappa = f2_validationKappa,
              recall = f2_validationRecall,
              precision = f2_validationPrec,
              F2 = f2_validationF2,
              F3 = f3score,
              F4 = f4score))
}

addSlash = function(string) {
  # Adds a trailing forward slash to the end of a string (ex path to a driectory) if it is not present
  if (substr(x = string,
             start = nchar(string),
             stop = nchar(string)) != "/") {
    string = paste(string, "/", sep = "")
  }
  return(string)
}

colfunc = colorRampPalette(c("red","goldenrod","forestgreen","royalblue","darkviolet"))

#######################
# Read in data and set up models
#######################

set.seed(27)  

# Read in data
ligTags = read.delim(file = './analysis/training/data_in/ligTags.tsv', sep = '\t', stringsAsFactors = F)
predFeats = read.delim(file = './analysis/training/data_in/predFeats.csv', sep = ',', stringsAsFactors = F)
bsResiDat = read.delim(file = './analysis/training/data_in/bsResiDat.tsv', sep = '\t', stringsAsFactors = F)

scaledFeats = predFeats  # Scale features between 0 & 1
for(i in 1:ncol(scaledFeats)){
  scaledFeats[,i] = (scaledFeats[,i] - min(scaledFeats[,i])) / (max(scaledFeats[,i]) - min(scaledFeats[,i]))
}

bwBSiteDists = distance(scaledFeats)
medPairwiseDist = median(bwBSiteDists[upper.tri(bwBSiteDists)])

clusLst = unique(bsResiDat$seqClust50)

# all(row.names(bsResiDat) == row.names(predFeats))

folds = 10

CVrepeats = 100

default_mtry = round(sqrt(ncol(predFeats)), 0)
default_ntree = 2000

# tune.grid = expand.grid(.mtry=default_mtry)
half_mtry = round(0.5*default_mtry,0)
tune.grid <- expand.grid(.mtry= c(-half_mtry:half_mtry) + default_mtry)

#########################
# Train and validate based on current ligand class 
# CV Validation from sampling
  # Revamped sampling to limit similar negative examples
  # Train with F2
  # mtry within [sqrt(feature count) +/- 50%]
#########################

dirInd = as.integer(dirInd) # Subdirectory from which script was called to indicate which ligTag column (& which ligand) to train classifier for

lig = ligTags[,dirInd]

# Define dataframes to hold model results
trainOut = as.data.frame(matrix(0, nrow = CVrepeats , ncol = 7))
row.names(trainOut) = 1:CVrepeats
colnames(trainOut) = c('mtry', 'kappa', 'recall', 'TP', 'TN', 'FP', 'FN')

clusBinding = rep(F, length(clusLst)) # Whether a cluster has any positve examples of binding with the current ligand/ligand class
for (j in (1:length(clusLst))){
  clusBinding[j] = any(lig[bsResiDat$seqClust50 == clusLst[j]])
}
# sum(clusBinding)

testCases = clusLst[clusBinding] # Clusters with any binding occurences

predictions = as.data.frame(matrix(nrow = length(row.names(bsResiDat)[bsResiDat$seqClust50 %in% testCases]), ncol = CVrepeats))
row.names(predictions) = row.names(bsResiDat)[bsResiDat$seqClust50 %in% testCases]

featImp = as.data.frame(matrix(0, nrow = ncol(predFeats), ncol = CVrepeats))
row.names(featImp) = colnames(predFeats)
colnames(featImp) = 1:CVrepeats

for (j in (1:CVrepeats)){

  foldClusIDs = createFolds(y = clusBinding, k = folds)
  
  trainDat = sampleDiverseSitesByLig(clusterIDs = bsResiDat$seqClust50, 
                                    testClust = outClust,
                                    featureSet = predFeats, 
                                    ligandTag = lig, 
                                    distThresh = medPairwiseDist, 
                                    scaledFeatureSet = scaledFeats)

  for(n in 1:folds){
    foldClusIDs[[n]] = as.integer((1:nrow(trainDat))[trainDat$clus %in% clusLst[foldClusIDs[[n]]]])
  }
  
  trainDat$clus <- NULL
  
  train.control = trainControl(index = foldClusIDs,
                               method = 'cv', 
                               number = folds,
                               sampling = 'down',
                               summaryFunction = f2)
  
  rfFit <- train(bound ~ .,
                 data = trainDat, 
                 method = "rf", 
                 trControl = train.control,
                 tuneGrid = tune.grid,
                 maximize = TRUE,
                 verbose = TRUE,
                 importance = TRUE, 
                 ntree = default_ntree,
                 metric = "F2")
  
  #rfFit$finalModel$votes
  
  trainKappa = Kappa(rfFit$finalModel$confusion[1:2,1:2])$Unweighted[1]
  trainRecall = 1 - rfFit$finalModel$confusion[2,3]
  trainAcc = (rfFit$finalModel$confusion[1,1] + rfFit$finalModel$confusion[2,2]) / sum(rfFit$finalModel$confusion)
  
  trainOut$mtry[j] = unname(rfFit$bestTune)[,1]
  trainOut$kappa[j] = trainKappa
  trainOut$recall[j] = trainRecall
  trainOut$TP[j] = rfFit$finalModel$confusion[2,2]
  trainOut$TN[j] = rfFit$finalModel$confusion[1,1]
  trainOut$FP[j] = rfFit$finalModel$confusion[1,2]
  trainOut$FN[j] = rfFit$finalModel$confusion[2,1]
  
  featImp[, j] = rfFit$finalModel$importance[,4]
  
  cat("train:\n\tRecall = ", trainRecall, "\n\tKappa = ", trainKappa,"\n\tAccuracy = ", trainAcc, '\n\tMtry = ', trainOut$mtry[j], '\n\n')
}

#########################
# output prelim stats and save training/validation performance files
#########################

ligColors = colfunc(ncol(ligTags))

OG_Dir = addSlash(OG_Dir)

cat("Ligand: ", colnames(ligTags)[dirInd], '\n\n')


# Save out files
write.table(x = trainOut, file = paste(OG_Dir,colnames(ligTags)[dirInd], '_training.csv', sep = ''), quote = F, sep = ',')
write.table(x = featImp, file = paste(OG_Dir,colnames(ligTags)[dirInd], '_features.csv', sep = ''), quote = F, sep = ',')


# Training
cat('Training performance\n')

print(apply(trainOut[,4:7], 2, mean))
print(apply(apply(trainOut[,4:7], 1, pCnt), 1, mean))

trainOut$f2 = 0
trainOut$prec = 0
for(i in 1:nrow(trainOut)){
  r = trainOut$recall[i]
  p = trainOut$TP[i] / (trainOut$TP[i] + trainOut$FP[i])
  trainOut$prec[i] = p
  trainOut$f2[i] = ((1+(2^2)) * p * r) / (2^2 * p + r)
  # trainOut$f3[i] = ((1+(3^2)) * p * r) / (3^2 * p + r)
  # trainOut$f4[i] = ((1+(4^2)) * p * r) / (4^2 * p + r)
}

pdf(file = paste(OG_Dir,colnames(ligTags)[dirInd], '_training_metrics.pdf', sep = ''),
    width = 4, # The width of the plot in inches
    height = 6) # The height of the plot in inches
boxplot(trainOut$kappa, trainOut$recall, trainOut$prec, trainOut$f2,
        ylim = c(0, 1), 
        names = c('Kappa', 'Recall', 'Prec.', 'F2'),
        main = paste(colnames(ligTags)[dirInd], 'Training', sep = ''),
        col = ligColors[dirInd])
dev.off()




