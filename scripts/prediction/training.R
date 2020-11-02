
tmp <- readline(prompt="Press any key to run through entire script:") 

library(randomForest)
library(reshape)
library(ggplot2)
library(caret)
library(MLmetrics)
library(philentropy)
library(vcd)
library(PRROC)


# indicate home dir for projet
homeDir = '/Users/dmattox/cbk/lec_gly_binding/'
setwd(homeDir)

##########################
# functions
###########################

# hush=function(code){ # Based on https://stackoverflow.com/questions/2723034/suppress-output-of-a-function
#   if ( .Platform$OS.type == 'unix' ) {
#     sink('/dev/null')
#   } else {
#     sink('NUL')
#   }
#   tmp = code
#   sink()
#   return(tmp)
# }


# sampleClustMembers <- function(clusterIDs, featureSet, ligandTag, testClust, n = 3){
#   # Sample n (default=3) random binding sites from each cluster besides the cluster withheld for testing (LOO, cluster index specified by testClust)
#   # Returns a training dataset with the features for each cluster representative, along with the last column indicating whether the ligand/ligand class of interest is contained in that binding site
#   
#   # Drop the excluded cluster
#   uniClusts = unique(clusterIDs)
#   uniClusts = uniClusts[ ! uniClusts %in% testClust]
#   
#   dat = as.data.frame(matrix(0, nrow = (length(uniClusts) * n), ncol = ncol(featureSet)))
#   colnames(dat) = colnames(predFeats)
#   dat$bound = F
#   
#   for (i in 1:length(uniClusts)) {
#     j = (i-1) * n + 1 # index for dat df to write subsampled cluster members to
#     inds =  (1:nrow(featureSet))[clusterIDs == uniClusts[i]] # all the row indices matching a specific cluster
#     inds = sample(x = inds, size = n, replace = T)
#     
#     dat[(j:(j+(n-1))),1:ncol(featureSet)] = featureSet[inds,]
#     dat$bound[(j:(j+(n-1)))] = ligandTag[inds]
#   }
#   return(dat)
# }
# 
# meanClustMembers <- function(clusterIDs, featureSet, ligandTag, testClust){
#   # Take the mean of each feature foor binding sites with the ligand & the binding sites w/o the ligand from each cluster besides the cluster withheld for testing (LOO, cluster index specified by testClust)
#   # Returns a training dataset with the features for each cluster representative, along with the last column indicating whether the ligand/ligand class of interest is contained in that binding site
#   
#   # Drop the excluded cluster
#   uniClusts = unique(clusterIDs)
#   uniClusts = uniClusts[ ! uniClusts %in% testClust]
#   
#   dat = as.data.frame(matrix(0, nrow = (length(uniClusts) * 2), ncol = ncol(featureSet)))
#   colnames(dat) = colnames(predFeats)
#   dat$bound = F
#   
#   dropTag = rep(F, nrow(dat))
#   
#   for (i in 1:length(uniClusts)) {
#     j = (i - 1) * 2 + 1 # index for dat df to write subsampled cluster members to
#     inds =  (1:nrow(featureSet))[clusterIDs == uniClusts[i]] # all the row indices matching a specific cluster
#     
#     if (any(ligandTag[inds])) {
#       dat[j, 1:ncol(featureSet)] = apply(X = featureSet[inds[ligandTag[inds]], ], MARGIN = 2, FUN = mean)
#     } else{
#       dropTag[j] = T
#     }
#     if (any(!ligandTag[inds])) {
#       dat[(j + 1), 1:ncol(featureSet)] = apply(X = featureSet[inds[!ligandTag[inds]], ], MARGIN = 2, FUN = mean)
#     } else{
#       dropTag[(j + 1)] = T
#     }
#     dat$bound[j] = T
#   }
#   dat = dat[!dropTag,]
#   return(dat)
# }
# 
# balancedClustMembers <- function(clusterIDs, featureSet, ligandTag, testClust, n = 2){
#   # Sample n (default=3) binding sites from the postive & negative cases (ligand present/absent) for each cluster besides the cluster withheld for testing (LOO, cluster index specified by testClust), giving a maximum of n*2 binding sites/cluster
#   # Returns a training dataset with the features for each cluster representative, along with the last column indicating whether the ligand/ligand class of interest is contained in that binding site
#   
#   # Drop the excluded cluster
#   uniClusts = unique(clusterIDs)
#   uniClusts = uniClusts[ ! uniClusts %in% testClust]
#   
#   dat = as.data.frame(matrix(0, nrow = (length(uniClusts) * n * 2), ncol = ncol(featureSet)))
#   colnames(dat) = colnames(predFeats)
#   dat$bound = F
#   
#   for (i in 1:length(uniClusts)) {
#     j = (i - 1) * n * 2 + 1 # index for dat df to write subsampled cluster members to
#     inds =  (1:nrow(featureSet))[clusterIDs == uniClusts[i]] # all the row indices matching a specific cluster
#     
#     if (any(ligandTag[inds])) {
#       posInds = sample(x = inds[ligandTag[inds]], size = n, replace = T)
#       dat[(j:(j+n-1)), 1:ncol(featureSet)] = featureSet[posInds, ]
#     }
#     if (any(!ligandTag[inds])) {
#       negInds = sample(x = inds[!ligandTag[inds]], size = n, replace = T)
#       dat[((j+n):(j+(n*2)-1)), 1:ncol(featureSet)] = featureSet[negInds, ]
#     }
#     dat$bound[(j:(j+n-1))] = T
#   }
#   dat = dat[dat$numBSresis_bin1 != 0,]
#   return(dat)
# }
# 
# meanAllClusts <- function(clusterIDs, featureSet, ligandTag){
#   # Take the mean of each feature foor binding sites with the ligand & the binding sites w/o the ligand from every cluster
#   # Returns a training dataset with the features for each cluster representative, along with the last column indicating whether the ligand/ligand class of interest is contained in that binding site
#   
#   uniClusts = unique(clusterIDs)
# 
#   dat = as.data.frame(matrix(0, nrow = (length(uniClusts) * 2), ncol = ncol(featureSet)))
#   colnames(dat) = colnames(predFeats)
#   dat$bound = F
#   
#   dropTag = rep(F, nrow(dat))
#   
#   for (i in 1:length(uniClusts)) {
#     j = (i - 1) * 2 + 1 # index for dat df to write subsampled cluster members to
#     inds =  (1:nrow(featureSet))[clusterIDs == uniClusts[i]] # all the row indices matching a specific cluster
#     
#     if (any(ligandTag[inds])) {
#       dat[j, 1:ncol(featureSet)] = apply(X = featureSet[inds[ligandTag[inds]], ], MARGIN = 2, FUN = mean)
#     } else{
#       dropTag[j] = T
#     }
#     if (any(!ligandTag[inds])) {
#       dat[(j + 1), 1:ncol(featureSet)] = apply(X = featureSet[inds[!ligandTag[inds]], ], MARGIN = 2, FUN = mean)
#     } else{
#       dropTag[(j + 1)] = T
#     }
#     dat$bound[j] = T
#   }
#   dat = dat[!dropTag,]
#   return(dat)
# }
# 
# meanPerClust <- function(clusterIDs, featureSet, ligandTag){
#   # Take the mean of each feature foor binding sites with the ligand & the binding sites w/o the ligand from every cluster
#   # Returns a training dataset with the features for each cluster representative, along with the last column indicating whether the ligand/ligand class of interest is contained in that binding site
#   
#   uniClusts = unique(clusterIDs)
#   
#   dat = as.data.frame(matrix(0, nrow = (length(uniClusts) * 2), ncol = ncol(featureSet)))
#   colnames(dat) = colnames(predFeats)
#   dat$bound = F
#   
#   dropTag = rep(F, nrow(dat))
#   
#   for (i in 1:length(uniClusts)) {
#     j = (i - 1) * 2 + 1 # index for dat df to write subsampled cluster members to
#     inds =  (1:nrow(featureSet))[clusterIDs == uniClusts[i]] # all the row indices matching a specific cluster
#     
#     if (any(ligandTag[inds])) {
#       dat[j, 1:ncol(featureSet)] = apply(X = featureSet[inds[ligandTag[inds]], ], MARGIN = 2, FUN = mean)
#     } else{
#       dropTag[j] = T
#     }
#     if (any(!ligandTag[inds])) {
#       dat[(j + 1), 1:ncol(featureSet)] = apply(X = featureSet[inds[!ligandTag[inds]], ], MARGIN = 2, FUN = mean)
#     } else{
#       dropTag[(j + 1)] = T
#     }
#     dat$bound[j] = T
#   }
#   dat = dat[!dropTag,]
#   return(dat)
# }


# 
# pullDiverseCases <- function(scaledData, startInds, thresh){
#   if (length(startInds) == 1){
#     
#     out = startInds
#     
#   } else if (length(startInds) == 2){
#     
#     if(suppressMessages(distance(scaledFeats[startInds,])) >= thresh){
#       out = startInds # Two binding sites are greater than the threshold distance from each other, keep both
#     } else{
#       out = sample(startInds, size = 1) # Two binding sites are within the threshold distance from each other, pick one at random
#     }
#     
#   } else {
#     
#     cnter = 1
#     out = rep(0, length(startInds)) # Hold the set of diverse indices
#     
#     while( length(startInds) >= 2 ){
#       
#       out[cnter] = sample(startInds, size = 1) # Sample an index randomly
#       cnter = cnter + 1 # increase the sample count
#       
#       startInds = startInds[! startInds %in% out] # Drop sampled indices from vector of remaining indices
#       
#       distMat = suppressMessages(distance(x = rbind(scaledData[out[out != 0],], scaledData[startInds,]))) # Find all pairwise distances b/w binding sites in cluster
#       distMat = distMat[1:sum(out != 0),-1*(1:sum(out != 0))] # Get the top n rows of the distance matrix dropping the first n columns, where n = the number of indices already sampled
#       
#       if (is.matrix(distMat)){
#         distMat = apply(X = distMat, MARGIN = 2, FUN = min) # For each of the remaining binding sites, take the minimum pairwise distance to any sampled binding site
#       }
#       
#       dropInds = startInds[distMat < thresh ]
#       startInds = startInds[! startInds %in% dropInds]
#       
#       if(length(startInds) == 1){
#         out[cnter] = startInds
#       }
#       
#     }
# 
#     out = out[out != 0]
#   }
#   return(out)
# }
# 
# sampleDiverseSitesByLig <- function(clusterIDs, testClust, featureSet, ligandTag, distThresh, scaledFeatureSet = featureSet){
#   # Sample diverse binding sites with ligand and without ligand [ligandTag] for each cluster in clusterIDs, except the cluster held out for LO(C)O validation indicated by testClust
#   # Samples binding sites from the specified feature set, calculates Euclidean distance b/w binding sites from scaledFeatureSet
#   # Binding sites sampled randomly if Euc. distance to any previously sampled binding sites is greater than distThresh (median pariwise distance between all binding sites)
#   
#     # Drop the excluded cluster
#     uniClusts = unique(clusterIDs)
#     uniClusts = uniClusts[ ! uniClusts %in% testClust]
# 
#     dat = as.data.frame(matrix(0, nrow = nrow(featureSet), ncol = ncol(featureSet)))
#     colnames(dat) = colnames(predFeats)
#     dat$bound = F
#     dat$clus = 0
#     
#     j = 1 # index for writing to returned dataframe (dat)
# 
#     for (i in 1:length(uniClusts)) {
#       inds =  (1:nrow(featureSet))[clusterIDs == uniClusts[i]] # all the row indices matching a specific cluster 
#       negInds = inds[! ligandTag[inds]] # Indices of binding sites w/o ligand
#       posInds = inds[ligandTag[inds]] # Indices of binding sites w/ ligand
#       
#       if (length(negInds > 0)){
#         negInds = pullDiverseCases(scaledData = scaledFeatureSet, startInds = negInds, thresh = distThresh)
#       }
#       if (length(posInds) > 0){
#         posInds = pullDiverseCases(scaledData = scaledFeatureSet, startInds = posInds, thresh = distThresh)
#       }
#       
#       inds = c(negInds, posInds)
#       
#       dat[(j:(j+length(inds) - 1)), (1:ncol(predFeats))] = predFeats[inds, ] # set feature values for representative binding sites
#       dat$bound[(j:(j+length(inds) - 1))] = ligandTag[inds] # Set bound variable
#       dat$clus[(j:(j+length(inds) - 1))] = uniClusts[i] # Set cluster ID
# 
#       j = j + length(inds)
#     }
#     dat = dat[-1*(j:nrow(dat)), ]
#     dat$bound = as.factor(dat$bound)
#     return(dat)
# }


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
  uniClusts = uniClusts[ ! uniClusts %in% testClust]
  
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
# 
# f4 <- function (data, lev = NULL, model = NULL, beta = 4) {
#   precision <- posPredValue(data$pred, data$obs, positive = "TRUE")
#   recall  <- sensitivity(data$pred, data$obs, postive = "TRUE")
#   f4_val <- ((1 + beta^2) * precision * recall) / (beta^2 * precision + recall)
#   names(f4_val) <- c("F4")
#   f4_val
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

##########################
# Set up models
###########################
# Set random seed
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

all(row.names(bsResiDat) == row.names(predFeats))

# Prep RF model specifications
default_mtry = floor(sqrt(ncol(predFeats))) # 209 --> 14


#########################
# Very 1st pass prediction - LO(C)O
#########################

featImp = as.data.frame(matrix(0, nrow = ncol(predFeats), ncol = ncol(ligTags)))
row.names(featImp) = colnames(predFeats)
colnames(featImp) = colnames(ligTags)

set.seed(27)
for(j in 1:ncol(ligTags)){
  
  results = as.data.frame(matrix(0, nrow = length(clusLst), ncol = 4))
  row.names(results) = clusLst
  colnames(results) = c('train_kappa', 'train_recall', 'test_kappa', 'test_recall')
  
  for(i in 1:length(clusLst)){
    # for(i in 1:10){
    
    splitOut = clusLst[i]
    
    # trainDat = sampleClustMembers(clusterIDs = bsResiDat$seqClust50, featureSet = predFeats, ligandTag = ligTags[,j], testClust =splitOut)
    trainDat = meanClustMembers(clusterIDs = bsResiDat$seqClust50, featureSet = predFeats, ligandTag = ligTags[,j], testClust =splitOut)
    # trainDat = balancedClustMembers(clusterIDs = bsResiDat$seqClust50, featureSet = predFeats, ligandTag = ligTags[,j], testClust =splitOut)
    train_binding = as.factor(trainDat$bound)
    trainDat$bound <- NULL
    
    testDat = predFeats[bsResiDat$seqClust50 %in% splitOut,]
    test_binding = factor(ligTags[,j][bsResiDat$seqClust50 %in% splitOut], levels = levels(train_binding))
    
    model = randomForest(x= trainDat, y = train_binding,
                         xtest = testDat, ytest = test_binding,
                         ntree = 1000,
                         mtry = default_mtry,
                         importance = T)
    
    featImp[row.names(model$importance),j] = featImp[row.names(model$importance),j] + model$importance[,4]
    
    results$train_kappa[i] = Kappa(model$confusion[1:2,1:2])$Unweighted[1]
    results$test_kappa[i] = Kappa(model$test$confusion[1:2,1:2])$Unweighted[1]
    
    # results$test_kappa[i] = model$test$votes[1,2]
    
    results$train_recall[i] = 1 - model$confusion[2,3]
    results$test_recall[i] = 1 - model$test$confusion[2,3]
    
    # results$test_TRUE_class.error[i] = model$test$votes[2,2]
  }
  
  featImp[,j] = featImp[,j] / length(clusLst)
  
  colnames(results) = paste(colnames(ligTags)[j], colnames(results), sep = '-')
  
  if(j == 1){
    allResults = results
  } else {
    allResults = cbind(allResults, results)
  }
}


# ex1 = siaResults$test_kappa[1:10]
# ex2 = siaResults$test_TRUE_class.error[1:10]
# boxplot(ex1, ex2)
# points(rep(1,10), ex1)
# points(rep(2,10), ex2)

mResults = melt(allResults, na.rm = T)
mResults$stage = 'train'
mResults$variable = as.character(mResults$variable)

mResults$stage[grepl('-test_', mResults$variable)] = 'validate'

mResults$ligand = gsub('(.*)-t.*', '\\1', mResults$variable)
mResults$variable = gsub('(.*)-(t.*)', '\\2', mResults$variable)

mKappa = mResults[grepl('kappa', mResults$variable), ]
mRecall = mResults[grepl('recall', mResults$variable), ]
mRecall$value = mRecall$value * 100
# 
# # Change from arror to accuracy
# siaResults$train_TRUE_class.error = 1-siaResults$train_TRUE_class.error
# siaResults$test_TRUE_class.error = 1-siaResults$test_TRUE_class.error
# 
# mSia = melt(data = siaResults, na.rm = T)
# mSia$stage = 'train'
# mSia$stage[grepl('^test_', mSia$variable)] = 'validate'
# mSia$variable = as.character(mSia$variable)
# mSia$variable[grep('kappa', mSia$variable)] = 'kappa'
# mSia$variable[grep('error', mSia$variable)] = 'PositiveClass_Accuracy'
# 
# ggplot(data = mSia, aes(x = variable, y = value, col = stage, fill = variable)) + 
#   geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.4) +
#   geom_boxplot(outlier.alpha = 0) + 
#   ylim(c(-0.2, 1)) +
#   scale_fill_manual(values = alpha(c('navajowhite','navajowhite'), 0.6), guide =F) + 
#   scale_color_manual(values = c('tomato','steelblue1')) +
#   labs(title = 'Leave one (Cluster) out - Sialic acid', x = "Performance metric", y = "Value") +
#   theme_dark(base_size = 22)

ggplot(data = mKappa, aes(x = ligand, y = value, col = stage)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.4) +
  geom_boxplot(outlier.alpha = 0) +
  ylim(c(-0.2, 1)) +
  scale_fill_manual(values = alpha(c('navajowhite','navajowhite'), 0.6), guide =F) +
  scale_color_manual(values = c('tomato','steelblue1')) +
  labs(title = 'Leave one (Cluster) out - Kappa values', x = "Ligand types", y = "Kappa") +
  theme_dark(base_size = 22) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))

ggplot(data = mRecall, aes(x = ligand, y = value, col = stage)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.4) +
  geom_boxplot(outlier.alpha = 0) +
  # ylim(c(-0.2, 100)) +
  scale_fill_manual(values = alpha(c('navajowhite','navajowhite'), 0.6), guide =F) +
  scale_color_manual(values = c('tomato','steelblue1')) +
  labs(title = 'Leave one (Cluster) out - Recall values', x = "Ligand types", y = "Recall (%)") +
  theme_dark(base_size = 22) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))


clusSizes = rep(0,length(clusLst))
for (i in 1:length(clusSizes)){
  clusSizes[i] = sum(bsResiDat$seqClust50 == clusLst[i])
}
# 
# tag = ! is.na(siaResults$test_TRUE_class.error)
# plot(clusSizes[tag], siaResults$test_kappa[tag])

ligCnts = rep(0,length(clusLst))
for (i in 1:length(clusLst)) {
  ligCnts[i] = any(ligTags$Sialic_Acid[bsResiDat$seqClust50 == clusLst[i]])
}

ligPercent = rep(0,length(clusLst))
for (i in 1:length(clusSizes)) {
  ligPercent[i] = mean(ligTags$Sialic_Acid[bsResiDat$seqClust50 == clusLst[i]])
}
hist(ligPercent[ligPercent != 0], xlab = 'Fraction of sites with sia in a cluster', main = 'how common is sialic acid when present in a cluster')
sum(ligPercent[ligPercent != 0] < 0.5)

siaPcnt = rep(0,length(clusLst))
for(i in 1:length(clusLst)){
  splitOut = clusLst[i]
  
  trainDat = sampleClustMembers(clusterIDs = bsResiDat$seqClust50, featureSet = predFeats, ligandTag = ligTags$Sialic_Acid, testClust =splitOut)
  
  siaPcnt[i] = mean(trainDat$bound)
}

boxplot(siaPcnt)
abline(h = sum(ligTags$Sialic_Acid)/length(ligTags$Sialic_Acid))

#########################
# 2nd pass prediction - 5x CV from means
#########################

folds = 5
rep = 3
tune.grid <- expand.grid(.mtry= c(-7:7) + default_mtry)

siaDat = meanAllClusts(clusterIDs = bsResiDat$seqClust50, featureSet = predFeats, ligandTag = ligTags$Sialic_Acid) # ~24.4% have binding interactions
siaDat$bound = as.factor(siaDat$bound)


set.seed(27)
inds = createMultiFolds(y = siaDat$bound, k = folds, times = rep)
# for (i in 1:length(inds)){
#   cat(sum(as.logical(siaDat$bound[inds[[i]]]))/length(inds[[i]]))
#   cat('\n')
# }

train.control = trainControl(index = inds,
                  method = 'repeatedcv', 
                  number = folds,
                  repeats = rep,
                  search = 'grid',
                  sampling = 'down'
                  )

rfFit <- train(bound ~ ., data = siaDat, 
               method = "rf", 
               trControl = train.control,
               tuneGrid = tune.grid, 
               maximize = TRUE,
               verbose = TRUE,
               importance = TRUE, 
               ntree = 1000)

rfFit
rfFit$finalModel
varImp(rfFit)

#########################
# 3rd pass prediction - 5x CV training from sampling + LO(C)O validation
#########################

folds = 5
reps = 3
tune.grid <- expand.grid(.mtry= c(-7:7) + default_mtry)

set.seed(27)  

# Loop through ligands with index i here

i=2 #sialic acid

lig = ligTags[,i]

# Define dataframes to hold model results
trainOut = as.data.frame(matrix(0, nrow = length(clusLst), ncol = 7))
row.names(trainOut) = clusLst
colnames(trainOut) = c('final_mtry', 'kappa', 'recall', 'TP', 'TN', 'FP', 'FN')

testOut = as.data.frame(matrix(0, nrow = length(clusLst), ncol = 4))
row.names(testOut) = clusLst
colnames(testOut) = c('TP', 'TN', 'FP', 'FN')

featImp = as.data.frame(matrix(0, nrow = ncol(predFeats), ncol = length(clusLst)))
row.names(featImp) = colnames(predFeats)
colnames(featImp) = clusLst

clusBinding = rep(F, length(clusLst)) # Whether a cluster has any positve examples of binding with the current ligand/ligand class
for (j in (1:length(clusLst))){
  clusBinding[j] = any(lig[bsResiDat$seqClust50 == clusLst[j]])
}
# sum(clusBinding)

testCases = clusLst[clusBinding] # Clusters with any binding occurences to iterativelty withold for validation in LO(C)O validation

predicitions = as.list(rep('', length(row.names(bsResiDat)[bsResiDat$seqClust50 %in% testCases])))
names(predicitions) = row.names(bsResiDat)[bsResiDat$seqClust50 %in% testCases]

for (j in (1:length(testCases))){
  
  outClust = testCases[j]
  cat("testing on clust #", outClust, '\n')
  
  trainingClusts = clusLst[! clusLst == outClust]
  trainingClustBinding = clusBinding[! clusLst == outClust]
  
  foldClusIDs = createMultiFolds(y = trainingClustBinding, k = folds, times = reps)
  
  repDatSampCnt = rep(0, reps)
  for (m in 1:reps){
     sampDat = sampleDiverseSitesByLig(clusterIDs = bsResiDat$seqClust50, 
                                       testClust = outClust,
                                       featureSet = predFeats, 
                                       ligandTag = lig, 
                                       distThresh = medPairwiseDist, 
                                       scaledFeatureSet = scaledFeats)
     
    repDatSampCnt[m] = nrow(sampDat)
    if (m == 1) {
      trainDat = sampDat
    } else{
      trainDat = rbind(trainDat, sampDat)
    }
    
    
    for(n in 1:folds){
      foldInd = ((m - 1) * folds) + n
      if (m == 1){
        foldClusIDs[[foldInd]] = (1:nrow(sampDat))[sampDat$clus %in% trainingClusts[foldClusIDs[[foldInd]]]]
      } else {
        foldClusIDs[[foldInd]] = as.integer(repDatSampCnt[m-1] + (1:nrow(sampDat))[sampDat$clus %in% trainingClusts[foldClusIDs[[foldInd]]]])
      }
      
    }
  }
  
  trainDat$clus <- NULL
  
  train.control = trainControl(index = foldClusIDs,
                               method = 'repeatedcv', 
                               number = folds,
                               repeats = reps,
                               search = 'grid',
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
                 ntree = 1500,
                 metric = "F2")
  
  trainKappa = Kappa(rfFit$finalModel$confusion[1:2,1:2])$Unweighted[1]
  trainRecall = 1 - rfFit$finalModel$confusion[2,3]
  trainAcc = (rfFit$finalModel$confusion[1,1] + rfFit$finalModel$confusion[2,2]) / sum(rfFit$finalModel$confusion)
  best_mtry = unname(rfFit$bestTune)[,1]
  
  trainTag = row.names(trainOut) == outClust
  trainOut$final_mtry[trainTag] = best_mtry
  trainOut$kappa[trainTag] = trainKappa
  trainOut$recall[trainTag] = trainRecall
  trainOut$TP[trainTag] = rfFit$finalModel$confusion[2,2]
  trainOut$TN[trainTag] = rfFit$finalModel$confusion[1,1]
  trainOut$FP[trainTag] = rfFit$finalModel$confusion[1,2]
  trainOut$FN[trainTag] = rfFit$finalModel$confusion[2,1]
  
  featImp[, j] = rfFit$finalModel$importance[,4]
  
  cat("train:\n\tRecall = ", trainRecall, "\n\tKappa = ", trainKappa,"\n\tAccuracy = ", trainAcc, '\n\tmtry = ', best_mtry, '\n')
  
  testDat = predFeats[bsResiDat$seqClust50 == outClust,]
  testObs = factor(lig[bsResiDat$seqClust50 == outClust], levels = levels(trainDat$bound))
  
  validate = predict(rfFit$finalModel, newdata = testDat, type = 'prob')
  
  for(m in 1:nrow(validate)) {
    bsName = row.names(validate)[m]
    predicitions[[bsName]] = validate[bsName,2]
  }
  
  TP = sum(validate[,2] >= 0.5 & testObs == "TRUE")
  TN = sum(validate[,2] < 0.5 & testObs == "FALSE")
  FN = sum(validate[,2] < 0.5 & testObs == "TRUE")
  FP = sum(validate[,2] >= 0.5 & testObs == "FALSE")
  
  randAcc = ((TN+FP)*(TN+FN) + (FN+TP)*(FP+TP)) / nrow(validate)^2
  testAcc = (TP+TN)/(TP+TN+FP+FN)
  testKappa = (testAcc - randAcc) / (1 - randAcc)
  testRecall = TP / (TP + FN)
  
  testTag = row.names(testOut) == outClust
  testOut$TP[testTag] = TP
  testOut$TN[testTag] = TN
  testOut$FN[testTag] = FN
  testOut$FP[testTag] = FP
  
  cat("test:\n\tRecall = ", testRecall, "\n\tKappa = ", testKappa,"\n\tAccuracy = ", testAcc, '\n__________________\n\n')
}

# featImp[, i] = featImp[, i] / length(testCases)

# for(k in 1:25){cat(sum(trainingClustBinding[foldClusIDs[[k]]])); cat('\n')}


#############
# Model trained to maximize Kappa score
#############
# kappaBased_featImp = featImp[,2]
# kappaBased_trainOut = trainOut[trainOut$final_mtry != 0, ]
# kappaBased_testOut = testOut[row.names(testOut) %in% as.character(testCases),]


kappaConfusion = kappaBased_trainOut[,4:7]
apply(kappaConfusion, MARGIN = 2, FUN = mean)
plot(density(apply(kappaConfusion, MARGIN = 1, FUN = sum)), main = 'Kappa trained - number sampled')
kappaConfusion_normalized = kappaConfusion
for(i in 1:nrow(kappaConfusion_normalized)){
  kappaConfusion_normalized[i,] = kappaConfusion[i,] / sum(kappaConfusion[i,])
}
apply(kappaConfusion_normalized, MARGIN = 2, FUN = mean)

kappaBased_trainOut$f2 = 0
kappaBased_trainOut$prec = 0
for(i in 1:nrow(kappaBased_trainOut)){
  r = kappaBased_trainOut$recall[i]
  p = kappaBased_trainOut$TP[i] / (kappaBased_trainOut$TP[i] + kappaBased_trainOut$FP[i])
  kappaBased_trainOut$prec[i] = p
  kappaBased_trainOut$f2[i] = ((1+(2^2)) * p * r) / (2^2 * p + r)
}

# Validation results
kappaTestCon = apply(X = kappaBased_testOut, MARGIN = 2, FUN = sum)

TP = kappaTestCon[[1]]
TN = kappaTestCon[[2]]
FP = kappaTestCon[[3]]
FN = kappaTestCon[[4]]

kappa_validationRecall = TP / ( TP + FN)
kappa_validationPrec = TP / (TP + FP)
kappa_validationF2 = ((1+(2^2)) * kappa_validationPrec * kappa_validationRecall) / (2^2 * kappa_validationPrec + kappa_validationRecall)
randAcc = ((TN+FP)*(TN+FN) + (FN+TP)*(FP+TP)) / sum(c(TP + TN + FP + FN))^2
testAcc = (TP+TN)/(TP+TN+FP+FN)
kappa_validationKappa = (testAcc - randAcc) / (1 - randAcc)

kappaBased_testOut$recall = kappaBased_testOut$TP / (kappaBased_testOut$TP + kappaBased_testOut$FN )
kappaBased_testOut$prec = kappaBased_testOut$TP / (kappaBased_testOut$TP + kappaBased_testOut$FP )
kappaBased_testOut$f2 = kappaBased_testOut$kappa = 0
for ( i in 1:nrow(kappaBased_testOut)){
  kappaBased_testOut$f2[i] = ((1+(2^2)) * kappaBased_testOut$prec[i] * kappaBased_testOut$recall[i]) / (2^2 * kappaBased_testOut$prec[i] + kappaBased_testOut$recall[i])
  
  TP = kappaBased_testOut$TP[i]
  TN = kappaBased_testOut$TN[i]
  FP = kappaBased_testOut$FP[i]
  FN = kappaBased_testOut$FN[i]
  
  randAcc = ((TN+FP)*(TN+FN) + (FN+TP)*(FP+TP)) / sum(c(TP + TN + FP + FN))^2
  testAcc = (TP+TN)/(TP+TN+FP+FN)
  kappaBased_testOut$kappa[i] = (testAcc - randAcc) / (1 - randAcc)
}
kappaBased_testOut$f2[is.na(kappaBased_testOut$f2)] = 0

kappaTestConNormed = kappaBased_testOut[,1:4]
for(i in 1:nrow(kappaTestConNormed)){
  kappaTestConNormed[i,] = kappaTestConNormed[i,] / sum(kappaTestConNormed[i,])
}


#############
# Model trained to maximize F2 score
#############
# f2Based_featImp = featImp[,2]
# f2Based_trainOut = trainOut[trainOut$final_mtry != 0, ]
# f2Based_testOut = testOut[row.names(testOut) %in% as.character(testCases),]


f2Confusion = f2Based_trainOut[,4:7]

apply(X = f2Confusion, MARGIN = 2, FUN = mean)
f2Confusion_norm = f2Confusion
for(i in 1:nrow(f2Confusion_norm)){
  f2Confusion_norm[i,] = f2Confusion[i,] / sum(f2Confusion[i,])
}
apply(f2Confusion_norm, MARGIN = 2, FUN = mean)

# calculate f2 scores
f2Based_trainOut$f2 = 0
f2Based_trainOut$precision = 0
for(i in 1:nrow(f2Based_trainOut)){
  r = f2Based_trainOut$recall[i]
  p = f2Based_trainOut$TP[i] / (f2Based_trainOut$TP[i] + f2Based_trainOut$FP[i])
  f2Based_trainOut$precision[i] = p
  f2Based_trainOut$f2[i] = ((1+(2^2)) * p * r) / (2^2 * p + r)
}

# Validation results
f2TestCon = apply(X = f2Based_testOut, MARGIN = 2, FUN = sum)

TP = f2TestCon[[1]]
TN = f2TestCon[[2]]
FP = f2TestCon[[3]]
FN = f2TestCon[[4]]

f2_validationRecall = TP / ( TP + FN)
f2_validationPrec = TP / (TP + FP)
f2_validationF2 = ((1+(2^2)) * f2_validationPrec * f2_validationRecall) / (2^2 * f2_validationPrec + f2_validationRecall)
randAcc = ((TN+FP)*(TN+FN) + (FN+TP)*(FP+TP)) / sum(c(TP + TN + FP + FN))^2
testAcc = (TP+TN)/(TP+TN+FP+FN)
f2_validationKappa = (testAcc - randAcc) / (1 - randAcc)


f2Based_testOut$recall = f2Based_testOut$TP / (f2Based_testOut$TP + f2Based_testOut$FN )
f2Based_testOut$prec = f2Based_testOut$TP / (f2Based_testOut$TP + f2Based_testOut$FP )
f2Based_testOut$f2 = f2Based_testOut$kappa = 0
for ( i in 1:nrow(f2Based_testOut)){
  f2Based_testOut$f2[i] = ((1+(2^2)) * f2Based_testOut$prec[i] * f2Based_testOut$recall[i]) / (2^2 * f2Based_testOut$prec[i] + f2Based_testOut$recall[i])
  
  TP = f2Based_testOut$TP[i]
  TN = f2Based_testOut$TN[i]
  FP = f2Based_testOut$FP[i]
  FN = f2Based_testOut$FN[i]
  
  randAcc = ((TN+FP)*(TN+FN) + (FN+TP)*(FP+TP)) / sum(c(TP + TN + FP + FN))^2
  testAcc = (TP+TN)/(TP+TN+FP+FN)
  f2Based_testOut$kappa[i] = (testAcc - randAcc) / (1 - randAcc)
}
f2Based_testOut$f2[is.na(f2Based_testOut$f2)] = 0

f2TestConNormed = f2Based_testOut[,1:4]
for(i in 1:nrow(f2TestConNormed)){
  f2TestConNormed[i,] = f2TestConNormed[i,] / sum(f2TestConNormed[i,])
}



################
# Plots
################
# Binding site smapling
plot(density(apply(f2Confusion, MARGIN = 1, FUN = sum)), main = 'Binding sites sampled across validations',
     xlab = 'Number of binding sites sampled',
     xlim = c(440, 500), ylim = c(0,0.06),
     col = 'darkorange1', lwd =3)
par(new = T)
plot(density(apply(kappaConfusion, MARGIN = 1, FUN = sum)),
     main = '', xlab = '', ylab ='', axes = F,
     xlim = c(440, 500), ylim = c(0,0.06),
     col = 'darkorchid1', lwd =3)
legend(x='topright', lty = 1, lwd = 3, col = c('darkorange1', 'darkorchid1'), legend = c('Kappa trained', ' F2 trained'))

# Training scores
par(mfrow=c(1,2))
boxplot(kappaBased_trainOut$kappa, kappaBased_trainOut$recall, kappaBased_trainOut$prec, kappaBased_trainOut$f2,
        ylim = c(0.7, 1),
        names = c('Kappa', 'Recall', 'Prec.', 'F2'),
        main = 'Kappa-based training performance',
        col = 'darkorange1')
boxplot(f2Based_trainOut$kappa, f2Based_trainOut$recall, f2Based_trainOut$prec, f2Based_trainOut$f2,
        ylim = c(0.7, 1), 
        names = c('Kappa', 'Recall', 'Prec.', 'F2'),
        main = 'F2-based training performance',
        col = 'darkorchid1')

# Validation performance by cluster
par(mfrow=c(1,2))
boxplot(kappaBased_testOut$kappa, kappaBased_testOut$recall, kappaBased_testOut$prec, kappaBased_testOut$f2,
        ylim = c(0, 1),
        names = c('Kappa', 'Recall', 'Prec.', 'F2'),
        main = 'Kappa-based validation performance',
        col = 'darkorange1')
boxplot(f2Based_testOut$kappa, f2Based_testOut$recall, f2Based_testOut$prec, f2Based_testOut$f2,
        ylim = c(0, 1), 
        names = c('Kappa', 'Recall', 'Prec.', 'F2'),
        main = 'F2-based validation performance',
        col = 'darkorchid1')

# Density distribution correlation matrices
par(mfrow=c(2,2))
plot(density(kappaTestConNormed$TN),
     main = 'True Neg', 
     xlab = '% of outcomes',
     xlim = c(0,1),
     col = 'darkorange1', lwd =3)
plot(density(kappaTestConNormed$FP),
     main = 'False Pos', 
     xlab = '% of outcomes',
     xlim = c(0,1),
     col = 'darkorange1', lwd =3)
plot(density(kappaTestConNormed$FN),
     main = 'False Neg', 
     xlab = '% of outcomes',
     xlim = c(0,1),
     col = 'darkorange1', lwd =3)
plot(density(kappaTestConNormed$TP),
     main = 'True Pos', 
     xlab = '% of outcomes',
     xlim = c(0,1),
     col = 'darkorange1', lwd =3)

par(mfrow=c(2,2))
plot(density(f2TestConNormed$TN),
     main = 'True Neg', 
     xlab = '% of outcomes',
     xlim = c(0,1),
     col = 'darkorchid1', lwd =3)
plot(density(f2TestConNormed$FP),
     main = 'False Pos', 
     xlab = '% of outcomes',
     xlim = c(0,1),
     col = 'darkorchid1', lwd =3)
plot(density(f2TestConNormed$FN),
     main = 'False Neg', 
     xlab = '% of outcomes',
     xlim = c(0,1),
     col = 'darkorchid1', lwd =3)
plot(density(f2TestConNormed$TP),
     main = 'True Pos', 
     xlab = '% of outcomes',
     xlim = c(0,1),
     col = 'darkorchid1', lwd =3)

# train vs test
par(mfrow=c(1,2))
plot(kappaBased_trainOut$f2, kappaBased_testOut$f2,
     main = 'Kappa-trained',
     xlim = c(0.9,1), ylim = c(0,1),
     xlab = 'Training F2 Score', ylab = 'Test F2 Score',
     col = 'darkorange1', lwd =3, cex = 1.5)
cor(kappaBased_trainOut$f2, kappaBased_testOut$f2)

plot(f2Based_trainOut$f2, f2Based_testOut$f2,
     main = 'F2-trained',
     xlim = c(0.9,1), ylim = c(0,1),
     xlab = 'Training F2 Score', ylab = 'Test F2 Score',
     col = 'darkorchid1', lwd =3, cex = 1.5)
cor(f2Based_trainOut$f2, f2Based_testOut$f2)

# size vs test
par(mfrow=c(1,2))
plot(apply(kappaBased_testOut[,1:4], 1, sum), kappaBased_testOut$f2,
     main = 'Kappa-trained',
     ylim = c(0,1),
     xlab = 'Number of binding sites in excluded (test) cluster', ylab = 'Test F2 Score',
     col = 'darkorange1', lwd =3, cex = 1.5)
plot(apply(f2Based_testOut[,1:4], 1, sum), f2Based_testOut$f2,
     main = 'F2-trained',
     ylim = c(0,1),
     xlab = 'Number of binding sites in excluded (test) cluster', ylab = 'Test F2 Score',
     col = 'darkorchid1', lwd =3, cex = 1.5)

# mtry vs train
par(mfrow=c(1,2))
plot(kappaBased_trainOut$final_mtry, kappaBased_trainOut$f2,
     main = 'Kappa-trained',
     xlim = c(7,21),
     ylim = c(0.85,1),
     xlab = 'mtry used in final model', ylab = 'Train F2 Score',
     col = 'darkorange1', lwd =3, cex = 1.5)
par(new = T)
plot(density(kappaBased_trainOut$final_mtry),
     xlim = c(7,21), axes = F, xlab = '', ylab = '', main = '')
legend(x = 'bottomright', legend = 'mtry density', col = 'black', lty = 1)

plot(f2Based_trainOut$final_mtry, f2Based_trainOut$f2,
     main = 'F2-trained',
     xlim = c(7,21),
     ylim = c(0.85,1),
     xlab = 'mtry used in final model', ylab = 'Train F2 Score',
     col = 'darkorchid1', lwd =3, cex = 1.5)
par(new = T)
plot(density(f2Based_trainOut$final_mtry),
     xlim = c(7,21), axes = F, xlab = '', ylab = '', main = '')
legend(x = 'bottomright', legend = 'mtry density', col = 'black', lty = 1)

# mtry vs test
par(mfrow=c(1,2))
plot(kappaBased_trainOut$final_mtry, kappaBased_testOut$f2,
     main = 'Kappa-trained',
     xlim = c(7,21),
     ylim = c(0,1),
     xlab = 'mtry used in final model', ylab = 'Train F2 Score',
     col = 'darkorange1', lwd =3, cex = 1.5)
plot(f2Based_trainOut$final_mtry, f2Based_testOut$f2,
     main = 'F2-trained',
     xlim = c(7,21),
     ylim = c(0,1),
     xlab = 'mtry used in final model', ylab = 'Train F2 Score',
     col = 'darkorchid1', lwd =3, cex = 1.5)

# mtry vs mtry
dev.off()
plot(jitter(kappaBased_trainOut$final_mtry, 0.1), jitter(f2Based_trainOut$final_mtry, 0.1),
     xlab = 'final mtry from kappa based model', ylab = 'final mtry from F2 based model',
     xlim = c(7,21), ylim = c(7,21))
cor(kappaBased_trainOut$final_mtry, f2Based_trainOut$final_mtry)

# Feature importance
plot(kappaBased_featImp, f2Based_featImp, xlab = 'Kappa trained feature importances', ylab = 'F2 trained feature importances', main = 'Feature importance (Mean Decrease in Gini)')
abline(a = 0, b=1)

plot(density(f2Based_featImp), xlab = 'Mean Decrease in Gini')

par(mar=c(4,10,4,2))
barplot(f2Based_featImp[order(f2Based_featImp, decreasing = T)][20:1], names.arg = row.names(featImp)[order(f2Based_featImp, decreasing = T)][20:1], 
        horiz = T, xlab = 'Mean Decrease in Gini Impurity', las=1)


################
# Get individual labels (AND probabilities)
################
outcomes = as.data.frame(matrix('', nrow = length(predicitions), ncol = 2))
row.names(outcomes) = names(predicitions)
colnames(outcomes) = c("Obs", "Pred")

outcomes$Obs = ligTags[row.names(outcomes), 2]
outcomes$Pred = as.numeric(predicitions[row.names(outcomes)])

# test = runif(n = nrow(outcomes), min = 0, max = 1)

sum(outcomes$Pred >= 0.5 & outcomes$Obs == T)
sum(testOut$TP)

tp = sum(outcomes$Pred >= 0.5 & outcomes$Obs == T)
fp = sum(outcomes$Pred >= 0.5 & outcomes$Obs == F)
fn = sum(outcomes$Pred < 0.5 & outcomes$Obs == T)

P =  tp / (tp + fp)
R = tp / (tp + fn)


pr = pr.curve(outcomes$Pred[outcomes$Obs == T], outcomes$Pred[outcomes$Obs == F], curve= T, rand.compute = T)

plot(pr$curve[,1:2], type = 'l', lwd = 3, col = 'darkorchid1',
     xlim = c(0,1), ylim = c(0,1),
     xlab = 'Recall', ylab = 'Precision', main = paste('PR Curve\nAUC = ', as.character(round(pr$auc.integral, digits = 5)), sep = ''))
abline(h = pr$rand$auc.integral, lty = 1, lwd = 2)
lines(x = c(R,R), y = c(-1,P), lty = 2)
lines(x = c(-1,R), y = c(P,P), lty = 2)
points(R,P, pch = 19)

################
# Investigate False Positives
################

lecIDs = unique(bsResiDat$uniparc)

iupacLigCnt = rep(0, length(lecIDs))

for (i in 1:length(lecIDs)){
  id = lecIDs[i]
  iupacLigCnt[i] = length(unique(bsResiDat$iupac[bsResiDat$uniparc == id]))
}

plot(density(iupacLigCnt), xlab = 'Number of different IUPAC ligands per UniProt ID', main = 'Density distribution of lectin promiscuity')

# Get ligand tags again
uniLigs = unique(bsResiDat$iupac)

parenCnt = bracCnt = manCnt = neuCnt = bracCnt = rep(0,length(uniLigs))
for (i in 1:length(uniLigs)){
  lig = uniLigs[i]
  parenCnt[i] = lengths(regmatches(lig, gregexpr("\\(", lig)))
  bracCnt[i] = lengths(regmatches(lig, gregexpr("\\[", lig)))
  manCnt[i] = lengths(regmatches(lig, gregexpr("Man", lig)))
  neuCnt[i] = lengths(regmatches(lig, gregexpr("NeuAc", lig)))
}

mTag = parenCnt == 0 & bracCnt == 0 # Monosaccharides
dTag = parenCnt == 1 & bracCnt == 0 # Disaccharides
tTag = (parenCnt == 2 & bracCnt == 0) | (parenCnt == 2 & bracCnt == 1) # Trisaccharides
qTag = (parenCnt == 3 & bracCnt == 0) | (parenCnt == 3 & bracCnt == 1) | (parenCnt == 1 & bracCnt == 2) # Tetrasaccharides
pTag = !(mTag | dTag | tTag | qTag) # 5+ sugars
bTag = bracCnt >= 1 # Branched glycans

manTag = manCnt > 3 # High mannose
uniLigs[manTag]
neuTag = neuCnt >= 1 # Has sialic acid
uniLigs[neuTag]
fucTag = grepl('^Fuc',uniLigs) # Has a terminal fucose
uniLigs[fucTag]

# Look at lectin specificity
lecSpec = as.data.frame(matrix(nrow = length(lecIDs), ncol = length(uniLigs)))
row.names(lecSpec) = lecIDs
colnames(lecSpec) = uniLigs
for (i in 1:length(lecIDs)){
  lecSpec[i,] = uniLigs %in% unique((bsResiDat$iupac[bsResiDat$uniparc == lecIDs[i]]))
}

hist(apply(lecSpec, 1, sum), xlim = c(0,15), breaks = 15,
     main = 'Count of unique bound ligands per lectin', xlab = 'Number of unique glycans bound',
     col = 'grey50')
sum(apply(lecSpec, 1, sum) == 1) / nrow(lecSpec)

structCnt = rep(0,length(lecIDs))
for(i in 1:length(lecIDs)){
  structCnt[i] = length(unique(bsResiDat$pdb[bsResiDat$uniparc == lecIDs[i]]))
}

plot(structCnt, apply(lecSpec, 1, sum),
     xlab = 'Number of structures avilable for a lectin', ylab = 'Number of observed unique ligands', main = 'Number of Ligands per Lectin vs. Number of Stuctures per Lectin',
     ylim = c(0,15),
     pch = 19, col = alpha('firebrick1',0.2))
abline(a = 0, b=1)

plot(density( (apply(lecSpec, 1, sum) / structCnt) ))
sum((apply(lecSpec, 1, sum) / structCnt) == 1)

all(row.names(lecSpec)[apply(lecSpec[,neuTag], 1, any)] == unique(bsResiDat$uniparc[bsResiDat$iupac %in% uniLigs[neuTag]])) # confirming ligTags match

siaIDs = unique(bsResiDat$uniparc[bsResiDat$iupac %in% uniLigs[neuTag]])

hist(apply(lecSpec[siaIDs,], 1, sum), breaks = 12,
     xlim = c(0,12),
     main = 'Count of unique bound ligands per sialic-acid-binding lectin', xlab = 'Number of unique glycans bound',
     col = 'grey50')
plot(structCnt[lecIDs %in% siaIDs], apply(lecSpec[siaIDs,], 1, sum),
     xlab = 'Number of structures avilable for a lectin', ylab = 'Number of observed unique ligands', main = 'Number of Ligands per Lectin vs. Number of Stuctures per Lectin\nSialic Acid Binding Lectins',
     ylim = c(0,15),
     pch = 19, col = alpha('firebrick1',0.2))
abline(a = 0, b=1)


sum(apply(lecSpec[siaIDs,!neuTag], 1, sum) != 0)
hist(apply(lecSpec[siaIDs,!neuTag], 1, sum), breaks = 12,
     xlim = c(0,12),
     main = 'Count of unique, bound, non-sialic-acid ligands per sialic-acid-binding lectin', xlab = 'Number of unique, non-sialic-acid glycans bound',
     col = 'grey50')

otherNonNeuLigs = apply(lecSpec[siaIDs,!neuTag], 2, sum)
otherNonNeuLigs = otherNonNeuLigs[otherNonNeuLigs != 0]
otherNonNeuLigs = otherNonNeuLigs[order(otherNonNeuLigs, decreasing = T)]

fpBS = row.names(outcomes)[outcomes$Pred >= 0.5 & outcomes$Obs == F] # Binding sites that ended up as false positives
fpLecs = unique(bsResiDat$uniparc[row.names(bsResiDat) %in% fpBS])
length(fpLecs)

sum(fpLecs %in% siaIDs)

hist(apply(lecSpec[fpLecs[fpLecs %in% siaIDs], !neuTag], 1, sum), breaks = 12,
     xlim = c(0,12),
     main = 'Count of unique, bound, non-sialic-acid ligands per FP sialic-acid-binding lectin', xlab = 'Number of unique, non-sialic-acid glycans bound',
     col = 'grey50')

hist(apply(lecSpec[fpLecs[ ! fpLecs %in% siaIDs], ], 1, sum), breaks = 12,
     xlim = c(0,12),
     main = 'Count of unique, bound ligands per other FP lectin', xlab = 'Number of unique glycans bound',
     col = 'grey50')

#########################
# 4th pass prediction - 5x CV training from sampling + LO(C)O validation
# Revamped sampling to limit similar negative examples
# Train with beta = [2,3,4]
# Lockdown mtry @ 14 (sqrt(feature count))
#########################

folds = 5
reps = 3

default_mtry = round(sqrt(ncol(predFeats)), 0)
default_ntree = 2000

tune.grid = expand.grid(.mtry=default_mtry)
# tune.grid <- expand.grid(.mtry= c(-7:7) + default_mtry)
# tune.grid = expand.grid(.mtry=21)



set.seed(27)  

# Loop through ligands with index i here

i=2 #sialic acid

lig = ligTags[,i]

# Define dataframes to hold model results
trainOut = as.data.frame(matrix(0, nrow = length(clusLst), ncol = 6))
row.names(trainOut) = clusLst
colnames(trainOut) = c('kappa', 'recall', 'TP', 'TN', 'FP', 'FN')

testOut = as.data.frame(matrix(0, nrow = length(clusLst), ncol = 4))
row.names(testOut) = clusLst
colnames(testOut) = c('TP', 'TN', 'FP', 'FN')

clusBinding = rep(F, length(clusLst)) # Whether a cluster has any positve examples of binding with the current ligand/ligand class
for (j in (1:length(clusLst))){
  clusBinding[j] = any(lig[bsResiDat$seqClust50 == clusLst[j]])
}
# sum(clusBinding)

testCases = clusLst[clusBinding] # Clusters with any binding occurences to iterativelty withold for validation in LO(C)O validation

predicitions = as.list(rep('', length(row.names(bsResiDat)[bsResiDat$seqClust50 %in% testCases])))
names(predicitions) = row.names(bsResiDat)[bsResiDat$seqClust50 %in% testCases]

featImp = as.data.frame(matrix(0, nrow = ncol(predFeats), ncol = length(testCases)))
row.names(featImp) = colnames(predFeats)
colnames(featImp) = testCases

for (j in (1:length(testCases))){
  
  outClust = testCases[j]
  cat("testing on clust #", outClust, '\n')
  
  trainingClusts = clusLst[! clusLst == outClust]
  trainingClustBinding = clusBinding[! clusLst == outClust]
  
  foldClusIDs = createMultiFolds(y = trainingClustBinding, k = folds, times = reps)
  
  repDatSampCnt = rep(0, reps)
  for (m in 1:reps){
    sampDat = sampleDiverseSitesByLig(clusterIDs = bsResiDat$seqClust50, 
                                      testClust = outClust,
                                      featureSet = predFeats, 
                                      ligandTag = lig, 
                                      distThresh = medPairwiseDist, 
                                      scaledFeatureSet = scaledFeats)
    
    repDatSampCnt[m] = nrow(sampDat)
    if (m == 1) {
      trainDat = sampDat
    } else{
      trainDat = rbind(trainDat, sampDat)
    }
    
    
    for(n in 1:folds){
      foldInd = ((m - 1) * folds) + n
      if (m == 1){
        foldClusIDs[[foldInd]] = (1:nrow(sampDat))[sampDat$clus %in% trainingClusts[foldClusIDs[[foldInd]]]]
      } else {
        foldClusIDs[[foldInd]] = as.integer(repDatSampCnt[m-1] + (1:nrow(sampDat))[sampDat$clus %in% trainingClusts[foldClusIDs[[foldInd]]]])
      }
      
    }
  }
  
  trainDat$clus <- NULL
  
  train.control = trainControl(index = foldClusIDs,
                               method = 'repeatedcv', 
                               number = folds,
                               repeats = reps,
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
  
  trainKappa = Kappa(rfFit$finalModel$confusion[1:2,1:2])$Unweighted[1]
  trainRecall = 1 - rfFit$finalModel$confusion[2,3]
  trainAcc = (rfFit$finalModel$confusion[1,1] + rfFit$finalModel$confusion[2,2]) / sum(rfFit$finalModel$confusion)

  trainTag = row.names(trainOut) == outClust
  trainOut$kappa[trainTag] = trainKappa
  trainOut$recall[trainTag] = trainRecall
  trainOut$TP[trainTag] = rfFit$finalModel$confusion[2,2]
  trainOut$TN[trainTag] = rfFit$finalModel$confusion[1,1]
  trainOut$FP[trainTag] = rfFit$finalModel$confusion[1,2]
  trainOut$FN[trainTag] = rfFit$finalModel$confusion[2,1]
  
  featImp[, j] = rfFit$finalModel$importance[,4]
  
  cat("train:\n\tRecall = ", trainRecall, "\n\tKappa = ", trainKappa,"\n\tAccuracy = ", trainAcc, '\n')
  
  testDat = predFeats[bsResiDat$seqClust50 == outClust,]
  testObs = factor(lig[bsResiDat$seqClust50 == outClust], levels = levels(trainDat$bound))
  
  validate = predict(rfFit$finalModel, newdata = testDat, type = 'prob')
  
  for(m in 1:nrow(validate)) {
    bsName = row.names(validate)[m]
    predicitions[[bsName]] = validate[bsName,2]
  }
  
  TP = sum(validate[,2] >= 0.5 & testObs == "TRUE")
  TN = sum(validate[,2] < 0.5 & testObs == "FALSE")
  FN = sum(validate[,2] < 0.5 & testObs == "TRUE")
  FP = sum(validate[,2] >= 0.5 & testObs == "FALSE")
  
  randAcc = ((TN+FP)*(TN+FN) + (FN+TP)*(FP+TP)) / nrow(validate)^2
  testAcc = (TP+TN)/(TP+TN+FP+FN)
  testKappa = (testAcc - randAcc) / (1 - randAcc)
  testRecall = TP / (TP + FN)
  
  testTag = row.names(testOut) == outClust
  testOut$TP[testTag] = TP
  testOut$TN[testTag] = TN
  testOut$FN[testTag] = FN
  testOut$FP[testTag] = FP
  
  cat("test:\n\tRecall = ", testRecall, "\n\tKappa = ", testKappa,"\n\tAccuracy = ", testAcc, '\n__________________\n\n')
}

#################
# F2 vs F3 vs F4
#################

# betaCols =  colorRampPalette(c("slateblue","cyan3"))(3)

# newf2_train = trainOut[as.character(testCases),]
# newf2_test = testOut[as.character(testCases),]
# newf2_feats = featImp[,as.character(testCases)]
# newf2_preds = predicitions

apply(newf2_train[,3:6], 2, mean)
apply(apply(newf2_train[,3:6], 1, pCnt), 1, mean)

newf2_train$f2 = newf2_train$f3 = newf2_train$f4 = 0
newf2_train$prec = 0
for(i in 1:nrow(f2Based_trainOut)){
  r = newf2_train$recall[i]
  p = newf2_train$TP[i] / (newf2_train$TP[i] + newf2_train$FP[i])
  newf2_train$prec[i] = p
  newf2_train$f2[i] = ((1+(2^2)) * p * r) / (2^2 * p + r)
  newf2_train$f3[i] = ((1+(3^2)) * p * r) / (3^2 * p + r)
  newf2_train$f4[i] = ((1+(4^2)) * p * r) / (4^2 * p + r)
}

par(mfrow=c(1,2))
boxplot(f2Based_trainOut$kappa, f2Based_trainOut$recall, f2Based_trainOut$prec, f2Based_trainOut$f2,
        ylim = c(0.7, 1), 
        names = c('Kappa', 'Recall', 'Prec.', 'F2'),
        main = 'OLD F2 trained',
        col = 'darkorchid1')
boxplot(newf2_train$kappa, newf2_train$recall, newf2_train$prec, newf2_train$f2,
        ylim = c(0.7, 1), 
        names = c('Kappa', 'Recall', 'Prec.', 'F2'),
        main = 'NEW F2 trained',
        col = 'slateblue')

f2TestMetrics = getKPRFb(newf2_test)

f2_outcomes = as.data.frame(matrix('', nrow = length(newf2_preds), ncol = 2))
row.names(f2_outcomes) = names(newf2_preds)
colnames(f2_outcomes) = c("Obs", "Pred")

f2_outcomes$Obs = ligTags[row.names(f2_outcomes), 2]
f2_outcomes$Pred = as.numeric(predicitions[row.names(f2_outcomes)])

pr = pr.curve(f2_outcomes$Pred[f2_outcomes$Obs == T], f2_outcomes$Pred[f2_outcomes$Obs == F], curve= T, rand.compute = T)

R = f2TestMetrics[['recall']]
p = f2TestMetrics[['precision']]


dev.off()
plot(pr$curve[,1:2], type = 'l', lwd = 3, col = betaCols[1],
     xlim = c(0,1), ylim = c(0,1),
     xlab = 'Recall', ylab = 'Precision', main = paste('PR Curve\nAUC = ', as.character(round(pr$auc.integral, digits = 5)), sep = ''))
abline(h = pr$rand$auc.integral, lty = 1, lwd = 2)
lines(x = c(R,R), y = c(-1,P), lty = 2)
lines(x = c(-1,R), y = c(P,P), lty = 2)
points(R,P, pch = 19)



# newf3_train = trainOut[as.character(testCases),]
# newf3_test = testOut[as.character(testCases),]
# newf3_feats = featImp
# newf3_preds = predicitions

apply(newf3_train[,3:6], 2, mean)
apply(apply(newf3_train[,3:6], 1, pCnt), 1, mean)

newf3_train$f2 = newf3_train$f3 = newf3_train$f4 = 0
newf3_train$prec = 0
for(i in 1:nrow(newf3_train)){
  r = newf3_train$recall[i]
  p = newf3_train$TP[i] / (newf3_train$TP[i] + newf3_train$FP[i])
  newf3_train$prec[i] = p
  newf3_train$f2[i] = ((1+(2^2)) * p * r) / (2^2 * p + r)
  newf3_train$f3[i] = ((1+(3^2)) * p * r) / (3^2 * p + r)
  newf3_train$f4[i] = ((1+(4^2)) * p * r) / (4^2 * p + r)
}

apply(newf3_test, 2, sum)

f3TestMetrics = getKPRFb(newf3_test)

f3_outcomes = as.data.frame(matrix('', nrow = length(newf3_preds), ncol = 2))
row.names(f3_outcomes) = names(newf3_preds)
colnames(f3_outcomes) = c("Obs", "Pred")

f3_outcomes$Obs = ligTags[row.names(f3_outcomes), 2]
f3_outcomes$Pred = as.numeric(predicitions[row.names(f3_outcomes)])

pr = pr.curve(f3_outcomes$Pred[f3_outcomes$Obs == T], f3_outcomes$Pred[f3_outcomes$Obs == F], curve= T, rand.compute = T)

R = f3TestMetrics[['recall']]
p = f3TestMetrics[['precision']]


dev.off()
plot(pr$curve[,1:2], type = 'l', lwd = 3, col = betaCols[2],
     xlim = c(0,1), ylim = c(0,1),
     xlab = 'Recall', ylab = 'Precision', main = paste('PR Curve\nAUC = ', as.character(round(pr$auc.integral, digits = 5)), sep = ''))
abline(h = pr$rand$auc.integral, lty = 1, lwd = 2)
lines(x = c(R,R), y = c(-1,P), lty = 2)
lines(x = c(-1,R), y = c(P,P), lty = 2)
points(R,P, pch = 19)

# newf4_train = trainOut[as.character(testCases),]
# newf4_test = testOut[as.character(testCases),]
# newf4_feats = featImp
# newf4_preds = predicitions

apply(newf4_train[,3:6], 2, mean)
apply(apply(newf4_train[,3:6], 1, pCnt), 1, mean)

newf4_train$f2 = newf4_train$f3 = newf4_train$f4 = 0
newf4_train$prec = 0
for(i in 1:nrow(newf4_train)){
  r = newf4_train$recall[i]
  p = newf4_train$TP[i] / (newf4_train$TP[i] + newf4_train$FP[i])
  newf4_train$prec[i] = p
  newf4_train$f2[i] = ((1+(2^2)) * p * r) / (2^2 * p + r)
  newf4_train$f3[i] = ((1+(3^2)) * p * r) / (3^2 * p + r)
  newf4_train$f4[i] = ((1+(4^2)) * p * r) / (4^2 * p + r)
}

apply(newf4_test, 2, sum)

f4TestMetrics = getKPRFb(newf4_test)


f4_outcomes = as.data.frame(matrix('', nrow = length(newf4_preds), ncol = 2))
row.names(f4_outcomes) = names(newf4_preds)
colnames(f4_outcomes) = c("Obs", "Pred")

f4_outcomes$Obs = ligTags[row.names(f4_outcomes), 2]
f4_outcomes$Pred = as.numeric(predicitions[row.names(f4_outcomes)])

pr = pr.curve(f4_outcomes$Pred[f4_outcomes$Obs == T], f4_outcomes$Pred[f4_outcomes$Obs == F], curve= T, rand.compute = T)

R = f4TestMetrics[['recall']]
p = f4TestMetrics[['precision']]


dev.off()
plot(pr$curve[,1:2], type = 'l', lwd = 3, col = betaCols[3],
     xlim = c(0,1), ylim = c(0,1),
     xlab = 'Recall', ylab = 'Precision', main = paste('PR Curve\nAUC = ', as.character(round(pr$auc.integral, digits = 5)), sep = ''))
abline(h = pr$rand$auc.integral, lty = 1, lwd = 2)
lines(x = c(R,R), y = c(-1,P), lty = 2)
lines(x = c(-1,R), y = c(P,P), lty = 2)
points(R,P, pch = 19)








par(mfrow=c(2,2))
boxplot(newf2_train$kappa, newf2_train$recall, newf2_train$prec, newf2_train$f2, newf2_train$f3, newf2_train$f4,
        ylim = c(0.7, 1), 
        names = c('Kappa', 'Recall', 'Prec.', 'F2', 'F3', 'F4'),
        main = 'F2 trained',
        col = betaCols[1])
boxplot(newf3_train$kappa, newf3_train$recall, newf3_train$prec, newf3_train$f2, newf3_train$f3, newf3_train$f4,
        ylim = c(0.7, 1), 
        names = c('Kappa', 'Recall', 'Prec.', 'F2', 'F3', 'F4'),
        main = 'F3 trained',
        col = betaCols[2])
boxplot(newf4_train$kappa, newf4_train$recall, newf4_train$prec, newf4_train$f2, newf4_train$f3, newf4_train$f4,
        ylim = c(0.7, 1), 
        names = c('Kappa', 'Recall', 'Prec.', 'F2', 'F3', 'F4'),
        main = 'F4 trained',
        col = betaCols[3])

##########
# Vanilla train & test validation
##########

trainOut = trainOut[as.character(testCases),]
testOut = testOut[as.character(testCases),]

# Training
apply(trainOut[,3:6], 2, mean)
apply(apply(trainOut[,3:6], 1, pCnt), 1, mean)

trainOut$f2 = 0
trainOut$prec = 0
for(i in 1:nrow(trainOut)){
  r = trainOut$recall[i]
  p = trainOut$TP[i] / (trainOut$TP[i] + trainOut$FP[i])
  trainOut$prec[i] = p
  trainOut$f2[i] = ((1+(2^2)) * p * r) / (2^2 * p + r)
  trainOut$f3[i] = ((1+(3^2)) * p * r) / (3^2 * p + r)
  trainOut$f4[i] = ((1+(4^2)) * p * r) / (4^2 * p + r)
}

boxplot(trainOut$kappa, trainOut$recall, trainOut$prec, trainOut$f2,
        ylim = c(0.7, 1), 
        names = c('Kappa', 'Recall', 'Prec.', 'F2'),
        main = 'New sampling, mtry = 21, F2 trained',
        col = 'seagreen')

# Validation
apply(testOut, 2, sum)

validationMetrics = getKPRFb(testOut)

outcomes = as.data.frame(matrix('', nrow = length(predicitions), ncol = 2))
row.names(outcomes) = names(predicitions)
colnames(outcomes) = c("Obs", "Pred")

outcomes$Obs = ligTags[row.names(outcomes), 2]
outcomes$Pred = as.numeric(predicitions[row.names(outcomes)])

pr = pr.curve(outcomes$Pred[outcomes$Obs == T], outcomes$Pred[outcomes$Obs == F], curve= T, rand.compute = T)

R = validationMetrics[['recall']]
P = validationMetrics[['precision']]


dev.off()
plot(pr$curve[,1:2], type = 'l', lwd = 3, col = 'slateblue',
     xlim = c(0,1), ylim = c(0,1),
     xlab = 'Recall', ylab = 'Precision', main = paste('PR Curve\nAUC = ', as.character(round(pr$auc.integral, digits = 5)), sep = ''))
abline(h = pr$rand$auc.integral, lty = 1, lwd = 2)
lines(x = c(R,R), y = c(-1,P), lty = 2)
lines(x = c(-1,R), y = c(P,P), lty = 2)
points(R,P, pch = 19)



#########################
# of the negatives, how many come from lectins that DO bind sialic acid in other structures?
#########################

all(row.names(lecSpec)[apply(lecSpec[,neuTag], 1, any)] == unique(bsResiDat$uniparc[bsResiDat$iupac %in% uniLigs[neuTag]])) # confirming ligTags match

siaIDs = unique(bsResiDat$uniparc[bsResiDat$iupac %in% uniLigs[neuTag]])

fpBS = row.names(outcomes)[outcomes$Pred >= 0.5 & outcomes$Obs == F] # Binding sites that ended up as false positives
fpLecs = unique(bsResiDat$uniparc[row.names(bsResiDat) %in% fpBS]) # unique lectin uniprot ids of lectins that show up as false positives
fpIDs = bsResiDat$uniparc[row.names(bsResiDat) %in% fpBS] # UniProt ID for FP binding site
length(fpLecs)

sum(fpLecs %in% siaIDs) # number of FP lectins that do bind sialic acid
sum(fpIDs %in% siaIDs) # Number of False Positives that come from lectins that DO bind sialic acid in other structures


sum(outcomes$Obs == F)
negBS = row.names(outcomes)[outcomes$Obs == F]
negIDs = bsResiDat$uniparc[row.names(bsResiDat) %in% negBS]
sum(negIDs %in% siaIDs) # Number of negative binding sites that come from lectins that DO bind sialic acid in other structures             

bound_fpBS = fpBS[fpIDs %in% siaIDs]
noFP_out = outcomes[! row.names(outcomes) %in% bound_fpBS,]
dropFPs = pr.curve(noFP_out$Pred[noFP_out$Obs == T], noFP_out$Pred[noFP_out$Obs == F], curve= T, rand.compute = T)

bound_NegBS = negBS[negIDs %in% siaIDs]
noFPTN_out = outcomes[! row.names(outcomes) %in% bound_NegBS,]
dropFPTNs = pr.curve(noFPTN_out$Pred[noFPTN_out$Obs == T], noFPTN_out$Pred[noFPTN_out$Obs == F], curve= T, rand.compute = T)


R=0.793
P=0.775
plot(dropFPs$curve[,1:2], type = 'l', lwd = 3, col = 'aquamarine3',
     xlim = c(0,1), ylim = c(0,1),
     xlab = 'Recall', ylab = 'Precision', main = paste('PR Curve\nAUC = ', as.character(round(dropFPs$auc.integral, digits = 5)), sep = ''))
abline(h = dropFPs$rand$auc.integral, lty = 1, lwd = 2)
lines(x = c(R,R), y = c(-1,P), lty = 2)
lines(x = c(-1,R), y = c(P,P), lty = 2)
points(R,P, pch = 19)

plot(dropFPTNs$curve[,1:2], type = 'l', lwd = 3, col = 'aquamarine3',
     xlim = c(0,1), ylim = c(0,1),
     xlab = 'Recall', ylab = 'Precision', main = paste('PR Curve\nAUC = ', as.character(round(dropFPTNs$auc.integral, digits = 5)), sep = ''))
abline(h = dropFPTNs$rand$auc.integral, lty = 1, lwd = 2)
lines(x = c(R,R), y = c(-1,P), lty = 2)
lines(x = c(-1,R), y = c(P,P), lty = 2)
points(R,P, pch = 19)


dev.off()
par(mfrow = c(2,1))
plot(density(outcomes[bound_NegBS, "Pred"]), main = 'Pred Scores for unbound sites from bound lectins',
     xlim = c(0,1), ylim = c(0,2),
     xlab = 'Prediciton score (% of votes)',
     lwd = 3,
     col = 'aquamarine3')

plot(density(outcomes[negBS, "Pred"]), main = 'Pred Scores for all unbound sites',
     xlim = c(0,1), ylim = c(0,2),
     xlab = 'Prediciton score (% of votes)',
     lwd = 3,
     col = 'aquamarine3')


R = validationMetrics[['recall']]
P = validationMetrics[['precision']]

dev.off()
par(mar =c(4,4,4,2)) # bottom, left, top, right
plot(pr$curve[,1:2], type = 'l', lwd = 4, col = 'slateblue',
     xlim = c(0,1), ylim = c(0,1),
     xlab = 'Recall', ylab = 'Precision', main = paste('PR Curve\nAUC = ', as.character(round(pr$auc.integral, digits = 5)), sep = ''))
abline(h = pr$rand$auc.integral, lty = 1, lwd = 2)
lines(x = c(R,R), y = c(-1,P), lty = 2)
lines(x = c(-1,R), y = c(P,P), lty = 2)
points(R,P, pch = 19, col = 'black', cex =1.3)
tag = pr$curve[,3] %in% outcomes[bound_NegBS, "Pred"]
rug((pr$curve[tag,1]), col = alpha('black',0.5))
rug((pr$curve[tag,2]), col = alpha('black',0.5), side = 2)


# 0.8398510 0.5895425 0.4395
nR = 0.8398510
nP = 0.5895425
points(nR,nP, pch = 19, col = 'red', cex =0.5)
lines(x = c(nR,nR), y = c(-1,nP), lty = 2, col ='red')

thresh = 0.4395
TP = sum(outcomes$Pred >= thresh & outcomes$Obs == "TRUE")
TN = sum(outcomes$Pred < thresh & outcomes$Obs == "FALSE")
FN = sum(outcomes$Pred < thresh & outcomes$Obs == "TRUE")
FP = sum(outcomes$Pred >= thresh & outcomes$Obs == "FALSE")


dev.off()
# par(mar =c(4,4,4,2)) # bottom, left, top, right
# plot(pr$curve[,1:2], type = 'l', lwd = 4, col = 'aquamarine3',
#      xlim = c(0,1), ylim = c(0,1),
#      xlab = 'Recall', ylab = 'Precision', main = paste('PR Curve\nAUC = ', as.character(round(pr$auc.integral, digits = 5)), sep = ''))
# abline(h = pr$rand$auc.integral, lty = 1, lwd = 2)
# lines(x = c(R,R), y = c(-1,P), lty = 2)
# lines(x = c(-1,R), y = c(P,P), lty = 2)
# points(R,P, pch = 19, col = 'black', cex =1.3)
# tag = pr$curve[,3] %in% outcomes[bound_fpBS, "Pred"]
# rug((pr$curve[tag,1]), col = alpha('black',0.5))
# rug((pr$curve[tag,2]), col = alpha('black',0.5), side = 2)

#
# bOutcomes = outcomes
# bFeatImp = featImp
#

meanImp = apply(featImp, 1, mean)
sdImp = apply(featImp, 1, sd)
lb = meanImp - sdImp
lb = lb[order(meanImp, decreasing = T)][20:1]
ub = meanImp + sdImp
ub = ub[order(meanImp, decreasing = T)][20:1]

par(mar=c(4,10,4,2))
p = barplot(meanImp[order(meanImp, decreasing = T)][20:1], names.arg = row.names(meanImp)[order(meanImp, decreasing = T)][20:1], xlim = c(0,max(ub)),
        horiz = T, xlab = 'Mean Decrease in Gini Impurity', las=1)

arrows(x0 = lb, y0 = p, x1 = ub, y1 = p, length = 0)

fCorr = cor(predFeats)

corrplot(corr = fCorr, order = 'hclust', addgrid.col = NA)

plot(density(fCorr[upper.tri(fCorr)]))

plot(apply(abs(fCorr), 2, mean), meanImp, xlab = 'mean absolute correlation of feature', ylab = 'mean feat imp')
plot(apply(abs(fCorr), 2, mean), sdImp, xlab = 'mean absolute correlation of feature', ylab = 'feat imp sd')

