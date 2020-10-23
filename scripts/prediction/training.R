
tmp <- readline(prompt="Press any key to run through entire script:") 

library(randomForest)
library(reshape)
library(ggplot2)
library(caret)
library(MLmetrics)
library(philentropy)
library(vcd)


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



pullDiverseCases <- function(scaledData, startInds, thresh){
  if (length(startInds) == 1){
    
    out = startInds
    
  } else if (length(startInds) == 2){
    
    if(suppressMessages(distance(scaledFeats[startInds,])) >= thresh){
      out = startInds # Two binding sites are greater than the threshold distance from each other, keep both
    } else{
      out = sample(startInds, size = 1) # Two binding sites are within the threshold distance from each other, pick one at random
    }
    
  } else {
    
    cnter = 1
    out = rep(0, length(startInds)) # Hold the set of diverse indices
    
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

    out = out[out != 0]
  }
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
      
      if (length(negInds > 0)){
        negInds = pullDiverseCases(scaledData = scaledFeatureSet, startInds = negInds, thresh = distThresh)
      }
      if (length(posInds) > 0){
        posInds = pullDiverseCases(scaledData = scaledFeatureSet, startInds = posInds, thresh = distThresh)
      }
      
      inds = c(negInds, posInds)
      
      dat[(j:(j+length(inds) - 1)), (1:ncol(predFeats))] = predFeats[inds, ] # set feature values for representative binding sites
      dat$bound[(j:(j+length(inds) - 1))] = ligandTag[inds] # Set bound variable
      dat$clus[(j:(j+length(inds) - 1))] = uniClusts[i] # Set cluster ID

      j = j + length(inds)
    }
    dat = dat[-1*(j:nrow(dat)), ]
    dat$bound = as.factor(dat$bound)
    return(dat)
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

# Define dataframes to hold model results
trainOut = as.data.frame(matrix(0, nrow = length(clusLst), ncol = 7))
row.names(trainOut) = clusLst
colnames(trainOut) = c('final_mtry', 'kappa', 'recall', 'TP', 'TN', 'FP', 'FN')

testOut = as.data.frame(matrix(0, nrow = length(clusLst), ncol = 4))
row.names(testOut) = clusLst
colnames(testOut) = c('TP', 'TN', 'FP', 'FN')

featImp = as.data.frame(matrix(0, nrow = ncol(predFeats), ncol = ncol(ligTags)))
row.names(featImp) = colnames(predFeats)
colnames(featImp) = colnames(ligTags)

# Loop through ligands with index i here

i=2 #sialic acid

lig = ligTags[,i]

clusBinding = rep(F, length(clusLst)) # Whether a cluster has any positve examples of binding with the current ligand/ligand class
for (j in (1:length(clusLst))){
  clusBinding[j] = any(lig[bsResiDat$seqClust50 == clusLst[j]])
}
# sum(clusBinding)

testCases = clusLst[clusBinding] # Clusters with any binding occurences to iterativelty withold for validation in LO(C)O validation

for (j in (1:2length(testCases))){
  
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
                               sampling = 'down')
  
  rfFit <- train(bound ~ .,
                 data = trainDat, 
                 method = "rf", 
                 trControl = train.control,
                 tuneGrid = tune.grid, 
                 maximize = TRUE,
                 verbose = TRUE,
                 importance = TRUE, 
                 ntree = 1500,
                 metric = "Kappa")
  
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
  
  featImp[, i] = featImp[, i] + rfFit$finalModel$importance[,4]
  
  cat("train:\n\tRecall = ", trainRecall, "\n\tKappa = ", trainKappa,"\n\tAccuracy = ", trainAcc, '\n\tmtry = ', best_mtry, '\n')
  
  testDat = predFeats[bsResiDat$seqClust50 == outClust,]
  testObs = factor(lig[bsResiDat$seqClust50 == outClust], levels = levels(trainDat$bound))
  
  validate = predict(rfFit$finalModel, newdata = testDat)
  
  TP = sum(validate == "TRUE" & testObs == "TRUE")
  TN = sum(validate == "FALSE" & testObs == "FALSE")
  FN = sum(validate == "FALSE" & testObs == "TRUE")
  FP = sum(validate == "TRUE" & testObs == "FALSE")
  
  randAcc = ((TN+FP)*(TN+FN) + (FN+TP)*(FP+TP)) / length(validate)^2
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

featImp[, i] = featImp[, i] / length(testCases)

# for(k in 1:25){cat(sum(trainingClustBinding[foldClusIDs[[k]]])); cat('\n')}






