
library(randomForest)
library(reshape)
library(ggplot2)
library(caret)
library(MLmetrics)

# library(vcd)


# indicate home dir for projet
homeDir = '/Users/dmattox/cbk/lec_gly_binding/'
setwd(homeDir)

##########################
# functions
###########################

sampleClustMembers <- function(clusterIDs, featureSet, ligandTag, testClust, n = 3){
  # Sample n (default=3) random binding sites from each cluster besides the cluster withheld for testing (LOO, cluster index specified by testClust)
  # Returns a training dataset with the features for each cluster representative, along with the last column indicating whether the ligand/ligand class of interest is contained in that binding site
  
  # Drop the excluded cluster
  uniClusts = unique(clusterIDs)
  uniClusts = uniClusts[ ! uniClusts %in% testClust]
  
  dat = as.data.frame(matrix(0, nrow = (length(uniClusts) * n), ncol = ncol(featureSet)))
  colnames(dat) = colnames(predFeats)
  dat$bound = F
  
  for (i in 1:length(uniClusts)) {
    j = (i-1) * n + 1 # index for dat df to write subsampled cluster members to
    inds =  (1:nrow(featureSet))[clusterIDs == uniClusts[i]] # all the row indices matching a specific cluster
    inds = sample(x = inds, size = n, replace = T)
    
    dat[(j:(j+(n-1))),1:ncol(featureSet)] = featureSet[inds,]
    dat$bound[(j:(j+(n-1)))] = ligandTag[inds]
  }
  return(dat)
}

meanClustMembers <- function(clusterIDs, featureSet, ligandTag, testClust){
  # Take the mean of each feature foor binding sites with the ligand & the binding sites w/o the ligand from each cluster besides the cluster withheld for testing (LOO, cluster index specified by testClust)
  # Returns a training dataset with the features for each cluster representative, along with the last column indicating whether the ligand/ligand class of interest is contained in that binding site
  
  # Drop the excluded cluster
  uniClusts = unique(clusterIDs)
  uniClusts = uniClusts[ ! uniClusts %in% testClust]
  
  dat = as.data.frame(matrix(0, nrow = (length(uniClusts) * 2), ncol = ncol(featureSet)))
  colnames(dat) = colnames(predFeats)
  dat$bound = F
  
  dropTag = rep(F, nrow(dat))
  
  for (i in 1:length(uniClusts)) {
    j = (i - 1) * 2 + 1 # index for dat df to write subsampled cluster members to
    inds =  (1:nrow(featureSet))[clusterIDs == uniClusts[i]] # all the row indices matching a specific cluster
    
    if (any(ligandTag[inds])) {
      dat[j, 1:ncol(featureSet)] = apply(X = featureSet[inds[ligandTag[inds]], ], MARGIN = 2, FUN = mean)
    } else{
      dropTag[j] = T
    }
    if (any(!ligandTag[inds])) {
      dat[(j + 1), 1:ncol(featureSet)] = apply(X = featureSet[inds[!ligandTag[inds]], ], MARGIN = 2, FUN = mean)
    } else{
      dropTag[(j + 1)] = T
    }
    dat$bound[j] = T
  }
  dat = dat[!dropTag,]
  return(dat)
}

balancedClustMembers <- function(clusterIDs, featureSet, ligandTag, testClust, n = 2){
  # Sample n (default=3) binding sites from the postive & negative cases (ligand present/absent) for each cluster besides the cluster withheld for testing (LOO, cluster index specified by testClust), giving a maximum of n*2 binding sites/cluster
  # Returns a training dataset with the features for each cluster representative, along with the last column indicating whether the ligand/ligand class of interest is contained in that binding site
  
  # Drop the excluded cluster
  uniClusts = unique(clusterIDs)
  uniClusts = uniClusts[ ! uniClusts %in% testClust]
  
  dat = as.data.frame(matrix(0, nrow = (length(uniClusts) * n * 2), ncol = ncol(featureSet)))
  colnames(dat) = colnames(predFeats)
  dat$bound = F
  
  for (i in 1:length(uniClusts)) {
    j = (i - 1) * n * 2 + 1 # index for dat df to write subsampled cluster members to
    inds =  (1:nrow(featureSet))[clusterIDs == uniClusts[i]] # all the row indices matching a specific cluster
    
    if (any(ligandTag[inds])) {
      posInds = sample(x = inds[ligandTag[inds]], size = n, replace = T)
      dat[(j:(j+n-1)), 1:ncol(featureSet)] = featureSet[posInds, ]
    }
    if (any(!ligandTag[inds])) {
      negInds = sample(x = inds[!ligandTag[inds]], size = n, replace = T)
      dat[((j+n):(j+(n*2)-1)), 1:ncol(featureSet)] = featureSet[negInds, ]
    }
    dat$bound[(j:(j+n-1))] = T
  }
  dat = dat[dat$numBSresis_bin1 != 0,]
  return(dat)
}

meanAllClusts <- function(clusterIDs, featureSet, ligandTag){
  # Take the mean of each feature foor binding sites with the ligand & the binding sites w/o the ligand from every cluster
  # Returns a training dataset with the features for each cluster representative, along with the last column indicating whether the ligand/ligand class of interest is contained in that binding site
  
  uniClusts = unique(clusterIDs)

  dat = as.data.frame(matrix(0, nrow = (length(uniClusts) * 2), ncol = ncol(featureSet)))
  colnames(dat) = colnames(predFeats)
  dat$bound = F
  
  dropTag = rep(F, nrow(dat))
  
  for (i in 1:length(uniClusts)) {
    j = (i - 1) * 2 + 1 # index for dat df to write subsampled cluster members to
    inds =  (1:nrow(featureSet))[clusterIDs == uniClusts[i]] # all the row indices matching a specific cluster
    
    if (any(ligandTag[inds])) {
      dat[j, 1:ncol(featureSet)] = apply(X = featureSet[inds[ligandTag[inds]], ], MARGIN = 2, FUN = mean)
    } else{
      dropTag[j] = T
    }
    if (any(!ligandTag[inds])) {
      dat[(j + 1), 1:ncol(featureSet)] = apply(X = featureSet[inds[!ligandTag[inds]], ], MARGIN = 2, FUN = mean)
    } else{
      dropTag[(j + 1)] = T
    }
    dat$bound[j] = T
  }
  dat = dat[!dropTag,]
  return(dat)
}

meanPerClust <- function(clusterIDs, featureSet, ligandTag){
  # Take the mean of each feature foor binding sites with the ligand & the binding sites w/o the ligand from every cluster
  # Returns a training dataset with the features for each cluster representative, along with the last column indicating whether the ligand/ligand class of interest is contained in that binding site
  
  uniClusts = unique(clusterIDs)
  
  dat = as.data.frame(matrix(0, nrow = (length(uniClusts) * 2), ncol = ncol(featureSet)))
  colnames(dat) = colnames(predFeats)
  dat$bound = F
  
  dropTag = rep(F, nrow(dat))
  
  for (i in 1:length(uniClusts)) {
    j = (i - 1) * 2 + 1 # index for dat df to write subsampled cluster members to
    inds =  (1:nrow(featureSet))[clusterIDs == uniClusts[i]] # all the row indices matching a specific cluster
    
    if (any(ligandTag[inds])) {
      dat[j, 1:ncol(featureSet)] = apply(X = featureSet[inds[ligandTag[inds]], ], MARGIN = 2, FUN = mean)
    } else{
      dropTag[j] = T
    }
    if (any(!ligandTag[inds])) {
      dat[(j + 1), 1:ncol(featureSet)] = apply(X = featureSet[inds[!ligandTag[inds]], ], MARGIN = 2, FUN = mean)
    } else{
      dropTag[(j + 1)] = T
    }
    dat$bound[j] = T
  }
  dat = dat[!dropTag,]
  return(dat)
}

##########################
# Set up models
###########################
# Set a random seed
set.seed(27)

# Read in data
ligTags = read.delim(file = './analysis/training/data_in/ligTags.tsv', sep = '\t', stringsAsFactors = F)
predFeats = read.delim(file = './analysis/training/data_in/predFeats.csv', sep = ',', stringsAsFactors = F)
bsResiDat = read.delim(file = './analysis/training/data_in/bsResiDat.tsv', sep = '\t', stringsAsFactors = F)

clusLst = unique(bsResiDat$seqClust50)

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
# 2nd pass prediction - 5x CV
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










