
library(randomForest)
library(reshape)
library(ggplot2)
library(vcd)

# library(caret)
# library(readxl)

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
    j = (i-1) * 2 + 1 # index for dat df to write subsampled cluster members to
    inds =  (1:nrow(featureSet))[clusterIDs == uniClusts[i]] # all the row indices matching a specific cluster
    inds = sample(x = inds, size = 3, replace = T)
    
    dat[(j:(j+2)),1:ncol(featureSet)] = featureSet[inds,]
    dat$bound[(j:(j+2))] = ligandTag[inds]
  }
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


###########################
# BS example with sialic acid


siaResults = as.data.frame(matrix(0, nrow = length(clusLst), ncol = 4))
row.names(siaResults) = clusLst
colnames(siaResults) = c('train_kappa', 'train_TRUE_class.error', 'test_kappa', 'test_TRUE_class.error')

set.seed(27)
# for(i in 1:length(clusLst)){
for(i in 1:10){
  splitOut = clusLst[1]
  
  trainDat = sampleClustMembers(clusterIDs = bsResiDat$seqClust50, featureSet = predFeats, ligandTag = ligTags$Sialic_Acid, testClust =splitOut)
  train_binding = as.factor(trainDat$bound)
  trainDat$bound <- NULL
  
  testDat = predFeats[bsResiDat$seqClust50 %in% splitOut,]
  test_binding = factor(ligTags$Sialic_Acid[bsResiDat$seqClust50 %in% splitOut], levels = levels(train_binding))
  
  model = randomForest(x= trainDat, y = train_binding,
                       xtest = testDat, ytest = test_binding,
                       ntree = 1000,
                       mtry = default_mtry,
                       importance = T)
  
  siaResults$train_kappa[i] = Kappa(model$confusion[1:2,1:2])$Unweighted[1]
  # siaResults$test_kappa[i] = Kappa(model$test$confusion[1:2,1:2])$Unweighted[1]
  
  siaResults$test_kappa[i] = model$test$votes[1,2]
  
  siaResults$train_TRUE_class.error[i] = model$confusion[2,3]
  # siaResults$test_TRUE_class.error[i] = model$test$confusion[2,3]
  
  siaResults$test_TRUE_class.error[i] = model$test$votes[2,2]
}

# Change from arror to accuracy
siaResults$train_TRUE_class.error = 1-siaResults$train_TRUE_class.error
siaResults$test_TRUE_class.error = 1-siaResults$test_TRUE_class.error

mSia = melt(data = siaResults, na.rm = T)
mSia$stage = 'train'
mSia$stage[grepl('^test_', mSia$variable)] = 'validate'
mSia$variable = as.character(mSia$variable)
mSia$variable[grep('kappa', mSia$variable)] = 'kappa'
mSia$variable[grep('error', mSia$variable)] = 'PositiveClass_Accuracy'

ggplot(data = mSia, aes(x = variable, y = value, col = stage, fill = variable)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.2) +
  geom_boxplot(outlier.alpha = 0) + 
  scale_fill_manual(values = alpha(c('navajowhite','navajowhite'), 0.75), guide =F) + 
  scale_color_manual(values = c('tomato','steelblue1')) +
  labs(title = 'Leave one (Cluster) out - Sialic acid', x = "Performance metric", y = "Value") +
  theme_dark(base_size = 22)

# clusSizes = rep(0,length(clusLst))
# for (i in 1:length(clusSizes)){
#   clusSizes[i] = sum(bsResiDat$seqClust50 == clusLst[i])
# }
# 
# tag = ! is.na(siaResults$test_TRUE_class.error)
# plot(clusSizes[tag], siaResults$test_kappa[tag])

ligCnts = rep(0,length(clusLst))
for (i in 1:length(clusSizes)) {
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

