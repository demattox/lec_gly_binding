

library(survey)
library(randomForest)
library(caret)
library(vcd)
library(scales)
library(philentropy)
library(vioplot)
library(pheatmap)
library(Cairo)


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
  # f3score = ((1+(3^2)) * f2_validationPrec * f2_validationRecall) / (3^2 * f2_validationPrec + f2_validationRecall)
  # f4score = ((1+(4^2)) * f2_validationPrec * f2_validationRecall) / (4^2 * f2_validationPrec + f2_validationRecall)
  randAcc = ((TN+FP)*(TN+FN) + (FN+TP)*(FP+TP)) / sum(c(TP + TN + FP + FN))^2
  testAcc = (TP+TN)/(TP+TN+FP+FN)
  f2_validationKappa = (testAcc - randAcc) / (1 - randAcc)
  return(list(kappa = f2_validationKappa,
              recall = f2_validationRecall,
              F2 = f2_validationF2,
              precision = f2_validationPrec))
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

perc.rank <- function(x) trunc(rank(x))/length(x)

###########################

homeDir = '/Users/dmattox/cbk/lec_gly_binding/'
setwd(homeDir)

ligTags = read.delim(file = './analysis/training/data_in/ligTags.tsv', sep = '\t', stringsAsFactors = F)
predFeats = read.delim(file = './analysis/training/data_in/predFeats.csv', sep = ',', stringsAsFactors = F)
bsResiDat = read.delim(file = './analysis/training/data_in/bsResiDat.tsv', sep = '\t', stringsAsFactors = F)

# load('./analysis/training/surveyObject.RData')

stats = read.delim(file = './analysis/training/weightedWMW_stats.tsv', sep = '\t', stringsAsFactors = F)

# sum(stats[grepl('_bin1', row.names(stats)),grepl('_adj', colnames(stats))] < 0.01)
# sum(stats[grepl('_bin2', row.names(stats)),grepl('_adj', colnames(stats))] < 0.01)
# sum(stats[grepl('_bin3', row.names(stats)),grepl('_adj', colnames(stats))] < 0.01)
# sum(stats[grepl('_bin4', row.names(stats)),grepl('_adj', colnames(stats))] < 0.01)
# 
# boxplot(unlist(abs(stats[grepl('_bin1', row.names(stats)),grepl('_effectSize', colnames(stats))])),
#   unlist(abs(stats[grepl('_bin2', row.names(stats)),grepl('_effectSize', colnames(stats))])),
#   unlist(abs(stats[grepl('_bin3', row.names(stats)),grepl('_effectSize', colnames(stats))])),
#   unlist(abs(stats[grepl('_bin4', row.names(stats)),grepl('_effectSize', colnames(stats))])),
#   notch = T, col = c('firebrick3', 'darkorange2', 'darkgoldenrod2', 'gold2'))

# Colors to features for plotting
featColors = rep('', ncol(predFeats))
resiFeats = colorRampPalette(c("plum1","tomato", "firebrick4"))(4)
pocketFeats = colorRampPalette(c('paleturquoise3', 'deepskyblue', 'mediumblue'))(4)

featColors[1:11] = 'forestgreen'

featColors[grep('^vol_4Ang$', colnames(predFeats)) : grep('^leftskew_10Ang$', colnames(predFeats))] = pocketFeats[2] # features within the d2Feats range
featColors[grepl('^vol_', colnames(predFeats)) | grepl('^pcntSurf_', colnames(predFeats))] = pocketFeats[1] # General pocket descriptors
featColors[grepl('^binnedD2', colnames(predFeats))] = pocketFeats[3] # PCs from the binned D2 measures
featColors[grepl('^zern', colnames(predFeats))] = pocketFeats[4] # PCs from the 3DZDs

featColors[grepl('^numBSresis', colnames(predFeats))] = resiFeats[1] # number of residues in binding site features
featColors[gsub('_bin\\d{1}', '', colnames(predFeats)) %in% c('H', 'B', 'E', 'G', 'T', 'S', 'X.')] = resiFeats[2] # secondary structure features
featColors[gsub('_bin\\d{1}', '', colnames(predFeats)) %in% c('nonpolar', 'polar', 'posCharge', 'negCharge', 'aromatic')] = resiFeats[3] # amino acid properties
featColors[grepl('^[[:upper:]]{3}_', colnames(predFeats)) | grepl('^CA$', colnames(predFeats))] = resiFeats[4] # amino acid identities

resiFeatTag = featColors %in% resiFeats
pocketFeatTag = featColors %in% pocketFeats

########################
## Lac/LacNAc - Descriptive
########################

###
# Weights
clusWeight = nrow(bsResiDat)/length(unique(bsResiDat$seqClust50)) # The proportion of the total weight allotted to each cluster
cWeights = rep(0, nrow(bsResiDat)) # weights based on cluster membership only (all binding sites with each cluster receives the same weight)
cuWeights = rep(0,nrow(bsResiDat)) # weights based on cluster membership and UniProt ID (all clusters get the same weight, and all unique lectins within each cluster also receive the same weight)

for (i in 1:length(unique(bsResiDat$seqClust50))){
  clus = unique(bsResiDat$seqClust50)[i]
  tag = bsResiDat$seqClust50 == clus
  cWeights[tag] = clusWeight/sum(tag) # Equally divide allotted weight between all binding sites within a given cluster
  
  lectinWeight = clusWeight / length(unique(bsResiDat$uniparc[tag])) # Weight allotted to each lectin
  for(j in 1:length(unique(bsResiDat$uniparc[tag]))){ # Further distribute weights within clusters based on uniprot ids
    uniTag = bsResiDat$uniparc[tag] == unique(bsResiDat$uniparc[tag])[j]
    cuWeights[tag][uniTag] = lectinWeight / sum(uniTag)
  }
}
round(sum(cWeights),5) == nrow(bsResiDat) # Sum of weights is equal to number of binding sites
round(sum(cuWeights),5) == nrow(bsResiDat) # Sum of weights is equal to number of binding sites

wmwFeats = cbind(predFeats,ligTags)

des <- svydesign(ids = ~1, data = wmwFeats, weights = 1/cuWeights)
###


ligNames = c('Lactose', 'N-Acetyllactosamine')
ligColors = c('goldenrod2', 'goldenrod2')

lacTags = ligTags[,c(4,11)]

###
# Ligand stats
all(row.names(lacTags) == row.names(bsResiDat))

# how many clusters are the ligands found in
length(unique(bsResiDat$seqClust50[lacTags[,1]]))
length(unique(bsResiDat$seqClust50[lacTags[,1]])) / length(unique(bsResiDat$seqClust50))
length(unique(bsResiDat$seqClust50[lacTags[,2]]))
length(unique(bsResiDat$seqClust50[lacTags[,2]])) / length(unique(bsResiDat$seqClust50))

# number of interactions
sum(lacTags[,1]) 
sum(lacTags[,2])

# cumulative weight of interactions
sum(cuWeights[lacTags[,1]])
sum(cuWeights[lacTags[,2]])

# shared clusters?
sum(unique(bsResiDat$seqClust50[lacTags[,2]]) %in% unique(bsResiDat$seqClust50[lacTags[,1]]))
###

####
## Limit to shared clusters
####

sharedClusts = unique(bsResiDat$seqClust50[lacTags[,2]])[unique(bsResiDat$seqClust50[lacTags[,2]]) %in% unique(bsResiDat$seqClust50[lacTags[,1]])]

# number of interactions
sum(lacTags[,1] & bsResiDat$seqClust50 %in% sharedClusts) 
sum(lacTags[,2] & bsResiDat$seqClust50 %in% sharedClusts) 

# cumulative weight of interactions
sum(cuWeights[lacTags[,1] & bsResiDat$seqClust50 %in% sharedClusts])
sum(cuWeights[lacTags[,2] & bsResiDat$seqClust50 %in% sharedClusts])

lacTags[,1] = lacTags[,1] & bsResiDat$seqClust50 %in% sharedClusts
lacTags[,2] = lacTags[,2] & bsResiDat$seqClust50 %in% sharedClusts


lacStats = as.data.frame(matrix(0, nrow = ncol(predFeats), ncol = (2*4)))
row.names(lacStats) = colnames(predFeats)
colnames(lacStats) = c(paste(colnames(lacTags), 'p', sep = '_'), paste(colnames(lacTags), '_FC', sep = '_'), paste(colnames(lacTags), 'effectSize', sep = '_'), paste(colnames(lacTags), 'adj', sep = '_'))

des_allLAC = subset(des, subset = (apply(lacTags, 1, any)))


for (i in 1:ncol(lacTags)){
  
  des_w = subset(des, subset = lacTags[,i]) # temporary design object holding all interactions with the ligand of interest
  des_wo = subset(des, subset = lacTags[,-i]) # same as above but OTHER the ligands of interest
  # scaled_des_w = subset(scaled_des, subset = ligTags[,i]) # 
  
  for(k in 1:ncol(predFeats)){ # for each feature k
    ligTest = svyranktest(formula = as.formula(paste(colnames(predFeats)[k], ' ~ ', colnames(lacTags)[i], sep = '')), # Wilcoxon–Mann–Whitney test, sample sizes can be small (~5% of 230 clusters ~= 10), no reason to assume distribution is normal as it likely isn't
                          design = des_allLAC, 
                          test = 'wilcoxon') 
    
    if (ligTest$p.value != 0){
      lacStats[k, grepl('_p$', colnames(lacStats))][i] = ligTest$p.value # Raw p-value
    } else{
      lacStats[k, grepl('_p$', colnames(lacStats))][i] = 5e-324 # If p-value is too small and R rounds to 0, reset as the smallest positive double as referenced by "?.Machine"
    }
    
    lacStats[k, grepl('_effectSize$', colnames(lacStats))][i] = ligTest$estimate # common language effect size (0.5 transformed to 0)
    
    lacStats[k, grepl('_FC$', colnames(lacStats))][i] = svymean(predFeats[lacTags[,i],k], des_w)[1] / svymean(predFeats[lacTags[,-i],k], des_wo)[1] # Fold change in weighted means
    
    # ligSpefic_feat_means[k,i] = svymean(scaledFeats[lacTags[,i],k], scaled_des_w)[1] # Get the weighted mean for each feature in interactions with each glycan of interest
  }
  
  lacStats[,grepl('_adj$', colnames(lacStats))][,i] = p.adjust(lacStats[grepl('_p$', colnames(lacStats))][,i], method = "BH") # Benjamini-Hochberg MHT correction (FDR)
}

all(round(lacStats$Gal.b1.4.Glc_adj,3) == round(lacStats$Gal.b1.4.GlcNAc_adj,3))

sum(lacStats[,grepl('_adj$', colnames(lacStats))] < 1e-16)
# superSigTag = lacStats[,grepl('_adj$', colnames(lacStats))] < 1e-16
# lacStats[,grepl('_adj$', colnames(lacStats))][superSigTag] <- 10**(-1*runif(sum(superSigTag), max = -log10(3e-19), min = -log10(1e-16))) # Sample from a log-uniform distribution


#####
## Volcano plots
#####

par(mfrow=c(1,2))

xLim = c(-0.5,0.5)
yLim = c(0,11)
for(i in 1:ncol(lacTags)){
  
  # yLim = c(0,max(-log10(0.1), -log10(min(lacStats[,grepl('_adj$', colnames(lacStats))][,i]))) + 1)
  
  tag = lacStats[,grepl('_adj$', colnames(lacStats))][,i] < 0.01
  
  # dev.off()
  plot(0,0,axes = F, main = '', xlab = '', ylab = '', pch = NA)
  bg = "seashell2"
  fg = "ivory"
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = bg)
  # abline(v = c(-1,-.5,0,.5,1), lwd = 6, col = fg)
  # abline(v = c(-1.25,-.75,-.25,.25,.75,1.25), lwd = 3, col = fg)
  # abline(h = c(0,1,2,3,4,5,6), lwd = 6, col = fg)
  # abline(h = c(0.5,1.5,2.5,3.5,4.5,5.5,6.5), lwd = 3, col = fg)
  
  abline(v = 0, lty=2, lwd = 4, col = 'white')
  
  par(new=T)
  
  plot(lacStats[,grepl('_effectSize$', colnames(lacStats))][,i], -log10(lacStats[,grepl('_adj$', colnames(lacStats))][,i]), # Plot all points w/ color @ alpha 0.5
       xlab = "Effect size", ylab = "-log10(FDR)", main = '',
       pch = 19, cex = 2, col = alpha(featColors, 0.33),
       cex.axis = 1.5, cex.lab = 1.5,
       xlim = xLim, ylim = yLim)
  
  abline(h= -log10(0.01), lwd = 4, col = 'white')
  abline(h = -log10(1e-16), lwd = 1.5, lty = 2, col = 'white')
  
  title(main = ligNames[i], col.main = ligColors[i], cex.main = 1.8, font.main  = 2)
  # mtext(ligNames[i], col = ligColors[i],
  #       side=3, adj=0, outer = T,
  #       line=1.2, cex=1, font=2)
  
  par(new=T)
  plot(lacStats[,grepl('_effectSize$', colnames(lacStats))][tag,i], -log10(lacStats[,grepl('_adj$', colnames(lacStats))][tag,i]), # Plot stat sig points again with alpha 1
       pch = 19, col = featColors[tag], cex = 2,
       axes = F, xlab = "", ylab = "", main = "",
       xlim = xLim, ylim = yLim)
  par(new=T)
  plot(lacStats[,grepl('_effectSize$', colnames(lacStats))][tag,i], -log10(lacStats[,grepl('_adj$', colnames(lacStats))][tag,i]), # Outline stat sig points in black to highlight them
       col = alpha('black',0.2), cex = 2.05,
       xlab = "", ylab = "", axes = F, main = "",
       xlim = xLim, ylim = yLim)
}

#comapre to global trends

plot(stats$Gal.b1.4.Glc_effectSize, lacStats$Gal.b1.4.Glc_effectSize,
     col = featColors, pch =19,
     main = 'Lactose',
     xlab = 'Effect Size vs background', ylab = 'Effect size vs LacNAc')
abline(v=0, h=0)
par(new=T)
plot(stats$Gal.b1.4.Glc_effectSize[lacStats$Gal.b1.4.Glc_adj <= 0.01], lacStats$Gal.b1.4.Glc_effectSize[lacStats$Gal.b1.4.Glc_adj <= 0.01],
     col = 'black', pch =19,
     cex = 1.7,
     main = '', xlab = '', ylab = '', axes = F)
par(new=T)
plot(stats$Gal.b1.4.Glc_effectSize[lacStats$Gal.b1.4.Glc_adj <= 0.01], lacStats$Gal.b1.4.Glc_effectSize[lacStats$Gal.b1.4.Glc_adj <= 0.01],
     col = featColors[lacStats$Gal.b1.4.GlcNAc_adj <= 0.01], pch =19,
     cex = 1.5,
     main = '', xlab = '', ylab = '', axes = F)
cor.test(stats$Gal.b1.4.Glc_effectSize, lacStats$Gal.b1.4.Glc_effectSize)




plot(stats$Gal.b1.4.GlcNAc_effectSize, lacStats$Gal.b1.4.GlcNAc_effectSize,
     col = alpha(featColors, 0.7), pch =19,
     main = 'N-Acetyllactosamine',
     xlab = 'Effect Size vs background', ylab = 'Effect size vs Lac')
abline(v=0, h=0)
par(new=T)
plot(stats$Gal.b1.4.GlcNAc_effectSize[lacStats$Gal.b1.4.GlcNAc_adj <= 0.01], lacStats$Gal.b1.4.GlcNAc_effectSize[lacStats$Gal.b1.4.GlcNAc_adj <= 0.01],
     col = 'black', pch =19,
     cex = 1.7,
     main = '', xlab = '', ylab = '', axes = F)
par(new=T)
plot(stats$Gal.b1.4.GlcNAc_effectSize[lacStats$Gal.b1.4.GlcNAc_adj <= 0.01], lacStats$Gal.b1.4.GlcNAc_effectSize[lacStats$Gal.b1.4.GlcNAc_adj <= 0.01],
     col = featColors[lacStats$Gal.b1.4.GlcNAc_adj <= 0.01], pch =19,
     cex = 1.5,
     main = '', xlab = '', ylab = '', axes = F)
cor.test(stats$Gal.b1.4.GlcNAc_effectSize, lacStats$Gal.b1.4.GlcNAc_effectSize)


########################
## Lac/LacNAc - Classification
########################

###
# Set up
###
scaledFeats = predFeats  # Scale features between 0 & 1
for(i in 1:ncol(scaledFeats)){
  scaledFeats[,i] = (scaledFeats[,i] - min(scaledFeats[,i])) / (max(scaledFeats[,i]) - min(scaledFeats[,i]))
}

bwBSiteDists = distance(scaledFeats)
medPairwiseDist = median(bwBSiteDists[upper.tri(bwBSiteDists)])

all(row.names(bsResiDat) == row.names(predFeats))
all(row.names(bsResiDat) == row.names(lacTags))

lac = apply(lacTags, 1, any)

predFeats = predFeats[lac,] # limit to only interactions with either Lac or LacNAc
bsResiDat = bsResiDat[lac,]
lacTags = lacTags[lac,]


folds = 5
# reps = 1

testReps = 10
LOCO_reps = 10

default_mtry = round(sqrt(ncol(predFeats)), 0)
default_ntree = 2000

# tune.grid = expand.grid(.mtry=default_mtry)
half_mtry = round(0.5*default_mtry,0)
tune.grid <- expand.grid(.mtry= c(-half_mtry:half_mtry) + default_mtry)


clusLst = unique(bsResiDat$seqClust50) # List of unique clusters with any binding for any of the ligands of interest

set.seed(27)

####################
# Run predictive RF - LacNAc
####################
# rm(trainResults, testResults, feats)


lig = lacTags[,2]


for (z in 1:LOCO_reps){
  # Define dataframes to hold model results
  trainOut = as.data.frame(matrix(0, nrow = length(clusLst), ncol = 7))
  row.names(trainOut) = clusLst
  colnames(trainOut) = c('mtry', 'kappa', 'recall', 'TP', 'TN', 'FP', 'FN')
  
  testOut = as.data.frame(matrix(0, nrow = testReps * length(clusLst), ncol = 4))
  for (j in 1:testReps){
    row.names(testOut)[(1:length(clusLst)) + (length(clusLst) * (j-1))] = paste(clusLst, j, sep = '_')
  }
  colnames(testOut) = c('TP', 'TN', 'FP', 'FN')
  
  clusBinding = rep(F, length(clusLst)) # Whether a cluster has any positive examples of binding with the current glcyan
  for (j in (1:length(clusLst))){
    clusBinding[j] = any(lig[bsResiDat$seqClust50 == clusLst[j]])
  }
  # sum(clusBinding)
  
  testCases = clusLst[clusBinding] # Clusters with any binding occurrences to iteratively withhold for validation in LO(C)O validation
  
  predictions = as.data.frame(matrix(nrow = length(row.names(bsResiDat)[bsResiDat$seqClust50 %in% testCases]), ncol = testReps))
  row.names(predictions) = row.names(bsResiDat)[bsResiDat$seqClust50 %in% testCases]
  
  featImp = as.data.frame(matrix(0, nrow = ncol(predFeats), ncol = length(testCases)))
  row.names(featImp) = colnames(predFeats)
  colnames(featImp) = testCases
  
  
  for (j in (1:length(testCases))){
    
    outClust = testCases[j]
    cat("testing on clust #", outClust, '[',as.character(round(100*j/length(testCases), 2)),'% done ]\n')
    
    trainingClusts = clusLst[! clusLst == outClust]
    trainingClustBinding = clusBinding[! clusLst == outClust]
    
    foldClusIDs = createFolds(y = trainingClustBinding, k = folds)
    
    trainDat = sampleDiverseSitesByLig(clusterIDs = bsResiDat$seqClust50, 
                                       testClust = outClust,
                                       featureSet = predFeats, 
                                       ligandTag = lig, 
                                       distThresh = medPairwiseDist, 
                                       scaledFeatureSet = scaledFeats)
    
    cat(sum(trainDat$bound == 'TRUE'), 'positve cases in training ( out of', nrow(trainDat), ')\n')
    
    for(n in 1:folds){
      foldClusIDs[[n]] = (1:nrow(trainDat))[trainDat$clus %in% trainingClusts[foldClusIDs[[n]]]]
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
    
    trainKappa = Kappa(rfFit$finalModel$confusion[1:2,1:2])$Unweighted[1]
    trainRecall = 1 - rfFit$finalModel$confusion[2,3]
    trainAcc = (rfFit$finalModel$confusion[1,1] + rfFit$finalModel$confusion[2,2]) / sum(rfFit$finalModel$confusion)
    prec = rfFit$finalModel$confusion[1,1] / (rfFit$finalModel$confusion[1,1] + rfFit$finalModel$confusion[2,1])
    
    trainTag = row.names(trainOut) == outClust
    trainOut$mtry[trainTag] = unname(rfFit$bestTune)[,1]
    trainOut$kappa[trainTag] = trainKappa
    trainOut$recall[trainTag] = trainRecall
    trainOut$TP[trainTag] = rfFit$finalModel$confusion[2,2]
    trainOut$TN[trainTag] = rfFit$finalModel$confusion[1,1]
    trainOut$FP[trainTag] = rfFit$finalModel$confusion[1,2]
    trainOut$FN[trainTag] = rfFit$finalModel$confusion[2,1]
    
    featImp[, j] = rfFit$finalModel$importance[,4]
    
    cat("train:\n\tRecall = ", trainRecall, "\n\tPrec = ", prec, "\n\tKappa = ", trainKappa,"\n\tAccuracy = ", trainAcc, '\n\tMtry = ', trainOut$mtry[trainTag], '\n\n')
    
    for(m in 1:testReps){
      inds =  (1:nrow(predFeats))[bsResiDat$seqClust50 == outClust] # all the row indices matching the validation cluster
      negInds = inds[! lig[inds]] # Indices of binding sites w/o ligand
      posInds = inds[lig[inds]] # Indices of binding sites w/ ligand
      
      outInds = pullDiverseCases(scaledData = scaledFeats, startInds = posInds, thresh = medPairwiseDist) # Diverse sample from examples of positive interactions
      if (length(negInds > 0)){
        outInds = pullDiverseCases(scaledData = scaledFeats, startInds = negInds, thresh = medPairwiseDist, prevSampled = outInds) # Diverse sample from examples of negative interactions
      }
      
      testDat = predFeats[outInds,]
      testObs = factor(lig[outInds], levels = levels(trainDat$bound))
      
      validate = predict(rfFit$finalModel, newdata = testDat, type = 'prob')
      
      for(n in 1:nrow(validate)) {
        bsName = row.names(validate)[n]
        predictions[bsName, m] = validate[bsName,2]
      }
      
      TP = sum(validate[,2] >= 0.5 & testObs == "TRUE")
      TN = sum(validate[,2] < 0.5 & testObs == "FALSE")
      FN = sum(validate[,2] < 0.5 & testObs == "TRUE")
      FP = sum(validate[,2] >= 0.5 & testObs == "FALSE")
      
      randAcc = ((TN+FP)*(TN+FN) + (FN+TP)*(FP+TP)) / nrow(validate)^2
      testAcc = (TP+TN)/(TP+TN+FP+FN)
      testKappa = (testAcc - randAcc) / (1 - randAcc)
      testRecall = TP / (TP + FN)
      testPrec = TP / (TP + FP)
      
      testTag = row.names(testOut) == paste(outClust, m, sep = '_')
      testOut$TP[testTag] = TP
      testOut$TN[testTag] = TN
      testOut$FN[testTag] = FN
      testOut$FP[testTag] = FP
      
      cat("test:\n\tRecall = ", testRecall, "\n\tPrecision = ", testPrec, "\n\tKappa = ", testKappa,"\n\tAccuracy = ", testAcc, '\n\n')
    }
    cat('__________________\n\n')
  }
  
  cat('\n\nLOCO round ', z, ' of ', LOCO_reps, ' done!\n\n\n\n')
  
  # Store training performance
  trainOut = trainOut[as.character(testCases),]
  trainOut$prec = trainOut$TP / (trainOut$TP + trainOut$FP)
  trainOut$acc = (trainOut$TP + trainOut$TN) / (trainOut$TP + trainOut$FP + trainOut$TN + trainOut$FN)
  
  trainOut$mode = 'pred'
  # trainOut$ligand = ligName
  if (! exists('trainResults')){
    trainResults = trainOut
  }else{
    trainResults = rbind(trainResults, trainOut)
  }

  
  
  # Store validation results
  testOut = testOut[gsub('_\\d*$', '', row.names(testOut)) %in% as.character(testCases),]
  
  outTmp = as.data.frame(matrix(0,nrow= 10, ncol = 5))
  colnames(outTmp) = c('kappa', 'recall', 'f2', 'prec', 'acc')
  
  for(k in 1:10){
    tag = grepl(paste('_', as.character(k), '$', sep = ''), row.names(testOut))
    repDat = testOut[tag,]
    outTmp[k,1:4] = as.numeric(getKPRFb(repDat)) # Kapppa, recall, f2, precision
    accTmp = apply(repDat, 2, sum)
    outTmp[k,5] = sum(accTmp[1:2]) / sum(accTmp)
  }
  outTmp$mode = 'pred'
  # outTmp$ligand = ligName
  outTmp$batch = z
  if (! exists('testResults')){
    testResults = outTmp
  }else{
    testResults = rbind(testResults, outTmp)
  }
  
  # Store feature importance data
  if (! exists('feats')){
    feats = featImp
  }else{
    feats = cbind(feats, featImp)
  }
}

# save(feats, testResults, trainResults, file = 'analysis/fine_specificity/lactose/train_test_featImp.Rdata')


####################
# Run random RF - LacNAc
####################

for (z in 1:LOCO_reps){
  
  # Shuffle binding labels
  lig = lig[sample(x = 1:nrow(lacTags), size = nrow(lacTags), replace = F)] # Shuffle the labels of the training data
  
  # Define dataframes to hold model results
  trainOut = as.data.frame(matrix(0, nrow = length(clusLst), ncol = 7))
  row.names(trainOut) = clusLst
  colnames(trainOut) = c('mtry', 'kappa', 'recall', 'TP', 'TN', 'FP', 'FN')
  
  testOut = as.data.frame(matrix(0, nrow = testReps * length(clusLst), ncol = 4))
  for (j in 1:testReps){
    row.names(testOut)[(1:length(clusLst)) + (length(clusLst) * (j-1))] = paste(clusLst, j, sep = '_')
  }
  colnames(testOut) = c('TP', 'TN', 'FP', 'FN')
  
  clusBinding = rep(F, length(clusLst)) # Whether a cluster has any positive examples of binding with the current glcyan
  for (j in (1:length(clusLst))){
    clusBinding[j] = any(lig[bsResiDat$seqClust50 == clusLst[j]])
  }
  # sum(clusBinding)
  
  testCases = clusLst[clusBinding] # Clusters with any binding occurrences to iteratively withhold for validation in LO(C)O validation
  
  predictions = as.data.frame(matrix(nrow = length(row.names(bsResiDat)[bsResiDat$seqClust50 %in% testCases]), ncol = testReps))
  row.names(predictions) = row.names(bsResiDat)[bsResiDat$seqClust50 %in% testCases]
  
  featImp = as.data.frame(matrix(0, nrow = ncol(predFeats), ncol = length(testCases)))
  row.names(featImp) = colnames(predFeats)
  colnames(featImp) = testCases
  
  
  for (j in (1:length(testCases))){
    
    outClust = testCases[j]
    cat("testing on clust #", outClust, '[',as.character(round(100*j/length(testCases), 2)),'% done ]\n')
    
    trainingClusts = clusLst[! clusLst == outClust]
    trainingClustBinding = clusBinding[! clusLst == outClust]
    
    foldClusIDs = createFolds(y = trainingClustBinding, k = folds)
    
    trainDat = sampleDiverseSitesByLig(clusterIDs = bsResiDat$seqClust50, 
                                       testClust = outClust,
                                       featureSet = predFeats, 
                                       ligandTag = lig, 
                                       distThresh = medPairwiseDist, 
                                       scaledFeatureSet = scaledFeats)
    
    cat(sum(trainDat$bound == 'TRUE'), 'positve cases in training ( out of', nrow(trainDat), ')\n')
    
    for(n in 1:folds){
      foldClusIDs[[n]] = (1:nrow(trainDat))[trainDat$clus %in% trainingClusts[foldClusIDs[[n]]]]
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
    
    trainKappa = Kappa(rfFit$finalModel$confusion[1:2,1:2])$Unweighted[1]
    trainRecall = 1 - rfFit$finalModel$confusion[2,3]
    trainAcc = (rfFit$finalModel$confusion[1,1] + rfFit$finalModel$confusion[2,2]) / sum(rfFit$finalModel$confusion)
    prec = rfFit$finalModel$confusion[1,1] / (rfFit$finalModel$confusion[1,1] + rfFit$finalModel$confusion[2,1])
    
    trainTag = row.names(trainOut) == outClust
    trainOut$mtry[trainTag] = unname(rfFit$bestTune)[,1]
    trainOut$kappa[trainTag] = trainKappa
    trainOut$recall[trainTag] = trainRecall
    trainOut$TP[trainTag] = rfFit$finalModel$confusion[2,2]
    trainOut$TN[trainTag] = rfFit$finalModel$confusion[1,1]
    trainOut$FP[trainTag] = rfFit$finalModel$confusion[1,2]
    trainOut$FN[trainTag] = rfFit$finalModel$confusion[2,1]
    
    featImp[, j] = rfFit$finalModel$importance[,4]
    
    cat("train:\n\tRecall = ", trainRecall, "\n\tPrec = ", prec, "\n\tKappa = ", trainKappa,"\n\tAccuracy = ", trainAcc, '\n\tMtry = ', trainOut$mtry[trainTag], '\n\n')
    
    for(m in 1:testReps){
      inds =  (1:nrow(predFeats))[bsResiDat$seqClust50 == outClust] # all the row indices matching the validation cluster
      negInds = inds[! lig[inds]] # Indices of binding sites w/o ligand
      posInds = inds[lig[inds]] # Indices of binding sites w/ ligand
      
      outInds = pullDiverseCases(scaledData = scaledFeats, startInds = posInds, thresh = medPairwiseDist) # Diverse sample from examples of positive interactions
      if (length(negInds > 0)){
        outInds = pullDiverseCases(scaledData = scaledFeats, startInds = negInds, thresh = medPairwiseDist, prevSampled = outInds) # Diverse sample from examples of negative interactions
      }
      
      testDat = predFeats[outInds,]
      testObs = factor(lig[outInds], levels = levels(trainDat$bound))
      
      validate = predict(rfFit$finalModel, newdata = testDat, type = 'prob')
      
      for(n in 1:nrow(validate)) {
        bsName = row.names(validate)[n]
        predictions[bsName, m] = validate[bsName,2]
      }
      
      TP = sum(validate[,2] >= 0.5 & testObs == "TRUE")
      TN = sum(validate[,2] < 0.5 & testObs == "FALSE")
      FN = sum(validate[,2] < 0.5 & testObs == "TRUE")
      FP = sum(validate[,2] >= 0.5 & testObs == "FALSE")
      
      randAcc = ((TN+FP)*(TN+FN) + (FN+TP)*(FP+TP)) / nrow(validate)^2
      testAcc = (TP+TN)/(TP+TN+FP+FN)
      testKappa = (testAcc - randAcc) / (1 - randAcc)
      testRecall = TP / (TP + FN)
      testPrec = TP / (TP + FP)
      
      testTag = row.names(testOut) == paste(outClust, m, sep = '_')
      testOut$TP[testTag] = TP
      testOut$TN[testTag] = TN
      testOut$FN[testTag] = FN
      testOut$FP[testTag] = FP
      
      cat("test:\n\tRecall = ", testRecall, "\n\tPrecision = ", testPrec, "\n\tKappa = ", testKappa,"\n\tAccuracy = ", testAcc, '\n\n')
    }
    cat('__________________\n\n')
  }
  
  cat('\n\nLOCO round ', z, ' of ', LOCO_reps, ' done!\n\n\n\n')
  
  # Store training performance
  trainOut = trainOut[as.character(testCases),]
  trainOut$prec = trainOut$TP / (trainOut$TP + trainOut$FP)
  trainOut$acc = (trainOut$TP + trainOut$TN) / (trainOut$TP + trainOut$FP + trainOut$TN + trainOut$FN)
  
  trainOut$mode = 'rand'
  # trainOut$ligand = ligName
  if (! exists('trainResults')){
    trainResults = trainOut
  }else{
    trainResults = rbind(trainResults, trainOut)
  }
  
  
  
  # Store validation results
  testOut = testOut[gsub('_\\d*$', '', row.names(testOut)) %in% as.character(testCases),]
  
  outTmp = as.data.frame(matrix(0,nrow= 10, ncol = 5))
  colnames(outTmp) = c('kappa', 'recall', 'f2', 'prec', 'acc')
  
  for(k in 1:10){
    tag = grepl(paste('_', as.character(k), '$', sep = ''), row.names(testOut))
    repDat = testOut[tag,]
    outTmp[k,1:4] = as.numeric(getKPRFb(repDat)) # Kapppa, recall, f2, precision
    accTmp = apply(repDat, 2, sum)
    outTmp[k,5] = sum(accTmp[1:2]) / sum(accTmp)
  }
  outTmp$mode = 'rand'
  # outTmp$ligand = ligName
  outTmp$batch = z
  if (! exists('testResults')){
    testResults = outTmp
  }else{
    testResults = rbind(testResults, outTmp)
  }

}

save(feats, testResults, trainResults, file = 'analysis/fine_specificity/lactose/train_test_featImp.Rdata')


##########################
# RF results
##########################

boxplot(trainResults$acc[trainResults$mode == 'pred'], trainResults$acc[trainResults$mode == 'rand'], 
        notch = T, ylim = c(0,1),names = c('Pred', 'Rand'),
        main = 'Training accuracy'
        )
boxplot(testResults$acc[testResults$mode == 'pred'], testResults$acc[testResults$mode == 'rand'], 
        notch = T, ylim = c(0,1),names = c('Pred', 'Rand'),
        main = 'Validation accuracy'
)

# Training
## P & R same plot with violin & box plots
xLim = c(0.5,1.5)
yLim = c(0,1)

par(cex.lab=1.5, mar = c(8.6, 5.1, 5.1, 4.1))
plot(0,0,col = 'white', xlab = '',ylab = '',type="n",
     axes=FALSE,ann=FALSE,
     xlim = xLim,
     ylim = c(0,1))
abline(h = c(0.2,0.4,0.6,0.8,1.0), lwd = 1, col = 'grey80')
abline(h = c(0.3,0.5,0.7,0.9), lwd = 1, lty = 2, col = 'grey80')
j=1
i = ((j-1) * 2) + 1
vioplot(trainResults$recall[trainResults$mode == 'pred'], at = i, side = 'left',
        col = ligColors[j],
        xlim = xLim, ylim = yLim,
        add = T,
        plotCentre = "line")
par(new = T)
boxplot(trainResults$recall[trainResults$mode == 'rand'], at = i-0.33, notch = T, outline = F,
        col = alpha(ligColors[j],0.2), border = 'grey50',
        axes = F, xlab = '', ylab = '', xlim = xLim, ylim = yLim)
vioplot(trainResults$prec[trainResults$mode == 'pred'], at = i, side = 'right',
        add = T,
        col = alpha(ligColors[j],0.8),
        xlim = xLim, ylim = yLim,
        plotCentre = "line")
par(new = T)
boxplot(trainResults$prec[trainResults$mode == 'rand'], at = i+0.33, notch = T, outline = F,
        col = alpha(ligColors[j],0.2), border = 'grey50',
        axes = F, xlab = '', ylab = '', xlim = xLim, ylim = yLim)

axis(side=1,at=1, labels = F)
axis(side=2,at=pretty(c(0,1)))
axis(side = 4, at = pretty(c(0,1)))  # Add second axis
mtext("Precision", side = 4, line = 3, cex = 2.5, col = alpha('black',0.85))  # Add second axis label
title(main = "N-Acetyllactosamine Classifier \n Training Precision & Recall vs random", xlab = "", ylab = "Recall", cex.lab = 2.5, cex.main = 2)
# labs(title = , x = "Ligands", y = "Recall") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))
text(x = seq.int(1),                   
     y = par("usr")[3] - 0.05,
     labels = ligNames[2],
     col = ligColors,
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 0,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.5,
     ## Increase label size.
     cex = 2)


# Validation
## P & R same plot with violin & box plots
xLim = c(0.5,1.5)
yLim = c(0,1)

par(cex.lab=1.5, mar = c(8.6, 5.1, 5.1, 4.1))
plot(0,0,col = 'white', xlab = '',ylab = '',type="n",
     axes=FALSE,ann=FALSE,
     xlim = xLim,
     ylim = c(0,1))
abline(h = c(0.2,0.4,0.6,0.8,1.0), lwd = 1, col = 'grey80')
abline(h = c(0.3,0.5,0.7,0.9), lwd = 1, lty = 2, col = 'grey80')
j=1
i = ((j-1) * 2) + 1
vioplot(testResults$recall[trainResults$mode == 'pred'], at = i, side = 'left',
        col = ligColors[j],
        xlim = xLim, ylim = yLim,
        add = T,
        plotCentre = "line")
par(new = T)
boxplot(testResults$recall[testResults$mode == 'rand'], at = i-0.33, notch = T, outline = F,
        col = alpha(ligColors[j],0.2), border = 'grey50',
        axes = F, xlab = '', ylab = '', xlim = xLim, ylim = yLim)
vioplot(testResults$prec[testResults$mode == 'pred'], at = i, side = 'right',
        add = T,
        col = alpha(ligColors[j],0.8),
        xlim = xLim, ylim = yLim,
        plotCentre = "line")
par(new = T)
boxplot(testResults$prec[testResults$mode == 'rand'], at = i+0.33, notch = T, outline = F,
        col = alpha(ligColors[j],0.2), border = 'grey50',
        axes = F, xlab = '', ylab = '', xlim = xLim, ylim = yLim)

axis(side=1,at=1, labels = F)
axis(side=2,at=pretty(c(0,1)))
axis(side = 4, at = pretty(c(0,1)))  # Add second axis
mtext("Precision", side = 4, line = 3, cex = 2.5, col = alpha('black',0.85))  # Add second axis label
title(main = "N-Acetyllactosamine Classifier \n Validation Precision & Recall vs random", xlab = "", ylab = "Recall", cex.lab = 2.5, cex.main = 2)
# labs(title = , x = "Ligands", y = "Recall") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))
text(x = seq.int(1),                   
     y = par("usr")[3] - 0.05,
     labels = ligNames[2],
     col = ligColors,
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 0,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.5,
     ## Increase label size.
     cex = 2)

###
# Feature importance
###
RESIpercents = apply(feats[resiFeatTag,], 2, perc.rank)
POCKpercents = apply(feats[pocketFeatTag,], 2, perc.rank)
PLIPpercents = apply(feats[1:11,], 2, perc.rank)

RESImeds = apply(RESIpercents, 1, median)
POCKmeds = apply(POCKpercents, 1, median)
PLIPmeds = apply(PLIPpercents, 1, median)

percentThresh = 0.75 # >= 75th percentile in feature importance in each feature class

par(mfrow=c(1,3))
plot(RESImeds[order(RESImeds, decreasing = T)], pch = 19, col = featColors[resiFeatTag][order(RESImeds, decreasing = T)], main = 'Residue feats')
abline(h = percentThresh)
plot(POCKmeds[order(POCKmeds, decreasing = T)], pch = 19, col = featColors[pocketFeatTag][order(POCKmeds, decreasing = T)], main = 'Pocket feats')
abline(h = percentThresh)
plot(PLIPmeds[order(PLIPmeds, decreasing = T)], pch = 19, col = featColors[1:11][order(PLIPmeds, decreasing = T)], main = 'PLIP feats')
abline(h = percentThresh)
 
dev.off()



stratAllFeatsImp = c(PLIPmeds, RESImeds, POCKmeds)
all(names(stratAllFeatsImp) == row.names(stats))

topImpfeats = stratAllFeatsImp >= percentThresh # Logical vector, indicates if feature passes threshold for stratified median feature importance percentiles

sigFeats = lacStats[,'Gal.b1.4.GlcNAc_adj'] < 0.01 # Logical vector, indicates if feature passes threshold for significance from WMW test, FDR of 1%

i=2
plot(abs(lacStats[,'Gal.b1.4.GlcNAc_effectSize']), stratAllFeatsImp,
     xlab = '|Effect size|', ylab = 'Stratified feat imp percentile', main = paste(ligNames[i], 'vs', ligNames[1]), col.main = ligColors[i],
     col = alpha(featColors, 0.4),
     xlim = c(0,0.4),
     ylim = c(0,1),
     pch = 19)
par(new=T)
plot(abs(lacStats[topImpfeats & sigFeats,'Gal.b1.4.GlcNAc_effectSize']), stratAllFeatsImp[topImpfeats & sigFeats],
     xlab = '', ylab = '', main = '', axes = F,
     xlim = c(0,0.4),
     ylim = c(0,1),
     col = alpha(featColors[topImpfeats & sigFeats], 1),
     pch = 19,
     cex = 1.5)
text(x = 0.35,
     y=0.05,
     labels = paste('R=', as.character(round(cor(abs(lacStats[,'Gal.b1.4.GlcNAc_effectSize']), stratAllFeatsImp), 3)), sep = ''),
     cex = 1.4)
abline(h = percentThresh, lty = 2)

##########################
# Heatmaps
##########################

# Set-up for 3A (altered from descriptive.R)

stats = stats[,colnames(lacStats)]

colnames(lacStats) = gsub('Gal.b1.4.Glc', 'Lac', colnames(lacStats))

stats = cbind(stats, lacStats)

hImpFeats = rbind(rep(F,length(topImpfeats)),rep(F,length(topImpfeats)),topImpfeats,topImpfeats)
hSigFeats = rbind(rep(F,length(sigFeats)),rep(F,length(sigFeats)),sigFeats,sigFeats)

allFeatCorrs = cor(stats[,grepl('_effectSize$', colnames(stats))], stats[,grepl('_effectSize$', colnames(stats))], method = 'pearson')
row.names(allFeatCorrs)  = gsub('_effectSize$', '', row.names(allFeatCorrs))

cWidth = 10
cHeight = 20 

breakLst = seq(-0.5,0.5,0.01)

annot <- data.frame(Feature_Type = rep("", nrow(stats)))
row.names(annot) = row.names(stats)

annot$Feature_Type[featColors == 'forestgreen'] <- 'PLIP interaction counts'

annot$Feature_Type[grep('^vol_4Ang$', row.names(stats)) : grep('^leftskew_10Ang$', row.names(stats))] <- 'D2 distribution features'
annot$Feature_Type[grepl('^vol_', row.names(stats)) | grepl('^pcntSurf_', row.names(stats))] = 'Pocket descriptors' # General pocket descriptors
annot$Feature_Type[grepl('^binnedD2', row.names(stats))] <- 'D2 Principal Components'
annot$Feature_Type[grepl('^zern', row.names(stats))] <- '3DZD Principal Components'

annot$Feature_Type[grepl('^numBSresis', row.names(stats))] <- 'Residue counts/bin'
annot$Feature_Type[gsub('_bin\\d{1}', '', row.names(stats)) %in% c('H', 'B', 'E', 'G', 'T', 'S', 'X.')] <- 'Residue sec struct.'
annot$Feature_Type[gsub('_bin\\d{1}', '', row.names(stats)) %in% c('nonpolar', 'polar', 'posCharge', 'negCharge', 'aromatic')] <- 'Amino acid property counts'
annot$Feature_Type[grepl('^[[:upper:]]{3}_', row.names(stats)) | grepl('^CA$', row.names(stats))] <- 'Residue identities'

annot$Feature_Type <- factor(annot$Feature_Type, levels = unique(annot$Feature_Type))


####
## Residue features only
####


resiAnnot = data.frame(Feature_Type = rep("", sum(resiFeatTag)))
row.names(resiAnnot) = row.names(stats)[resiFeatTag]

row.names(resiAnnot) = gsub('^H_bin', 'a-helix_bin', row.names(resiAnnot))
row.names(resiAnnot) = gsub('^B_bin', 'b-bridge_bin', row.names(resiAnnot))
row.names(resiAnnot) = gsub('^E_bin', 'b-strand_bin', row.names(resiAnnot))
row.names(resiAnnot) = gsub('^G_bin', '3/10-helix_bin', row.names(resiAnnot))
row.names(resiAnnot) = gsub('^T_bin', 'turn_bin', row.names(resiAnnot))
row.names(resiAnnot) = gsub('^S_bin', 'bend_bin', row.names(resiAnnot))
row.names(resiAnnot) = gsub('^X._bin', 'loop_bin', row.names(resiAnnot))

row.names(resiAnnot) = gsub('^CA$', 'Ca2+', row.names(resiAnnot))


resiAnnot$Feature_Type = as.character(annot$Feature_Type[resiFeatTag])
resiAnnot$Feature_Type <- factor(resiAnnot$Feature_Type, levels = unique(resiAnnot$Feature_Type))

resiAnnot$Bin = rep("", sum(resiFeatTag))
resiAnnot$Bin[grepl('_bin1$', row.names(resiAnnot))] = 'Bin 1'
resiAnnot$Bin[grepl('_bin2$', row.names(resiAnnot))] = 'Bin 2'
resiAnnot$Bin[grepl('_bin3$', row.names(resiAnnot))] = 'Bin 3'
resiAnnot$Bin[grepl('_bin4$', row.names(resiAnnot))] = 'Bin 4'
resiAnnot$Bin[grepl('^Ca2\\+$', row.names(resiAnnot))] = 'Bin 1'
resiAnnot$Bin <- factor(resiAnnot$Bin, levels = unique(resiAnnot$Bin))

resiFeat_stats = stats
row.names(resiFeat_stats)[resiFeatTag] = row.names(resiAnnot)

Feature_Type <- unique(featColors[resiFeatTag])
names(Feature_Type) <- levels(resiAnnot$Feature_Type)

Bin <- c('firebrick3', 'darkorange2', 'darkgoldenrod2', 'gold2')
names(Bin) <- levels(resiAnnot$Bin)

resiAnnot_cols <- list(Feature_Type = Feature_Type, Bin = Bin)

resiFeat_stats = t(resiFeat_stats[resiFeatTag,grepl('_effectSize$', colnames(resiFeat_stats))])
row.names(resiFeat_stats) = gsub('_effectSize$', '', row.names(resiFeat_stats))


CairoPDF(file = paste('./analysis/fine_specificity/lactose/', 
                      'resiFeats_heatmap',
                      '.pdf', sep = ''),
         width = 28,
         height = 7)
pheatmap(resiFeat_stats,
         color = colorRampPalette(c("royalblue1", "ivory", "gold1"))(length(breakLst)),
         border_color = 'white',
         cellwidth = cWidth,
         cellheight = cHeight,
         cluster_rows = F,
         # clustering_distance_rows = 'correlation',
         display_numbers = ifelse((hImpFeats & hSigFeats)[,resiFeatTag] , "\u2022", ""),
         fontsize_number = 20, number_color = 'black',
         labels_row = c('Lac vs all', 'LacNAc vs all', 'Lac vs LacNAc', 'LacNAc vs Lac'),
         annotation_col = resiAnnot,
         annotation_colors = resiAnnot_cols,
         # legend_breaks = c(1,5),
         main = '',
         gaps_row = c(2),
         breaks = breakLst,
         show_colnames = T,
         fontsize_col = 7,
         angle_col = 45)
dev.off()

####
## Pocket features only
####

pockAnnot = data.frame(Feature_Type = rep("", sum(pocketFeatTag)))
row.names(pockAnnot) = row.names(stats)[pocketFeatTag]

pockAnnot$Feature_Type = as.character(annot$Feature_Type[pocketFeatTag])
pockAnnot$Feature_Type <- factor(pockAnnot$Feature_Type, levels = unique(pockAnnot$Feature_Type))

pockAnnot$Threshold = rep("", sum(pocketFeatTag))
pockAnnot$Threshold[grepl('_4Ang$', row.names(pockAnnot))] = '4 Ang thresh'
pockAnnot$Threshold[grepl('_6Ang$', row.names(pockAnnot))] = '6 Ang thresh'
pockAnnot$Threshold[grepl('_8Ang$', row.names(pockAnnot))] = '8 Ang thresh'
pockAnnot$Threshold[grepl('_10Ang$', row.names(pockAnnot))] = '10 Ang thresh'
pockAnnot$Threshold[pockAnnot$Threshold == ""] = 'NA'
pockAnnot$Threshold <- factor(pockAnnot$Threshold, levels = unique(pockAnnot$Threshold))



Feature_Type <- unique(featColors[pocketFeatTag])
names(Feature_Type) <- levels(pockAnnot$Feature_Type)

Threshold <- c('firebrick3', 'darkorange2', 'darkgoldenrod2', 'gold2', 'grey75')
names(Threshold) <- levels(pockAnnot$Threshold)

pockAnnot_cols <- list(Feature_Type = Feature_Type, Threshold = Threshold)


pocketFeat_stats = t(stats[pocketFeatTag,grepl('_effectSize$', colnames(stats))])
row.names(pocketFeat_stats) = gsub('_effectSize$', '', row.names(pocketFeat_stats))

CairoPDF(file = paste('./analysis/fine_specificity/lactose/', 
                      'pocketFeats_heatmap',
                      '.pdf', sep = ''),
         width = 24,
         height = 7)
pheatmap(pocketFeat_stats,
         color = colorRampPalette(c("royalblue1", "ivory", "gold1"))(length(breakLst)),
         border_color = 'white',
         cellwidth = cWidth + 3,
         cellheight = cHeight,
         cluster_rows = F,
         # clustering_distance_rows = 'correlation',
         display_numbers = ifelse((hImpFeats & hSigFeats)[,pocketFeatTag] , "\u2022", ""),
         fontsize_number = 20, number_color = 'black',
         labels_row = c('Lac vs all', 'LacNAc vs all', 'Lac vs LacNAc', 'LacNAc vs Lac'),
         annotation_col = pockAnnot,
         annotation_colors = pockAnnot_cols,
         # legend_breaks = c(1,5),
         main = '',
         gaps_row = c(2),
         breaks = breakLst,
         show_colnames = T,
         fontsize_col = 9,
         angle_col = 45)
dev.off()


####
## PLIP interaction counts
####

plipFeat_stats = t(stats[1:11,grepl('_effectSize$', colnames(stats))])
row.names(plipFeat_stats) = gsub('_effectSize$', '', row.names(plipFeat_stats))


CairoPDF(file = paste('./analysis/fine_specificity/lactose/',
                      'plipFeats_heatmap',
                      '.pdf', sep = ''),
         width = 24,
         height = 7)
pheatmap(plipFeat_stats,
         color = colorRampPalette(c("royalblue1", "ivory", "gold1"))(length(breakLst)),
         border_color = 'ivory',
         cellwidth = cWidth+3,
         cellheight = cHeight, 
         cluster_rows = F,
         display_numbers = ifelse((hImpFeats & hSigFeats)[,1:11] , "\u2022", ""),
         fontsize_number = 20, number_color = 'black',
         labels_row = c('Lac vs all', 'LacNAc vs all', 'Lac vs LacNAc', 'LacNAc vs Lac'),
         main = '',
         gaps_row = c(2),
         breaks = breakLst,
         show_colnames = T,
         fontsize_col = 9,
         angle_col = 45)
dev.off()

