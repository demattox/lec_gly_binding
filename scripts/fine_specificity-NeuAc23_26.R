

library(survey)
# library(randomForest)
# library(caret)
# library(vcd)
library(scales)
library(colorspace)
library(philentropy)
library(vioplot)
library(pheatmap)
library(Cairo)
library(protr)
library(doParallel)
library(umap)

#######################
# Functions
#######################
# 
# pullDiverseCases <- function(scaledData, startInds, thresh, prevSampled = NA){
#   
#   if (any(is.na(prevSampled))) {
#     # if no inds being passed in that were sampled from the positive class
#     if (length(startInds) == 1) {
#       out = startInds
#     } else if (length(startInds) == 2) {
#       if (suppressMessages(distance(scaledFeats[startInds,])) >= thresh) {
#         out = startInds # Two binding sites are greater than the threshold distance from each other, keep both
#       } else{
#         out = sample(startInds, size = 1) # Two binding sites are within the threshold distance from each other, pick one at random
#       }
#       
#     } else {
#       cnter = 1
#       out = rep(0, length(startInds)) # Hold the set of diverse indices
#     }
#   } else{
#     # if inds are passed in from the positive class
#     
#     out = c(prevSampled, rep(0, length(startInds)))
#     cnter = length(prevSampled) + 1
#     
#     if (length(startInds) == 1) {
#       # One new neg site to compare to previously smapled positive sites
#       
#       distMat = suppressMessages(distance(x = rbind(scaledData[out[out != 0], ], scaledData[startInds, ]))) # Find all pairwise distances b/w binding sites in cluster
#       if (is.matrix(distMat)){ # if more than one previously sampled binding site
#         distMat = distMat[1:sum(out != 0), -1 * (1:sum(out != 0))] # Get the top n rows of the distance matrix dropping the first n columns, where n = the number of indices already sampled
#       }
#       if (!any(distMat < thresh)) {
#         out[cnter] = startInds
#       }
#     }
#   }
#   
#   
#   if (any(out == 0)){
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
#   }
#   
#   out = out[out != 0]
#   return(out)
# }
# 
# sampleDiverseSitesByLig <- function(clusterIDs, testClust, featureSet, ligandTag, distThresh, scaledFeatureSet = featureSet){
#   # Sample diverse binding sites with ligand and without ligand [ligandTag] for each cluster in clusterIDs, except the cluster held out for LO(C)O validation indicated by testClust
#   # Samples binding sites from the specified feature set, calculates Euclidean distance b/w binding sites from scaledFeatureSet
#   # Binding sites sampled randomly if Euc. distance to any previously sampled binding sites is greater than distThresh (median pariwise distance between all binding sites)
#   
#   # Drop the excluded cluster
#   uniClusts = unique(clusterIDs)
#   uniClusts = uniClusts[ ! uniClusts %in% testClust]
#   
#   dat = as.data.frame(matrix(0, nrow = nrow(featureSet), ncol = ncol(featureSet)))
#   colnames(dat) = colnames(predFeats)
#   dat$bound = F
#   dat$clus = 0
#   
#   j = 1 # index for writing to returned dataframe (dat)
#   
#   for (i in 1:length(uniClusts)) {
#     inds =  (1:nrow(featureSet))[clusterIDs == uniClusts[i]] # all the row indices matching a specific cluster 
#     negInds = inds[! ligandTag[inds]] # Indices of binding sites w/o ligand
#     posInds = inds[ligandTag[inds]] # Indices of binding sites w/ ligand
#     
#     if (length(posInds) > 0){
#       outInds = pullDiverseCases(scaledData = scaledFeatureSet, startInds = posInds, thresh = distThresh)
#     } else {
#       outInds = NA
#     }
#     if (length(negInds > 0)){
#       outInds = pullDiverseCases(scaledData = scaledFeatureSet, startInds = negInds, thresh = distThresh, prevSampled = outInds)
#     }
#     
#     dat[(j:(j+length(outInds) - 1)), (1:ncol(predFeats))] = predFeats[outInds, ] # set feature values for representative binding sites
#     dat$bound[(j:(j+length(outInds) - 1))] = ligandTag[outInds] # Set bound variable
#     dat$clus[(j:(j+length(outInds) - 1))] = uniClusts[i] # Set cluster ID
#     
#     j = j + length(outInds)
#   }
#   dat = dat[-1*(j:nrow(dat)), ]
#   dat$bound = as.factor(dat$bound)
#   return(dat)
# }

# f2 <- function (data, lev = NULL, model = NULL, beta = 2) {
#   precision <- posPredValue(data$pred, data$obs, positive = "TRUE")
#   recall  <- sensitivity(data$pred, data$obs, postive = "TRUE")
#   f2_val <- ((1 + beta^2) * precision * recall) / (beta^2 * precision + recall)
#   names(f2_val) <- c("F2")
#   f2_val
# }

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

# getKPRFb <- function(conMatDF){
#   # conMatDF has rows of different confusion matrices, with the columns ordered as TP, TN, FP, FN
#   # sums each column and finds performance metrics
#   f2TestCon = apply(X = conMatDF, MARGIN = 2, FUN = sum)
#   
#   TP = f2TestCon[[1]]
#   TN = f2TestCon[[2]]
#   FP = f2TestCon[[3]]
#   FN = f2TestCon[[4]]
#   
#   f2_validationRecall = TP / ( TP + FN)
#   f2_validationPrec = TP / (TP + FP)
#   f2_validationF2 = ((1+(2^2)) * f2_validationPrec * f2_validationRecall) / (2^2 * f2_validationPrec + f2_validationRecall)
#   # f3score = ((1+(3^2)) * f2_validationPrec * f2_validationRecall) / (3^2 * f2_validationPrec + f2_validationRecall)
#   # f4score = ((1+(4^2)) * f2_validationPrec * f2_validationRecall) / (4^2 * f2_validationPrec + f2_validationRecall)
#   randAcc = ((TN+FP)*(TN+FN) + (FN+TP)*(FP+TP)) / sum(c(TP + TN + FP + FN))^2
#   testAcc = (TP+TN)/(TP+TN+FP+FN)
#   f2_validationKappa = (testAcc - randAcc) / (1 - randAcc)
#   return(list(kappa = f2_validationKappa,
#               recall = f2_validationRecall,
#               F2 = f2_validationF2,
#               precision = f2_validationPrec))
# }

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

# perc.rank <- function(x) trunc(rank(x))/length(x)

###########################

homeDir = '/Users/dmattox/cbk/lec_gly_binding/'
setwd(homeDir)

ligTags = read.delim(file = './analysis/training/data_in/ligTags.tsv', sep = '\t', stringsAsFactors = F)
predFeats = read.delim(file = './analysis/training/data_in/predFeats.csv', sep = ',', stringsAsFactors = F)
bsResiDat = read.delim(file = './analysis/training/data_in/bsResiDat.tsv', sep = '\t', stringsAsFactors = F)

# load('./analysis/training/surveyObject.RData')

stats = read.delim(file = './analysis/training/weightedWMW_stats.tsv', sep = '\t', stringsAsFactors = F)


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
## Sialic acid glys
########################

uniLigs = unique(bsResiDat$iupac)

cpl50 = rep(0, length(uniLigs))
for (i in 1:length(uniLigs)){
  cpl50[i] = length(unique(bsResiDat$seqClust50[bsResiDat$iupac == uniLigs[i]]))
}

ligSort50 = uniLigs[order(cpl50, decreasing = T)]
cpl50 = cpl50[order(cpl50, decreasing = T)]


parenCnt = bracCnt = manCnt = neuCnt = bracCnt = rep(0,length(ligSort50))
for (i in 1:length(ligSort50)){
  lig = ligSort50[i]
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

manTag = manCnt >= 3 & grepl('^Man',ligSort50) # High mannose
neuTag = grepl('^NeuAc',ligSort50) & !mTag # Has terminal sialic acid and is not a monosacc.
fucTag = grepl('^Fuc',ligSort50) & !mTag # Has a terminal fucose and is not a monosacc.

# ligSort50[mTag]
# ligSort50[dTag]
# ligSort50[tTag]
# ligSort50[qTag]
# ligSort50[pTag]
#  
# ligSort50[bTag]

ligSort50[grepl('Kdo',ligSort50)]
ligSort50[grepl('Kdn',ligSort50)]


# NeuGc-containing glycans
ngnaTag = grepl('NeuGc',ligSort50) # Has terminal sialic acid and is not a monosacc.
ligSort50[ngnaTag]
cpl50[ngnaTag]


# a2-3 vs a2-6 matched glycans
neu23Tag = grepl('^NeuAc\\(a2-3\\)',ligSort50)
ligSort50[neu23Tag]
neu26Tag = grepl('^NeuAc\\(a2-6\\)',ligSort50)
ligSort50[neu26Tag]

sum(gsub('^NeuAc\\(a2-3\\)', 'NeuAc(a2-6)', ligSort50[neu23Tag]) %in% ligSort50[neu26Tag])

# matched23 = ligSort50[neu23Tag][gsub('^NeuAc\\(a2-3\\)', 'NeuAc(a2-6)', ligSort50[neu23Tag]) %in% ligSort50[neu26Tag]]
# matched26 = ligSort50[neu26Tag][gsub('^NeuAc\\(a2-6\\)', 'NeuAc(a2-3)', ligSort50[neu26Tag]) %in% ligSort50[neu23Tag]]
# 
# matched26 = matched26[c(3,1,2)]

# for (i in 1:3){
#   cat(matched23[i], '\n\t')
#   cat(cpl50[ligSort50 == matched23[i]], '\t')
#   cat(sum(bsResiDat$iupac == matched23[i]),'\n\n')
# }
# for(i in 1:3){
#   cat(matched26[i], '\n\t')
#   cat(cpl50[ligSort50 == matched26[i]], '\t')
#   cat(sum(bsResiDat$iupac == matched26[i]),'\n\n')
# }
# 
# 
# for (i in 1:3){
#   cat(matched23[i],'\t')
#   cat(sum(unique(bsResiDat$seqClust50[bsResiDat$iupac == matched23[i]]) %in% unique(bsResiDat$seqClust50[bsResiDat$iupac == matched26[i]])), '\n2-3: ')
#   cat(unique(bsResiDat$seqClust50[bsResiDat$iupac == matched23[i]])[order(unique(bsResiDat$seqClust50[bsResiDat$iupac == matched23[i]]))], '\n2-6: ')
#   cat(unique(bsResiDat$seqClust50[bsResiDat$iupac == matched26[i]])[order(unique(bsResiDat$seqClust50[bsResiDat$iupac == matched26[i]]))], '\n\n')
# }
# 
# length(unique(bsResiDat$seqClust50[bsResiDat$iupac %in% matched23]))
# 
# length(unique(bsResiDat$seqClust50[bsResiDat$iupac %in% matched26]))
# 
# sum(unique(bsResiDat$seqClust50[bsResiDat$iupac %in% matched23]) %in% unique(bsResiDat$seqClust50[bsResiDat$iupac %in% matched26]))
# 

tag23 = bsResiDat$iupac %in% ligSort50[neu23Tag]
tag26 = bsResiDat$iupac %in% ligSort50[neu26Tag]

neuClusts = unique(bsResiDat$seqClust50[tag23 | tag26])
neuClusts = neuClusts[order(neuClusts)]
for (i in 1:length(neuClusts)){
  clu = neuClusts[i]
  cat(clu, '\t\t')

  cat('\t', unique(bsResiDat$origine[bsResiDat$seqClust50 == clu]),'\n')
  cat('', unique(bsResiDat$espece[bsResiDat$seqClust50 == clu]),'\n')
}


########################
## Neua2-3/a2-6 - Descriptive
########################
tag = tag23 | tag26
bsResiDat = bsResiDat[tag,]
predFeats = predFeats[tag,]
tag23 = tag23[tag]
tag26 = tag26[tag]


scaledFeats = predFeats  # Scale features between 0 & 1
zeroCol = rep(F, ncol(predFeats))
for(i in 1:ncol(scaledFeats)){
  if (all(scaledFeats[,i] == 0)){
    zeroCol[i] = T
  } else{
    scaledFeats[,i] = (scaledFeats[,i] - min(scaledFeats[,i])) / (max(scaledFeats[,i]) - min(scaledFeats[,i]))
  }
}
scaledFeats = scaledFeats[,!zeroCol]

neu.umap = umap(scaledFeats)

summary(neu.umap$layout)

# UMAP plot
plot(x = neu.umap$layout[,1], y = neu.umap$layout[,2],
     pch = 19, #col = alpha(colors[labels.int], 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = "UMAP vis for 3' & 6' interactions",
     xlim = c(-7.8,7), ylim = c(-10, 7.8))


# Color by 3' vs 6'
plot(x = neu.umap$layout[tag23,1], y = neu.umap$layout[tag23,2],
     pch = 19, cex = 1.5,
     col = alpha('dodgerblue1', 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = "UMAP vis for 3' & 6' interactions",
     xlim = c(-7.8,7), ylim = c(-10, 7.8))
par(new = T)
plot(x = neu.umap$layout[tag26,1], y = neu.umap$layout[tag26,2],
    #xlim = c(-6,7), ylim = c(-6,7),
    pch = 15, cex = 1.5,
    # col = alpha(ligColors[2], 0.6),
    col = alpha('firebrick2', 0.6),
    xlab = '', ylab = '', main = '',
    xlim = c(-7.8,7), ylim = c(-10, 7.8))
legend(x = 'bottomright', legend = c("3' SA", "6' SA"), col = c('dodgerblue1', 'firebrick2'), pch = c(19, 15))


# Color by species of origin
specList = unique(bsResiDat$origine)
specCols = colfunc(length(specList))

plot(x = neu.umap$layout[bsResiDat$origine == specList[1],1], y = neu.umap$layout[bsResiDat$origine == specList[1],2],
     #xlim = c(-6,7), ylim = c(-6,7),
     pch = 19, cex = 1.5,
     col = alpha(specCols[1], 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = "UMAP vis for 3' & 6' interactions",
     xlim = c(-7.8,7), ylim = c(-10, 7.8))
for (i in 2:length(specList)){
  par(new = T)
  plot(x = neu.umap$layout[bsResiDat$origine == specList[i],1], y = neu.umap$layout[bsResiDat$origine == specList[i],2],
       #xlim = c(-6,7), ylim = c(-6,7),
       pch = 19, cex = 1.5,
       # col = alpha(ligColors[2], 0.6),
       col = alpha(specCols[i], 0.6),
       xlab = '', ylab = '', main = '', axes = F,
       xlim = c(-7.8,7), ylim = c(-10, 7.8))
}

legend(x = 'bottomright', legend = specList, col = specCols, pch = 19)

# par(new = T)
# plot(x = neu.umap$layout[bsResiDat$seqClust50 %in% c(25,27,28,29),1], y = neu.umap$layout[bsResiDat$seqClust50 %in% c(25,27,28,29),2],
#      col = 'black',
#      xlim = c(-7.4,10), ylim = c(-10.2, 4))

# Plot by cluster ID
clusLst = unique(bsResiDat$seqClust50)
clusLst = clusLst[order(clusLst)]
clusCols = colfunc(length(clusLst))
set.seed(1)
clusCols = clusCols[sample(1:length(clusCols), size = length(clusCols), replace = F)]

plot(x = neu.umap$layout[bsResiDat$seqClust50 == clusLst[1],1], y = neu.umap$layout[bsResiDat$seqClust50 == clusLst[1],2],
     pch = 19, cex = 1.5,
     col = alpha(clusCols[1], 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = "UMAP vis for 3' & 6' interactions",
     xlim = c(-7.8,7), ylim = c(-10, 7.8))
for (i in 2:length(clusLst)){
  par(new = T)
  plot(x = neu.umap$layout[bsResiDat$seqClust50 == clusLst[i],1], y = neu.umap$layout[bsResiDat$seqClust50 == clusLst[i],2],
       pch = 19, cex = 1.5,
       # col = alpha(ligColors[2], 0.6),
       col = alpha(clusCols[i], 0.6),
       xlab = '', ylab = '', main = '', axes = F,
       xlim = c(-7.8,7), ylim = c(-10, 7.8))
}

# par(new = T)
# plot(x = neu.umap$layout[bsResiDat$seqClust50 %in% c(25,27,28,29),1], y = neu.umap$layout[bsResiDat$seqClust50 %in% c(25,27,28,29),2],
#      col = alpha('black',0.4), cex = 1.5,
#      xlab = '', ylab = '', axes = F,
#      xlim = c(-7.4,10), ylim = c(-10.2, 4))
par(new = T)
plot(x = neu.umap$layout[bsResiDat$seqClust50 == 47,1], y = neu.umap$layout[bsResiDat$seqClust50 == 47,2],
     col = alpha('black',0.4), cex = 1.5,
     xlab = '', ylab = '', axes = F,
     xlim = c(-7.8,7), ylim = c(-10, 7.8))



plot(0, axes = F, type = 'n', xlab = '',ylab='')
legend(x = 'center', legend = clusLst, col = clusCols, pch = 19, bty = 'n')


############################
## HA only
############################

unique(bsResiDat$seqClust50[grepl('flu', bsResiDat$espece)])

tag = bsResiDat$seqClust50 %in% c(25,27,28,29)
bsResiDat = bsResiDat[tag,]
predFeats = predFeats[tag,]
tag23 = tag23[tag]
tag26 = tag26[tag]


scaledFeats = predFeats  # Scale features between 0 & 1
zeroCol = rep(F, ncol(predFeats))
for(i in 1:ncol(scaledFeats)){
  if ((max(scaledFeats[,i]) - min(scaledFeats[,i]) == 0)){
    zeroCol[i] = T
  } else{
    scaledFeats[,i] = (scaledFeats[,i] - min(scaledFeats[,i])) / (max(scaledFeats[,i]) - min(scaledFeats[,i]))
  }
}
scaledFeats = scaledFeats[,!zeroCol]

neu.umap = umap(scaledFeats)

summary(neu.umap$layout)

xLim = c(-4.4,4.1)
yLim = c(-4, 7)

# UMAP plot
plot(x = neu.umap$layout[,1], y = neu.umap$layout[,2],
     pch = 19, #col = alpha(colors[labels.int], 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = "UMAP vis for 3' & 6' interactions",
     xlim = xLim, ylim = yLim)


# Color by 3' vs 6'
plot(x = neu.umap$layout[tag23,1], y = neu.umap$layout[tag23,2],
     pch = 19, cex = 1.5,
     col = alpha('dodgerblue1', 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = "UMAP vis for 3' & 6' interactions",
     xlim = xLim, ylim = yLim)
par(new = T)
plot(x = neu.umap$layout[tag26,1], y = neu.umap$layout[tag26,2],
     pch = 19, cex = 1.5,
     # col = alpha(ligColors[2], 0.6),
     col = alpha('firebrick2', 0.6),
     xlab = '', ylab = '', main = '',
     xlim = xLim, ylim = yLim)
legend(x = 'bottomright', legend = c("3' SA", "6' SA"), col = c('dodgerblue1', 'firebrick2'), pch = c(19, 19))

# Map UniProt IDs to flu type/subtype
bsResiDat$strain = ''
unique(bsResiDat$uniparc)

bsResiDat$strain[bsResiDat$uniparc == 'Q82500'] = 'Type A - H1N1'
bsResiDat$strain[bsResiDat$uniparc == 'P03452'] = 'Type A - H1N1'
bsResiDat$strain[bsResiDat$uniparc == 'C7C6F1'] = 'Type A - H1N1'
bsResiDat$strain[bsResiDat$uniparc == 'C3W5S1'] = 'Type A - H1N1'

bsResiDat$strain[bsResiDat$uniparc == 'P03437'] = 'Type A - H3N2'
bsResiDat$strain[bsResiDat$uniparc == 'P03442'] = 'Type A - H3N8'

bsResiDat$strain[bsResiDat$uniparc == 'Q6DQ34'] = 'Type A - H5N1'
bsResiDat$strain[bsResiDat$uniparc == 'Q207Z6'] = 'Type A - H5N1'
bsResiDat$strain[bsResiDat$uniparc == 'Q5EP31'] = 'Type A - H5N1'

bsResiDat$strain[bsResiDat$uniparc == 'B7NYS1'] = 'Type A - H7N2'
bsResiDat$strain[bsResiDat$uniparc == 'G0KQM4'] = 'Type A - H7N3'
bsResiDat$strain[bsResiDat$uniparc == 'Q6GYW3'] = 'Type A - H7N3'
bsResiDat$strain[bsResiDat$uniparc == 'M4YV75'] = 'Type A - H7N9'

bsResiDat$strain[bsResiDat$uniparc == 'Q0A448'] = 'Type A - H10N7'

bsResiDat$strain[bsResiDat$uniparc == 'P03462'] = 'Type B'


# Plot UMAP by strain
strainLst = unique(bsResiDat$strain)
strainLst = strainLst[order(strainLst)]
strainLst = strainLst[c(2:8, 1, 9)]

# tmpCols = colorRampPalette(c("darkblue",'deepskyblue',"mediumaquamarine","cyan2"))(5)
tmpCols = rainbow_hcl(6)
strainCols = c(tmpCols[1], rep(tmpCols[2],2), tmpCols[3], rep(tmpCols[4],3), tmpCols[5], tmpCols[6])
strainPnts = c(19, 19, 15, 19, 19, 15, 18, 19, 19)


plot(x = neu.umap$layout[bsResiDat$strain == strainLst[1],1], y = neu.umap$layout[bsResiDat$strain == strainLst[1],2],
     pch = strainPnts[1],
     cex = 1.5,
     col = alpha(strainCols[1], 0.7),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = "UMAP vis for 3' & 6' interactions",
     xlim = xLim, ylim = yLim)
for (i in 2:length(strainCols)){
  par(new = T)
  plot(x = neu.umap$layout[bsResiDat$strain == strainLst[i],1], y = neu.umap$layout[bsResiDat$strain == strainLst[i],2],
       pch = strainPnts[i],
       cex = 1.5,
       col = alpha(strainCols[i], 0.7),
       xlab = '', ylab = '', main = "", axes = F,
       xlim = xLim, ylim = yLim)
}

legend(x = 'bottomright', legend = strainLst, col = strainCols, pch = strainPnts)

###
# Weight by type/HA subtype
bsResiDat$type = ''
bsResiDat$type[grepl('Type A - H1N',bsResiDat$strain)] = 'H1'
bsResiDat$type[grepl('Type A - H3N',bsResiDat$strain)] = 'H3'
bsResiDat$type[grepl('Type A - H5N',bsResiDat$strain)] = 'H5'
bsResiDat$type[grepl('Type A - H7N',bsResiDat$strain)] = 'H7'
bsResiDat$type[grepl('Type A - H10N',bsResiDat$strain)] = 'H10'
bsResiDat$type[grepl('Type B',bsResiDat$strain)] = 'B'

clusWeight = nrow(bsResiDat)/length(unique(bsResiDat$type)) # The proportion of the total weight allotted to influenza type/subtype
cWeights = rep(0, nrow(bsResiDat)) # weights based on subtype membership only (all binding sites with each subtype receives the same weight)
cuWeights = rep(0,nrow(bsResiDat)) # weights based on subtype membership and UniProt ID (all subtypes get the same weight, and all unique lectins within each subtype also receive the same weight)

for (i in 1:length(unique(bsResiDat$type))){
  clus = unique(bsResiDat$type)[i]
  tag = bsResiDat$type == clus
  cWeights[tag] = clusWeight/sum(tag) # Equally divide allotted weight between all binding sites within a given cluster
  
  lectinWeight = clusWeight / length(unique(bsResiDat$uniparc[tag])) # Weight allotted to each lectin
  for(j in 1:length(unique(bsResiDat$uniparc[tag]))){ # Further distribute weights within clusters based on uniprot ids
    uniTag = bsResiDat$uniparc[tag] == unique(bsResiDat$uniparc[tag])[j]
    cuWeights[tag][uniTag] = lectinWeight / sum(uniTag)
  }
}
round(sum(cWeights),5) == nrow(bsResiDat) # Sum of weights is equal to number of binding sites
round(sum(cuWeights),5) == nrow(bsResiDat) # Sum of weights is equal to number of binding sites




strainWts = (cuWeights - min(cuWeights))/(max(cuWeights)-min(cuWeights)) * (3 - 1) + 1 # Scale interaction weights between 1 and 3
plot(x = neu.umap$layout[bsResiDat$strain == strainLst[1],1], y = neu.umap$layout[bsResiDat$strain == strainLst[1],2],
     pch = strainPnts[1],
     cex = strainWts[bsResiDat$strain == strainLst[1]],
     col = alpha(strainCols[1], 0.7),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = "UMAP vis for 3' & 6' interactions",
     xlim = xLim, ylim = yLim)
# text(x = mean(neu.umap$layout[bsResiDat$strain == strainLst[1],1]), y = mean(neu.umap$layout[bsResiDat$strain == strainLst[1],2]),
#      labels = strainLst[1],
#      col = strainCols[1])
for (i in 2:length(strainLst)){
  par(new = T)
  plot(x = neu.umap$layout[bsResiDat$strain == strainLst[i],1], y = neu.umap$layout[bsResiDat$strain == strainLst[i],2],
       pch = strainPnts[i],
       cex = strainWts[bsResiDat$strain == strainLst[i]],
       col = alpha(strainCols[i], 0.7),
       xlab = '', ylab = '', main = "", axes = F,
       xlim = xLim, ylim = yLim)
  # text(x = mean(neu.umap$layout[bsResiDat$strain == strainLst[i],1]), y = mean(neu.umap$layout[bsResiDat$strain == strainLst[i],2]),
  #      labels = strainLst[i],
  #      col = strainCols[i])
}
legend(x = 'bottomright', legend = strainLst, col = strainCols, pch = strainPnts, pt.cex = 2)
typeLst = unique(bsResiDat$type)
typeLst = typeLst[order(typeLst)]
typeLst = typeLst[c(2,4:6,3,1)]

for(i in 1:length(typeLst)){
  tag = bsResiDat$type == typeLst[i]
  text(x = mean(neu.umap$layout[tag,1]), y = mean(neu.umap$layout[tag,2]),
       labels = typeLst[i],
       col = unique(strainCols)[i],
       cex = 2.5,
       pos = 1)
}




# Color by subtype & ligand
plot(x = neu.umap$layout[bsResiDat$strain == strainLst[1] & tag26,1], y = neu.umap$layout[bsResiDat$strain == strainLst[1] & tag26,2],
     pch = strainPnts[1],
     cex = 2.5,
     col = alpha('firebrick1', 0.7),
     xlab = '', ylab = '', main = "", axes = F,
     xlim = xLim, ylim = yLim)
par(new=T)
plot(x = neu.umap$layout[bsResiDat$strain == strainLst[1] & tag23,1], y = neu.umap$layout[bsResiDat$strain == strainLst[1] & tag23,2],
     pch = strainPnts[1],
     cex = 2.5,
     col = alpha('dodgerblue2', 0.7),
     xlab = '', ylab = '', main = "", axes = F,
     xlim = xLim, ylim = yLim)
par(new=T)
plot(x = neu.umap$layout[bsResiDat$strain == strainLst[1],1], y = neu.umap$layout[bsResiDat$strain == strainLst[1],2],
     pch = strainPnts[1],
     cex = 1.5,
     col = alpha(strainCols[1], 0.9),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = "UMAP vis. for 3' & 6' SA interactions with HA",
     xlim = xLim, ylim = yLim)
# text(x = mean(neu.umap$layout[bsResiDat$strain == strainLst[1],1]), y = mean(neu.umap$layout[bsResiDat$strain == strainLst[1],2]),
#      labels = strainLst[1],
#      col = strainCols[1])
for (i in 2:length(strainLst)){
  par(new = T)
  plot(x = neu.umap$layout[bsResiDat$strain == strainLst[i] & tag26,1], y = neu.umap$layout[bsResiDat$strain == strainLst[i] & tag26,2],
       pch = strainPnts[i],
       cex = 2.5,
       col = alpha('firebrick1', 0.7),
       xlab = '', ylab = '', main = "", axes = F,
       xlim = xLim, ylim = yLim)
  par(new=T)
  plot(x = neu.umap$layout[bsResiDat$strain == strainLst[i] & tag23,1], y = neu.umap$layout[bsResiDat$strain == strainLst[i] & tag23,2],
       pch = strainPnts[i],
       cex = 2.5,
       col = alpha('dodgerblue2', 0.7),
       xlab = '', ylab = '', main = "", axes = F,
       xlim = xLim, ylim = yLim)
  par(new=T)
  plot(x = neu.umap$layout[bsResiDat$strain == strainLst[i],1], y = neu.umap$layout[bsResiDat$strain == strainLst[i],2],
       pch = strainPnts[i],
       cex = 1.5,
       col = alpha(strainCols[i], 0.9),
       xlab = '', ylab = '', main = "", axes = F,
       xlim = xLim, ylim = yLim)
  # text(x = mean(neu.umap$layout[bsResiDat$strain == strainLst[i],1]), y = mean(neu.umap$layout[bsResiDat$strain == strainLst[i],2]),
  #      labels = strainLst[i],
  #      col = strainCols[i])
}
legend(x = 'bottomright', legend = c(strainLst, '', "6' SA", "3' SA"), col = c(strainCols, 'white', 'firebrick1', 'dodgerblue2'), pch = c(strainPnts, rep(1,3)), pt.cex = 2, bty = 'n')

# 900 x 770

######
# HA Only WMW
######
bsResiDat$type[tag23]
bsResiDat$type[tag26]
for(i in 1:length(unique(bsResiDat$type))){
  cat(typeLst[i],'\n\t')
  cat(sum(tag26 & bsResiDat$type == typeLst[i]), '\t', sum(tag23 & bsResiDat$type == typeLst[i]), '\n\t')
  cat(sum(cuWeights[tag26 & bsResiDat$type == typeLst[i]]), '\t', sum(cuWeights[tag23 & bsResiDat$type == typeLst[i]]), '\n')
}

sum(tag26)
sum(cuWeights[tag26])
sum(tag23)
sum(cuWeights[tag23])


neuTags = cbind(tag26, tag23)

ligNames = c("6' NeuAc", "3' NeuAc")
ligColors = c('firebrick1', 'dodgerblue1')

colnames(neuTags) = c('a26_NeuAc', 'a23_NeuAc')

neuStats = as.data.frame(matrix(0, nrow = ncol(predFeats), ncol = (ncol(neuTags)*4)))
row.names(neuStats) = colnames(predFeats)
colnames(neuStats) = c(paste(colnames(neuTags), 'p', sep = '_'), paste(colnames(neuTags), 'FC', sep = '_'), paste(colnames(neuTags), 'effectSize', sep = '_'), paste(colnames(neuTags), 'adj', sep = '_'))

ligSpefic_feat_means = as.data.frame(matrix(0,nrow = ncol(predFeats), ncol = ncol(neuTags)))
colnames(ligSpefic_feat_means) = colnames(neuTags)
row.names(ligSpefic_feat_means) = colnames(predFeats)

scaledFeats = predFeats  # Scale features between 0 & 1
zeroCol = rep(F, ncol(predFeats))
for(i in 1:ncol(scaledFeats)){
  if ((max(scaledFeats[,i]) - min(scaledFeats[,i]) == 0)){
    zeroCol[i] = T
  } else{
    scaledFeats[,i] = (scaledFeats[,i] - min(scaledFeats[,i])) / (max(scaledFeats[,i]) - min(scaledFeats[,i]))
  }
}
scaledFeats[,zeroCol] = 0.5


wmwFeats = cbind(predFeats,neuTags)
des <- svydesign(ids = ~1, data = wmwFeats, weights = 1/cuWeights) # Build survey object with binding indicators for glycans of interest

scaledWMW_feats = cbind(scaledFeats,neuTags)
scaled_des <- svydesign(ids = ~1, data = scaledWMW_feats, weights = 1/cuWeights)

for(i in 1:ncol(neuTags)){
  
  des_allNEU= subset(des, subset = apply(neuTags, 1, any))
  
  des_w = subset(des, subset = neuTags[,i]) # temporary design object holding all interactions with the ligand of interest
  des_wo = subset(des, subset = neuTags[,(-1*i)]) # same as above but OTHER the ligands of interest
  scaled_des_w = subset(scaled_des, subset = neuTags[,i]) #
  
  for(k in 1:ncol(predFeats)){ # for each feature k
    ligTest = svyranktest(formula = as.formula(paste(colnames(predFeats)[k], ' ~ ', colnames(neuTags)[i], sep = '')), # Wilcoxon–Mann–Whitney test, sample sizes can be small (~5% of 230 clusters ~= 10), no reason to assume distribution is normal as it likely isn't
                          design = des_allNEU, 
                          test = 'wilcoxon') 
    
    if (ligTest$p.value != 0){
      neuStats[k, grepl('_p$', colnames(neuStats))][i] = ligTest$p.value # Raw p-value
    } else{
      neuStats[k, grepl('_p$', colnames(neuStats))][i] = 5e-324 # If p-value is too small and R rounds to 0, reset as the smallest positive double as referenced by "?.Machine"
    }
    
    neuStats[k, grepl('_effectSize$', colnames(neuStats))][i] = ligTest$estimate # common language effect size (0.5 transformed to 0)
    
    neuStats[k, grepl('_FC$', colnames(neuStats))][i] = svymean(predFeats[neuTags[,i],k], des_w)[1] / svymean(predFeats[neuTags[,(-1*i)],k], des_wo)[1] # Fold change in weighted means
    
    ligSpefic_feat_means[k,i] = svymean(scaledFeats[neuTags[,i],k], scaled_des_w)[1] # Get the weighted mean for each feature in interactions with each glycan of interest
  }
  
  neuStats[,grepl('_adj$', colnames(neuStats))][,i] = p.adjust(neuStats[grepl('_p$', colnames(neuStats))][,i], method = "BH") # Benjamini-Hochberg MHT correction (FDR)
}


# neuStats = neuStats[,grepl('^a26_NeuAc', colnames(neuStats))]

sum(neuStats[,grepl('_adj$', colnames(neuStats))] < 1e-16)

# superSigTag = neuStats[,grepl('_adj$', colnames(neuStats))] < 1e-16
# neuStats[,grepl('_adj$', colnames(neuStats))][superSigTag] <- 10**(-1*runif(sum(superSigTag), max = -log10(3e-19), min = -log10(1e-16))) # Sample from a log-uniform distribution


#####
# HA Only Volcano Plots
#####

xLim = c(-0.5,0.5)
yLim = c(0,(-log10(min(neuStats[,grepl('_adj$', colnames(neuStats))])) + 1))
# for(i in 1:ncol(tag26)){
for(i in 1){
  # yLim = c(0,max(-log10(0.1), -log10(min(neuStats[,grepl('_adj$', colnames(neuStats))][,i]))) + 1)
  
  tag = neuStats[,grepl('_adj$', colnames(neuStats))][,i] < 0.01
  
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
  
  plot(neuStats[,grepl('_effectSize$', colnames(neuStats))][,i], -log10(neuStats[,grepl('_adj$', colnames(neuStats))][,i]), # Plot all points w/ color @ alpha 0.5
       xlab = "Effect size", ylab = "-log10(FDR)", main = '',
       pch = 19, cex = 2, col = alpha(featColors, 0.33),
       cex.axis = 1.5, cex.lab = 1.5,
       xlim = xLim, ylim = yLim)
  
  abline(h= -log10(0.01), lwd = 4, col = 'white')
  # abline(h = -log10(1e-16), lwd = 1.5, lty = 2, col = 'white')
  
  title(main = ligNames[i], col.main = ligColors[i], cex.main = 1.8, font.main  = 2)
  # mtext(ligNames[i], col = ligColors[i],
  #       side=3, adj=0, outer = T,
  #       line=1.2, cex=1, font=2)
  
  par(new=T)
  plot(neuStats[,grepl('_effectSize$', colnames(neuStats))][tag,i], -log10(neuStats[,grepl('_adj$', colnames(neuStats))][tag,i]), # Plot stat sig points again with alpha 1
       pch = 19, col = featColors[tag], cex = 2,
       axes = F, xlab = "", ylab = "", main = "",
       xlim = xLim, ylim = yLim)
  par(new=T)
  plot(neuStats[,grepl('_effectSize$', colnames(neuStats))][tag,i], -log10(neuStats[,grepl('_adj$', colnames(neuStats))][tag,i]), # Outline stat sig points in black to highlight them
       col = alpha('black',0.2), cex = 2.05,
       xlab = "", ylab = "", axes = F, main = "",
       xlim = xLim, ylim = yLim)
}



#####
# HA Only Heatmaps
#####

neuStats = neuStats[,!grepl('^a23', colnames(neuStats))]

stats = stats[,(grepl('NeuAc', colnames(stats)) | grepl('Sialic', colnames(stats)))]

stats = cbind(stats, neuStats)

sigFeats = neuStats[,grepl('_adj$', colnames(neuStats))] < 0.01

hSigFeats = rbind(rep(F,length(sigFeats)),rep(F,length(sigFeats)),rep(F,length(sigFeats)),sigFeats)

allFeatCorrs = cor(stats[,grepl('_effectSize$', colnames(stats))], stats[,grepl('_effectSize$', colnames(stats))], method = 'pearson')
row.names(allFeatCorrs)  = gsub('_effectSize$', '', row.names(allFeatCorrs))
row.names(allFeatCorrs)[3:4] = c("3' SL", "6' vs 3' NeuAc")
colnames(allFeatCorrs) = row.names(allFeatCorrs)

corrplot::corrplot(allFeatCorrs)

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

dev.off()

CairoPDF(file = paste('./analysis/fine_specificity/HA_23-26/', 
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
         display_numbers = ifelse((hSigFeats)[,resiFeatTag] , "\u2022", ""),
         fontsize_number = 20, number_color = 'black',
         labels_row = c('Sialic acid v all',
                        'NeuAc vs all',
                        'SL vs all',
                        
                        "6' vs 3' NeuAc (HA)"),
         annotation_col = resiAnnot,
         annotation_colors = resiAnnot_cols,
         # legend_breaks = c(1,5),
         main = '',
         gaps_row = c(3),
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

CairoPDF(file = paste('./analysis/fine_specificity/HA_23-26/', 
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
         display_numbers = ifelse((hSigFeats)[,pocketFeatTag] , "\u2022", ""),
         fontsize_number = 20, number_color = 'black',
         labels_row = c('Sialic acid v all',
                        'NeuAc vs all',
                        'SL vs all',
                        
                        "6' vs 3' NeuAc (HA)"),
         annotation_col = pockAnnot,
         annotation_colors = pockAnnot_cols,
         # legend_breaks = c(1,5),
         main = '',
         gaps_row = c(3),
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


CairoPDF(file = paste('./analysis/fine_specificity/HA_23-26/',
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
         display_numbers = ifelse((hSigFeats)[,1:11] , "\u2022", ""),
         fontsize_number = 20, number_color = 'black',
         labels_row = c('Sialic acid v all',
                        'NeuAc vs all',
                        'SL vs all',
                        
                        "6' vs 3' NeuAc (HA)"),
         main = '',
         gaps_row = c(3),
         breaks = breakLst,
         show_colnames = T,
         fontsize_col = 9,
         angle_col = 45)
dev.off()


#####
# HA Only significant features
#####

row.names(neuStats)[neuStats$a26_NeuAc_adj < 0.01 & neuStats$a26_NeuAc_effectSize > 0]
upFeats = neuStats[neuStats$a26_NeuAc_adj < 0.01 & neuStats$a26_NeuAc_effectSize > 0,]
upFeats$y = -log10(upFeats$a26_NeuAc_adj)
upFeats = upFeats[order(upFeats$y,decreasing = T),]
upFeats

row.names(neuStats)[neuStats$a26_NeuAc_adj < 0.01 & neuStats$a26_NeuAc_effectSize < 0]
downFeats = neuStats[neuStats$a26_NeuAc_adj < 0.01 & neuStats$a26_NeuAc_effectSize < 0,]
downFeats$y = -log10(downFeats$a26_NeuAc_adj)
downFeats = downFeats[order(downFeats$y,decreasing = T),]
downFeats


#####
# HA Only UMAP with sig feats
#####


scaledFeats = predFeats  # Scale features between 0 & 1
zeroCol = rep(F, ncol(predFeats))
for(i in 1:ncol(scaledFeats)){
  if ((max(scaledFeats[,i]) - min(scaledFeats[,i]) == 0)){
    zeroCol[i] = T
  } else{
    scaledFeats[,i] = (scaledFeats[,i] - min(scaledFeats[,i])) / (max(scaledFeats[,i]) - min(scaledFeats[,i]))
  }
}

neu.umap = umap(scaledFeats[,sigFeats])

summary(neu.umap$layout)

xLim = c(-3.6,2.9)
yLim = c(-4, 3.5)

# UMAP plot
plot(x = neu.umap$layout[,1], y = neu.umap$layout[,2],
     pch = 19, #col = alpha(colors[labels.int], 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = "UMAP vis for 3' & 6' interactions",
     xlim = xLim, ylim = yLim)


# Color by 3' vs 6'
plot(x = neu.umap$layout[tag23,1], y = neu.umap$layout[tag23,2],
     pch = 19, cex = 1.5,
     col = alpha('dodgerblue1', 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = "UMAP vis for 3' & 6' interactions",
     xlim = xLim, ylim = yLim)
par(new = T)
plot(x = neu.umap$layout[tag26,1], y = neu.umap$layout[tag26,2],
     pch = 19, cex = 1.5,
     # col = alpha(ligColors[2], 0.6),
     col = alpha('firebrick2', 0.6),
     xlab = '', ylab = '', main = '',
     xlim = xLim, ylim = yLim)
legend(x = 'topright', legend = c("3' SA", "6' SA"), col = c('dodgerblue1', 'firebrick2'), pch = c(19, 19))


# Color by subtype & ligand
plot(x = neu.umap$layout[bsResiDat$strain == strainLst[1] & tag26,1], y = neu.umap$layout[bsResiDat$strain == strainLst[1] & tag26,2],
     pch = strainPnts[1],
     cex = 2.5,
     col = alpha('firebrick2', 0.7),
     xlab = '', ylab = '', main = "", axes = F,
     xlim = xLim, ylim = yLim)
par(new=T)
plot(x = neu.umap$layout[bsResiDat$strain == strainLst[1] & tag23,1], y = neu.umap$layout[bsResiDat$strain == strainLst[1] & tag23,2],
     pch = strainPnts[1],
     cex = 2.5,
     col = alpha('dodgerblue1', 0.7),
     xlab = '', ylab = '', main = "", axes = F,
     xlim = xLim, ylim = yLim)
par(new=T)
plot(x = neu.umap$layout[bsResiDat$strain == strainLst[1],1], y = neu.umap$layout[bsResiDat$strain == strainLst[1],2],
     pch = strainPnts[1],
     cex = 1.5,
     col = alpha(strainCols[1], 0.9),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = "HA UMAP vis. for 3' & 6' interactions with sig. features",
     xlim = xLim, ylim = yLim)
# text(x = mean(neu.umap$layout[bsResiDat$strain == strainLst[1],1]), y = mean(neu.umap$layout[bsResiDat$strain == strainLst[1],2]),
#      labels = strainLst[1],
#      col = strainCols[1])
for (i in 2:length(strainLst)){
  par(new = T)
  plot(x = neu.umap$layout[bsResiDat$strain == strainLst[i] & tag26,1], y = neu.umap$layout[bsResiDat$strain == strainLst[i] & tag26,2],
       pch = strainPnts[i],
       cex = 2.5,
       col = alpha('firebrick2', 0.7),
       xlab = '', ylab = '', main = "", axes = F,
       xlim = xLim, ylim = yLim)
  par(new=T)
  plot(x = neu.umap$layout[bsResiDat$strain == strainLst[i] & tag23,1], y = neu.umap$layout[bsResiDat$strain == strainLst[i] & tag23,2],
       pch = strainPnts[i],
       cex = 2.5,
       col = alpha('dodgerblue1', 0.7),
       xlab = '', ylab = '', main = "", axes = F,
       xlim = xLim, ylim = yLim)
  par(new=T)
  plot(x = neu.umap$layout[bsResiDat$strain == strainLst[i],1], y = neu.umap$layout[bsResiDat$strain == strainLst[i],2],
       pch = strainPnts[i],
       cex = 1.5,
       col = alpha(strainCols[i], 0.9),
       xlab = '', ylab = '', main = "", axes = F,
       xlim = xLim, ylim = yLim)
  # text(x = mean(neu.umap$layout[bsResiDat$strain == strainLst[i],1]), y = mean(neu.umap$layout[bsResiDat$strain == strainLst[i],2]),
  #      labels = strainLst[i],
  #      col = strainCols[i])
}
legend(x = 'topright', legend = c(strainLst, '', "6' SA", "3' SA"), col = c(strainCols, 'white', 'firebrick2', 'dodgerblue1'), pch = c(strainPnts, rep(1,3)), pt.cex = 2, bty = 'n')




# Color by 3' vs 6' and look for example interactions
plot(x = neu.umap$layout[tag23,1], y = neu.umap$layout[tag23,2],
     pch = 19, cex = 1.5,
     col = alpha('dodgerblue1', 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = "UMAP vis for 3' & 6' interactions \nand representative interactions",
     xlim = xLim, ylim = yLim)
par(new = T)
plot(x = neu.umap$layout[tag26,1], y = neu.umap$layout[tag26,2],
     pch = 19, cex = 1.5,
     # col = alpha(ligColors[2], 0.6),
     col = alpha('firebrick2', 0.6),
     xlab = '', ylab = '', main = '', axes = F,
     xlim = xLim, ylim = yLim)
legend(x = 'topleft', legend = c("3' SA", "6' SA"), col = c('dodgerblue1', 'firebrick2'), pch = c(19, 19))


for (i in 1:ncol(neuTags)){
  cat(colnames(neuTags)[i], '\n')
  boundDat = rbind(ligSpefic_feat_means[sigFeats,i], scaledFeats[neuTags[,i],sigFeats])
  dists = distance(boundDat, method = 'euclidean', use.row.names = T)
  
  cat(row.names(boundDat)[which.min(dists[1,2:nrow(boundDat)]) +1],'\t')
  cat(min(dists[1,2:nrow(boundDat)]),'\n\n')
}

# neu26.umap = predict(neu.umap, t(ligSpefic_feat_means[sigFeats,]))

par(new=T)
plot(x = neu.umap$layout["3UBE_GAL:M:1",1], y = neu.umap$layout["3UBE_GAL:M:1",2],
     cex = 3,
     lwd = 3,
     col = alpha('firebrick2', 1),
     xlab = '', ylab = '', main = '', axes = F,
     xlim = xLim, ylim = yLim)
text(x = neu.umap$layout["3UBE_GAL:M:1",1], y = neu.umap$layout["3UBE_GAL:M:1",2],
     labels = "3UBE_GAL:M:1", pos = 3, offset = 1, col = 'firebrick2')
# par(new=T)
# plot(x = neu26.umap["a26_NeuAc",1], y = neu26.umap["a26_NeuAc",2],
#      cex =1,
#      pch=19,
#      col = alpha('firebrick2', 1),
#      xlab = '', ylab = '', main = '', axes = F,
#      xlim = xLim, ylim = yLim)




par(new=T)
plot(x = neu.umap$layout["1RVX_NAG:R:1",1], y = neu.umap$layout["1RVX_NAG:R:1",2],
     cex = 3,
     lwd = 3,
     col = alpha('dodgerblue1',1),
     xlab = '', ylab = '', main = '', axes = F,
     xlim = xLim, ylim = yLim)
text(x = neu.umap$layout["1RVX_NAG:R:1",1], y = neu.umap$layout["1RVX_NAG:R:1",2],
     labels = "1RVX_NAG:R:1", pos = 1, offset = 1, col = 'dodgerblue1')
# par(new=T)
# plot(x = neu26.umap["a23_NeuAc",1], y = neu26.umap["a23_NeuAc",2],
#      cex =1,
#      pch=19,
#      col = alpha('dodgerblue1', 1),
#      xlab = '', ylab = '', main = '', axes = F,
#      xlim = xLim, ylim = yLim)

#####
# HA Only UMAP with sequence
#####

allSeqs = readFASTA('data/structures/holo/seqs/nrSeqs/allSeqs.fst')
allSeqs = allSeqs[names(allSeqs) %in% bsResiDat$pdb]

globPsimmat = parSeqSim(allSeqs, cores = 4, type = "global", submat = "BLOSUM62")
row.names(globPsimmat) = names(allSeqs)
colnames(globPsimmat) = names(allSeqs)

plot(density(globPsimmat))

ha.umap = umap(globPsimmat)

summary(ha.umap$layout)

xLim = c(-3.8,7)
yLim = c(-7.5, 10.1)

row.names(ha.umap$layout)

ha23 = row.names(ha.umap$layout) %in% bsResiDat$pdb[tag23]
ha26 = row.names(ha.umap$layout) %in% bsResiDat$pdb[tag26]

haStrains = rep('',nrow(globPsimmat))
for (i in 1:nrow(globPsimmat)){
  haStrains[i] = unique(bsResiDat$strain[bsResiDat$pdb == row.names(globPsimmat)[i]])
}

haTypes = rep('',nrow(globPsimmat))
for (i in 1:nrow(globPsimmat)){
  haTypes[i] = unique(bsResiDat$type[bsResiDat$pdb == colnames(globPsimmat)[i]])
}

# UMAP plot
plot(x = ha.umap$layout[,1], y = ha.umap$layout[,2],
     pch = 19, #col = alpha(colors[labels.int], 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = "UMAP vis for HA sequence similarities",
     xlim = xLim, ylim = yLim)

# Color by ligand
plot(x = ha.umap$layout[ha23,1], y = ha.umap$layout[ha23,2],
     pch = 19, cex = 1.5,
     col = alpha('dodgerblue1', 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = "UMAP vis for HA sequence similarities",
     xlim = xLim, ylim = yLim)
par(new = T)
plot(x = ha.umap$layout[ha26,1], y = ha.umap$layout[ha26,2],
     pch = 19, cex = 1.5,
     # col = alpha(ligColors[2], 0.6),
     col = alpha('firebrick2', 0.6),
     xlab = '', ylab = '', main = '',
     xlim = xLim, ylim = yLim)
legend(x = 'topright', legend = c("3' SA", "6' SA"), col = c('dodgerblue1', 'firebrick2'), pch = c(19, 19))


# Color by subtype & ligand
plot(x = ha.umap$layout[haStrains == strainLst[1] & ha26,1], y = ha.umap$layout[haStrains == strainLst[1] & ha26,2],
     pch = strainPnts[1],
     cex = 2.5,
     col = alpha('firebrick2', 0.7),
     xlab = '', ylab = '', main = "", axes = F,
     xlim = xLim, ylim = yLim)
par(new=T)
plot(x = ha.umap$layout[haStrains == strainLst[1] & ha23,1], y = ha.umap$layout[haStrains == strainLst[1] & ha23,2],
     pch = strainPnts[1],
     cex = 2.5,
     col = alpha('dodgerblue1', 0.7),
     xlab = '', ylab = '', main = "", axes = F,
     xlim = xLim, ylim = yLim)
par(new=T)
plot(x = ha.umap$layout[haStrains == strainLst[1],1], y = ha.umap$layout[haStrains == strainLst[1],2],
     pch = strainPnts[1],
     cex = 1.5,
     col = alpha(strainCols[1], 0.9),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = "HA UMAP vis. for HA sequence similarities",
     xlim = xLim, ylim = yLim)
for (i in 2:length(strainLst)){
  par(new = T)
  plot(x = ha.umap$layout[haStrains == strainLst[i] & ha26,1], y = ha.umap$layout[haStrains == strainLst[i] & ha26,2],
       pch = strainPnts[i],
       cex = 2.5,
       col = alpha('firebrick2', 0.7),
       xlab = '', ylab = '', main = "", axes = F,
       xlim = xLim, ylim = yLim)
  par(new=T)
  plot(x = ha.umap$layout[haStrains == strainLst[i] & ha23,1], y = ha.umap$layout[haStrains == strainLst[i] & ha23,2],
       pch = strainPnts[i],
       cex = 2.5,
       col = alpha('dodgerblue1', 0.7),
       xlab = '', ylab = '', main = "", axes = F,
       xlim = xLim, ylim = yLim)
  par(new=T)
  plot(x = ha.umap$layout[haStrains == strainLst[i],1], y = ha.umap$layout[haStrains == strainLst[i],2],
       pch = strainPnts[i],
       cex = 1.5,
       col = alpha(strainCols[i], 0.9),
       xlab = '', ylab = '', main = "", axes = F,
       xlim = xLim, ylim = yLim)
}
legend(x = 'topright', legend = c(strainLst, '', "6' SA", "3' SA"), col = c(strainCols, 'white', 'firebrick2', 'dodgerblue1'), pch = c(strainPnts, rep(1,3)), pt.cex = 2, bty = 'n')

colnames(globPsimmat) = haStrains
corrplot::corrplot(globPsimmat, order = 'hclust')


globDistMat = globPsimmat*-1 + 1

ha.mds = cmdscale(globDistMat, eig = T, k=2)

plot(ha.mds$points[,1], ha.mds$points[,2])

summary(ha.mds$points)

plot(x = ha.mds$points[ha23,1], y = ha.mds$points[ha23,2],
     pch = 19, cex = 1.5,
     col = alpha('dodgerblue1', 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = "UMAP vis for HA sequence similarities",
     xlim = c(-0.5,0.3), ylim = c(-0.75, 0.15))
par(new = T)
plot(x = ha.mds$points[ha26,1], y = ha.mds$points[ha26,2],
     pch = 19, cex = 1.5,
     # col = alpha(ligColors[2], 0.6),
     col = alpha('firebrick2', 0.6),
     xlab = '', ylab = '', main = '', axes = F,
     xlim = c(-0.5,0.3), ylim = c(-0.75, 0.15))
legend(x = 'topright', legend = c("3' SA", "6' SA"), col = c('dodgerblue1', 'firebrick2'), pch = c(19, 19))


plot(x = ha.mds$points[haStrains == strainLst[1] & ha26,1], y = ha.mds$points[haStrains == strainLst[1] & ha26,2],
     pch = strainPnts[1],
     cex = 2.7,
     col = alpha('firebrick2', 0.7),
     xlab = '', ylab = '', main = "", axes = F,
     xlim = c(-0.5,0.3), ylim = c(-0.75, 0.15))
par(new=T)
plot(x = ha.mds$points[haStrains == strainLst[1] & ha23,1], y = ha.mds$points[haStrains == strainLst[1] & ha23,2],
     pch = strainPnts[1],
     cex = 2.5,
     col = alpha('dodgerblue1', 0.7),
     xlab = '', ylab = '', main = "", axes = F,
     xlim = c(-0.5,0.3), ylim = c(-0.75, 0.15))
par(new=T)
plot(x = ha.mds$points[haStrains == strainLst[1],1], y = ha.mds$points[haStrains == strainLst[1],2],
     pch = strainPnts[1],
     cex = 1.5,
     col = alpha(strainCols[1], 0.9),
     xlab = 'Dim 1', ylab = 'Dim 2', main = "MDS of HA sequence distances",
     xlim = c(-0.5,0.3), ylim = c(-0.75, 0.15))
for (i in 2:length(strainLst)){
  par(new = T)
  plot(x = ha.mds$points[haStrains == strainLst[i] & ha26,1], y = ha.mds$points[haStrains == strainLst[i] & ha26,2],
       pch = strainPnts[i],
       cex = 2.7,
       col = alpha('firebrick2', 0.7),
       xlab = '', ylab = '', main = "", axes = F,
       xlim = c(-0.5,0.3), ylim = c(-0.75, 0.15))
  par(new=T)
  plot(x = ha.mds$points[haStrains == strainLst[i] & ha23,1], y = ha.mds$points[haStrains == strainLst[i] & ha23,2],
       pch = strainPnts[i],
       cex = 2.5,
       col = alpha('dodgerblue1', 0.7),
       xlab = '', ylab = '', main = "", axes = F,
       xlim = c(-0.5,0.3), ylim = c(-0.75, 0.15))
  par(new=T)
  plot(x = ha.mds$points[haStrains == strainLst[i],1], y = ha.mds$points[haStrains == strainLst[i],2],
       pch = strainPnts[i],
       cex = 1.5,
       col = alpha(strainCols[i], 0.9),
       xlab = '', ylab = '', main = "", axes = F,
       xlim = c(-0.5,0.3), ylim = c(-0.75, 0.15))
}
legend(x = 'bottomright', legend = c(strainLst, '', "6' SA", "3' SA"), col = c(strainCols, 'white', 'firebrick2', 'dodgerblue1'), pch = c(strainPnts, rep(1,3)), pt.cex = 2, bty = 'n')

for(i in 1:length(typeLst)){
  tag = haTypes == typeLst[i]
  text(x = mean(ha.mds$points[tag,1]), y = mean(ha.mds$points[tag,2]),
       labels = typeLst[i],
       col = unique(strainCols)[i],
       cex = 2.5,
       pos = 3)
}



# Try MDS with feature representations

scaledFeats = predFeats  # Scale features between 0 & 1
zeroCol = rep(F, ncol(predFeats))
for(i in 1:ncol(scaledFeats)){
  if ((max(scaledFeats[,i]) - min(scaledFeats[,i]) == 0)){
    zeroCol[i] = T
  } else{
    scaledFeats[,i] = (scaledFeats[,i] - min(scaledFeats[,i])) / (max(scaledFeats[,i]) - min(scaledFeats[,i]))
  }
}

tmp = cor(t(scaledFeats[,sigFeats]))
row.names(tmp) = paste(bsResiDat$type, bsResiDat$pdb)
colnames(tmp) = rep('', ncol(tmp))
colnames(tmp)[tag23] = "3'"
colnames(tmp)[tag26] = "6'"

# dev.off()
# corrplot::corrplot(tmp, order = 'hclust', tl.cex = 1)

# Flip direction of correlation and scale between 0,1 for use as "distance" matrix in MDS
tmp = (tmp*-1) + 1
tmp = (tmp - min(tmp)) / (max(tmp) - min(tmp))
plot(density(tmp))

neu.mds = cmdscale(tmp, eig = T, k=2) # MDS of euclidean dist matrix from scaled sig features for each HA interaction with 3'/6' SA

summary(neu.mds$points)

plot(neu.mds$points)

plot(x = neu.mds$points[ha23,1], y = neu.mds$points[ha23,2],
     pch = 19, cex = 1.5,
     col = alpha('dodgerblue1', 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = "UMAP vis for HA sequence similarities",
     xlim = c(-0.4,0.6), ylim = c(-.3, 0.65))
par(new = T)
plot(x = neu.mds$points[ha26,1], y = neu.mds$points[ha26,2],
     pch = 19, cex = 1.5,
     # col = alpha(ligColors[2], 0.6),
     col = alpha('firebrick2', 0.6),
     xlab = '', ylab = '', main = '', axes = F,
     xlim = c(-0.4,0.6), ylim = c(-.3, 0.65))
legend(x = 'bottomright', legend = c("3' SA", "6' SA"), col = c('dodgerblue1', 'firebrick2'), pch = c(19, 19))


plot(x = neu.mds$points[haStrains == strainLst[1] & ha26,1], y = neu.mds$points[haStrains == strainLst[1] & ha26,2],
     pch = strainPnts[1],
     cex = 2.5,
     col = alpha('firebrick2', 0.7),
     xlab = '', ylab = '', main = "", axes = F,
     xlim = c(-0.4,0.6), ylim = c(-.3, 0.65))
par(new=T)
plot(x = neu.mds$points[haStrains == strainLst[1] & ha23,1], y = neu.mds$points[haStrains == strainLst[1] & ha23,2],
     pch = strainPnts[1],
     cex = 2.5,
     col = alpha('dodgerblue1', 0.7),
     xlab = '', ylab = '', main = "", axes = F,
     xlim = c(-0.4,0.6), ylim = c(-.3, 0.65))
par(new=T)
plot(x = neu.mds$points[haStrains == strainLst[1],1], y = neu.mds$points[haStrains == strainLst[1],2],
     pch = strainPnts[1],
     cex = 1.5,
     col = alpha(strainCols[1], 0.9),
     xlab = 'Dim 1', ylab = 'Dim 2', main = "MDS of HA interaction distances (Pearson Cor 28 sig feats)",
     xlim = c(-0.4,0.6), ylim = c(-.3, 0.65))
for (i in 2:length(strainLst)){
  par(new = T)
  plot(x = neu.mds$points[haStrains == strainLst[i] & ha26,1], y = neu.mds$points[haStrains == strainLst[i] & ha26,2],
       pch = strainPnts[i],
       cex = 2.5,
       col = alpha('firebrick2', 0.7),
       xlab = '', ylab = '', main = "", axes = F,
       xlim = c(-0.4,0.6), ylim = c(-.3, 0.65))
  par(new=T)
  plot(x = neu.mds$points[haStrains == strainLst[i] & ha23,1], y = neu.mds$points[haStrains == strainLst[i] & ha23,2],
       pch = strainPnts[i],
       cex = 2.5,
       col = alpha('dodgerblue1', 0.7),
       xlab = '', ylab = '', main = "", axes = F,
       xlim = c(-0.4,0.6), ylim = c(-.3, 0.65))
  par(new=T)
  plot(x = neu.mds$points[haStrains == strainLst[i],1], y = neu.mds$points[haStrains == strainLst[i],2],
       pch = strainPnts[i],
       cex = 1.5,
       col = alpha(strainCols[i], 0.9),
       xlab = '', ylab = '', main = "", axes = F,
       xlim = c(-0.4,0.6), ylim = c(-.3, 0.65))
}
legend(x = 'topright', legend = c(strainLst, '', "6' SA", "3' SA"), col = c(strainCols, 'white', 'firebrick2', 'dodgerblue1'), pch = c(strainPnts, rep(1,3)), pt.cex = 2, bty = 'n')


######
# Sequence vs features
######

# Cluster by sequence
corrplot::corrplot(globPsimmat, order = 'hclust')


annot <- data.frame(Ligand = rep("", nrow(globPsimmat)))
row.names(annot) = row.names(globPsimmat)

colnames(globPsimmat) = row.names(globPsimmat)

annot$Ligand[ha23] = "3' NeuAc"
annot$Ligand[ha26] = "6' NeuAc"

annot$Ligand <- factor(annot$Ligand, levels = unique(annot$Ligand))

annot$Type = haTypes

annot$Type <- factor(annot$Type, levels = typeLst)

Ligand = c('dodgerblue2', 'firebrick1')
names(Ligand) = levels(annot$Ligand)

Type = unique(strainCols)
names(Type) = levels(annot$Type)

annot_cols <- list(Ligand = Ligand, Type = Type)


cWidth = 8
cHeight = 8 

breakLst = seq(0,1,0.01)

# Upper triangular compatible
dev.off()
pdf(file = "analysis/fine_specificity/HA_23-26/wholeSeqSim.pdf",
    height = 9, width = 10)
pheatmap(globPsimmat,
         color = colorRampPalette(c("ivory", "darkslateblue"))(length(breakLst)),
         border_color = 'ivory',
         breaks = breakLst,
         cutree_cols = 6,
         cutree_rows = 6,
         cellwidth = cWidth, cellheight = cHeight,
         labels_col = haStrains,
         annotation_col = annot,
         annotation_colors = annot_cols, 
         clustering_method = 'complete')
 dev.off()

# Lower triangular compatible
pdf(file = "analysis/fine_specificity/HA_23-26/wholeSeqSim_lower.pdf",
    height = 9, width = 10)
pheatmap(globPsimmat,
         color = colorRampPalette(c("ivory", "darkslateblue"))(length(breakLst)),
         border_color = 'ivory',
         breaks = breakLst,
         cutree_cols = 6,
         cutree_rows = 6,
         cellwidth = cWidth, cellheight = cHeight,
         labels_col = paste0(gsub('^Type ', '', haStrains), rep(' (', length(allSeqs)), names(allSeqs), rep(')', length(allSeqs))),
         annotation_row = annot,
         annotation_colors = annot_cols, 
         clustering_method = 'complete')
dev.off()

# Cluster interactions by binding site sequences
bsSeqs = vector(mode = 'list', length = nrow(bsResiDat))
names(bsSeqs) = row.names(bsResiDat)

for (i in 1:nrow(bsResiDat)){
  bsSeqs[[i]] = gsub('\\|', '', bsResiDat$bsiteSequence[i])
}

globBSmat = parSeqSim(bsSeqs, cores = 4, type = "global", submat = "BLOSUM62")
locBSmat = parSeqSim(bsSeqs, cores = 4, type = "local", submat = "BLOSUM62")

colnames(globBSmat) = bsResiDat$type
row.names(globBSmat) = names(bsSeqs)

colnames(locBSmat) = bsResiDat$type
row.names(locBSmat) = names(bsSeqs)

corrplot::corrplot(globBSmat, order = 'hclust', method = 'square')
corrplot::corrplot(locBSmat, order = 'hclust')

colnames(globBSmat) = names(bsSeqs)
colnames(locBSmat) = names(bsSeqs)

breakLst = seq(0,1,0.01)

annot <- data.frame(Ligand = rep("", nrow(globBSmat)))
row.names(annot) = row.names(globBSmat)

annot$Ligand[tag23] = "3' NeuAc"
annot$Ligand[tag26] = "6' NeuAc"

annot$Ligand <- factor(annot$Ligand, levels = unique(annot$Ligand))

annot$Type = bsResiDat$type

annot$Type <- factor(annot$Type, levels = typeLst)

Ligand = c('dodgerblue2', 'firebrick1')
names(Ligand) = levels(annot$Ligand)

Type = unique(strainCols)
names(Type) = levels(annot$Type)

annot_cols <- list(Ligand = Ligand, Type = Type)

pheatmap(globBSmat,
         color = colorRampPalette(c("ivory", "darkslateblue"))(length(breakLst)),
         border_color = 'ivory',
         breaks = breakLst,
         cutree_cols = 7,
         cutree_rows = 7,
         labels_row = bsResiDat$pdb,
         labels_col = bsResiDat$strain,
         annotation_col = annot,
         annotation_colors = annot_cols,
         clustering_method = 'complete')


pheatmap(locBSmat,
         color = colorRampPalette(c("ivory", "darkslateblue"))(length(breakLst)),
         border_color = 'ivory',
         breaks = breakLst,
         cutree_cols = 5,
         cutree_rows = 5,
         labels_row = bsResiDat$pdb,
         labels_col = bsResiDat$strain,
         annotation_col = annot,
         annotation_colors = annot_cols,
         clustering_method = 'complete')


# Cluster interactions by features

scaledFeats = predFeats  # Scale features between 0 & 1
zeroCol = rep(F, ncol(predFeats))
for(i in 1:ncol(scaledFeats)){
  if ((max(scaledFeats[,i]) - min(scaledFeats[,i]) == 0)){
    zeroCol[i] = T
  } else{
    scaledFeats[,i] = (scaledFeats[,i] - min(scaledFeats[,i])) / (max(scaledFeats[,i]) - min(scaledFeats[,i]))
  }
}
cor28Feats = cor(t(scaledFeats[,sigFeats]))



annot <- data.frame(Ligand = rep("", nrow(cor28Feats)))
row.names(annot) = row.names(cor28Feats)

annot$Ligand[tag23] = "3' NeuAc"
annot$Ligand[tag26] = "6' NeuAc"

annot$Ligand <- factor(annot$Ligand, levels = unique(annot$Ligand))

annot$Type = bsResiDat$type

annot$Type <- factor(annot$Type, levels = typeLst)

Ligand = c('dodgerblue2', 'firebrick1')
names(Ligand) = levels(annot$Ligand)

Type = unique(strainCols)
names(Type) = levels(annot$Type)

annot_cols <- list(Ligand = Ligand, Type = Type)

breakLst = seq(-1,1,0.01)

cWidth = 7.5
cHeight = 7.5

dev.off()
pdf(file = "analysis/fine_specificity/HA_23-26/28FeatCor.pdf",
    height =20, width = 20)
pheatmap(cor28Feats,
         color = colorRampPalette(c("gold2", "ivory", "darkslateblue"))(length(breakLst)),
         border_color = 'ivory',
         na_col = 'white',
         breaks = breakLst,
         cutree_cols = 8,
         cutree_rows = 8,
         labels_row = paste0(gsub('^Type ', '', bsResiDat$strain), rep(' (', nrow(bsResiDat)), bsResiDat$pdb, rep(')', nrow(bsResiDat))),
         labels_col = bsResiDat$pdb,
         cellwidth = cWidth, cellheight = cHeight,
         annotation_col = annot,
         annotation_colors = annot_cols, fontsize_row = 10)
         # clustering_method = 'complete')
dev.off()
# 1350 x 1350

# Cluster structures by features

scaledStructFeats = as.data.frame(matrix(0, nrow = length(unique(bsResiDat$pdb)), ncol = sum(sigFeats))) # get feature means for each structure
row.names(scaledStructFeats) = row.names(globPsimmat)
colnames(scaledStructFeats) = colnames(predFeats)[sigFeats]

for (i in 1:nrow(scaledStructFeats)){
  pdb = row.names(scaledStructFeats)[i]
  scaledStructFeats[i,] = apply(predFeats[bsResiDat$pdb == pdb, sigFeats], 2, mean)
}

corrplot::corrplot(cor(t(scaledStructFeats)))

# Scale features between 0 & 1
for(i in 1:ncol(scaledStructFeats)){
  scaledStructFeats[,i] = (scaledStructFeats[,i] - min(scaledStructFeats[,i])) / (max(scaledStructFeats[,i]) - min(scaledStructFeats[,i]))
}

cor28FeatsStruct = cor(t(scaledStructFeats))
corrplot::corrplot(cor28FeatsStruct)

breakLst = seq(-1,1,0.01)

annot <- data.frame(Ligand = rep("", nrow(cor28FeatsStruct)))
row.names(annot) = row.names(cor28FeatsStruct)

annot$Ligand[ha23] = "3' NeuAc"
annot$Ligand[ha26] = "6' NeuAc"

annot$Ligand <- factor(annot$Ligand, levels = unique(annot$Ligand))

annot$Type = haTypes

annot$Type <- factor(annot$Type, levels = typeLst)

Ligand = c('dodgerblue2', 'firebrick1')
names(Ligand) = levels(annot$Ligand)

Type = unique(strainCols)
names(Type) = levels(annot$Type)

annot_cols <- list(Ligand = Ligand, Type = Type)

pheatmap(cor28FeatsStruct,
         color = colorRampPalette(c("ivory", "darkslateblue"))(length(breakLst)),
         border_color = 'ivory',
         breaks = breakLst,
         cutree_cols = 6,
         cutree_rows = 6,
         # labels_row = bsResiDat$pdb,
         labels_col = haStrains,
         annotation_col = annot,
         annotation_colors = annot_cols,
         clustering_method = "complete")



############################
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

# wmwFeats = cbind(predFeats,ligTags)
# 
# des <- svydesign(ids = ~1, data = wmwFeats, weights = 1/cuWeights)
###

# tag23 = matrix(F,nrow = nrow(bsResiDat), ncol = 3)
# row.names(tag23) = row.names(bsResiDat)
# tag26 = tag23
# colnames(tag23) = matched23
# colnames(tag26) = matched26
# 
# for (i in 1:3){
#   tag23[,i] = bsResiDat$iupac == matched23[i]
#   tag26[,i] = bsResiDat$iupac == matched26[i]
# }


tag23 = bsResiDat$iupac %in% ligSort50[neu23Tag]
tag26 = bsResiDat$iupac %in% ligSort50[neu26Tag]


###
# Ligand stats
all(row.names(tag23) == row.names(bsResiDat))
all(row.names(tag26) == row.names(bsResiDat))

# how many clusters are the ligands found in
length(unique(bsResiDat$seqClust50[tag23]))
length(unique(bsResiDat$seqClust50[tag23])) / length(unique(bsResiDat$seqClust50))

length(unique(bsResiDat$seqClust50[tag26]))
length(unique(bsResiDat$seqClust50[tag26])) / length(unique(bsResiDat$seqClust50))

# number of interactions
sum(tag23)
sum(tag26)

# cumulative weight of interactions
for (i in 1:ncol(tag23)){
  cat(colnames(tag23)[i], '\t')
  cat(sum(cuWeights[tag23[,i]]), '\n')
}
sum(cuWeights[tag23])

for (i in 1:ncol(tag26)){
  cat(colnames(tag26)[i], '\t')
  cat(sum(cuWeights[tag26[,i]]), '\n')
}
sum(cuWeights[tag26])


# shared clusters?
sum(unique(bsResiDat$seqClust50[tag23]) %in% unique(bsResiDat$seqClust50[tag26]))
unique(bsResiDat$seqClust50[tag23])[unique(bsResiDat$seqClust50[tag23]) %in% unique(bsResiDat$seqClust50[tag26])]

####
## Limit to shared clusters
####
# 
# for(i in 1:ncol(tag23)){
#   cat(colnames(tag23)[i],'\n')
#   
#   sharedClusts = unique(bsResiDat$seqClust50[tag23[,i]])[unique(bsResiDat$seqClust50[tag23[,i]]) %in% unique(bsResiDat$seqClust50[tag26[,i]])]
#   cat(sharedClusts[order(sharedClusts)],'\n')
#   
#   # number of interactions
#   cat('int cnts\n\t')
#   cat('a2-3: ', sum(tag23[,i] & bsResiDat$seqClust50 %in% sharedClusts), '\n\t')
#   cat('a2-6: ', sum(tag26[,i] & bsResiDat$seqClust50 %in% sharedClusts) , '\n')
#   
#   # cumulative weight of interactions
#   cat('int wts\n\t')
#   cat('a2-3: ', sum(cuWeights[tag23[,i] & bsResiDat$seqClust50 %in% sharedClusts]), '\n\t')
#   cat('a2-6: ', sum(cuWeights[tag26[,i] & bsResiDat$seqClust50 %in% sharedClusts]), '\n\n')
#   
#   # Limit tags to only interactions in shared clusters for each pair
#   tag23[,i] = tag23[,i] & bsResiDat$seqClust50 %in% sharedClusts
#   tag26[,i] = tag26[,i] & bsResiDat$seqClust50 %in% sharedClusts
# }

###
#re-weight within shared clusters
# sharedClusts = unique(bsResiDat$seqClust50[apply(tag23,1,any)])[unique(bsResiDat$seqClust50[apply(tag23,1,any)]) %in% unique(bsResiDat$seqClust50[apply(tag26,1,any)])]

# tag = bsResiDat$seqClust50 %in% sharedClusts & (apply(tag23, 1, any) | apply(tag26, 1, any))
tag = tag23 | tag26

bsResiDat = bsResiDat[tag,]
predFeats = predFeats[tag,]
tag23 = tag23[tag]
tag26 = tag26[tag]


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



# for(i in 1:ncol(tag23)){
#   cat(colnames(tag23)[i],'\n')
#   
#   sharedClusts = unique(bsResiDat$seqClust50[tag23[,i]])[unique(bsResiDat$seqClust50[tag23[,i]]) %in% unique(bsResiDat$seqClust50[tag26[,i]])]
#   cat(sharedClusts[order(sharedClusts)],'\n')
#   
#   # number of interactions
#   cat('int cnts\n\t')
#   cat('a2-3: ', sum(tag23[,i] & bsResiDat$seqClust50 %in% sharedClusts), '\n\t')
#   cat('a2-6: ', sum(tag26[,i] & bsResiDat$seqClust50 %in% sharedClusts) , '\n')
#   
#   # cumulative weight of interactions
#   cat('int wts\n\t')
#   cat('a2-3: ', sum(cuWeights[tag23[,i] & bsResiDat$seqClust50 %in% sharedClusts]), '\n\t')
#   cat('a2-6: ', sum(cuWeights[tag26[,i] & bsResiDat$seqClust50 %in% sharedClusts]), '\n\n')
#   
#   # Limit tags to only interactions in shared clusters for each pair
#   tag23[,i] = tag23[,i] & bsResiDat$seqClust50 %in% sharedClusts
#   tag26[,i] = tag26[,i] & bsResiDat$seqClust50 %in% sharedClusts
# }


# how many clusters are the ligands found in
length(unique(bsResiDat$seqClust50[tag23]))

length(unique(bsResiDat$seqClust50[tag26]))

sum(cuWeights[tag23])
sum(cuWeights[tag26])



# tag23 = apply(tag23, 1, any)
# tag26 = apply(tag26, 1, any)

# neuTags = cbind(tag26, tag23)

# number of interactions
cat('int cnts\n\t')
cat('a2-3: ', sum(tag23), '\n\t')
cat('a2-6: ', sum(tag26) , '\n')

# cumulative weight of interactions
cat('int wts\n\t')
cat('a2-3: ', sum(cuWeights[tag23]), '\n\t')
cat('a2-6: ', sum(cuWeights[tag26]), '\n\n')

###



# # Drop ligands with <=1 interactions in shared clusters?
# apply(tag23, 2, sum)
# apply(tag26, 2, sum)
# tag23 = tag23[,apply(tag23, 2, sum) > 1]
# tag26 = tag26[,apply(tag26, 2, sum) > 1]
# # tag26 = tag26[,apply(tag26, 2, any)] 

neuTags = cbind(tag26, tag23)


#################
# Weighted WMW - Matched NeuAc/Gc pairs one v one (shared clusters)
#################

ligNames = c("6' vs 3' NeuAc", "6' vs 3' sialyated Gal")
ligColors = rep('darkviolet', ncol(neuTags))

# Reformat lignames to be formula-friendly
colnames(neuTags) = gsub('\\(', '.', colnames(neuTags))
colnames(neuTags) = gsub('\\)', '.', colnames(neuTags))
colnames(neuTags) = gsub('\\)', '.', colnames(neuTags))
colnames(neuTags) = gsub('-', '', colnames(neuTags))

# colnames(tag26) = colnames(neuTags)[1:2]
# colnames(tag23) = colnames(neuTags)[3:4]



neuStats = as.data.frame(matrix(0, nrow = ncol(predFeats), ncol = (ncol(neuTags)*4)))
row.names(neuStats) = colnames(predFeats)
colnames(neuStats) = c(paste(colnames(neuTags), 'p', sep = '_'), paste(colnames(neuTags), 'FC', sep = '_'), paste(colnames(neuTags), 'effectSize', sep = '_'), paste(colnames(neuTags), 'adj', sep = '_'))


wmwFeats = cbind(predFeats,neuTags)
des <- svydesign(ids = ~1, data = wmwFeats, weights = 1/cuWeights) # Build survey object with binding indicators for glycans of interest

# for(i in 1:ncol(tag26)){
i = 1

des_allNEU= subset(des, subset = apply(neuTags, 1, any)) # limit to interactions with current NeuGc glycan or matched NeuAc counter-part

des_w = subset(des, subset = neuTags[,i]) # temporary design object holding all interactions with the ligand of interest
des_wo = subset(des, subset = neuTags[,(-1*i)]) # same as above but OTHER the ligands of interest
# scaled_des_w = subset(scaled_des, subset = ligTags[,i]) # 

for(k in 1:ncol(predFeats)){ # for each feature k
  ligTest = svyranktest(formula = as.formula(paste(colnames(predFeats)[k], ' ~ ', colnames(neuTags)[i], sep = '')), # Wilcoxon–Mann–Whitney test, sample sizes can be small (~5% of 230 clusters ~= 10), no reason to assume distribution is normal as it likely isn't
                        design = des_allNEU, 
                        test = 'wilcoxon') 
  
  if (ligTest$p.value != 0){
    neuStats[k, grepl('_p$', colnames(neuStats))][i] = ligTest$p.value # Raw p-value
  } else{
    neuStats[k, grepl('_p$', colnames(neuStats))][i] = 5e-324 # If p-value is too small and R rounds to 0, reset as the smallest positive double as referenced by "?.Machine"
  }
  
  neuStats[k, grepl('_effectSize$', colnames(neuStats))][i] = ligTest$estimate # common language effect size (0.5 transformed to 0)
  
  neuStats[k, grepl('_FC$', colnames(neuStats))][i] = svymean(predFeats[neuTags[,i],k], des_w)[1] / svymean(predFeats[neuTags[,(-1*i)],k], des_wo)[1] # Fold change in weighted means
  
  # ligSpefic_feat_means[k,i] = svymean(scaledFeats[neuTags[,i],k], scaled_des_w)[1] # Get the weighted mean for each feature in interactions with each glycan of interest
}

neuStats[,grepl('_adj$', colnames(neuStats))][,i] = p.adjust(neuStats[grepl('_p$', colnames(neuStats))][,i], method = "BH") # Benjamini-Hochberg MHT correction (FDR)
# }


sum(neuStats[,grepl('_adj$', colnames(neuStats))] < 1e-16)

superSigTag = neuStats[,grepl('_adj$', colnames(neuStats))] < 1e-16
neuStats[,grepl('_adj$', colnames(neuStats))][superSigTag] <- 10**(-1*runif(sum(superSigTag), max = -log10(3e-19), min = -log10(1e-16))) # Sample from a log-uniform distribution


neuStats

#####
## Volcano plots
#####

# par(mfrow=c(1,2))
xLim = c(-0.53,0.53)
yLim = c(0,(-log10(min(neuStats[,grepl('_adj$', colnames(neuStats))])) + 1))
# for(i in 1:ncol(tag26)){
for(i in 1){
  # yLim = c(0,max(-log10(0.1), -log10(min(neuStats[,grepl('_adj$', colnames(neuStats))][,i]))) + 1)
  
  tag = neuStats[,grepl('_adj$', colnames(neuStats))][,i] < 0.01
  
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
  
  plot(neuStats[,grepl('_effectSize$', colnames(neuStats))][,i], -log10(neuStats[,grepl('_adj$', colnames(neuStats))][,i]), # Plot all points w/ color @ alpha 0.5
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
  plot(neuStats[,grepl('_effectSize$', colnames(neuStats))][tag,i], -log10(neuStats[,grepl('_adj$', colnames(neuStats))][tag,i]), # Plot stat sig points again with alpha 1
       pch = 19, col = featColors[tag], cex = 2,
       axes = F, xlab = "", ylab = "", main = "",
       xlim = xLim, ylim = yLim)
  par(new=T)
  plot(neuStats[,grepl('_effectSize$', colnames(neuStats))][tag,i], -log10(neuStats[,grepl('_adj$', colnames(neuStats))][tag,i]), # Outline stat sig points in black to highlight them
       col = alpha('black',0.2), cex = 2.05,
       xlab = "", ylab = "", axes = F, main = "",
       xlim = xLim, ylim = yLim)
}


cor(neuStats[,grepl('_effectSize$', colnames(neuStats))])
dev.off()
plot(neuStats[,grepl('_effectSize$', colnames(neuStats))], col = featColors, pch = 19,
     xlab = ligNames[1], ylab = ligNames[2])
abline(v= 0, h =0)

cor(neuStats$NeuAc.a26.Gal.b14.GlcNAc_effectSize, stats$NeuAc.a2.3.Gal.b1.4.Glc_effectSize)
plot(stats$NeuAc.a2.3.Gal.b1.4.Glc_effectSize, neuStats$NeuAc.a26.Gal.b14.GlcNAc_effectSize)
##########################
# Heatmaps
##########################

# Set-up for 3A (altered from descriptive.R)

neuStats = neuStats[,!grepl('^a23', colnames(neuStats))]

stats = stats[,(grepl('NeuAc', colnames(stats)) | grepl('Sialic', colnames(stats)))]

stats = cbind(stats, neuStats)

sigFeats = neuStats[,grepl('_adj$', colnames(neuStats))] < 0.01

hSigFeats = rbind(rep(F,length(sigFeats)),rep(F,length(sigFeats)),rep(F,length(sigFeats)),sigFeats)

allFeatCorrs = cor(stats[,grepl('_effectSize$', colnames(stats))], stats[,grepl('_effectSize$', colnames(stats))], method = 'pearson')
row.names(allFeatCorrs)  = gsub('_effectSize$', '', row.names(allFeatCorrs))
row.names(allFeatCorrs)[3:4] = c("3' SL", "6' vs 3' NeuAc")
colnames(allFeatCorrs) = row.names(allFeatCorrs)

corrplot::corrplot(allFeatCorrs)

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


CairoPDF(file = paste('./analysis/fine_specificity/23-26/', 
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
         display_numbers = ifelse((hSigFeats)[,resiFeatTag] , "\u2022", ""),
         fontsize_number = 20, number_color = 'black',
         labels_row = c('Sialic acid v all',
                        'NeuAc vs all',
                        'SL vs all',
                        
                        "6' vs 3' NeuAc"),
         annotation_col = resiAnnot,
         annotation_colors = resiAnnot_cols,
         # legend_breaks = c(1,5),
         main = '',
         gaps_row = c(3),
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

CairoPDF(file = paste('./analysis/fine_specificity/23-26/', 
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
         display_numbers = ifelse((hSigFeats)[,pocketFeatTag] , "\u2022", ""),
         fontsize_number = 20, number_color = 'black',
         labels_row = c('Sialic acid v all',
                        'NeuAc vs all',
                        'SL vs all',
                        
                        "6' vs 3' NeuAc"),
         annotation_col = pockAnnot,
         annotation_colors = pockAnnot_cols,
         # legend_breaks = c(1,5),
         main = '',
         gaps_row = c(3),
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


CairoPDF(file = paste('./analysis/fine_specificity/23-26/',
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
         display_numbers = ifelse((hSigFeats)[,1:11] , "\u2022", ""),
         fontsize_number = 20, number_color = 'black',
         labels_row = c('Sialic acid v all',
                        'NeuAc vs all',
                        'SL vs all',
                        
                        "6' vs 3' NeuAc"),
         main = '',
         gaps_row = c(3),
         breaks = breakLst,
         show_colnames = T,
         fontsize_col = 9,
         angle_col = 45)
dev.off()



#####################

cat('cluster\t\ta2-3\t\ta2-6\n')
for(i in 1:length(unique(bsResiDat$seqClust50))){
  clu = unique(bsResiDat$seqClust50)[order(unique(bsResiDat$seqClust50))][i]
  cat(clu, '\t\t')
  
  cat(sum(cuWeights[tag23 & bsResiDat$seqClust50 == clu]), '\t\t')
  cat(sum(cuWeights[tag26 & bsResiDat$seqClust50 == clu]), '\n')
  cat('\t', unique(bsResiDat$origine[bsResiDat$seqClust50 == clu]),'\n')
  cat('\t', unique(bsResiDat$espece[bsResiDat$seqClust50 == clu]),'\n\n')
}









