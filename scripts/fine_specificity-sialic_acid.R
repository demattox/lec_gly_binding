

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

ligSort50[mTag]
ligSort50[dTag]
ligSort50[tTag]
ligSort50[qTag]
ligSort50[pTag]
 
ligSort50[bTag]

ligSort50[manTag]
for(i in 1:length(ligSort50[manTag])){
  cat(ligSort50[manTag][i], '\n')
}
ligSort50[neuTag]
for(i in 1:length(ligSort50[neuTag])){
  cat(ligSort50[neuTag][i], '\n')
}
ligSort50[fucTag]
for(i in 1:length(ligSort50[fucTag])){
  cat(ligSort50[fucTag][i], '\n')
}

ligSort50[grepl('Kdo',ligSort50)]
ligSort50[grepl('Kdn',ligSort50)]


# NeuGc-containing glycans
ngnaTag = grepl('NeuGc',ligSort50) # Has terminal sialic acid and is not a monosacc.
ligSort50[ngnaTag]
cpl50[ngnaTag]

gcTags = as.data.frame(matrix(F, nrow = nrow(bsResiDat), ncol = sum(ngnaTag)))
row.names(gcTags) = row.names(bsResiDat)
colnames(gcTags) = ligSort50[ngnaTag]

for(i in 1:length(ligSort50[ngnaTag])){
  tag = bsResiDat$iupac == ligSort50[ngnaTag][i]
  gcTags[,i] = tag
  cat(ligSort50[ngnaTag][i], '\t')
  cat(sum(tag), '\t')
  gcClusts = unique(bsResiDat$seqClust50[tag])
  cat(gcClusts, '\n')
  
  neugly = gsub('NeuGc', 'NeuAc', ligSort50[ngnaTag][i])
  cat(neugly, '\t', sum(bsResiDat$iupac == neugly), '\t', length(unique(bsResiDat$seqClust50[bsResiDat$iupac == neugly])), '\t', unique(bsResiDat$seqClust50[bsResiDat$iupac == neugly])[order(unique(bsResiDat$seqClust50[bsResiDat$iupac == neugly]))])
  cat('\n\n')
}

# a2-3 vs a2-6 matched glycans
neu23Tag = grepl('^NeuAc\\(a2-3\\)',ligSort50)
ligSort50[neu23Tag]
neu26Tag = grepl('^NeuAc\\(a2-6\\)',ligSort50)
ligSort50[neu26Tag]

sum(gsub('^NeuAc\\(a2-3\\)', 'NeuAc(a2-6)', ligSort50[neu23Tag]) %in% ligSort50[neu26Tag])
ligSort50[neu23Tag][gsub('^NeuAc\\(a2-3\\)', 'NeuAc(a2-6)', ligSort50[neu23Tag]) %in% ligSort50[neu26Tag]]


########################
## NeuAc/Gc Monosaccharides - Descriptive
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

# wmwFeats = cbind(predFeats,ligTags)
# 
# des <- svydesign(ids = ~1, data = wmwFeats, weights = 1/cuWeights)
###


ligNames = c('Terminal NeuAc', 'Terminal NeuGc')
ligColors = c('darkviolet', 'deepskyblue ')

neuTags = cbind(gcTags, ligTags[,8], bsResiDat$iupac == 'NeuAc(a2-3)Gal(b1-3)GalNAc', ligTags[,9], bsResiDat$iupac == 'NeuAc(a2-3)Gal(b1-3)GlcNAc')
colnames(neuTags)[5:8] = c('NeuAc',
                           'NeuAc(a2-3)Gal(b1-3)GalNAc',
                           'NeuAc(a2-3)Gal(b1-3)Glc',
                           'NeuAc(a2-3)Gal(b1-3)GlcNAc')

###
# Ligand stats
all(row.names(neuTags) == row.names(bsResiDat))

# how many clusters are the ligands found in
length(unique(bsResiDat$seqClust50[apply(neuTags[,1:4], 1, any)]))
length(unique(bsResiDat$seqClust50[apply(neuTags[,1:4], 1, any)])) / length(unique(bsResiDat$seqClust50))

length(unique(bsResiDat$seqClust50[apply(neuTags[,5:8], 1, any)]))
length(unique(bsResiDat$seqClust50[apply(neuTags[,5:8], 1, any)])) / length(unique(bsResiDat$seqClust50))

# number of interactions
apply(neuTags, 2, sum)
sum(apply(neuTags, 2, sum)[1:4])
sum(apply(neuTags, 2, sum)[5:8])

# cumulative weight of interactions
for (i in 1:ncol(neuTags)){
  cat(colnames(neuTags)[i], '\t')
  cat(sum(cuWeights[neuTags[,i]]), '\n')
}
sum(cuWeights[apply(neuTags[,1:4], 1, any)])
sum(cuWeights[apply(neuTags[,5:8], 1, any)])


# shared clusters?
sum(unique(bsResiDat$seqClust50[apply(neuTags[,1:4], 1, any)]) %in% unique(bsResiDat$seqClust50[apply(neuTags[,5:8], 1, any)]))

#################
# # Weighted WMW - NeuAc vs NeuGc monosaccharides
# #################
# 
# wmwFeats = cbind(predFeats,neuTags)
# des <- svydesign(ids = ~1, data = wmwFeats, weights = 1/cuWeights) # Build survey object with 
# 
# neuTags = neuTags[,c(1,5)]
# 
# neuStats = as.data.frame(matrix(0, nrow = ncol(predFeats), ncol = (2*4)))
# row.names(neuStats) = colnames(predFeats)
# colnames(neuStats) = c(paste(colnames(neuTags), 'p', sep = '_'), paste(colnames(neuTags), 'FC', sep = '_'), paste(colnames(neuTags), 'effectSize', sep = '_'), paste(colnames(neuTags), 'adj', sep = '_'))
# 
# des_allNEU= subset(des, subset = (apply(neuTags, 1, any)))
# 
# 
# for (i in 1:ncol(neuTags)){
#   
#   des_w = subset(des, subset = neuTags[,i]) # temporary design object holding all interactions with the ligand of interest
#   des_wo = subset(des, subset = neuTags[,-i]) # same as above but OTHER the ligands of interest
#   # scaled_des_w = subset(scaled_des, subset = ligTags[,i]) # 
#   
#   for(k in 1:ncol(predFeats)){ # for each feature k
#     ligTest = svyranktest(formula = as.formula(paste(colnames(predFeats)[k], ' ~ ', colnames(neuTags)[i], sep = '')), # Wilcoxon–Mann–Whitney test, sample sizes can be small (~5% of 230 clusters ~= 10), no reason to assume distribution is normal as it likely isn't
#                           design = des_allNEU, 
#                           test = 'wilcoxon') 
#     
#     if (ligTest$p.value != 0){
#       neuStats[k, grepl('_p$', colnames(neuStats))][i] = ligTest$p.value # Raw p-value
#     } else{
#       neuStats[k, grepl('_p$', colnames(neuStats))][i] = 5e-324 # If p-value is too small and R rounds to 0, reset as the smallest positive double as referenced by "?.Machine"
#     }
#     
#     neuStats[k, grepl('_effectSize$', colnames(neuStats))][i] = ligTest$estimate # common language effect size (0.5 transformed to 0)
#     
#     neuStats[k, grepl('_FC$', colnames(neuStats))][i] = svymean(predFeats[neuTags[,i],k], des_w)[1] / svymean(predFeats[neuTags[,-i],k], des_wo)[1] # Fold change in weighted means
#     
#     # ligSpefic_feat_means[k,i] = svymean(scaledFeats[neuTags[,i],k], scaled_des_w)[1] # Get the weighted mean for each feature in interactions with each glycan of interest
#   }
#   
#   neuStats[,grepl('_adj$', colnames(neuStats))][,i] = p.adjust(neuStats[grepl('_p$', colnames(neuStats))][,i], method = "BH") # Benjamini-Hochberg MHT correction (FDR)
# }
# 
# all(round(neuStats$NeuAc_adj,3) == round(neuStats$NeuGc_adj,3))
# 
# # superSigTag = neuStats[,grepl('_adj$', colnames(neuStats))] < 1e-16
# # neuStats[,grepl('_adj$', colnames(neuStats))][superSigTag] <- 10**(-1*runif(sum(superSigTag), max = -log10(3e-19), min = -log10(1e-16))) # Sample from a log-uniform distribution
# 
# 
# #####
# ## Volcano plots
# #####
# 
# par(mfrow=c(1,2))
# 
# xLim = c(-0.5,0.5)
# yLim = c(0,(-log10(min(neuStats[,grepl('_adj$', colnames(neuStats))])) + 1))
# for(i in 1:ncol(neuTags)){
#   
#   # yLim = c(0,max(-log10(0.1), -log10(min(neuStats[,grepl('_adj$', colnames(neuStats))][,i]))) + 1)
#   
#   tag = neuStats[,grepl('_adj$', colnames(neuStats))][,i] < 0.01
#   
#   # dev.off()
#   plot(0,0,axes = F, main = '', xlab = '', ylab = '', pch = NA)
#   bg = "seashell2"
#   fg = "ivory"
#   rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = bg)
#   # abline(v = c(-1,-.5,0,.5,1), lwd = 6, col = fg)
#   # abline(v = c(-1.25,-.75,-.25,.25,.75,1.25), lwd = 3, col = fg)
#   # abline(h = c(0,1,2,3,4,5,6), lwd = 6, col = fg)
#   # abline(h = c(0.5,1.5,2.5,3.5,4.5,5.5,6.5), lwd = 3, col = fg)
#   
#   abline(v = 0, lty=2, lwd = 4, col = 'white')
#   
#   par(new=T)
#   
#   plot(neuStats[,grepl('_effectSize$', colnames(neuStats))][,i], -log10(neuStats[,grepl('_adj$', colnames(neuStats))][,i]), # Plot all points w/ color @ alpha 0.5
#        xlab = "Effect size", ylab = "-log10(FDR)", main = '',
#        pch = 19, cex = 2, col = alpha(featColors, 0.33),
#        cex.axis = 1.5, cex.lab = 1.5,
#        xlim = xLim, ylim = yLim)
#   
#   abline(h= -log10(0.01), lwd = 4, col = 'white')
#   abline(h = -log10(1e-16), lwd = 1.5, lty = 2, col = 'white')
#   
#   title(main = ligNames[i], col.main = ligColors[i], cex.main = 1.8, font.main  = 2)
#   # mtext(ligNames[i], col = ligColors[i],
#   #       side=3, adj=0, outer = T,
#   #       line=1.2, cex=1, font=2)
#   
#   par(new=T)
#   plot(neuStats[,grepl('_effectSize$', colnames(neuStats))][tag,i], -log10(neuStats[,grepl('_adj$', colnames(neuStats))][tag,i]), # Plot stat sig points again with alpha 1
#        pch = 19, col = featColors[tag], cex = 2,
#        axes = F, xlab = "", ylab = "", main = "",
#        xlim = xLim, ylim = yLim)
#   par(new=T)
#   plot(neuStats[,grepl('_effectSize$', colnames(neuStats))][tag,i], -log10(neuStats[,grepl('_adj$', colnames(neuStats))][tag,i]), # Outline stat sig points in black to highlight them
#        col = alpha('black',0.2), cex = 2.05,
#        xlab = "", ylab = "", axes = F, main = "",
#        xlim = xLim, ylim = yLim)
# }

#compare to global trends

plot(stats$NeuAc_effectSize, neuStats$NeuAc_effectSize,
     col = featColors, pch =19,
     main = 'NeuAc',
     xlab = 'Effect Size vs background', ylab = 'Effect Size vs NeuGc')
abline(v=0, h=0)
par(new=T)
plot(stats$NeuAc_effectSize[neuStats$NeuAc_adj <= 0.01], neuStats$NeuAc_effectSize[neuStats$NeuAc_adj <= 0.01],
     col = 'black', pch =19,
     cex = 1.7,
     main = '', xlab = '', ylab = '', axes = F)
par(new=T)
plot(stats$NeuAc_effectSize[neuStats$NeuAc_adj <= 0.01], neuStats$NeuAc_effectSize[neuStats$NeuAc_adj <= 0.01],
     col = featColors[neuStats$NeuAc_adj <= 0.01], pch =19,
     cex = 1.5,
     main = '', xlab = '', ylab = '', axes = F)
cor.test(stats$NeuAc_effectSize, neuStats$NeuAc_effectSize,)

#################
# Weighted WMW - NeuAc vs NeuGc groups
#################

neuTags$Term_NeuGc = apply(neuTags[1:4],1, any)
neuTags$Term_NeuAc = apply(neuTags[5:8],1, any)
neuTags = neuTags[,c(10,9)]

wmwFeats = cbind(predFeats,neuTags)
des <- svydesign(ids = ~1, data = wmwFeats, weights = 1/cuWeights) # Build survey object with binding indicators for glycans of interest

neuStats = as.data.frame(matrix(0, nrow = ncol(predFeats), ncol = (2*4)))
row.names(neuStats) = colnames(predFeats)
colnames(neuStats) = c(paste(colnames(neuTags), 'p', sep = '_'), paste(colnames(neuTags), 'FC', sep = '_'), paste(colnames(neuTags), 'effectSize', sep = '_'), paste(colnames(neuTags), 'adj', sep = '_'))

des_allNEU= subset(des, subset = (apply(neuTags, 1, any)))


for (i in 1:ncol(neuTags)){
  
  des_w = subset(des, subset = neuTags[,i]) # temporary design object holding all interactions with the ligand of interest
  des_wo = subset(des, subset = neuTags[,-i]) # same as above but OTHER the ligands of interest
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
    
    neuStats[k, grepl('_FC$', colnames(neuStats))][i] = svymean(predFeats[neuTags[,i],k], des_w)[1] / svymean(predFeats[neuTags[,-i],k], des_wo)[1] # Fold change in weighted means
    
    # ligSpefic_feat_means[k,i] = svymean(scaledFeats[neuTags[,i],k], scaled_des_w)[1] # Get the weighted mean for each feature in interactions with each glycan of interest
  }
  
  neuStats[,grepl('_adj$', colnames(neuStats))][,i] = p.adjust(neuStats[grepl('_p$', colnames(neuStats))][,i], method = "BH") # Benjamini-Hochberg MHT correction (FDR)
}

all(round(neuStats$Term_NeuAc_adj,3) == round(neuStats$Term_NeuGc_adj,3))

superSigTag = neuStats[,grepl('_adj$', colnames(neuStats))] < 1e-16
neuStats[,grepl('_adj$', colnames(neuStats))][superSigTag] <- 10**(-1*runif(sum(superSigTag), max = -log10(3e-19), min = -log10(1e-16))) # Sample from a log-uniform distribution


#####
## Volcano plots
#####

par(mfrow=c(1,2))

xLim = c(-0.5,0.5)
yLim = c(0,(-log10(min(neuStats[,grepl('_adj$', colnames(neuStats))])) + 1))
for(i in 1:ncol(neuTags)){
  
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

##########################
# Heatmaps
##########################

# Set-up for 3A (altered from descriptive.R)

stats = stats[,(grepl('NeuAc', colnames(stats)) | grepl('Sialic', colnames(stats)))]


stats = cbind(stats, neuStats)

hImpFeats = rbind(rep(F,length(topImpfeats)),rep(F,length(topImpfeats)),rep(F,length(topImpfeats)),topImpfeats,topImpfeats)
hSigFeats = rbind(rep(F,length(sigFeats)),rep(F,length(sigFeats)),rep(F,length(topImpfeats)),sigFeats,sigFeats)

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
         # labels_row = c('Lac vs all', 'LacNAc vs all', 'Lac vs LacNAc', 'LacNAc vs Lac'),
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

