
library(reshape)
library(ggplot2)
library(philentropy)
library(pheatmap)
library(PRROC)
library(vioplot)

###############
# Functions
###############
colfunc = colorRampPalette(c("red","goldenrod","forestgreen","royalblue","darkviolet"))

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


###############
# read in files
###############
homeDir = '/Users/dmattox/cbk/lec_gly_binding/'
setwd(homeDir)

ligTags = read.delim(file = './analysis/training/data_in/ligTags.tsv', sep = '\t', stringsAsFactors = F)
predFeats = read.delim(file = './analysis/training/data_in/predFeats.csv', sep = ',', stringsAsFactors = F)
bsResiDat = read.delim(file = './analysis/training/data_in/bsResiDat.tsv', sep = '\t', stringsAsFactors = F)

ligColors = colfunc(ncol(ligTags))

# inDir = './analysis/training/train_and_validate/seqID50/'
# inDir = './analysis/training/train_and_validate/seqID80/'
# inDir = './analysis/training/train_and_validate/beta3id50/'

predDirs = './analysis/training/train_and_validate/nr_CVLOCO_reps/pred/'
randDirs = './analysis/training/train_and_validate/nr_CVLOCO_reps/rand/'











# inFiles = dir(inDir)
# inTrain = inFiles[grepl('training.csv', inFiles)]
# inTest = inFiles[grepl('testing.csv', inFiles)]
# inOutcomes = inFiles[grepl('outcomes.csv', inFiles)]
# inFeats = inFiles[grepl('features.csv', inFiles)]
# 
# newOrd = match(colnames(ligTags), gsub('(.*)_outcomes.csv', '\\1', inOutcomes))
# inTrain = inTrain[newOrd]
# inTest = inTest[newOrd]
# inOutcomes = inOutcomes[newOrd]
# inFeats = inFeats[newOrd]

###############
# Look at lectin specificity
###############
lecIDs = unique(bsResiDat$uniparc)

iupacLigCnt = rep(0, length(lecIDs))

for (i in 1:length(lecIDs)){
  id = lecIDs[i]
  iupacLigCnt[i] = length(unique(bsResiDat$iupac[bsResiDat$uniparc == id]))
}

# Get ligand tags again
uniLigs = unique(bsResiDat$iupac)

# parenCnt = bracCnt = manCnt = neuCnt = bracCnt = rep(0,length(uniLigs))
# for (i in 1:length(uniLigs)){
#   lig = uniLigs[i]
#   parenCnt[i] = lengths(regmatches(lig, gregexpr("\\(", lig)))
#   bracCnt[i] = lengths(regmatches(lig, gregexpr("\\[", lig)))
#   manCnt[i] = lengths(regmatches(lig, gregexpr("Man", lig)))
#   neuCnt[i] = lengths(regmatches(lig, gregexpr("NeuAc", lig)))
# }
# 
# mTag = parenCnt == 0 & bracCnt == 0 # Monosaccharides
# dTag = parenCnt == 1 & bracCnt == 0 # Disaccharides
# tTag = (parenCnt == 2 & bracCnt == 0) | (parenCnt == 2 & bracCnt == 1) # Trisaccharides
# qTag = (parenCnt == 3 & bracCnt == 0) | (parenCnt == 3 & bracCnt == 1) | (parenCnt == 1 & bracCnt == 2) # Tetrasaccharides
# pTag = !(mTag | dTag | tTag | qTag) # 5+ sugars
# bTag = bracCnt >= 1 # Branched glycans
# 
# manTag = manCnt > 3 # High mannose
# # uniLigs[manTag]
# neuTag = neuCnt >= 1 # Has sialic acid
# # uniLigs[neuTag]
# fucTag = grepl('^Fuc',uniLigs) # Has a terminal fucose
# # uniLigs[fucTag]

# Look at lectin specificity
lecSpec = as.data.frame(matrix(nrow = length(lecIDs), ncol = length(uniLigs)))
row.names(lecSpec) = lecIDs
colnames(lecSpec) = uniLigs
for (i in 1:length(lecIDs)){
  lecSpec[i,] = uniLigs %in% unique((bsResiDat$iupac[bsResiDat$uniparc == lecIDs[i]]))
}


###############
# Training performance
###############

# sampSizes = rep(0, ncol(ligTags))

for(i in 1:length(dir(predDirs))){
  lig = colnames(ligTags)[as.numeric(dir(predDirs)[i])]
  cat(lig,'\n')
  curDir = paste(predDirs, dir(predDirs)[i], sep = '')
  subDirs = dir(curDir)
  for (j in 1:length(subDirs)){
    cat(j,'\n')
    inFiles = dir(paste(curDir, subDirs[j], sep = '/'))
    
    readFile = inFiles[grepl('training.csv', inFiles)] # change here to read other files
    
    tmp = read.delim(file = paste(curDir, subDirs[j], readFile, sep = '/'), header = T, sep = ',', stringsAsFactors = F)
    tmp$mode = 'pred'
    if (! exists('train')){
      train = tmp
    }else{
      train = rbind(train, tmp)
    }
  }
}

for(i in 1:length(dir(randDirs))){
  lig = colnames(ligTags)[i]
  cat(lig,'\n')
  curDir = paste(randDirs, dir(randDirs)[i], sep = '')
  subDirs = dir(curDir)
  for (j in 1:length(subDirs)){
    cat(j,'\n')
    inFiles = dir(paste(curDir, subDirs[j], sep = '/'))
    
    readFile = inFiles[grepl('training.csv', inFiles)] # change here to read other files
    
    tmp = read.delim(file = paste(curDir, subDirs[j], trainFile, sep = '/'), header = T, sep = ',', stringsAsFactors = F)
    tmp$mode = 'rand'
    if (! exists('train')){
      train = tmp
    }else{
      train = rbind(train, tmp)
    }
  }
}
  


train$f2 = 0
train$prec = 0
for(i in 1:nrow(train)){
  r = train$recall[i]
  p = train$TP[i] / (train$TP[i] + train$FP[i])
  train$prec[i] = p
  train$f2[i] = ((1+(2^2)) * p * r) / (2^2 * p + r)
  # trainOut$f3[i] = ((1+(3^2)) * p * r) / (3^2 * p + r)
  # trainOut$f4[i] = ((1+(4^2)) * p * r) / (4^2 * p + r)
}

mTrain = melt(train, id.vars = 'ligand', measure.vars = c('kappa', 'recall', 'prec'))
colnames(mTrain) = c('ligand', 'metric', 'value')



# ggplot(data = mTrain, aes(x = metric, y = value, col = ligand)) +
#   geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.4) +
#   geom_boxplot(outlier.alpha = 0) +
#   ylim(c(0.5, 1)) +
#   scale_fill_manual(values = alpha(c('navajowhite','navajowhite'), 0.6), guide =F) +
#   scale_color_manual(values = ligColors) +
#   labs(title = '5x CV with LOCO validation - Training Performance', x = "Metric type", y = "Metric value") +
#   theme_dark(base_size = 22) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))
# 
# ggplot(data = mTrain, aes(x = ligand, y = value, col = metric)) +
#   geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.4) +
#   geom_boxplot(outlier.alpha = 0) +
#   ylim(c(0.5, 1)) +
#   scale_fill_manual(values = alpha(c('navajowhite','navajowhite'), 0.6), guide =F) +
#   scale_color_manual(values = c('darkgreen', 'firebrick1', 'darkorchid1', 'dodgerblue')) +
#   labs(title = '5x CV with LOCO validation - Training Performance', x = "Ligand for prediciton", y = "Metric value") +
#   theme_dark(base_size = 22) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))


# Light colors
ggplot(data = mTrain, aes(x = metric, y = value, col = ligand, fill = metric)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.4) +
  geom_boxplot(outlier.alpha = 0) +
  ylim(c(0.5, 1)) +
  scale_fill_manual(values = alpha(rep('snow3',length(unique(mTrain$metric))), 0.6), guide =F) +
  scale_color_manual(values = ligColors) +
  labs(title = '5x CV with LOCO validation - Training Performance', x = "Metric type", y = "Metric value") +
  theme_light(base_size = 22) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))



ggplot(data = mTrain, aes(x = ligand, y = value, col = metric, fill = metric)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.4) +
  geom_boxplot(outlier.alpha = 0) +
  ylim(c(0.5, 1)) +
  scale_fill_manual(values = alpha(rep('snow3',length(unique(mTrain$metric))), 0.6), guide =F) +
  scale_color_manual(values = c('darkgreen', 'firebrick1', 'darkorchid1', 'dodgerblue')) +
  labs(title = '5x CV with LOCO validation - Training Performance', x = "Ligand for prediciton", y = "Metric value") +
  theme_light(base_size = 22) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))


###############
# Validation performance
###############
test = as.data.frame(matrix(0,nrow= ncol(ligTags)*10, ncol = 5))
colnames(test) = c('ligand', 'kappa', 'recall', 'f2', 'prec')
# colnames(test) = c('kappa', 'recall', 'f2', 'precision')

for(i in 1:length(inTest)){
  lig = gsub('(.*)_testing.csv', '\\1', inTest)[i]

  tmp = read.delim(file = paste(inDir, inTest[i], sep = ''), header = T, sep = ',', stringsAsFactors = F)
  
  for(j in 1:10){
    tag = grepl(paste('_', as.character(j), '$', sep = ''), row.names(tmp))
    repDat = tmp[tag,]
    test[(j+(10*(i-1))),] = c(lig, as.numeric(getKPRFb(repDat))) # Kapppa, recall, f2, precision
  }
  # R =  repDat$TP/(repDat$TP + repDat$FN)
  # R = as.data.frame(cbind(rep(lig, nrow(repDat)), R))
  # colnames(R) = c('ligand', 'recall')
  # if (i ==1){
  #   recall = R
  # } else {
  #   recall = rbind(recall, R)
  # }
  # test[lig,'ligand'] = lig
}

test[,2:5] = apply(test[,2:5], 2, as.numeric)
test$ligand = factor(test$ligand, levels = levels(mTrain$ligand))

mTest= melt(test, id.vars = 'ligand', measure.vars = c('kappa', 'recall', 'prec'))
colnames(mTest) = c('ligand', 'metric', 'value')


ggplot(data = mTest, aes(x = metric, y = value, col = ligand, fill = metric)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.4) +
  geom_boxplot(outlier.alpha = 0) +
  ylim(c(-0.15, 1)) +
  scale_fill_manual(values = alpha(rep('snow3',length(unique(mTest$metric))), 0.6), guide =F) +
  scale_color_manual(values = ligColors) +
  labs(title = '5x CV with LOCO validation - Validation Performance', x = "Metric type", y = "Metric value") +
  theme_light(base_size = 22) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))

ggplot(data = mTest, aes(x = ligand, y = value, col = metric, fill = metric)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.4) +
  geom_boxplot(outlier.alpha = 0) +
  ylim(c(0, 1)) +
  scale_fill_manual(values = alpha(rep('snow3',length(unique(mTest$metric))), 0.6), guide =F) +
  scale_color_manual(values = c('darkgreen', 'firebrick1', 'darkorchid1', 'dodgerblue')) +
  labs(title = '5x CV with LOCO validation - Validation Performance', x = "Ligand for prediciton", y = "Metric value") +
  theme_light(base_size = 22) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))

# recall$recall = as.numeric(as.character(recall$recall))

# ggplot(data = recall, aes(x = ligand, y = recall, col = ligand, fill = ligand)) +
#   geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.4) +
#   geom_boxplot(outlier.alpha = 0) +
#   ylim(c(0, 1)) +
#   scale_fill_manual(values = alpha(rep('snow3',ncol(ligTags)), 0.6), guide =F) +
#   scale_color_manual(values = ligColors, guide = F) +
#   labs(title = 'Validation Recall for each excluded cluster', x = "Ligand for prediciton", y = "Cluster-specific Recall value") +
#   theme_light(base_size = 22) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))




breakLst = seq(0,1,0.01)
mycol = colorRampPalette(c("ivory", "cornflowerblue", "navy"))(n = length(breakLst))

testMeds = as.data.frame(matrix(0,nrow  = length(unique(test$ligand)), ncol = 4))
row.names(testMeds) = gsub('(.*)_training.csv', '\\1', inTrain)
colnames(testMeds) = colnames(test)[-1]

for(i in 1:length(unique(test$ligand))){
  tag = grepl(unique(test$ligand)[i], test$ligand)
  testMeds[i,] = apply(test[tag,2:5], 2, median)
}

testMeds$f2 <- NULL

pheatmap(testMeds, color =  mycol,
         cluster_cols = F, cluster_rows = F,
         display_numbers = T, number_color = 'darkgoldenrod', fontsize_number = 28,
         breaks = breakLst, main = 'Median validation metrics for each ligand class', 
         gaps_row = c(3), gaps_col = c(1))



plot(sampSizes, testMeds$kappa, pch = 19, cex = 1.5,
     ylim = c(0, 0.5),
     xlim = c(0,550),
     xlab = 'Average number of samples sites for training',
     ylab = 'Median Validation Kappa',
     cex.lab=1.5, cex.axis=1.5)

cor.test(sampSizes, testMeds$kappa)


# PR curves
par(mfrow = c(3,5))
for(i in 1:ncol(ligTags)){
  lig = gsub('(.*)_outcomes.csv', '\\1', inOutcomes)[i]
  outcomes = read.delim(file = paste(inDir, inOutcomes[i], sep = ''), header = T, sep = ',', stringsAsFactors = F)
  
  boundLigs = unique(bsResiDat$iupac[ligTags[,i]]) # UniProt IDs of all lectins that do have any examples of binding
  boundIDs = unique(bsResiDat$uniparc[bsResiDat$iupac %in% boundLigs])
  fpBS = row.names(outcomes)[outcomes$Pred >= 0.5 & outcomes$Obs == F] # Binding sites that ended up as false positives
  fpLecs = unique(bsResiDat$uniparc[row.names(bsResiDat) %in% fpBS]) # unique lectin uniprot ids of lectins that show up as false positives
  fpIDs = bsResiDat$uniparc[row.names(bsResiDat) %in% fpBS] # UniProt ID for FP binding site
  
  negBS = row.names(outcomes)[outcomes$Obs == F]
  negIDs = bsResiDat$uniparc[row.names(bsResiDat) %in% negBS] # UniProt IDs of lectins with negative binding sites that DO bind ligand in other structures
  bound_NegBS = negBS[negIDs %in% boundIDs] # Negative binding sites from 

  pr = pr.curve(outcomes$Pred[outcomes$Obs == T], outcomes$Pred[outcomes$Obs == F], curve= T, rand.compute = T)
  
  R = test$recall[i]
  P = test$precision[i]
  
  plot(pr$curve[,1:2], type = 'l', lwd = 4, col = ligColors[i],
       xlim = c(0,1), ylim = c(0,1),
       xlab = 'Recall', ylab = 'Precision', main = paste('PR Curve - ', lig, '\nAUC = ', as.character(round(pr$auc.integral, digits = 5)), sep = ''))
  abline(h = pr$rand$auc.integral, lty = 1, lwd = 2)
  lines(x = c(R,R), y = c(-1,P), lty = 2)
  lines(x = c(-1,R), y = c(P,P), lty = 2)
  points(R,P, pch = 19)
  tag = pr$curve[,3] %in% outcomes[bound_NegBS, "Pred"]
  rug((pr$curve[tag,1][pr$curve[tag,1] > R]), col = alpha('black',0.7))
  rug((pr$curve[tag,1][pr$curve[tag,1] <= R]), col = alpha('red',0.7))
  rug((pr$curve[tag,2][pr$curve[tag,2] < P]), col = alpha('black',0.7), side = 2)
  rug((pr$curve[tag,2][pr$curve[tag,2] >= P]), col = alpha('red',0.7), side = 2)
  text(x= 0.7, y = 0.9, labels = paste(round(100 * sum(fpIDs %in% boundIDs)/length(fpIDs),2), '% of ', length(fpIDs), ' FPs', sep = ''), col = 'red', cex = 1.2)
}

# Feature importance
for(i in 1:length(inFeats)){
  lig = gsub('(.*)_features.csv', '\\1', inFeats)[i]
  tmp = read.delim(file = paste(inDir, inFeats[i], sep = ''), header = T, sep = ',', stringsAsFactors = F)
  tmp = apply(tmp, 1, mean)
  if (i == 1){
    feats = tmp
  } else {
    feats = cbind(feats, tmp)
    colnames(feats)[i] = lig
  }
}
colnames(feats)[1] = gsub('(.*)_features.csv', '\\1', inFeats)[1]
feats = as.data.frame(feats)

neuFeats = c('negCharge_bin1','ASP_bin1', 'majorMaxima_8Ang')

neuInds = (1:nrow(feats))[row.names(feats)[order(feats$Sialic_Acid, decreasing = T)] %in% neuFeats]
row.names(feats)[order(feats$Sialic_Acid, decreasing = T)][neuInds]
100 - (neuInds/nrow(feats))*100

neuInds = (1:nrow(feats))[row.names(feats)[order(feats$NeuAc, decreasing = T)] %in% neuFeats]
row.names(feats)[order(feats$NeuAc, decreasing = T)][neuInds]
100 - (neuInds/nrow(feats))*100

neuInds = (1:nrow(feats))[row.names(feats)[order(feats$NeuAc.a2.3.Gal.b1.4.Glc, decreasing = T)] %in% neuFeats]
row.names(feats)[order(feats$NeuAc.a2.3.Gal.b1.4.Glc, decreasing = T)][neuInds]
100 - (neuInds/nrow(feats))*100


fucFeats = c('aromatic_bin2','binnedD2_PC7', 'zern_PC4')

fucInds = (1:nrow(feats))[row.names(feats)[order(feats$Fuc, decreasing = T)] %in% fucFeats]
row.names(feats)[order(feats$Fuc, decreasing = T)][fucInds]
100 - (fucInds/nrow(feats))*100

fucInds = (1:nrow(feats))[row.names(feats)[order(feats$Terminal_Fucose, decreasing = T)] %in% fucFeats]
row.names(feats)[order(feats$Terminal_Fucose, decreasing = T)][fucInds]
100 - (fucInds/nrow(feats))*100


manFeats = c('hydrophobics', 'nonpolar_bin2', 'aromatic_bin2', 'var_4Ang', 'med_4Ang', 'q3_4Ang', 'var_6Ang', 'med_6Ang', 'q3_6Ang', 'binnedD2_PC5')

manInds = (1:nrow(feats))[row.names(feats)[order(feats$Man, decreasing = T)] %in% manFeats]
row.names(feats)[order(feats$Man, decreasing = T)][manInds]
100 - (manInds/nrow(feats))*100
for(i in 1:length(manFeats)){
  cat(row.names(feats)[order(feats$Man, decreasing = T)][manInds][i], ' ', (100 - (manInds/nrow(feats))*100)[i], '\n')
}

manInds = (1:nrow(feats))[row.names(feats)[order(feats$High_Mannose, decreasing = T)] %in% manFeats]
row.names(feats)[order(feats$High_Mannose, decreasing = T)][manInds]
100 - (manInds/nrow(feats))*100
for(i in 1:length(manFeats)){
  cat(row.names(feats)[order(feats$High_Mannose, decreasing = T)][manInds][i], ' ', (100 - (manInds/nrow(feats))*100)[i], '\n')
}

manInds = (1:nrow(feats))[row.names(feats)[order(feats$Man.a1.2.Man, decreasing = T)] %in% manFeats]
row.names(feats)[order(feats$Man.a1.2.Man, decreasing = T)][manInds]
100 - (manInds/nrow(feats))*100
for(i in 1:length(manFeats)){
  cat(row.names(feats)[order(feats$Man.a1.2.Man, decreasing = T)][manInds][i], ' ', (100 - (manInds/nrow(feats))*100)[i], '\n')
}
