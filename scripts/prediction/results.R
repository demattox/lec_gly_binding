

library(reshape)
library(ggplot2)
library(philentropy)

library(pheatmap)

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

inDir = './analysis/training/train_and_validate/seqID50/'

inFiles = dir(inDir)
inTrain = inFiles[grepl('training.csv', inFiles)]
inTest = inFiles[grepl('testing.csv', inFiles)]
inOutcomes = inFiles[grepl('outcomes.csv', inFiles)]
inFeats = inFiles[grepl('features.csv', inFiles)]

newOrd = match(colnames(ligTags), gsub('(.*)_outcomes.csv', '\\1', inOutcomes))
inTrain = inTrain[newOrd]
inTest = inTest[newOrd]
inOutcomes = inOutcomes[newOrd]
inFeats = inFeats[newOrd]

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
# uniLigs[manTag]
neuTag = neuCnt >= 1 # Has sialic acid
# uniLigs[neuTag]
fucTag = grepl('^Fuc',uniLigs) # Has a terminal fucose
# uniLigs[fucTag]

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

sampSizes = rep(0, ncol(ligTags))

for(i in 1:length(inTrain)){
  lig = gsub('(.*)_training.csv', '\\1', inTrain)[i]
  tmp = read.delim(file = paste(inDir, inTrain[i], sep = ''), header = T, sep = ',', stringsAsFactors = F)
  sampSizes[i] = mean(apply(tmp[3:6], 1, sum)) 
  tmp = cbind(ligand = rep(lig, nrow(tmp)), tmp[,(-3:-6)])
  if (i == 1){
    train = tmp
  } else {

    train = rbind(train, tmp)
  }
}

mTrain = melt(train, id.vars = 'ligand')
colnames(mTrain) = c('ligand', 'metric', 'value')



ggplot(data = mTrain, aes(x = metric, y = value, col = ligand)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.4) +
  geom_boxplot(outlier.alpha = 0) +
  ylim(c(0.5, 1)) +
  scale_fill_manual(values = alpha(c('navajowhite','navajowhite'), 0.6), guide =F) +
  scale_color_manual(values = ligColors) +
  labs(title = '5x CV with LOCO validation - Training Performance', x = "Metric type", y = "Metric value") +
  theme_dark(base_size = 22) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))

ggplot(data = mTrain, aes(x = ligand, y = value, col = metric)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.4) +
  geom_boxplot(outlier.alpha = 0) +
  ylim(c(0.5, 1)) +
  scale_fill_manual(values = alpha(c('navajowhite','navajowhite'), 0.6), guide =F) +
  scale_color_manual(values = c('darkgreen', 'firebrick1', 'darkorchid1', 'dodgerblue')) +
  labs(title = '5x CV with LOCO validation - Training Performance', x = "Ligand for prediciton", y = "Metric value") +
  theme_dark(base_size = 22) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))


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
test = as.data.frame(matrix(0,nrow= ncol(ligTags), ncol = 4))
row.names(test) = colnames(ligTags)
# colnames(test) = c('ligand', 'kappa', 'recall', 'f2', 'precision')
colnames(test) = c('kappa', 'recall', 'f2', 'precision')

for(i in 1:length(inTest)){
  lig = gsub('(.*)_testing.csv', '\\1', inTest)[i]

  tmp = read.delim(file = paste(inDir, inTest[i], sep = ''), header = T, sep = ',', stringsAsFactors = F)
  
  R =  tmp$TP/(tmp$TP + tmp$FN)
  R = as.data.frame(cbind(rep(lig, nrow(tmp)), R))
  colnames(R) = c('ligand', 'recall')
  if (i ==1){
    recall = R
  } else {
    recall = rbind(recall, R)
  }
  # test[lig,'ligand'] = lig
  test[lig,] = as.numeric(getKPRFb(tmp)) # Kapppa, recall, f2, precision
  
}
# mTest= melt(test, id.vars = 'ligand')

breakLst = seq(0,1,0.01)
mycol = colorRampPalette(c("ivory", "cornflowerblue", "navy"))(n = length(breakLst))


pheatmap(test, color =  mycol,
         cluster_cols = F, cluster_rows = F,
         display_numbers = T, number_color = 'red', fontsize_number = 14,
         breaks = breakLst, main = 'Aggregate validation metrics for each ligand class', 
         gaps_row = c(3), gaps_col = c(1))

plot(sampSizes, test$kappa, pch = 19, cex = 1.5,
     ylim = c(-0.15, 0.4),
     xlim = c(0,500),
     xlab = 'Average number of samples sites for training',
     ylab = 'Validation Kappa (aggregate)',
     cex.lab=1.5, cex.axis=1.5)

cor.test(sampSizes, test$kappa)


recall$recall = as.numeric(as.character(recall$recall))

ggplot(data = recall, aes(x = ligand, y = recall, col = ligand, fill = ligand)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.4) +
  geom_boxplot(outlier.alpha = 0) +
  ylim(c(0, 1)) +
  scale_fill_manual(values = alpha(rep('snow3',ncol(ligTags)), 0.6), guide =F) +
  scale_color_manual(values = ligColors, guide = F) +
  labs(title = 'Validation Recall for each excluded cluster', x = "Ligand for prediciton", y = "Cluster-specific Recall value") +
  theme_light(base_size = 22) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))

# PR curves
par(mfrow = c(3,5))
for(i in 1:ncol(ligTags)){
  lig = gsub('(.*)_outcomes.csv', '\\1', inOutcomes)[i]
  outcomes = read.delim(file = paste(inDir, inOutcomes[i], sep = ''), header = T, sep = ',', stringsAsFactors = F)
  
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
}



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
