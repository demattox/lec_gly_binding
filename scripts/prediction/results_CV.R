
library(reshape)
library(ggplot2)
library(philentropy)
library(pheatmap)
library(PRROC)
library(vioplot)

###############
# Functions
###############
colfunc <- colorRampPalette(c("red","goldenrod","forestgreen","royalblue","darkviolet"))

getKPRFb <- function(conMatDF){
  # conMatDF has rows of different confusion matrices, with the columns ordered as TP, TN, FP, FN
  # sums each column and finds performance metrics
  f2TestCon <- apply(X = conMatDF, MARGIN = 2, FUN = sum)
  
  TP <- f2TestCon[[1]]
  TN <- f2TestCon[[2]]
  FP <- f2TestCon[[3]]
  FN <- f2TestCon[[4]]
  
  f2_validationRecall <- TP / ( TP + FN)
  f2_validationPrec <- TP / (TP + FP)
  f2_validationF2 <- ((1+(2^2)) * f2_validationPrec * f2_validationRecall) / (2^2 * f2_validationPrec + f2_validationRecall)
  # f3score <- ((1+(3^2)) * f2_validationPrec * f2_validationRecall) / (3^2 * f2_validationPrec + f2_validationRecall)
  # f4score <- ((1+(4^2)) * f2_validationPrec * f2_validationRecall) / (4^2 * f2_validationPrec + f2_validationRecall)
  randAcc <- ((TN+FP)*(TN+FN) + (FN+TP)*(FP+TP)) / sum(c(TP + TN + FP + FN))^2
  testAcc <- (TP+TN)/(TP+TN+FP+FN)
  f2_validationKappa <- (testAcc - randAcc) / (1 - randAcc)
  return(list(kappa = f2_validationKappa,
              recall = f2_validationRecall,
              F2 = f2_validationF2,
              precision = f2_validationPrec))
}


###############
# read in files
###############
homeDir <- '/Users/dmattox/cbk/lec_gly_binding/'
setwd(homeDir)

ligTags <- read.delim(file = './analysis/training/data_in/ligTags.tsv', sep = '\t', stringsAsFactors = F)
predFeats <- read.delim(file = './analysis/training/data_in/predFeats.csv', sep = ',', stringsAsFactors = F)
bsResiDat <- read.delim(file = './analysis/training/data_in/bsResiDat.tsv', sep = '\t', stringsAsFactors = F)

ligColors <- colfunc(ncol(ligTags))

# inDir <- './analysis/training/train_and_validate/seqID50/'
# inDir <- './analysis/training/train_and_validate/seqID80/'
# inDir <- './analysis/training/train_and_validate/beta3id50/'

#inDir <- './analysis/training/train_and_validate/nr_TrainVal/seqID50/'
# inDir <- './analysis/training/train_and_validate/nr_TrainVal/seqID80/'

inDir <- './analysis/training/train_and_validate/nr_10xCVonly/pred/'
inRand <- './analysis/training/train_and_validate/nr_10xCVonly/rand/'

# Read in results
inFiles <- dir(inDir)
inTrain <- inFiles[grepl('training.csv', inFiles)]
#inTest <- inFiles[grepl('testing.csv', inFiles)]
#inOutcomes <- inFiles[grepl('outcomes.csv', inFiles)]
inFeats <- inFiles[grepl('features.csv', inFiles)]

newOrd <- match(colnames(ligTags), gsub('(.*)_features.csv', '\\1', inFeats))
inTrain <- inTrain[newOrd]
#inTest <- inTest[newOrd]
#inOutcomes <- inOutcomes[newOrd]
inFeats <- inFeats[newOrd]


# Read in results from classifier trained on shuffled labels
randFiles <- dir(inRand)
randTrain <- randFiles[grepl('training.csv', randFiles)]
#randTest <- randFiles[grepl('testing.csv', randFiles)]
#randOutcomes <- randFiles[grepl('outcomes.csv', randFiles)]
randFeats <- randFiles[grepl('features.csv', randFiles)]

newOrd <- match(colnames(ligTags), gsub('(.*)_features.csv', '\\1', randFeats))
randTrain <- randTrain[newOrd]
#randTest <- randTest[newOrd]
#randOutcomes <- randOutcomes[newOrd]
randFeats <- randFeats[newOrd]

###############
# Look at lectin specificity
###############
lecIDs <- unique(bsResiDat$uniparc)

iupacLigCnt <- rep(0, length(lecIDs))

for (i in 1:length(lecIDs)){
  id <- lecIDs[i]
  iupacLigCnt[i] <- length(unique(bsResiDat$iupac[bsResiDat$uniparc == id]))
}

# Get ligand tags again
uniLigs <- unique(bsResiDat$iupac)

# parenCnt <- bracCnt <- manCnt <- neuCnt <- bracCnt <- rep(0,length(uniLigs))
# for (i in 1:length(uniLigs)){
#   lig <- uniLigs[i]
#   parenCnt[i] <- lengths(regmatches(lig, gregexpr("\\(", lig)))
#   bracCnt[i] <- lengths(regmatches(lig, gregexpr("\\[", lig)))
#   manCnt[i] <- lengths(regmatches(lig, gregexpr("Man", lig)))
#   neuCnt[i] <- lengths(regmatches(lig, gregexpr("NeuAc", lig)))
# }
# 
# mTag <- parenCnt == 0 & bracCnt == 0 # Monosaccharides
# dTag <- parenCnt == 1 & bracCnt == 0 # Disaccharides
# tTag <- (parenCnt == 2 & bracCnt == 0) | (parenCnt == 2 & bracCnt == 1) # Trisaccharides
# qTag <- (parenCnt == 3 & bracCnt == 0) | (parenCnt == 3 & bracCnt == 1) | (parenCnt == 1 & bracCnt == 2) # Tetrasaccharides
# pTag <- !(mTag | dTag | tTag | qTag) # 5+ sugars
# bTag <- bracCnt >= 1 # Branched glycans
# 
# manTag <- manCnt > 3 # High mannose
# # uniLigs[manTag]
# neuTag <- neuCnt >= 1 # Has sialic acid
# # uniLigs[neuTag]
# fucTag <- grepl('^Fuc',uniLigs) # Has a terminal fucose
# # uniLigs[fucTag]

# Look at lectin specificity
lecSpec <- as.data.frame(matrix(nrow = length(lecIDs), ncol = length(uniLigs)))
row.names(lecSpec) <- lecIDs
colnames(lecSpec) <- uniLigs
for (i in 1:length(lecIDs)){
  lecSpec[i,] <- uniLigs %in% unique((bsResiDat$iupac[bsResiDat$uniparc == lecIDs[i]]))
}


###############
# Training performance
###############

sampSizes <- rep(0, ncol(ligTags))

for(i in 1:length(inTrain)){
  lig <- gsub('(.*)_training.csv', '\\1', inTrain)[i]
  tmp <- read.delim(file = paste(inDir, inTrain[i], sep = ''), header = T, sep = ',', stringsAsFactors = F)
  sampSizes[i] <- mean(apply(tmp[4:7], 1, sum)) 
  tmp <- cbind(ligand = rep(lig, nrow(tmp)), tmp[,c(-1)])
  if (i == 1){
    train <- tmp
  } else {
    train <- rbind(train, tmp)
  }
}
train$f2 <- 0
train$prec <- 0
for(i in 1:nrow(train)){
  r <- train$recall[i]
  p <- train$TP[i] / (train$TP[i] + train$FP[i])
  train$prec[i] <- p
  train$f2[i] <- ((1+(2^2)) * p * r) / (2^2 * p + r)
  # trainOut$f3[i] <- ((1+(3^2)) * p * r) / (3^2 * p + r)
  # trainOut$f4[i] <- ((1+(4^2)) * p * r) / (4^2 * p + r)
}

mTrain <- melt(train, id.vars = 'ligand', measure.vars = c('kappa', 'recall', 'prec'))
colnames(mTrain) <- c('ligand', 'metric', 'value')
# Light colors
ggplot(data = mTrain, aes(x = metric, y = value, col = ligand, fill = metric)) +
  geom_hline(yintercept = 0) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.4) +
  geom_boxplot(outlier.alpha = 0) +
  ylim(c(-1, 1)) +
  scale_fill_manual(values = alpha(rep('snow3',length(unique(mTrain$metric))), 0.6), guide =F) +
  scale_color_manual(values = ligColors) +
  labs(title = '10x CV validation performance', x = "Metric type", y = "Metric value") +
  theme_light(base_size = 22) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))



for(i in 1:length(randTrain)){
  lig <- gsub('(.*)_training.csv', '\\1', randTrain)[i]
  tmp <- read.delim(file = paste(inRand, randTrain[i], sep = ''), header = T, sep = ',', stringsAsFactors = F)
  sampSizes[i] <- mean(apply(tmp[4:7], 1, sum)) 
  tmp <- cbind(ligand = rep(lig, nrow(tmp)), tmp[,c(-1)])
  if (i == 1){
    trainRand <- tmp
  } else {
    trainRand <- rbind(trainRand, tmp)
  }
}
trainRand$f2 <- 0
trainRand$prec <- 0
for(i in 1:nrow(trainRand)){
  r <- trainRand$recall[i]
  p <- trainRand$TP[i] / (trainRand$TP[i] + trainRand$FP[i])
  trainRand$prec[i] <- p
  trainRand$f2[i] <- ((1+(2^2)) * p * r) / (2^2 * p + r)
  # trainRandOut$f3[i] <- ((1+(3^2)) * p * r) / (3^2 * p + r)
  # trainRandOut$f4[i] <- ((1+(4^2)) * p * r) / (4^2 * p + r)
}

par(cex.lab=1.5)
plot(0,0,col = 'white', xlab = '',ylab = '',type="n",
     axes=FALSE,ann=FALSE,
     xlim = c(0,16),
     ylim = c(-1,1))
i<-1 
vioplot(train$kappa[train$ligand == colnames(ligTags)[i]], at = i, side = 'left',
        add=T,
        col = ligColors[i],
        xlim = c(0,16), ylim = c(-1,1),
        plotCentre = "line")
vioplot(trainRand$kappa[trainRand$ligand == colnames(ligTags)[i]], at = i, side = 'right',
        add = T,
        col = alpha(ligColors[i],0.33),
        xlim = c(0,16), ylim = c(-1,1),
        plotCentre = "line")
for(i in 2:ncol(ligTags)){
  vioplot(train$kappa[train$ligand == colnames(ligTags)[i]], at = i, side = 'left',
          col = ligColors[i],
          xlim = c(0,16), ylim = c(-1,1),
          add = T,
          plotCentre = "line")
  vioplot(trainRand$kappa[trainRand$ligand == colnames(ligTags)[i]], at = i, side = 'right',
          add = T,
          col = alpha(ligColors[i],0.3),
          xlim = c(0,16), ylim = c(-1,1),
          plotCentre = "line")
}
abline(h=0, lty = 2)
axis(side=1,at=1:15, labels = rep('', 15))
axis(side=2,at=seq.int(-1,1,0.5))
title(xlab = "Ligands", ylab = "Kappa", main = "10x CV - Kappa compared to random classifiers")




par(cex.lab=1.5)
plot(0,0,col = 'white', xlab = '',ylab = '',type="n",
     axes=FALSE,ann=FALSE,
     xlim = c(0,16),
     ylim = c(0,1))
i<-1 
vioplot(train$prec[train$ligand == colnames(ligTags)[i]], at = i, side = 'left',
        add=T,
        col = ligColors[i],
        xlim = c(0,16), ylim = c(0,1),
        plotCentre = "line")
vioplot(trainRand$prec[trainRand$ligand == colnames(ligTags)[i]], at = i, side = 'right',
        add = T,
        col = alpha(ligColors[i],0.33),
        xlim = c(0,16), ylim = c(0,1),
        plotCentre = "line")
for(i in 2:ncol(ligTags)){
  vioplot(train$prec[train$ligand == colnames(ligTags)[i]], at = i, side = 'left',
          col = ligColors[i],
          xlim = c(0,16), ylim = c(0,1),
          add = T,
          plotCentre = "line")
  vioplot(trainRand$prec[trainRand$ligand == colnames(ligTags)[i]], at = i, side = 'right',
          add = T,
          col = alpha(ligColors[i],0.3),
          xlim = c(0,16), ylim = c(0,1),
          plotCentre = "line")
}
axis(side=1,at=1:15, labels = rep('', 15))
axis(side=2,at=seq.int(0,1,0.25))
title(xlab = "Ligands", ylab = "Precision", main = "10x CV - Precision compared to random classifiers")



par(cex.lab=1.5)
plot(0,0,col = 'white', xlab = '',ylab = '',type="n",
     axes=FALSE,ann=FALSE,
     xlim = c(0,16),
     ylim = c(0,1))
i<-1 
vioplot(train$recall[train$ligand == colnames(ligTags)[i]], at = i, side = 'left',
        add=T,
        col = ligColors[i],
        xlim = c(0,16), ylim = c(0,1),
        plotCentre = "line")
vioplot(trainRand$recall[trainRand$ligand == colnames(ligTags)[i]], at = i, side = 'right',
        add = T,
        col = alpha(ligColors[i],0.33),
        xlim = c(0,16), ylim = c(0,1),
        plotCentre = "line")
for(i in 2:ncol(ligTags)){
  vioplot(train$recall[trainRand$ligand == colnames(ligTags)[i]], at = i, side = 'left',
          col = ligColors[i],
          xlim = c(0,16), ylim = c(0,1),
          add = T,
          plotCentre = "line")
  vioplot(trainRand$recall[trainRand$ligand == colnames(ligTags)[i]], at = i, side = 'right',
          add = T,
          col = alpha(ligColors[i],0.3),
          xlim = c(0,16), ylim = c(0,1),
          plotCentre = "line")
}
axis(side=1,at=1:15, labels = rep('', 15))
axis(side=2,at=seq.int(0,1,0.25))
title(xlab = "Ligands", ylab = "Recall", main = "10x CV - Recall compared to random classifiers")


## P & R same plot
xLim = c(0,30)
yLim = c(0,1)

par(cex.lab=1.5, mar = c(5, 4, 4, 4) + 0.3)
plot(0,0,col = 'white', xlab = '',ylab = '',type="n",
     axes=FALSE,ann=FALSE,
     xlim = xLim,
     ylim = c(0,1))
i<-1 
vioplot(train$recall[train$ligand == colnames(ligTags)[i]], at = i, side = 'left',
        add=T,
        col = ligColors[i],
        xlim = xLim, ylim = yLim,
        plotCentre = "line")
par(new = T)
boxplot(trainRand$recall[trainRand$ligand == colnames(ligTags)[i]], at = i-0.33, notch = T, outline = F,
        col = alpha(ligColors[i],0.2), border = 'grey50',
        axes = F, xlab = '', ylab = '', xlim = xLim, ylim = yLim)
vioplot(train$prec[trainRand$ligand == colnames(ligTags)[i]], at = i, side = 'right',
        add = T,
        col = alpha(ligColors[i],0.8),
        xlim = xLim, ylim = yLim,
        plotCentre = "line")
par(new = T)
boxplot(trainRand$prec[trainRand$ligand == colnames(ligTags)[i]], at = i+0.33, notch = T, outline = F,
        col = alpha(ligColors[i],0.2), border = 'grey50',
        axes = F, xlab = '', ylab = '', xlim = xLim, ylim = yLim)
for(j in 2:ncol(ligTags)){
  i = ((j-1) * 2) + 1
  vioplot(train$recall[trainRand$ligand == colnames(ligTags)[j]], at = i, side = 'left',
          col = ligColors[j],
          xlim = xLim, ylim = yLim,
          add = T,
          plotCentre = "line")
  par(new = T)
  boxplot(trainRand$recall[trainRand$ligand == colnames(ligTags)[j]], at = i-0.33, notch = T, outline = F,
          col = alpha(ligColors[j],0.2), border = 'grey50',
          axes = F, xlab = '', ylab = '', xlim = xLim, ylim = yLim)
  vioplot(train$prec[trainRand$ligand == colnames(ligTags)[j]], at = i, side = 'right',
          add = T,
          col = alpha(ligColors[j],0.8),
          xlim = xLim, ylim = yLim,
          plotCentre = "line")
  par(new = T)
  boxplot(trainRand$prec[trainRand$ligand == colnames(ligTags)[j]], at = i+0.33, notch = T, outline = F,
          col = alpha(ligColors[j],0.2), border = 'grey50',
          axes = F, xlab = '', ylab = '', xlim = xLim, ylim = yLim)
}
axis(side=1,at=seq.int(1,29,2), labels = rep('', 15))
axis(side=2,at=pretty(c(0,1)))
axis(side = 4, at = pretty(c(0,1)))      # Add second axis
abline(h = 0.5, lwd = 0.75, lty = 2)
mtext("Precision", side = 4, line = 3, cex = 1.5, col = 'black')             # Add second axis label
title(xlab = "Ligands", ylab = "Recall", main = "10x CV - Precision & Recall vs random")

















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

for(i in 1:length(randFeats)){
  lig = gsub('(.*)_features.csv', '\\1', randFeats)[i]
  tmp = read.delim(file = paste(inRand, randFeats[i], sep = ''), header = T, sep = ',', stringsAsFactors = F)
  tmp = apply(tmp, 1, mean)
  if (i == 1){
    featsRand = tmp
  } else {
    featsRand = cbind(featsRand, tmp)
    colnames(featsRand)[i] = lig
  }
}
colnames(featsRand)[1] = gsub('(.*)_features.csv', '\\1', randFeats)[1]
featsRand = as.data.frame(featsRand)

