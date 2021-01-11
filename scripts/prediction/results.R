
library(reshape)
library(ggplot2)
library(philentropy)
library(pheatmap)
library(PRROC)
library(vioplot)
library(VennDiagram)

# devtools::install_github("copenhagencenterforglycomics/ggsugar")
# library(ggsugar)


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

####
ligNames = colnames(ligTags)
ligNames = gsub('_', ' ', ligNames)

ligNames[12] = expression(bold(paste("2", alpha, "-Mannobiose", sep = '')))

ligNames[1] = expression(bold("Terminal NeuAc Group"))
ligNames[2] = expression(bold("High Mannose Group"))
ligNames[3] = expression(bold("Terminal Fuc Group"))
ligNames[4] = expression(bold("Lactose"))
ligNames[5] = expression(bold("Galactose"))
ligNames[6] = expression(bold("Mannose"))
ligNames[7] = expression(bold("N-Acetyl Galactosamine"))
ligNames[8] = expression(bold("N-Acetyl Neuraminic Acid"))
ligNames[9] = expression(bold("3'-Siayllactose"))
ligNames[10] = expression(bold("Glucose"))
ligNames[11] = expression(bold("N-Acetyl Lactosamine"))
ligNames[13] = expression(bold("N-Acetyl Glucosamine"))
ligNames[14] = expression(bold("Fucose"))
ligNames[15] = expression(bold("TF Antigen"))


ligColors = rep('', ncol(ligTags))
ligColors[1] = 'purple2' # Sialic Acid
ligColors[8] = 'darkviolet' # NeuAc monosacc.
ligColors[9] = 'purple2' # 3' Sialyllactose

ligColors[2] = 'forestgreen' # High mannose
ligColors[6] = 'darkgreen' # Mannose monosacc.
ligColors[12] = 'forestgreen' # 2alpha mannobiose

ligColors[3] = 'red1' # Terminal Fuc
ligColors[14] = 'firebrick3' # Fuc monosacc.

ligColors[4] = 'goldenrod2' # Lactose
ligColors[5] = 'darkgoldenrod3' # Gal monosacc.
ligColors[7] = 'darkgoldenrod3' # GalNAc (Tn antigen)
ligColors[11] = 'goldenrod2' # N-Acetyllactosamine (LacNAc)
ligColors[15] = 'goldenrod2' # TF antigen

ligColors[10] = 'mediumblue' # Glc monosacc.
ligColors[13] = 'royalblue2' # GlcNAc
####

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

for(i in 1:length(dir(predDirs))){
  lig = colnames(ligTags)[as.numeric(dir(predDirs)[i])]
  cat(lig,'\n')
  curDir = paste(predDirs, dir(predDirs)[i], sep = '')
  subDirs = dir(curDir)
  for (j in 1:length(subDirs)){
    cat(j,'\t')
    inFiles = dir(paste(curDir, subDirs[j], sep = '/'))
    
    readFile = inFiles[grepl('training.csv', inFiles)] # change here to read other files
    
    tmp = read.delim(file = paste(curDir, subDirs[j], readFile, sep = '/'), header = T, sep = ',', stringsAsFactors = F)
    tmp$mode = 'pred'
    tmp$ligand = lig
    if (!exists('trainDat')){
      trainDat = tmp
    } else{
      trainDat = rbind(trainDat, tmp)
    }
  }
  cat('\n')
}


for(i in 1:length(dir(randDirs))){
  lig = colnames(ligTags)[i]
  cat(lig,'\n')
  curDir = paste(randDirs, dir(randDirs)[i], sep = '')
  subDirs = dir(curDir)
  for (j in 1:length(subDirs)){
    cat(j,'\t')
    inFiles = dir(paste(curDir, subDirs[j], sep = '/'))
    
    readFile = inFiles[grepl('training.csv', inFiles)] # change here to read other files
    
    tmp = read.delim(file = paste(curDir, subDirs[j], readFile, sep = '/'), header = T, sep = ',', stringsAsFactors = F)
    tmp$mode = 'rand'
    tmp$ligand = lig
    if (! exists('trainDat')){
      trainDat = tmp
    }else{
      trainDat = rbind(trainDat, tmp)
    }
  }
  cat('\n')
}
  


trainDat$f2 = 0
trainDat$prec = 0
for(i in 1:nrow(trainDat)){
  r = trainDat$recall[i]
  p = trainDat$TP[i] / (trainDat$TP[i] + trainDat$FP[i])
  trainDat$prec[i] = p
  trainDat$f2[i] = ((1+(2^2)) * p * r) / (2^2 * p + r)
  # trainDatOut$f3[i] = ((1+(3^2)) * p * r) / (3^2 * p + r)
  # trainDatOut$f4[i] = ((1+(4^2)) * p * r) / (4^2 * p + r)
}

mTrain = melt(trainDat, id.vars = c('ligand', 'mode'), measure.vars = c('kappa', 'recall', 'prec'))
colnames(mTrain) = c('ligand', 'classifier', 'metric', 'value')

# Light colors all box plots
# ggplot(data = mTrain, aes(x = metric, y = value, col = ligand, fill = classifier)) +
#   geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.1) +
#   geom_boxplot(outlier.alpha = 0) +
#   ylim(c(-1, 1)) +
#   scale_fill_manual(values = c(alpha('snow3', 0.6), alpha('black',0.8))) +
#   scale_color_manual(values = ligColors) +
#   labs(title = '5x CV with LOCO validation - Training Performance', x = "Metric type", y = "Metric value") +
#   theme_light(base_size = 22) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))


## P & R same plot with violin & box plots
xLim = c(0,30)
yLim = c(0,1)

par(cex.lab=1.5, mar = c(8.6, 5.1, 5.1, 4.1))
plot(0,0,col = 'white', xlab = '',ylab = '',type="n",
     axes=FALSE,ann=FALSE,
     xlim = xLim,
     ylim = c(0,1))
for(j in 1:ncol(ligTags)){
  i = ((j-1) * 2) + 1
  vioplot(trainDat$recall[(trainDat$ligand == colnames(ligTags)[j]) & (trainDat$mode == 'pred')], at = i, side = 'left',
          col = ligColors[j],
          xlim = xLim, ylim = yLim,
          add = T,
          plotCentre = "line")
  par(new = T)
  boxplot(trainDat$recall[(trainDat$ligand == colnames(ligTags)[j]) & (trainDat$mode == 'rand')], at = i-0.33, notch = T, outline = F,
          col = alpha(ligColors[j],0.2), border = 'grey50',
          axes = F, xlab = '', ylab = '', xlim = xLim, ylim = yLim)
  vioplot(trainDat$prec[(trainDat$ligand == colnames(ligTags)[j]) & (trainDat$mode == 'pred')], at = i, side = 'right',
          add = T,
          col = alpha(ligColors[j],0.8),
          xlim = xLim, ylim = yLim,
          plotCentre = "line")
  par(new = T)
  boxplot(trainDat$prec[(trainDat$ligand == colnames(ligTags)[j]) & (trainDat$mode == 'rand')], at = i+0.33, notch = T, outline = F,
          col = alpha(ligColors[j],0.2), border = 'grey50',
          axes = F, xlab = '', ylab = '', xlim = xLim, ylim = yLim)
}
axis(side=1,at=seq.int(1,29,2), labels = F)
axis(side=2,at=pretty(c(0,1)))
axis(side = 4, at = pretty(c(0,1)))  # Add second axis
# abline(h = 0.5, lwd = 0.75, lty = 2)
abline(v = 6, lwd = 0.75, lty=1)
mtext("Precision", side = 4, line = 3, cex = 2, col = alpha('black',0.85))  # Add second axis label
title(main = "5x CV + LOCO - TRAINING\nPrecision & Recall vs random", xlab = "", ylab = "Recall", cex.lab = 2)
# labs(title = , x = "Ligands", y = "Recall") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))
text(x = seq.int(1,29,2),                   
     y = par("usr")[3] - 0.05,
     labels = ligNames,
     col = ligColors,
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965,
     ## Increase label size.
     cex = 1.2)


###############
# Validation performance
###############
# rm(testDat)
for(i in 1:length(dir(predDirs))){
  lig = colnames(ligTags)[as.numeric(dir(predDirs)[i])]
  cat(lig,'\n')
  curDir = paste(predDirs, dir(predDirs)[i], sep = '')
  subDirs = dir(curDir)
  for (j in 1:length(subDirs)){
    cat(j,'\t')
    inFiles = dir(paste(curDir, subDirs[j], sep = '/'))
    
    readFile = inFiles[grepl('testing.csv', inFiles)] # change here to read other files
    
    tmp = read.delim(file = paste(curDir, subDirs[j], readFile, sep = '/'), header = T, sep = ',', stringsAsFactors = F)
    
    outTmp = as.data.frame(matrix(0,nrow= 10, ncol = 4))
    colnames(outTmp) = c('kappa', 'recall', 'f2', 'prec')
    
    for(k in 1:10){
      tag = grepl(paste('_', as.character(k), '$', sep = ''), row.names(tmp))
      repDat = tmp[tag,]
      outTmp[k,] = as.numeric(getKPRFb(repDat)) # Kapppa, recall, f2, precision
    }
    
    outTmp$mode = 'pred'
    outTmp$ligand = lig
    outTmp$batch = as.numeric(subDirs[j])
    if (! exists('testDat')){
      testDat = outTmp
    }else{
      testDat = rbind(testDat, outTmp)
    }
  }
  cat('\n')
}

for(i in 1:length(dir(randDirs))){
  lig = colnames(ligTags)[i]
  cat(lig,'\n')
  curDir = paste(randDirs, dir(randDirs)[i], sep = '')
  subDirs = dir(curDir)
  for (j in 1:length(subDirs)){
    cat(j,'\t')
    inFiles = dir(paste(curDir, subDirs[j], sep = '/'))
    
    readFile = inFiles[grepl('testing.csv', inFiles)] # change here to read other files
    
    tmp = read.delim(file = paste(curDir, subDirs[j], readFile, sep = '/'), header = T, sep = ',', stringsAsFactors = F)
    
    outTmp = as.data.frame(matrix(0,nrow= 10, ncol = 4))
    colnames(outTmp) = c('kappa', 'recall', 'f2', 'prec')
    
    for(k in 1:10){
      tag = grepl(paste('_', as.character(k), '$', sep = ''), row.names(tmp))
      repDat = tmp[tag,]
      outTmp[k,] = as.numeric(getKPRFb(repDat)) # Kapppa, recall, f2, precision
    }
    
    outTmp$mode = 'rand'
    outTmp$ligand = lig
    outTmp$batch = as.numeric(subDirs[j])
    if (! exists('testDat')){
      testDat = outTmp
    }else{
      testDat = rbind(testDat, outTmp)
    }
  }
  cat('\n')
}


mTest= melt(testDat[testDat$mode == 'pred',], id.vars = 'ligand', measure.vars = c('kappa', 'recall', 'prec'))
colnames(mTest) = c('ligand', 'metric', 'value')


# ggplot(data = mTest, aes(x = metric, y = value, col = ligand, fill = metric)) +
#   geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.4) +
#   geom_boxplot(outlier.alpha = 0) +
#   ylim(c(-0.15, 1)) +
#   scale_fill_manual(values = alpha(rep('snow3',length(unique(mTest$metric))), 0.6), guide =F) +
#   scale_color_manual(values = ligColors) +
#   labs(title = '5x CV with LOCO validation - Validation Performance', x = "Metric type", y = "Metric value") +
#   theme_light(base_size = 22) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))

# Kappa vs random
par(mar = c(8.6, 5.1, 5.1, 4.1), # change the margins
    lwd = 2, # increase the line thickness
    cex.axis = 1.2) # increase default axis label size
plot(0,0,col = 'white', xlab = '',ylab = '',type="n",
     axes=FALSE,ann=FALSE,
     xlim = c(0,16),
     ylim = c(-1,1))
i<-1 
vioplot(testDat$kappa[(testDat$ligand == colnames(ligTags)[i]) & (testDat$mode == 'pred')], at = i, side = 'left',
        add=T,
        col = ligColors[i],
        xlim = c(0,16), ylim = c(-1,1),
        plotCentre = "line")
vioplot(testDat$kappa[(testDat$ligand == colnames(ligTags)[i]) & (testDat$mode == 'rand')], at = i, side = 'right',
        add = T,
        col = alpha(ligColors[i],0.33),
        xlim = c(0,16), ylim = c(-1,1),
        plotCentre = "line")
for(i in 2:ncol(ligTags)){
  vioplot(testDat$kappa[(testDat$ligand == colnames(ligTags)[i]) & (testDat$mode == 'pred')], at = i, side = 'left',
          col = ligColors[i],
          xlim = c(0,16), ylim = c(-1,1),
          add = T,
          plotCentre = "line")
  vioplot(testDat$kappa[(testDat$ligand == colnames(ligTags)[i]) & (testDat$mode == 'rand')], at = i, side = 'right',
          add = T,
          col = alpha(ligColors[i],0.3),
          xlim = c(0,16), ylim = c(-1,1),
          plotCentre = "line")
}
abline(h=0, lty = 2)
axis(side=1,at=1:15, labels = F)
axis(side=2,at=seq.int(-1,1,0.5))
title(xlab = "", ylab = "Kappa", main = "5x CV + LOCO\nKappa compared to random classifiers")
text(x = 1:15,
     y = par("usr")[3] - 0.08,
     labels = ligNames,
     col = ligColors,
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965,
     ## Increase label size.
     cex = 1.2)


# Check training sample sizes
trainDat$sampSizes = apply(trainDat[,4:7], MARGIN = 1, FUN = sum)

cor.test(trainDat$sampSizes[trainDat$mode == 'pred'], trainDat$kappa[trainDat$mode == 'pred'])

plot(0,0, type = 'n',
     xlim = c(0,175), ylim = c(-1, 1),
     xlab = '', ylab = '')
for (i in 1:ncol(ligTags)){
  par(new = T)
  plot(trainDat$sampSizes[trainDat$ligand == colnames(ligTags)[i] & trainDat$mode == 'pred'], trainDat$kappa[trainDat$ligand == colnames(ligTags)[i] & trainDat$mode == 'pred'],
       col = alpha(ligColors[i], 0.3), pch = 19,
       xlim = c(0,175), ylim = c(-1, 1),
       axes = F, xlab = '', ylab ='')
}
title(xlab = 'Number of samples used for training', ylab = 'Training Kappa', cex.lab = 1.5 )
text(x = 110, y = -0.8, labels = 'Pearson corr: 0.53 (p<0.001)', cex = 1.3)

CVperf_melt = melt(trainDat[trainDat$mode == 'pred',], id.vars = c("ligand", "sampSizes", "kappa", "f2"))


p = ggplot(data = CVperf_melt, aes(x = sampSizes, y = kappa, col = ligand))

p + geom_point() + scale_color_manual(values = alpha(ligColors, 0.3))

# p + geom_sugar(sugar='galnac')


# ggplot(data = mTest, aes(x = metric, y = value, col = ligand, fill = metric)) +
#   geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.4) +
#   geom_boxplot(outlier.alpha = 0) +
#   ylim(c(-0.15, 1)) +
#   scale_fill_manual(values = alpha(rep('snow3',length(unique(mTest$metric))), 0.6), guide =F) +
#   scale_color_manual(values = ligColors) +
#   labs(title = '5x CV with LOCO validation - Validation Performance', x = "Metric type", y = "Metric value") +
#   theme_light(base_size = 22) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))



breakLst = seq(0,1,0.01)
mycol = colorRampPalette(c("ivory", "cornflowerblue", "navy"))(n = length(breakLst))

# testMeds = as.data.frame(matrix(0,nrow  = length(unique(test$ligand)), ncol = 4))
# row.names(testMeds) = gsub('(.*)_training.csv', '\\1', inTrain)
# colnames(testMeds) = colnames(test)[-1]
# 
# for(i in 1:length(unique(test$ligand))){
#   tag = grepl(unique(test$ligand)[i], test$ligand)
#   testMeds[i,] = apply(test[tag,2:5], 2, median)
# }
# 
# testMeds$f2 <- NULL
# 
# pheatmap(testMeds, color =  mycol,
#          cluster_cols = F, cluster_rows = F,
#          display_numbers = T, number_color = 'darkgoldenrod', fontsize_number = 28,
#          breaks = breakLst, main = 'Median validation metrics for each ligand class', 
#          gaps_row = c(3), gaps_col = c(1))
# 
# 
# 
# plot(sampSizes, testMeds$kappa, pch = 19, cex = 1.5,
#      ylim = c(0, 0.5),
#      xlim = c(0,550),
#      xlab = 'Average number of samples sites for training',
#      ylab = 'Median Validation Kappa',
#      cex.lab=1.5, cex.axis=1.5)
# 
# cor.test(sampSizes, testMeds$kappa)

#########################
## Figure 3
#########################
## P & R same plot with violin & box plots
xLim = c(-1,30)
yLim = c(0,1)

par(mar = c(8.6, 5.1, 5.1, 4.1), # change the margins
    lwd = 2, # increase the line thickness
    cex.axis = 1.2 # increase default axis label size
)
# par(cex.lab=1.5, mar = c(5, 4, 4, 4) + 0.3)
plot(0,0,col = 'white', xlab = '',ylab = '',type="n",
     axes=FALSE,ann=FALSE,
     xlim = xLim,
     ylim = c(0,1))
for(j in 1:ncol(ligTags)){
  i = ((j-1) * 2) + 1
  vioplot(testDat$recall[(testDat$ligand == colnames(ligTags)[j]) & (testDat$mode == 'pred')], at = i, side = 'left',
          col = ligColors[j],
          xlim = xLim, ylim = yLim,
          add = T,
          plotCentre = "line")
  par(new = T)
  boxplot(testDat$recall[(testDat$ligand == colnames(ligTags)[j]) & (testDat$mode == 'rand')], at = i-0.33, notch = T, outline = F,
          col = alpha(ligColors[j],0.2), border = 'grey50',
          axes = F, xlab = '', ylab = '', xlim = xLim, ylim = yLim)
  vioplot(testDat$prec[(testDat$ligand == colnames(ligTags)[j]) & (testDat$mode == 'pred')], at = i, side = 'right',
          add = T,
          col = alpha(ligColors[j],0.8),
          xlim = xLim, ylim = yLim,
          plotCentre = "line")
  par(new = T)
  boxplot(testDat$prec[(testDat$ligand == colnames(ligTags)[j]) & (testDat$mode == 'rand')], at = i+0.33, notch = T, outline = F,
          col = alpha(ligColors[j],0.2), border = 'grey50',
          axes = F, xlab = '', ylab = '', xlim = xLim, ylim = yLim)
}

par(new = T) 
plot(0,0, col = 'white', type = 'n',
     xlim = xLim,
     ylim = c(1.3, 4),
     axes = F, xlab = '', ylab ='')
lines(x = c(1,30), y = c(0,0), lwd = 2)
for(j in 1:ncol(ligTags)){
  i = ((j-1) * 2) + 1
  par(new = T)
  boxplot(log10(trainDat$sampSizes[trainDat$ligand == colnames(ligTags)[j] & trainDat$mode == 'pred']),
          at = i,
          col = ligColors[j], border = ligColors[j],
          boxwex = 2.5,
          xlim = xLim,
          ylim = c(1.3, 4),
          axes = F, xlab = '', ylab ='')
}
logTicks = c(20,50,100,150)
axis(side=2,at=log10(logTicks), labels = logTicks, las = 2, pos = 0.3)


par(new = T)
plot(0,0,col = 'white', xlab = '',ylab = '',type="n",
     axes=FALSE,ann=FALSE,
     xlim = xLim,
     ylim = c(0,1))
axis(side=1,at=seq.int(1,29,2), labels = F)
axis(side=2,at=pretty(c(0,1)), las = 2)
axis(side = 4, at = pretty(c(0,1)), las = 2)  # Add second axis
# abline(v=6, lwd = 0.75)
mtext("Precision", side = 4, line = 3, cex = 2, col = alpha('black',0.85))  # Add second axis label
title(xlab = "", ylab = "Recall", main = "5x CV Random Forest - LO(C)O Validation performance\nPrecision & Recall vs random", cex.main = 1.5, cex.lab = 2)
text(x = seq.int(1,29,2),
     y = par("usr")[3] - 0.05,
     labels = ligNames,
     col = ligColors,
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965,
     ## Increase label size.
     cex = 1.2)







#######################
# PR curves
#######################

rawLigNames = rep('',length(ligNames))
rawLigNames[12] = "2 alpha-Mannobiose"

rawLigNames[1] = "Terminal NeuAc Group"
rawLigNames[2] = "High Mannose Group"
rawLigNames[3] = "Terminal Fuc Group"
rawLigNames[4] = "Lactose"
rawLigNames[5] = "Galactose"
rawLigNames[6] = "Mannose"
rawLigNames[7] = "N-Acetyl Galactosamine"
rawLigNames[8] = "N-Acetyl Neuraminic Acid"
rawLigNames[9] = "3'-Siayllactose"
rawLigNames[10] = "Glucose"
rawLigNames[11] = "N-Acetyl Lactosamine"
rawLigNames[13] = "N-Acetyl Glucosamine"
rawLigNames[14] = "Fucose"
rawLigNames[15] = "TF Antigen"


par(mfrow = c(3,5))

ligDirs = dir(predDirs)  
ligDirs  = ligDirs[order(as.numeric(ligDirs))]

allPredAUCs = as.data.frame(matrix(0, nrow = 100, ncol = ncol(ligTags)))
colnames(allPredAUCs) = colnames(ligTags)

allRandAUCs = allPredAUCs

for(i in 9:length(ligDirs)){
  lig = colnames(ligTags)[as.numeric(ligDirs[i])]
  cat(lig,'\n')
  
  curDir = paste(predDirs, ligDirs[i], sep = '')
  ranDir = gsub(pattern = 'pred', replacement = 'rand', x = curDir)
  
  subDirs = dir(curDir)
  
  # predAUC = randAUC = rep(0,100)
  
  for (j in 1:length(subDirs)){
    cat(j, '\t')
    #read in pred
    inFiles = dir(paste(curDir, subDirs[j], sep = '/'))
    
    readFile = inFiles[grepl('outcomes.csv', inFiles)] # change here to read other files
    
    tmp = read.delim(file = paste(curDir, subDirs[j], readFile, sep = '/'), header = T, sep = ',', stringsAsFactors = F)
    
    if (! exists('outcomes')){
      outcomes = tmp
    }else{
      outcomes = cbind(outcomes, tmp)
    }
    
    #read in rand
    inFiles = dir(paste(ranDir, subDirs[j], sep = '/'))
    
    readFile = inFiles[grepl('outcomes.csv', inFiles)] # change here to read other files
    
    tmp = read.delim(file = paste(ranDir, subDirs[j], readFile, sep = '/'), header = T, sep = ',', stringsAsFactors = F)
    
    extraRows = as.data.frame(matrix(NA, nrow = sum(!row.names(bsResiDat) %in% row.names(tmp)), ncol = 10))
    row.names(extraRows) = row.names(bsResiDat)[!row.names(bsResiDat) %in% row.names(tmp)]
    
    tmp = rbind(tmp, extraRows)
    tmp = tmp[row.names(bsResiDat),]
    
    if (! exists('Rand_outcomes')){
      Rand_outcomes = tmp

    }else{
      Rand_outcomes = cbind(Rand_outcomes, tmp)
    }
  }
  cat('\n')
  
  Rand_outcomes = Rand_outcomes[! apply(is.na(Rand_outcomes), 1, all), ] # drop rows that don't have any values
  Rand_outcomes = Rand_outcomes[row.names(outcomes),]
  
  obs = ligTags[row.names(outcomes), lig]

  al = 0.1
  wid = 2
  
  n=1
  pr = pr.curve(outcomes[(obs == T) & (!is.na(outcomes[,n])), n], outcomes[(obs == F) & (!is.na(outcomes[,n])), n], curve= T, rand.compute = T)
  allPredAUCs[n,i] = pr$auc.integral
  plot(pr$curve[,1:2], type = 'l', lwd = wid, col = alpha(ligColors[i],al),
       xlim = c(0,1), ylim = c(0,1),
       xlab = 'Recall', ylab = 'Precision', cex.lab = 1.5)
  
  rand_pr = pr.curve(Rand_outcomes[(obs == T) & (!is.na(Rand_outcomes[,n])), n], Rand_outcomes[(obs == F) & (!is.na(Rand_outcomes[,n])), n], curve= T, rand.compute = T)
  lines(rand_pr$curve[,1:2], lwd = wid, col = alpha('black',al))
  allRandAUCs[n,i] = rand_pr$auc.integral
  
  for(n in (2:ncol(outcomes))){
    pr = pr.curve(outcomes[(obs == T) & (!is.na(outcomes[,n])), n], outcomes[(obs == F) & (!is.na(outcomes[,n])), n], curve= T, rand.compute = T)
    allPredAUCs[n,i] = pr$auc.integral
    lines(pr$curve[,1:2], lwd = wid, col = alpha(ligColors[i],al))
    
    rand_pr = pr.curve(Rand_outcomes[(obs == T) & (!is.na(Rand_outcomes[,n])), n], Rand_outcomes[(obs == F) & (!is.na(Rand_outcomes[,n])), n], curve= T, rand.compute = T)
    lines(rand_pr$curve[,1:2], lwd = wid, col = alpha('black',al))
    allRandAUCs[n,i] = rand_pr$auc.integral
  }

  if (i != 12){
    title(main = bquote(atop("PR Curve" ~ - ~ bold(.(rawLigNames[i])),
                             "Mean AUC:" ~ .(round(mean(allPredAUCs[,i]),2)) ~ '(random:' ~ .(round(mean(allRandAUCs[,i]),2)) ~ ')')))
  }else{
    title(main = bquote(atop("PR Curve" ~ - ~ bold( "2" ~ alpha ~ "-Mannobiose"),
                             "Mean AUC:" ~ .(round(mean(allPredAUCs[,i]),2)) ~ '(random:' ~ .(round(mean(allRandAUCs[,i]),2)) ~ ')')))
  }

  rm(outcomes, Rand_outcomes)
  
}









# for(i in 1:ncol(ligTags)){
#   lig = gsub('(.*)_outcomes.csv', '\\1', inOutcomes)[i]
#   outcomes = read.delim(file = paste(inDir, inOutcomes[i], sep = ''), header = T, sep = ',', stringsAsFactors = F)
#   
#   boundLigs = unique(bsResiDat$iupac[ligTags[,i]]) # UniProt IDs of all lectins that do have any examples of binding
#   boundIDs = unique(bsResiDat$uniparc[bsResiDat$iupac %in% boundLigs])
#   fpBS = row.names(outcomes)[outcomes$Pred >= 0.5 & outcomes$Obs == F] # Binding sites that ended up as false positives
#   fpLecs = unique(bsResiDat$uniparc[row.names(bsResiDat) %in% fpBS]) # unique lectin uniprot ids of lectins that show up as false positives
#   fpIDs = bsResiDat$uniparc[row.names(bsResiDat) %in% fpBS] # UniProt ID for FP binding site
#   
#   negBS = row.names(outcomes)[outcomes$Obs == F]
#   negIDs = bsResiDat$uniparc[row.names(bsResiDat) %in% negBS] # UniProt IDs of lectins with negative binding sites that DO bind ligand in other structures
#   bound_NegBS = negBS[negIDs %in% boundIDs] # Negative binding sites from 
# 
#   pr = pr.curve(outcomes$Pred[outcomes$Obs == T], outcomes$Pred[outcomes$Obs == F], curve= T, rand.compute = T)
#   
#   R = test$recall[i]
#   P = test$precision[i]
#   
#   plot(pr$curve[,1:2], type = 'l', lwd = 4, col = ligColors[i],
#        xlim = c(0,1), ylim = c(0,1),
#        xlab = 'Recall', ylab = 'Precision', main = paste('PR Curve - ', lig, '\nAUC = ', as.character(round(pr$auc.integral, digits = 5)), sep = ''))
#   abline(h = pr$rand$auc.integral, lty = 1, lwd = 2)
#   lines(x = c(R,R), y = c(-1,P), lty = 2)
#   lines(x = c(-1,R), y = c(P,P), lty = 2)
#   points(R,P, pch = 19)
#   tag = pr$curve[,3] %in% outcomes[bound_NegBS, "Pred"]
#   rug((pr$curve[tag,1][pr$curve[tag,1] > R]), col = alpha('black',0.7))
#   rug((pr$curve[tag,1][pr$curve[tag,1] <= R]), col = alpha('red',0.7))
#   rug((pr$curve[tag,2][pr$curve[tag,2] < P]), col = alpha('black',0.7), side = 2)
#   rug((pr$curve[tag,2][pr$curve[tag,2] >= P]), col = alpha('red',0.7), side = 2)
#   text(x= 0.7, y = 0.9, labels = paste(round(100 * sum(fpIDs %in% boundIDs)/length(fpIDs),2), '% of ', length(fpIDs), ' FPs', sep = ''), col = 'red', cex = 1.2)
# }



#######################
# Feature importance
#######################

# rm(feats)
for(i in 1:length(dir(predDirs))){
  lig = colnames(ligTags)[as.numeric(dir(predDirs)[i])]
  cat(lig,'\n')
  curDir = paste(predDirs, dir(predDirs)[i], sep = '')
  subDirs = dir(curDir)
  for (j in 1:length(subDirs)){
    cat(j,'\t')
    inFiles = dir(paste(curDir, subDirs[j], sep = '/'))
    
    readFile = inFiles[grepl('features.csv', inFiles)] # change here to read other files
    
    tmp = read.delim(file = paste(curDir, subDirs[j], readFile, sep = '/'), header = T, sep = ',', stringsAsFactors = F)
    colnames(tmp) = rep(lig,ncol(tmp))
    if (! exists('feats')){
      feats = tmp
    }else{
      feats = cbind(feats, tmp)
    }
  }
  cat('\n')
}

featColors = rep('', nrow(feats))
resiFeats = colorRampPalette(c("plum1","tomato", "firebrick4"))(4)
pocketFeats = colorRampPalette(c('paleturquoise3', 'deepskyblue', 'mediumblue'))(4)

featColors[1:11] = 'forestgreen'

featColors[grep('^vol_4Ang$', row.names(feats)) : grep('^leftskew_10Ang$', row.names(feats))] = pocketFeats[2] # features within the d2Feats range
featColors[grepl('^vol_', row.names(feats)) | grepl('^pcntSurf_', row.names(feats))] = pocketFeats[1] # General pocket descriptors
featColors[grepl('^binnedD2', row.names(feats))] = pocketFeats[3] # PCs from the binned D2 measures
featColors[grepl('^zern', row.names(feats))] = pocketFeats[4] # PCs from the 3DZDs

featColors[grepl('^numBSresis', row.names(feats))] = resiFeats[1] # number of residues in binding site features
featColors[gsub('_bin\\d{1}', '', row.names(feats)) %in% c('H', 'B', 'E', 'G', 'T', 'S', 'X.')] = resiFeats[2] # secondary structure features
featColors[gsub('_bin\\d{1}', '', row.names(feats)) %in% c('nonpolar', 'polar', 'posCharge', 'negCharge', 'aromatic')] = resiFeats[3] # amino acid properties
featColors[grepl('^[[:upper:]]{3}_', row.names(feats)) | grepl('^CA$', row.names(feats))] = resiFeats[4] # amino acid identities

resiFeatTag = featColors %in% resiFeats
pocketFeatTag = featColors %in% pocketFeats

perc.rank <- function(x) trunc(rank(x))/length(x)

featPercentiles = apply(feats, 2, perc.rank)
RESIpercents = apply(feats[resiFeatTag,], 2, perc.rank)
POCKpercents = apply(feats[pocketFeatTag,], 2, perc.rank)
PLIPpercents = apply(feats[1:11,], 2, perc.rank)


medFeatPercentiles = as.data.frame(matrix(0, nrow = nrow(featPercentiles), ncol = length(unique(colnames(featPercentiles)))))
row.names(medFeatPercentiles) = row.names(featPercentiles)
colnames(medFeatPercentiles) = colnames(ligTags)

RESImeds = medFeatPercentiles[resiFeatTag,]
POCKmeds = medFeatPercentiles[pocketFeatTag,]
PLIPmeds = medFeatPercentiles[1:11,]


for (i in 1:length(unique(colnames(featPercentiles)))){
  lig = unique(colnames(featPercentiles))[i]
  medFeatPercentiles[,i] = apply(featPercentiles[,grepl(lig, colnames(featPercentiles))], 1, median)
  RESImeds[,i] = apply(RESIpercents[,grepl(lig, colnames(RESIpercents))], 1, median)
  POCKmeds[,i] = apply(POCKpercents[,grepl(lig, colnames(POCKpercents))], 1, median)
  PLIPmeds[,i] = apply(PLIPpercents[,grepl(lig, colnames(PLIPpercents))], 1, median)
  
}

plot(medFeatPercentiles$Fuc[order(medFeatPercentiles$Fuc, decreasing = T)], pch = 19, col = featColors[order(medFeatPercentiles$Fuc, decreasing = T)])
plot(RESImeds$Fuc[order(RESImeds$Fuc, decreasing = T)], pch = 19, col = featColors[resiFeatTag][order(RESImeds$Fuc, decreasing = T)])
plot(POCKmeds$Fuc[order(POCKmeds$Fuc, decreasing = T)], pch = 19, col = featColors[pocketFeatTag][order(POCKmeds$Fuc, decreasing = T)])
plot(PLIPmeds$Fuc[order(PLIPmeds$Fuc, decreasing = T)], pch = 19, col = featColors[1:11][order(PLIPmeds$Fuc, decreasing = T)])


par(mfrow = c(3,5))
for (i in 1:ncol(medFeatPercentiles)){
  plot(medFeatPercentiles[order(medFeatPercentiles[, i], decreasing = T), i],
       pch = 19,
       col = featColors[order(medFeatPercentiles[,i], decreasing = T)],
       ylab = 'Median importance percentile')
  title(main = ligNames[i], col.main = ligColors[i])
}

par(mfrow = c(3,5))
for (i in 1:ncol(RESImeds)){
  plot(RESImeds[order(RESImeds[, i], decreasing = T), i],
       pch = 19,
       col = featColors[resiFeatTag][order(RESImeds[,i], decreasing = T)],
       ylab = 'Median importance percentile',
       ylim = c(0,1))
  abline(h = 0.75)
  title(main = ligNames[i], col.main = ligColors[i])
}

par(mfrow = c(3,5))
for (i in 1:ncol(POCKmeds)){
  plot(POCKmeds[order(POCKmeds[, i], decreasing = T), i],
       pch = 19,
       col = featColors[pocketFeatTag][order(POCKmeds[,i], decreasing = T)],
       ylab = 'Median importance percentile',
       ylim = c(0,1))
  abline(h = 0.75)
  title(main = ligNames[i], col.main = ligColors[i])
}

par(mfrow = c(3,5))
for (i in 1:ncol(PLIPmeds)){
  plot(PLIPmeds[order(PLIPmeds[, i], decreasing = T), i],
       pch = 19,
       col = featColors[1:11][order(PLIPmeds[,i], decreasing = T)],
       ylab = 'Median importance percentile',
       ylim = c(0,1))
  abline(h = 0.75)
  title(main = ligNames[i], col.main = ligColors[i])
}

#######################
# Shared features by sig & feat importance
#######################
# Re-do and filter by median feature importance percentile (>= nth percentile)

stats = read.delim(file = './analysis/training/weightedWMW_stats.tsv', sep = '\t', stringsAsFactors = F)
# all(round(stats,4) == round(stats_weighted,4))

percentThresh = 0.75 # >= 75th percentile in feature importance in each feature class

all(gsub('_effectSize','',colnames(stats[,grepl('_effectSize$', colnames(stats))])) == colnames(medFeatPercentiles))

stratAllFeatsImp = rbind(PLIPmeds, RESImeds, POCKmeds)
all(row.names(stratAllFeatsImp) == row.names(stats))

topImpfeats = as.data.frame(stratAllFeatsImp >= percentThresh) # Logical dataframe, indicates if feature passes threshold for stratified median feature importance percentiles
all(colnames(medFeatPercentiles) == colnames(topImpfeats))


sigFeats = as.data.frame(stats[,grepl('_adj$', colnames(stats))] < 0.01) # Logical dataframe, indicates if feature passes threshold for significance from WMW test, FDR of 1%
colnames(sigFeats) = gsub('_adj$', '', colnames(sigFeats))
all(colnames(sigFeats) == colnames(topImpfeats))

# par(mfrow = c(3,5))
# for (i in 1:ncol(medFeatPercentiles)){
#   lig = colnames(medFeatPercentiles)[i]
#   
#   plot(abs(stats[,grepl('_effectSize$', colnames(stats))][,i]), medFeatPercentiles[,i],
#        xlab = '|Effect size|', ylab = 'Overall feat imp percentile', main = lig,
#        col = alpha(featColors, 0.4),
#        xlim = c(0,0.4),
#        ylim = c(0,1),
#        pch = 19)
#   par(new=T)
#   plot(abs(stats[,grepl('_effectSize$', colnames(stats))][topImpfeats[,i] & sigFeats[,i],i]), medFeatPercentiles[topImpfeats[,i] & sigFeats[,i],i],
#        xlab = '', ylab = '', main = '',
#        xlim = c(0,0.4),
#        ylim = c(0,1),
#        col = alpha(featColors[topImpfeats[,i] & sigFeats[,i]], 1),
#        pch = 19,
#        cex = 1.5)
#   text(x = 0.35,
#        y=0.05,
#        labels = paste('R=', as.character(round(cor(abs(stats[,grepl('_effectSize$', colnames(stats))][,i]), medFeatPercentiles[,i]), 3)), sep = ''))
# }

# Plot with STRATIFIED median feature importance percentiles on the y-axis
par(mfrow = c(3,5))
for (i in 1:ncol(medFeatPercentiles)){
  lig = colnames(medFeatPercentiles)[i]
  
  plot(abs(stats[,grepl('_effectSize$', colnames(stats))][,i]), stratAllFeatsImp[,i],
       xlab = '|Effect size|', ylab = 'Stratified feat imp percentile', main = ligNames[i], col.main = ligColors[i],
       col = alpha(featColors, 0.4),
       xlim = c(0,0.4),
       ylim = c(0,1),
       pch = 19)
  par(new=T)
  plot(abs(stats[,grepl('_effectSize$', colnames(stats))][topImpfeats[,i] & sigFeats[,i],i]), stratAllFeatsImp[topImpfeats[,i] & sigFeats[,i],i],
       xlab = '', ylab = '', main = '',
       xlim = c(0,0.4),
       ylim = c(0,1),
       col = alpha(featColors[topImpfeats[,i] & sigFeats[,i]], 1),
       pch = 19,
       cex = 1.5)
  text(x = 0.35,
       y=0.05,
       labels = paste('R=', as.character(round(cor(abs(stats[,grepl('_effectSize$', colnames(stats))][,i]), stratAllFeatsImp[,i]), 3)), sep = ''),
       cex = 1.4)
  abline(h = percentThresh, lty = 2)
}

# Sialic acid features
siaBindingFeats_pos = list(row.names(stats)[sigFeats$Sialic_Acid & stats$Sialic_Acid_effectSize > 0 & topImpfeats$Sialic_Acid],
                           row.names(stats)[sigFeats$NeuAc & stats$NeuAc_effectSize > 0 & topImpfeats$NeuAc],
                           row.names(stats)[sigFeats$NeuAc.a2.3.Gal.b1.4.Glc & stats$NeuAc.a2.3.Gal.b1.4.Glc_effectSize > 0 & topImpfeats$NeuAc.a2.3.Gal.b1.4.Glc])

siaBindingFeats_neg = list(row.names(stats)[sigFeats$Sialic_Acid & stats$Sialic_Acid_effectSize < 0 & topImpfeats$Sialic_Acid],
                           row.names(stats)[sigFeats$NeuAc & stats$NeuAc_effectSize < 0 & topImpfeats$NeuAc],
                           row.names(stats)[sigFeats$NeuAc.a2.3.Gal.b1.4.Glc & stats$NeuAc.a2.3.Gal.b1.4.Glc_effectSize < 0 & topImpfeats$NeuAc.a2.3.Gal.b1.4.Glc])

intersect(intersect(siaBindingFeats_pos[[1]], siaBindingFeats_pos[[2]]), siaBindingFeats_pos[[3]])
intersect(intersect(siaBindingFeats_neg[[1]], siaBindingFeats_neg[[2]]), siaBindingFeats_neg[[3]])

venn.diagram(
  x = siaBindingFeats_pos,
  category.names = c("Sialic acid" , "NeuAc" , "NeuAc(a2-3)Gal(b1-4)Glc"),
  filename = './manuscript/figures/subplots/neuac_feats_up.png',
  output=F,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#DAA520FF'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#DAA520FF',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#DAA520FF'),
  rotation = 1
)
venn.diagram(
  x = siaBindingFeats_neg,
  category.names = c("Sialic acid" , "NeuAc" , "NeuAc(a2-3)Gal(b1-4)Glc"),
  filename = './manuscript/figures/subplots/neuac_feats_down.png',
  output=F,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#DAA520FF'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#DAA520FF',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#DAA520FF'),
  rotation = 1
)


# Fucose binding
fucBindingFeats_pos = list(row.names(stats)[stats$Fuc_effectSize > 0 & sigFeats$Fuc & topImpfeats$Fuc],
                           row.names(stats)[stats$Terminal_Fucose_effectSize > 0 & sigFeats$Terminal_Fucose & topImpfeats$Terminal_Fucose])

fucBindingFeats_neg = list(row.names(stats)[stats$Fuc_effectSize < 0 & sigFeats$Fuc & topImpfeats$Fuc],
                           row.names(stats)[stats$Terminal_Fucose_effectSize < 0 & sigFeats$Terminal_Fucose & topImpfeats$Terminal_Fucose])

intersect(fucBindingFeats_pos[[1]], fucBindingFeats_pos[[2]])
intersect(fucBindingFeats_neg[[1]], fucBindingFeats_neg[[2]])

dev.off()
draw.pairwise.venn(area1 = length(fucBindingFeats_pos[[1]]),
                   area2 = length(fucBindingFeats_pos[[2]]),
                   cross.area = length(intersect(fucBindingFeats_pos[[1]], fucBindingFeats_pos[[2]])),
                   euler.d = T, scaled = T, ext.text = F, cex = 2,
                   category = c('Fucose (monosacc.)', 'Terminal fucose'), cat.cex = 2,
                   cat.pos = c(340,20),
                   cat.dist = c(0.04, 0.05),
                   fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
                   col = c("#440154ff", '#21908dff'),
                   cat.col = c("#440154ff", '#21908dff'),
                   cat.fontfamily = "sans",
                   fontfamily = "sans"
)
dev.off()
draw.pairwise.venn(area1 = length(fucBindingFeats_neg[[1]]),
                   area2 = length(fucBindingFeats_neg[[2]]),
                   cross.area = length(intersect(fucBindingFeats_neg[[1]], fucBindingFeats_neg[[2]])),
                   euler.d = T, scaled = T, ext.text = F, cex = 2,
                   category = c('Fucose (monosacc.)', 'Terminal fucose'), cat.cex = 2,
                   cat.pos = c(340,30),
                   fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
                   col = c("#440154ff", '#21908dff'),
                   cat.col = c("#440154ff", '#21908dff'),
                   cat.fontfamily = "sans",
                   fontfamily = "sans"
)

# Mannose binding
manBindingFeats_pos = list(row.names(stats)[stats$Man_effectSize > 0 & sigFeats$Man & topImpfeats$Man],
                           row.names(stats)[stats$High_Mannose_effectSize > 0 & sigFeats$High_Mannose & topImpfeats$High_Mannose],
                           row.names(stats)[stats$Man.a1.2.Man_effectSize > 0 & sigFeats$Man.a1.2.Man & topImpfeats$Man.a1.2.Man])

manBindingFeats_neg = list(row.names(stats)[stats$Man_effectSize < 0 & sigFeats$Man & topImpfeats$Man],
                           row.names(stats)[stats$High_Mannose_effectSize < 0 & sigFeats$High_Mannose & topImpfeats$High_Mannose],
                           row.names(stats)[stats$Man.a1.2.Man_effectSize < 0 & sigFeats$Man.a1.2.Man & topImpfeats$Man.a1.2.Man])

intersect(intersect(manBindingFeats_pos[[1]], manBindingFeats_pos[[2]]), manBindingFeats_pos[[3]])
intersect(intersect(manBindingFeats_neg[[1]], manBindingFeats_neg[[2]]), manBindingFeats_neg[[3]])

dev.off()
venn.diagram(
  x = manBindingFeats_pos,
  category.names = c("Man" , "High_mannose" , "Man(a1-2)Man"),
  filename = './manuscript/figures/subplots/man_feats_up.png',
  output=F,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#DAA520FF'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#DAA520FF',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#DAA520FF'),
  rotation = 1
)


venn.diagram(
  x = manBindingFeats_neg,
  category.names = c("Man" , "High_mannose" , "Man(a1-2)Man"),
  filename = './manuscript/figures/subplots/man_feats_down.png',
  output=F,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#DAA520FF'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#DAA520FF',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#DAA520FF'),
  rotation = 1
)

# Galactose binding
galBindingFeats_pos = list(row.names(stats)[stats$Gal.b1.4.Glc_effectSize > 0 & sigFeats$Gal.b1.4.Glc & topImpfeats$Gal.b1.4.Glc],
                           row.names(stats)[stats$Gal_effectSize > 0 & sigFeats$Gal & topImpfeats$Gal],
                           row.names(stats)[stats$GalNAc_effectSize > 0 & sigFeats$GalNAc & topImpfeats$GalNAc],
                           row.names(stats)[stats$Gal.b1.4.GlcNAc_effectSize > 0 & sigFeats$Gal.b1.4.GlcNAc & topImpfeats$Gal.b1.4.GlcNAc])

galBindingFeats_neg = list(row.names(stats)[stats$Gal.b1.4.Glc_effectSize < 0 & sigFeats$Gal.b1.4.Glc & topImpfeats$Gal.b1.4.Glc],
                           row.names(stats)[stats$Gal_effectSize < 0 & sigFeats$Gal & topImpfeats$Gal],
                           row.names(stats)[stats$GalNAc_effectSize < 0 & sigFeats$Gal.b1.4.GlcNAc & topImpfeats$Gal.b1.4.GlcNAc],
                           row.names(stats)[stats$Gal.b1.4.GlcNAc_effectSize < 0 & sigFeats$Gal.b1.4.GlcNAc & topImpfeats$Gal.b1.4.GlcNAc])


intersect(intersect(intersect(galBindingFeats_pos[[1]], galBindingFeats_pos[[2]]), galBindingFeats_pos[[3]]), galBindingFeats_pos[[4]])
intersect(intersect(intersect(galBindingFeats_neg[[1]], galBindingFeats_neg[[2]]), galBindingFeats_neg[[3]]), galBindingFeats_neg[[4]])

dev.off()
draw.quad.venn(area1 = length(galBindingFeats_pos[[1]]),
               area2 = length(galBindingFeats_pos[[2]]),
               area3 = length(galBindingFeats_pos[[3]]),
               area4 = length(galBindingFeats_pos[[4]]),
               
               n12 = length(intersect(galBindingFeats_pos[[1]], galBindingFeats_pos[[2]])),
               n13 = length(intersect(galBindingFeats_pos[[1]], galBindingFeats_pos[[3]])),
               n14 = length(intersect(galBindingFeats_pos[[1]], galBindingFeats_pos[[4]])),
               n23 = length(intersect(galBindingFeats_pos[[2]], galBindingFeats_pos[[3]])),
               n24 = length(intersect(galBindingFeats_pos[[2]], galBindingFeats_pos[[4]])),
               n34 = length(intersect(galBindingFeats_pos[[3]], galBindingFeats_pos[[4]])),
               
               n123 = length(intersect(intersect(galBindingFeats_pos[[1]], galBindingFeats_pos[[2]]), galBindingFeats_pos[[3]])),
               n124 = length(intersect(intersect(galBindingFeats_pos[[1]], galBindingFeats_pos[[2]]), galBindingFeats_pos[[4]])),
               n134 = length(intersect(intersect(galBindingFeats_pos[[1]], galBindingFeats_pos[[3]]), galBindingFeats_pos[[4]])),
               n234 = length(intersect(intersect(galBindingFeats_pos[[2]], galBindingFeats_pos[[3]]), galBindingFeats_pos[[4]])),
               
               n1234 = length(intersect(intersect(intersect(galBindingFeats_pos[[1]], galBindingFeats_pos[[2]]), galBindingFeats_pos[[3]]), galBindingFeats_pos[[4]])),
               
               category = c('Lac', 'Gal', 'GalNAc', 'LacNAc'), cat.cex = 2,
               cat.fontfamily = "sans",
               fontfamily = "sans",
               euler.d = T, scaled = T, ext.text = F, cex = 2,
               
               col=c("#440154ff", '#21908dff', '#DAA520FF', 'red3'),
               fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#DAA520FF',0.3), alpha('red3',0.3)),
               cat.col = c("#440154ff", '#21908dff', '#DAA520FF', 'red3')
)


dev.off()
draw.quad.venn(area1 = length(galBindingFeats_neg[[1]]),
               area2 = length(galBindingFeats_neg[[2]]),
               area3 = length(galBindingFeats_neg[[3]]),
               area4 = length(galBindingFeats_neg[[4]]),
               
               n12 = length(intersect(galBindingFeats_neg[[1]], galBindingFeats_neg[[2]])),
               n13 = length(intersect(galBindingFeats_neg[[1]], galBindingFeats_neg[[3]])),
               n14 = length(intersect(galBindingFeats_neg[[1]], galBindingFeats_neg[[4]])),
               n23 = length(intersect(galBindingFeats_neg[[2]], galBindingFeats_neg[[3]])),
               n24 = length(intersect(galBindingFeats_neg[[2]], galBindingFeats_neg[[4]])),
               n34 = length(intersect(galBindingFeats_neg[[3]], galBindingFeats_neg[[4]])),
               
               n123 = length(intersect(intersect(galBindingFeats_neg[[1]], galBindingFeats_neg[[2]]), galBindingFeats_neg[[3]])),
               n124 = length(intersect(intersect(galBindingFeats_neg[[1]], galBindingFeats_neg[[2]]), galBindingFeats_neg[[4]])),
               n134 = length(intersect(intersect(galBindingFeats_neg[[1]], galBindingFeats_neg[[3]]), galBindingFeats_neg[[4]])),
               n234 = length(intersect(intersect(galBindingFeats_neg[[2]], galBindingFeats_neg[[3]]), galBindingFeats_neg[[4]])),
               
               n1234 = length(intersect(intersect(intersect(galBindingFeats_neg[[1]], galBindingFeats_neg[[2]]), galBindingFeats_neg[[3]]), galBindingFeats_neg[[4]])),
               
               category = c('Lac', 'Gal', 'GalNAc', 'LacNAc'), cat.cex = 2,
               cat.fontfamily = "sans",
               fontfamily = "sans",
               euler.d = T, scaled = T, ext.text = F, cex = 2,
               
               col=c("#440154ff", '#21908dff', '#DAA520FF', 'red3'),
               fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#DAA520FF',0.3), alpha('red3',0.3)),
               cat.col = c("#440154ff", '#21908dff', '#DAA520FF', 'red3')
)


# Glucose binding
glcBindingFeats_pos = list(row.names(stats)[stats$Glc_effectSize > 0 & sigFeats$Glc & topImpfeats$Glc],
                           row.names(stats)[stats$GlcNAc_effectSize > 0 & sigFeats$GlcNAc & topImpfeats$GlcNAc])

glcBindingFeats_neg = list(row.names(stats)[stats$Glc_effectSize < 0 & sigFeats$Glc & topImpfeats$Glc],
                           row.names(stats)[stats$GlcNAc_effectSize < 0 & sigFeats$GlcNAc & topImpfeats$GlcNAc])

intersect(glcBindingFeats_pos[[1]], glcBindingFeats_pos[[2]])
intersect(glcBindingFeats_neg[[1]], glcBindingFeats_neg[[2]])

dev.off()
draw.pairwise.venn(area1 = length(glcBindingFeats_pos[[1]]),
                   area2 = length(glcBindingFeats_pos[[2]]),
                   cross.area = length(intersect(glcBindingFeats_pos[[1]], glcBindingFeats_pos[[2]])),
                   euler.d = T, scaled = T, ext.text = F, cex = 2,
                   category = c('Glucose', 'GlcNAc'), cat.cex = 2,
                   cat.pos = c(340,20),
                   cat.dist = c(0.04, 0.05),
                   fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
                   col = c("#440154ff", '#21908dff'),
                   cat.col = c("#440154ff", '#21908dff'),
                   cat.fontfamily = "sans",
                   fontfamily = "sans"
)
dev.off()
draw.pairwise.venn(area1 = length(glcBindingFeats_neg[[1]]),
                   area2 = length(glcBindingFeats_neg[[2]]),
                   cross.area = length(intersect(glcBindingFeats_neg[[1]], glcBindingFeats_neg[[2]])),
                   euler.d = T, scaled = T, ext.text = F, cex = 2,
                   category = c('Glucose', 'GlcNAc'), cat.cex = 2,
                   cat.pos = c(340,30),
                   fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
                   col = c("#440154ff", '#21908dff'),
                   cat.col = c("#440154ff", '#21908dff'),
                   cat.fontfamily = "sans",
                   fontfamily = "sans"
)

##################################
###### Check across all ligands
##################################

groupedLigNames = ligNames[c(1,8,9,
                             3,14,
                             2,6,12,
                             4,5,7,11,15,
                             10,13)]
groupedLigCols = ligColors[c(1,8,9,
                             3,14,
                             2,6,12,
                             4,5,7,11,15,
                             10,13)]

all(row.names(stats) == row.names(sigFeats))
all(row.names(stats) == row.names(topImpfeats))

minVal = 0.1
weakVal = 0.3
strongVal = 1

allSigFeats = as.data.frame(matrix(0, nrow = 15, ncol = nrow(stats)))
colnames(allSigFeats) = row.names(stats)
row.names(allSigFeats) = colnames(ligTags)[c(1,8,9,
                                             3,14,
                                             2,6,12,
                                             4,5,7,11,15,
                                             10,13)]

for (i in 1:nrow(allSigFeats)){
  lig = row.names(allSigFeats)[i]
  ligPattern = paste0('^', lig)
  
  up = row.names(stats)[stats[,grepl(paste0(ligPattern, '_effectSize$'), colnames(stats))] > 0 &
                          sigFeats[,grepl(paste0(ligPattern, '$'), colnames(sigFeats))] &
                          topImpfeats[,grepl(paste0(ligPattern, '$'), colnames(topImpfeats))]]
  
  down = row.names(stats)[stats[,grepl(paste0(ligPattern, '_effectSize$'), colnames(stats))] < 0 &
                            sigFeats[,grepl(paste0(ligPattern, '$'), colnames(sigFeats))] &
                            topImpfeats[,grepl(paste0(ligPattern, '$'), colnames(topImpfeats))]]
  
  weakUp = row.names(stats)[stats[,grepl(paste0(ligPattern, '_effectSize$'), colnames(stats))] > 0 &
                            (sigFeats[,grepl(paste0(ligPattern, '$'), colnames(sigFeats))] | topImpfeats[,grepl(paste0(ligPattern, '$'), colnames(topImpfeats))])]
  
  weakDown = row.names(stats)[stats[,grepl(paste0(ligPattern, '_effectSize$'), colnames(stats))] < 0 &
                            (sigFeats[,grepl(paste0(ligPattern, '$'), colnames(sigFeats))] | topImpfeats[,grepl(paste0(ligPattern, '$'), colnames(topImpfeats))])]
  
  allSigFeats[i,stats[,grepl(paste0(ligPattern, '_effectSize$'), colnames(stats))] > 0] = minVal
  allSigFeats[i,stats[,grepl(paste0(ligPattern, '_effectSize$'), colnames(stats))] < 0] = -1 * minVal
  
  allSigFeats[i,weakUp] = weakVal
  allSigFeats[i,weakDown] = -1 * weakVal
  
  allSigFeats[i,up] = strongVal
  allSigFeats[i,down] = -1 * strongVal
}

breakLst = c(-1*c(strongVal,weakVal+.1,minVal+.1), c(minVal-.05,weakVal-.1,strongVal-.1))
cols = c(alpha('royalblue', c(strongVal,weakVal,minVal)), alpha('gold2', c(minVal,weakVal,strongVal)))

pheatmap(allSigFeats, cluster_rows = F, cluster_cols = F, # all features
         color = cols, breaks = breakLst,
         gaps_row = c(3,5,7,13))

allSigFeats = allSigFeats[,apply(X = abs(allSigFeats) == 1, MARGIN = 2, FUN = any)]



grep('bin2', colnames(allSigFeats))

pheatmap(allSigFeats, cluster_rows = F, cluster_cols = F, # Sig for at least on ligand
         color = cols, breaks = breakLst,
         gaps_row = c(3,5,8,13))



allSigFeats = allSigFeats[,!colnames(allSigFeats) %in% row.names(stats)[featColors == pocketFeats[2]] & (!grepl('^vol_',colnames(allSigFeats)))]



pheatmap(allSigFeats, cluster_rows = F, cluster_cols = F, # Drop extra distribution descriptors
         color = cols, breaks = breakLst,
         gaps_row = c(3,5,8,13),
         labels_row = groupedLigNames)

allSigFeats = allSigFeats[c(1,3,
                            5,
                            6,7,8,
                            9,10,
                            14,15),]
groupedLigNames = groupedLigNames[c(1,3,5,6,7,8,9,10,14,15)]
groupedLigCols = groupedLigCols[c(1,3,5,6,7,8,9,10,14,15)]

pheatmap(allSigFeats, cluster_rows = F, cluster_cols = F, # Dropped low performing models
         color = cols, breaks = breakLst,
         # gaps_row = c(3,5,8,13),
         labels_row = groupedLigNames)


allSigFeats = allSigFeats[,apply(X = abs(allSigFeats) == 1, MARGIN = 2, FUN = any)]

grep('bin2', colnames(allSigFeats))
grep('bin4', colnames(allSigFeats))
grep('zern', colnames(allSigFeats))

colnames(allSigFeats) = gsub('^H_bin', 'a-helix_bin', colnames(allSigFeats))
colnames(allSigFeats) = gsub('^B_bin', 'b-bridge_bin', colnames(allSigFeats))
colnames(allSigFeats) = gsub('^E_bin', 'b-strand_bin', colnames(allSigFeats))
colnames(allSigFeats) = gsub('^G_bin', '3/10-helix_bin', colnames(allSigFeats))
colnames(allSigFeats) = gsub('^T_bin', 'turn_bin', colnames(allSigFeats))
colnames(allSigFeats) = gsub('^S_bin', 'bend_bin', colnames(allSigFeats))
colnames(allSigFeats) = gsub('^X._bin', 'loop_bin', colnames(allSigFeats))

cols = colorRampPalette(c("royalblue","lightblue2","white","lightgoldenrod1","goldenrod2"))(6)
breakLst = c(-1*c(strongVal,weakVal+.1,minVal+.1), c(minVal-.1,weakVal-.1,strongVal-.1))

# bk = allSigFeats
allSigFeats = bk

## Manually cluster groups of columns
# CLuster PLIP columns
clus = hclust(d = dist(t(allSigFeats[,1:5]), method = 'euc'), method = 'average')
allSigFeats = allSigFeats[,c(clus$order, 6:ncol(allSigFeats))]

# CLuster resi cnt columns
strt = 6
stp = 9
clus = hclust(d = dist(t(allSigFeats[,strt:stp]), method = 'euc'), method = 'average')
allSigFeats = allSigFeats[,c(1:(strt-1), ((strt-1)+clus$order), (stp+1):ncol(allSigFeats))]


# CLuster resi bin1 columns
strt = 10
stp = 26
clus = hclust(d = dist(t(allSigFeats[,strt:stp]), method = 'euc'), method = 'average')
allSigFeats = allSigFeats[,c(1:(strt-1), ((strt-1)+clus$order), (stp+1):ncol(allSigFeats))]


# CLuster resi bin2 columns
strt = 27
stp = 32
clus = hclust(d = dist(t(allSigFeats[,strt:stp]), method = 'euc'), method = 'average')
allSigFeats = allSigFeats[,c(1:(strt-1), ((strt-1)+clus$order), (stp+1):ncol(allSigFeats))]


# CLuster resi bin3 columns
strt = 33
stp = 41
clus = hclust(d = dist(t(allSigFeats[,strt:stp]), method = 'euc'), method = 'average')
allSigFeats = allSigFeats[,c(1:(strt-1), ((strt-1)+clus$order), (stp+1):ncol(allSigFeats))]

# CLuster resi bin4 columns
strt = 42
stp = 52
clus = hclust(d = dist(t(allSigFeats[,strt:stp]), method = 'euc'), method = 'average')
allSigFeats = allSigFeats[,c(1:(strt-1), ((strt-1)+clus$order), (stp+1):ncol(allSigFeats))]

# Cluster pocket columns
strt = 54
stp = 58
clus = hclust(d = dist(t(allSigFeats[,strt:stp]), method = 'euc'), method = 'average')
allSigFeats = allSigFeats[,c(1:(strt-1), ((strt-1)+clus$order), (stp+1):ncol(allSigFeats))]

pheatmap(allSigFeats, cluster_rows = F, cluster_cols = F, # Drop empty columns after pruning rows
         color = cols,
         breaks = breakLst,
         border_color = 'black',
         cellwidth = 14,
         cellheight = 20,
         gaps_row = c(2,3,6,8),
         gaps_col = c(5,9,26,32,41,52,53,58),
         labels_row = groupedLigNames,
         fontsize_col = 10,
         angle_col = 45)

## Cluster ligand groups by similarity
fig4_dists = as.data.frame(matrix(0,nrow = 5, ncol = ncol(allSigFeats)))
row.names(fig4_dists) = c('NeuAc', 'Fuc', 'Man', 'Gal', 'Glc')

fig4_dists[1,] = apply(allSigFeats[1:2,], 2, mean)
fig4_dists[2,] = allSigFeats[3,]
fig4_dists[3,] = apply(allSigFeats[4:6,], 2, mean)
fig4_dists[4,] = apply(allSigFeats[7:8,], 2, mean)
fig4_dists[5,] = apply(allSigFeats[9:10,], 2, mean)

clus = hclust(d = dist(fig4_dists, method = 'euc'), method = 'average')
clus$labels[clus$order]

allSigFeats = allSigFeats[c(3, # Fuc
                            1:2, # NeuAc
                            7:8, # Gal
                            4:6, # Man
                            9:10),] # Glc
groupedLigNames = groupedLigNames[c(3, # Fuc
                                    1:2, # NeuAc
                                    7:8, # Gal
                                    4:6, # Man
                                    9:10)] # Glc

## Cluster rows w/in Man group
clus = hclust(d = dist(allSigFeats[6:8,], method = 'euc'), method = 'average')
clus$labels[clus$order]

allSigFeats = allSigFeats[c(1:5,
                            7,6,8, # reordered rows with mannose
                            9:10),]
groupedLigNames = groupedLigNames[c(1:5,
                     7,6,8,
                     9:10)]

## Add column annotation
annot <- data.frame(Feature_Type = rep("", ncol(allSigFeats)))
row.names(annot) = colnames(allSigFeats)

names(featColors) = row.names(stats)

names(featColors) = gsub('^H_bin', 'a-helix_bin', names(featColors))
names(featColors) = gsub('^B_bin', 'b-bridge_bin', names(featColors))
names(featColors) = gsub('^E_bin', 'b-strand_bin', names(featColors))
names(featColors) = gsub('^G_bin', '3/10-helix_bin', names(featColors))
names(featColors) = gsub('^T_bin', 'turn_bin', names(featColors))
names(featColors) = gsub('^S_bin', 'bend_bin', names(featColors))
names(featColors) = gsub('^X._bin', 'loop_bin', names(featColors))

fig4_feat_colors = featColors[colnames(allSigFeats)]

annot$Feature_Type[fig4_feat_colors == 'forestgreen'] <- 'PLIP interaction counts'

annot$Feature_Type[grepl('^vol_', colnames(allSigFeats)) | grepl('^pcntSurf_', colnames(allSigFeats))] = 'Pocket descriptors' # General pocket descriptors
annot$Feature_Type[grepl('^binnedD2', colnames(allSigFeats))] <- 'D2 Principal Components'
annot$Feature_Type[grepl('^zern', colnames(allSigFeats))] <- '3DZD Principal Components'

annot$Feature_Type[grepl('^numBSresis', colnames(allSigFeats))] <- 'Residue counts/bin'
annot$Feature_Type[gsub('_bin\\d{1}', '', colnames(allSigFeats)) %in% c('a-helix', 'b-bridge', 'b-strand', '3/10-helix', 'turn', 'bend', 'loop')] <- 'Residue sec struct.'
annot$Feature_Type[gsub('_bin\\d{1}', '', colnames(allSigFeats)) %in% c('nonpolar', 'polar', 'posCharge', 'negCharge', 'aromatic')] <- 'Amino acid property counts'
annot$Feature_Type[grepl('^[[:upper:]]{3}_', colnames(allSigFeats)) | grepl('^CA$', colnames(allSigFeats))] <- 'Residue identities'

annot$Feature_Type = as.character(annot$Feature_Type)
annot$Feature_Type <- factor(annot$Feature_Type, levels = unique(annot$Feature_Type))

annot$Bin = "NA"
annot$Bin[grepl('_bin1$', row.names(annot))] = 'Bin 1'
annot$Bin[grepl('_bin2$', row.names(annot))] = 'Bin 2'
annot$Bin[grepl('_bin3$', row.names(annot))] = 'Bin 3'
annot$Bin[grepl('_bin4$', row.names(annot))] = 'Bin 4'

annot$Bin <- factor(annot$Bin, levels = c('Bin 1', 'Bin 2', 'Bin 3', 'Bin 4', 'NA'))

all(colnames(allSigFeats) == row.names(annot)) 

Feature_Type <- unique(fig4_feat_colors)
names(Feature_Type) <- levels(annot$Feature_Type)

Bin <- c('firebrick3', 'darkorange2', 'darkgoldenrod2', 'gold2','grey75')
names(Bin) <- levels(annot$Bin)

r_annot = data.frame(Gly_Type = rep("", nrow(allSigFeats)))
row.names(r_annot) = row.names(allSigFeats)

r_annot$Gly_Type[1] = 'Fucose'
r_annot$Gly_Type[2:3] = 'Sialic Acid'
r_annot$Gly_Type[4:5] = 'Galactose'
r_annot$Gly_Type[6:8] = 'Mannose'
r_annot$Gly_Type[9:10] = 'Glucose'

r_annot$Gly_Type = factor(r_annot$Gly_Type, levels = unique(r_annot$Gly_Type))

Gly_Type = c('red1', 'purple2', 'goldenrod2', 'forestgreen', 'royalblue2')
names(Gly_Type) = levels(r_annot$Gly_Type)

annot_cols <- list(Feature_Type = Feature_Type, Bin = Bin, Gly_Type = Gly_Type)

###
# FIGURE 4
###
pdf(file = paste('./manuscript/figures/subplots/', 
                 'important_feats',
                 '.pdf', sep = ''),
    width = 18,
    height = 7)
pheatmap(allSigFeats, cluster_rows = F, cluster_cols = F, # Drop empty columns after pruning rows
         color = cols,
         breaks = breakLst,
         border_color = 'black',
         cellwidth = 14,
         cellheight = 20,
         gaps_row = c(1,3,5,8),
         gaps_col = c(5,9,26,32,41,52,53,58),
         annotation_col = annot,
         annotation_row = r_annot,
         annotation_colors = annot_cols,
         labels_row = groupedLigNames,
         fontsize_col = 10,
         angle_col = 45)
dev.off()

plot(0,0, axes = F, xlab = '', ylab = '', col = 'white')
legend(x = 'center',
       pch = 22,
       pt.cex = 2,  
       col = c('white',
                rep('black',3),
                'white','white',
               rep('black',3)),
       pt.bg = c('white',
               cols[1:3],
               'white','white',
               cols[4:6]),
       legend = c(expression(bold('Depleted Features')),
                  'High Imp RF & Stat Sig',
                  'High Imp RF XOR Stat Sig',
                  'Low Imp RF & Non Sig',
                  '',
                  expression(bold('Enriched Features')),
                  'Low Imp RF & Non Sig',
                  'High Imp RF XOR Stat Sig',
                  'High Imp RF & Stat Sig'))




###
## Fig 0 barplots
###
scaledFeats = predFeats
for(i in 1:ncol(scaledFeats)){
  scaledFeats[,i] = (scaledFeats[,i] - min(scaledFeats[,i])) / (max(scaledFeats[,i]) - min(scaledFeats[,i]))
}

fig0_feats = c('hbonds','numBSresis_bin3', 'E_bin4', 'negCharge_bin4', 'CA', 'pcntSurf_4Ang', 'var_4Ang', 'binnedD2_PC3', 'zern_PC5')

fig0_means = as.data.frame(matrix(0, nrow = length(fig0_feats), ncol = 6))
row.names(fig0_means) = fig0_feats
colnames(fig0_means) = c('lac', 'no_lac',
                         'man', 'no_man',
                         'fuc', 'no_fuc')

fig0_means$lac = apply(scaledFeats[ligTags$Gal.b1.4.Glc, fig0_feats], 2, mean)
fig0_means$no_lac = apply(scaledFeats[! ligTags$Gal.b1.4.Glc, fig0_feats], 2, mean)

fig0_means$man = apply(scaledFeats[ligTags$High_Mannose, fig0_feats], 2, mean)
fig0_means$no_man = apply(scaledFeats[! ligTags$High_Mannose, fig0_feats], 2, mean)

fig0_means$fuc = apply(scaledFeats[ligTags$Fuc, fig0_feats], 2, mean)
fig0_means$no_fuc = apply(scaledFeats[! ligTags$Fuc, fig0_feats], 2, mean)

barplot(as.matrix(t(fig0_means[,1:2])), beside = T, col = c('gold2', 'grey60'), ylim = c(0,1))

barplot(as.matrix(t(fig0_means[,3:4])), beside = T, col = c('forestgreen', 'grey60'), ylim = c(0,1))

barplot(as.matrix(t(fig0_means[,5:6])), beside = T, col = c('red2', 'grey60'), ylim = c(0,.7))


