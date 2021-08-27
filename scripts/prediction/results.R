
library(reshape)
library(ggplot2)
library(philentropy)
library(pheatmap)
library(PRROC)
library(vioplot)
library(VennDiagram)
library(survey)
library(Cairo)

library(umap)

# devtools::install_github("copenhagencenterforglycomics/ggsugar")
# library(ggsugar)
# library(V8)


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

load('./analysis/training/surveyObject.RData')

####
ligNames = colnames(ligTags)
ligNames = gsub('_', ' ', ligNames)

ligNames[13] = expression(bold(paste("2", alpha, "-Mannobiose", sep = '')))

ligNames[1] = expression(bold("Terminal NeuAc Group"))
ligNames[2] = expression(bold("High Mannose Group"))
ligNames[3] = expression(bold("Terminal Fuc Group"))
ligNames[4] = expression(bold("Lactose"))
ligNames[5] = expression(bold("Galactose"))
ligNames[6] = expression(bold("Mannose"))
ligNames[7] = expression(bold("N-Acetylgalactosamine"))
ligNames[8] = expression(bold("N-Acetylneuraminic Acid"))
ligNames[9] = expression(bold("3'-Siayllactose"))
ligNames[10] = expression(bold("Glucose"))
ligNames[11] = expression(bold("N-Acetylglucosamine"))
ligNames[12] = expression(bold("N-Acetyllactosamine"))
ligNames[14] = expression(bold("Fucose"))
ligNames[15] = expression(bold("TF Antigen"))


ligColors = rep('', ncol(ligTags))
ligColors[1] = 'purple2' # Sialic Acid
ligColors[8] = 'darkviolet' # NeuAc monosacc.
ligColors[9] = 'purple2' # 3' Sialyllactose

ligColors[2] = 'forestgreen' # High mannose
ligColors[6] = 'darkgreen' # Mannose monosacc.
ligColors[13] = 'forestgreen' # 2alpha mannobiose

ligColors[3] = 'red1' # Terminal Fuc
ligColors[14] = 'firebrick3' # Fuc monosacc.

ligColors[4] = 'goldenrod2' # Lactose
ligColors[5] = 'darkgoldenrod3' # Gal monosacc.
ligColors[7] = 'darkgoldenrod3' # GalNAc (Tn antigen)
ligColors[12] = 'goldenrod2' # N-Acetyllactosamine (LacNAc)
ligColors[15] = 'goldenrod2' # TF antigen

ligColors[10] = 'mediumblue' # Glc monosacc.
ligColors[11] = 'royalblue2' # GlcNAc
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


scaledFeats = predFeats
for(i in 1:ncol(scaledFeats)){
  scaledFeats[,i] = (scaledFeats[,i] - min(scaledFeats[,i])) / (max(scaledFeats[,i]) - min(scaledFeats[,i]))
}

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
  if(j == 2){ # Swap Term Fuc and High Man
    i = 5
  }else if(j ==3){
    i = 3
  }
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
     labels = ligNames[c(1,3,2,4:15)],
     col = ligColors[c(1,3,2,4:15)],
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965,
     ## Increase label size.
     cex = 1.2)

# Copied @ 1400 x 1000


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
     ylim = c(0,1))
i<-1 
vioplot(testDat$f2[(testDat$ligand == colnames(ligTags)[i]) & (testDat$mode == 'pred')], at = i, side = 'left',
        add=T,
        col = ligColors[i],
        xlim = c(0,16), ylim = c(-1,1),
        plotCentre = "line")
vioplot(testDat$f2[(testDat$ligand == colnames(ligTags)[i]) & (testDat$mode == 'rand')], at = i, side = 'right',
        add = T,
        col = alpha(ligColors[i],0.33),
        xlim = c(0,16), ylim = c(-1,1),
        plotCentre = "line")
for(i in 2:ncol(ligTags)){
  vioplot(testDat$f2[(testDat$ligand == colnames(ligTags)[i]) & (testDat$mode == 'pred')], at = i, side = 'left',
          col = ligColors[i],
          xlim = c(0,16), ylim = c(-1,1),
          add = T,
          plotCentre = "line")
  vioplot(testDat$f2[(testDat$ligand == colnames(ligTags)[i]) & (testDat$mode == 'rand')], at = i, side = 'right',
          add = T,
          col = alpha(ligColors[i],0.3),
          xlim = c(0,16), ylim = c(-1,1),
          plotCentre = "line")
}
axis(side=1,at=1:15, labels = F)
axis(side=2,at=seq.int(-1,1,0.5))
title(xlab = "", ylab = "F2 Score", main = "5x CV + LOCO\nValidation F2 Scores compared to random classifiers")
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

cor.test(trainDat$sampSizes[trainDat$mode == 'pred'], trainDat$f2[trainDat$mode == 'pred'])


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
text(x = 110, y = -0.8, labels = 'Pearson corr: 0.55 (p<0.001)', cex = 1.3)

CVperf_melt = melt(trainDat[trainDat$mode == 'pred',], id.vars = c("ligand", "sampSizes", "kappa", "f2"))

# CVperf_melt = CVperf_melt[sample(x = (1:nrow(CVperf_melt)), size = nrow(CVperf_melt), replace = F),]

CVperf_ligMeans = as.data.frame(matrix(0, nrow = 15, ncol = 4))
colnames(CVperf_ligMeans) = colnames(CVperf_melt)[1:4]

CVperf_ligMeans$ligand = unique(CVperf_melt$ligand)
for(i in 1:nrow(CVperf_ligMeans)){
  gly = CVperf_ligMeans$ligand[i]
  CVperf_ligMeans$sampSizes[i] = mean(CVperf_melt$sampSizes[CVperf_melt$ligand == gly])
  CVperf_ligMeans$kappa[i] = mean(CVperf_melt$kappa[CVperf_melt$ligand == gly])
  CVperf_ligMeans$f2[i] = mean(CVperf_melt$f2[CVperf_melt$ligand == gly])
}




CVperf_melt$ligand = factor(CVperf_melt$ligand)

levels(CVperf_melt$ligand)

meltLigCols = rep('', ncol(ligTags))
meltLigCols[14] = 'purple2' # Sialic Acid
meltLigCols[12] = 'darkviolet' # NeuAc monosacc.
meltLigCols[13] = 'purple2' # 3' Sialyllactose

meltLigCols[9] = 'forestgreen' # High mannose
meltLigCols[10] = 'darkgreen' # Mannose monosacc.
meltLigCols[11] = 'forestgreen' # 2alpha mannobiose

meltLigCols[15] = 'red1' # Terminal Fuc
meltLigCols[1] = 'firebrick3' # Fuc monosacc.

meltLigCols[4] = 'goldenrod2' # Lactose
meltLigCols[2] = 'darkgoldenrod3' # Gal monosacc.
meltLigCols[6] = 'darkgoldenrod3' # GalNAc (Tn antigen)
meltLigCols[5] = 'goldenrod2' # N-Acetyllactosamine (LacNAc)
meltLigCols[3] = 'goldenrod2' # TF antigen

meltLigCols[7] = 'mediumblue' # Glc monosacc.
meltLigCols[8] = 'royalblue2' # GlcNAc

p = ggplot(data = CVperf_melt, aes(x = sampSizes, y = f2, col = ligand))

p = p +
  geom_point(size = 5) +
  theme_linedraw(base_size = 22) +
  scale_color_manual(values = alpha(meltLigCols, 0.3), guide = "none") +
  labs(x = "Number of interactions for training", y = "Training F2 score") +
  xlim(0, 125) + ylim(0,1)

p



CVperf_ligMeans$ligand = factor(CVperf_ligMeans$ligand, levels(CVperf_melt$ligand))
levels(CVperf_ligMeans$ligand)
row.names(CVperf_ligMeans) = CVperf_ligMeans$ligand
CVperf_ligMeans = CVperf_ligMeans[levels(CVperf_ligMeans$ligand),]

# mapping = rep(0,15)
# for (i in 1:15){
#   gly = levels(CVperf_ligMeans$ligand)[i]
#   gly = paste0('^', gly, '$')
#   mapping[i] = grep(gly, CVperf_ligMeans$ligand)
# }
# 
# meltLigCols = meltLigCols[c()]

cor.test(CVperf_melt$sampSizes, CVperf_melt$f2)

p + geom_label(data = CVperf_ligMeans, aes(x=sampSizes, y = f2, label = ligand), colour = meltLigCols) +
  annotate("text", label = 'Pearson corr: 0.39 (p<0.001)', x = 100, y = 0.2, size = 6, colour = "black")


p + 
  geom_text(data = CVperf_ligMeans, aes(x=sampSizes, y = f2, label = ligand), colour = meltLigCols, size = 8.02, family = "mono", fontface = 'bold') +
  geom_text(data = CVperf_ligMeans, aes(x=sampSizes, y = f2, label = ligand), colour = 'black', size = 7.98, family = "mono") +
  annotate("text", label = 'Pearson corr: 0.39 (p<0.001)', x = 100, y = 0.2, size = 6, colour = "black")
  

# p + geom_sugar(data = CVperf_ligMeans, aes(x=sampSizes, y = f2),
#                sugar = rep('n-linked-sialyl-biantennary', 15),
#                # sugar=c('fuc',
#                #         'gal',
#                #         # 'gal(b1-3)galnac',
#                #         'gal',
#                #
#                #         # 'gal(b1-4)glc',
#                #         'gal',
#                #
#                #         # 'gal(b1-4)glcnac',
#                #         'gal',
#                #
#                #         'galnac',
#                #         'glc',
#                #         'glcnac',
#                #         'high_man',
#                #         'man',
#                #         # 'man(a1-2)man',
#                #         'man',
#                #
#                #         'neuac',
#                #         # 'neuac(a2-3)gal(b1-4)glc',
#                #         'neuac',
#                #
#                #         'n-linked-sialyl-biantennary',
#                #         'fuc'),
#                size = 4)



# ggplot(data = mTest, aes(x = metric, y = value, col = ligand, fill = metric)) +
#   geom_point(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.4) +
#   geom_boxplot(outlier.alpha = 0) +
#   ylim(c(-0.15, 1)) +
#   scale_fill_manual(values = alpha(rep('snow3',length(unique(mTest$metric))), 0.6), guide =F) +
#   scale_color_manual(values = ligColors) +
#   labs(title = '5x CV with LOCO validation - Validation Performance', x = "Metric type", y = "Metric value") +
#   theme_light(base_size = 22) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))



# p <- ggplot(data.frame(x=(1:3),y=rep(1,3),sugar=c('o-fuc','STn','n-linked-sialyl-biantennary'))) + 
#   geom_sugar(aes(x,y,sugar=sugar),size=10,align="centre") + 
#   theme_minimal()
# p


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



mean(testDat$recall[testDat$mode == 'pred'])
mean(testDat$prec[testDat$mode == 'pred'])

mean(testDat$recall[testDat$mode == 'rand'])
mean(testDat$prec[testDat$mode == 'rand'])

mean(trainDat$recall[trainDat$mode == 'pred'])
mean(trainDat$prec[trainDat$mode == 'pred'])

mean(trainDat$recall[trainDat$mode == 'rand'])
mean(trainDat$prec[trainDat$mode == 'rand'])


summPerf = as.data.frame(matrix(0,nrow = 15, ncol = 4))
row.names(summPerf) = colnames(ligTags)
colnames(summPerf) = c('recall', 'precision', 'rand_recall', 'rand_prec')

for (i in 1:nrow(summPerf)){
  lig = row.names(summPerf)[i]
  mets = c(median(testDat$recall[testDat$mode == 'pred' & testDat$ligand == lig]),
           median(testDat$prec[testDat$mode == 'pred' & testDat$ligand == lig]),
           median(testDat$recall[testDat$mode == 'rand' & testDat$ligand == lig]),
           median(testDat$prec[testDat$mode == 'rand' & testDat$ligand == lig]))
  summPerf[i,] = mets
}

#########################
## Figure 2
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


###
###
# Alt Fig 2 with sample size boxplots dropped out into a separate panel
###

# Panel A -- P & R same plot with violin & box plots
xLim = c(0,30)
yLim = c(0.2,1)

par(mar = c(8.6, 5.1, 5.1, 4.1), # change the margins
    lwd = 2, # increase the line thickness
    cex.axis = 1.2 # increase default axis label size
)
# par(cex.lab=1.5, mar = c(5, 4, 4, 4) + 0.3)
plot(0,0,col = 'white', xlab = '',ylab = '',type="n",
     axes=FALSE,ann=FALSE,
     xlim = xLim,
     ylim = yLim)
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
plot(0,0,col = 'white', xlab = '',ylab = '',type="n",
     axes=FALSE,ann=FALSE,
     xlim = xLim,
     ylim = yLim)
axis(side=1,at=seq.int(1,29,2), labels = F)
axis(side=2,at=pretty(c(0.2,1)), las = 2)
axis(side = 4, at = pretty(c(0.2,1)), las = 2)  # Add second axis
# abline(v=6, lwd = 0.75)
mtext("Precision", side = 4, line = 3, cex = 2, col = alpha('black',0.85))  # Add second axis label
title(xlab = "", ylab = "Recall", main = "Random Forest - LO(C)O Validation performance\nPrecision & Recall vs random", cex.main = 1.5, cex.lab = 2)
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

# copy to clipboard @ 1800x800

# Panel B - Boxplots of number of binding interaction used to train models

plot(0,0, col = 'white', type = 'n',
     xlim = xLim,
     ylim = c(1.3, 4),
     axes = F, xlab = '', ylab = '')
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
axis(side=2,at=log10(logTicks), labels = logTicks, las = 2, pos = 0)
axis(side=3,at=seq.int(1,29,2), labels = F, pos = 2.2)


###
###
# Alt Fig 2 with sample size boxplots dropped out into a supp fig panel
## Diag labels
###

dev.off()
pdf(file = paste('./manuscript/figures/subplots/', 
                 'RF-violins',
                 '.pdf', sep = ''),
    width = 20,
    height = 12)#  P & R same plot with violin & box plots
xLim = c(0,30)
yLim = c(0.2,1)

par(mar = c(16.1, 9.1, 5.1, 4.1), # change the margins
    lwd = 2, # increase the line thickness
    cex.axis = 1.2 # increase default axis label size
)
# par(cex.lab=1.5, mar = c(5, 4, 4, 4) + 0.3)
plot(0,0,col = 'white', xlab = '',ylab = '',type="n",
     axes=FALSE,ann=FALSE,
     xlim = xLim,
     ylim = yLim)
abline(h = c(0.2,0.4,0.6,0.8,1.0), lwd = 1, col = 'grey80')
abline(h = c(0.3,0.5,0.7,0.9), lwd = 1, lty = 2, col = 'grey80')
for(j in 1:ncol(ligTags)){
  i = ((j-1) * 2) + 1
  if(j == 2){ # Swap Term Fuc and High Man
    i = 5
  }else if(j ==3){
    i = 3
  }
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
plot(0,0,col = 'white', xlab = '',ylab = '',type="n", # Add new plot to simplify labels/axes
     axes=FALSE,ann=FALSE,
     xlim = xLim,
     ylim = yLim)
axis(side=1,at=seq.int(1,29,2), labels = F)
axis(side=2,at=pretty(c(0.2,1)), las = 2)
axis(side = 4, at = pretty(c(0.2,1)), las = 2)  # Add second axis
# abline(v=6, lwd = 0.75)
mtext("Precision", side = 4, line = 3, cex = 2, col = alpha('black',0.85))  # Add second axis label
title(xlab = "", ylab = "Recall", main = "Random Forest Validation - Precision & Recall (Violin plots) \n Versus Random Classifier (Boxplots)", cex.main = 2.5, cex.lab = 2)
text(x = seq.int(1,29,2),
     y = par("usr")[3] - 0.05,
     labels = ligNames[c(1,3,2,4:15)],
     col = ligColors[c(1,3,2,4:15)],
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965,
     ## Increase label size.
     cex = 1.8)
dev.off()


# Panel B, (supp fig 4) - Boxplots of number of binding interaction used to train models
par(mar = c(16.1, 9.1, 5.1, 4.1), # change the margins
    lwd = 2, # increase the line thickness
    cex.axis = 1.2 # increase default axis label size
)
plot(0,0, col = 'white', type = 'n',
     xlim = xLim,
     ylim = c(1.3, 4),
     axes = F, xlab = '', ylab = '')
lines(x = c(1,30), y = c(0,0), lwd = 2)
logTicks = c(20,50,100,150)

for (i in 1:length(logTicks)){
  segments(x0=0.25,y0= log10(logTicks)[i], x1 = 30, y1 = log10(logTicks)[i], lwd = 1, col = 'grey80')
}


for(j in 1:ncol(ligTags)){
  i = ((j-1) * 2) + 1
  if(j == 2){ # Swap Term Fuc and High Man
    i = 5
  }else if(j ==3){
    i = 3
  }
  par(new = T)
  boxplot(log10(trainDat$sampSizes[trainDat$ligand == colnames(ligTags)[j] & trainDat$mode == 'pred']),
          at = i,
          col = ligColors[j], border = ligColors[j],
          boxwex = 2.5,
          xlim = xLim,
          ylim = c(1.3, 4),
          axes = F, xlab = '', ylab ='')
}
axis(side=2,at=log10(logTicks), labels = logTicks, las = 2, pos = 0)
axis(side=1,at=seq.int(1,29,2), labels = F, pos = 1.2)
text(x = seq.int(1,29,2)+.75,
     y = par("usr")[3] - 0.05,
     labels = ligNames[c(1,3,2,4:15)],
     col = ligColors[c(1,3,2,4:15)],
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 35,
     ## Adjust the labels to almost 100% right-justified.
     adj = 1.1,
     ## Increase label size.
     cex = 1.2)

# copied at 1200 x 800

#######################
# PR curves
#######################

# rawLigNames = rep('',length(ligNames))
# rawLigNames[13] = "2 alpha-Mannobiose"
# 
# rawLigNames[1] = "Terminal NeuAc Group"
# rawLigNames[2] = "High Mannose Group"
# rawLigNames[3] = "Terminal Fuc Group"
# rawLigNames[4] = "Lactose"
# rawLigNames[5] = "Galactose"
# rawLigNames[6] = "Mannose"
# rawLigNames[7] = "N-Acetylnalactosamine"
# rawLigNames[8] = "N-Acetylneuraminic Acid"
# rawLigNames[9] = "3'-Siayllactose"
# rawLigNames[10] = "Glucose"
# rawLigNames[12] = "N-Acetyllactosamine"
# rawLigNames[11] = "N-Acetylglucosamine"
# rawLigNames[14] = "Fucose"
# rawLigNames[15] = "TF Antigen"
# 
# 
# par(mfrow = c(3,5))
# 
# ligDirs = dir(predDirs)  
# ligDirs  = ligDirs[order(as.numeric(ligDirs))]
# 
# allPredAUCs = as.data.frame(matrix(0, nrow = 100, ncol = ncol(ligTags)))
# colnames(allPredAUCs) = colnames(ligTags)
# 
# allRandAUCs = allPredAUCs
# 
# for(i in 9:length(ligDirs)){
#   lig = colnames(ligTags)[as.numeric(ligDirs[i])]
#   cat(lig,'\n')
#   
#   curDir = paste(predDirs, ligDirs[i], sep = '')
#   ranDir = gsub(pattern = 'pred', replacement = 'rand', x = curDir)
#   
#   subDirs = dir(curDir)
#   
#   # predAUC = randAUC = rep(0,100)
#   
#   for (j in 1:length(subDirs)){
#     cat(j, '\t')
#     #read in pred
#     inFiles = dir(paste(curDir, subDirs[j], sep = '/'))
#     
#     readFile = inFiles[grepl('outcomes.csv', inFiles)] # change here to read other files
#     
#     tmp = read.delim(file = paste(curDir, subDirs[j], readFile, sep = '/'), header = T, sep = ',', stringsAsFactors = F)
#     
#     if (! exists('outcomes')){
#       outcomes = tmp
#     }else{
#       outcomes = cbind(outcomes, tmp)
#     }
#     
#     #read in rand
#     inFiles = dir(paste(ranDir, subDirs[j], sep = '/'))
#     
#     readFile = inFiles[grepl('outcomes.csv', inFiles)] # change here to read other files
#     
#     tmp = read.delim(file = paste(ranDir, subDirs[j], readFile, sep = '/'), header = T, sep = ',', stringsAsFactors = F)
#     
#     extraRows = as.data.frame(matrix(NA, nrow = sum(!row.names(bsResiDat) %in% row.names(tmp)), ncol = 10))
#     row.names(extraRows) = row.names(bsResiDat)[!row.names(bsResiDat) %in% row.names(tmp)]
#     
#     tmp = rbind(tmp, extraRows)
#     tmp = tmp[row.names(bsResiDat),]
#     
#     if (! exists('Rand_outcomes')){
#       Rand_outcomes = tmp
# 
#     }else{
#       Rand_outcomes = cbind(Rand_outcomes, tmp)
#     }
#   }
#   cat('\n')
#   
#   Rand_outcomes = Rand_outcomes[! apply(is.na(Rand_outcomes), 1, all), ] # drop rows that don't have any values
#   Rand_outcomes = Rand_outcomes[row.names(outcomes),]
#   
#   obs = ligTags[row.names(outcomes), lig]
# 
#   al = 0.1
#   wid = 2
#   
#   n=1
#   pr = pr.curve(outcomes[(obs == T) & (!is.na(outcomes[,n])), n], outcomes[(obs == F) & (!is.na(outcomes[,n])), n], curve= T, rand.compute = T)
#   allPredAUCs[n,i] = pr$auc.integral
#   plot(pr$curve[,1:2], type = 'l', lwd = wid, col = alpha(ligColors[i],al),
#        xlim = c(0,1), ylim = c(0,1),
#        xlab = 'Recall', ylab = 'Precision', cex.lab = 1.5)
#   
#   rand_pr = pr.curve(Rand_outcomes[(obs == T) & (!is.na(Rand_outcomes[,n])), n], Rand_outcomes[(obs == F) & (!is.na(Rand_outcomes[,n])), n], curve= T, rand.compute = T)
#   lines(rand_pr$curve[,1:2], lwd = wid, col = alpha('black',al))
#   allRandAUCs[n,i] = rand_pr$auc.integral
#   
#   for(n in (2:ncol(outcomes))){
#     pr = pr.curve(outcomes[(obs == T) & (!is.na(outcomes[,n])), n], outcomes[(obs == F) & (!is.na(outcomes[,n])), n], curve= T, rand.compute = T)
#     allPredAUCs[n,i] = pr$auc.integral
#     lines(pr$curve[,1:2], lwd = wid, col = alpha(ligColors[i],al))
#     
#     rand_pr = pr.curve(Rand_outcomes[(obs == T) & (!is.na(Rand_outcomes[,n])), n], Rand_outcomes[(obs == F) & (!is.na(Rand_outcomes[,n])), n], curve= T, rand.compute = T)
#     lines(rand_pr$curve[,1:2], lwd = wid, col = alpha('black',al))
#     allRandAUCs[n,i] = rand_pr$auc.integral
#   }
# 
#   if (i != 12){
#     title(main = bquote(atop("PR Curve" ~ - ~ bold(.(rawLigNames[i])),
#                              "Mean AUC:" ~ .(round(mean(allPredAUCs[,i]),2)) ~ '(random:' ~ .(round(mean(allRandAUCs[,i]),2)) ~ ')')))
#   }else{
#     title(main = bquote(atop("PR Curve" ~ - ~ bold( "2" ~ alpha ~ "-Mannobiose"),
#                              "Mean AUC:" ~ .(round(mean(allPredAUCs[,i]),2)) ~ '(random:' ~ .(round(mean(allRandAUCs[,i]),2)) ~ ')')))
#   }
# 
#   rm(outcomes, Rand_outcomes)
#   
# }
# 
# 







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

dev.off()
plot(medFeatPercentiles$Fuc[order(medFeatPercentiles$Fuc, decreasing = T)], pch = 19, col = featColors[order(medFeatPercentiles$Fuc, decreasing = T)])
plot(RESImeds$Fuc[order(RESImeds$Fuc, decreasing = T)], pch = 19, col = featColors[resiFeatTag][order(RESImeds$Fuc, decreasing = T)])
plot(POCKmeds$Fuc[order(POCKmeds$Fuc, decreasing = T)], pch = 19, col = featColors[pocketFeatTag][order(POCKmeds$Fuc, decreasing = T)])
plot(PLIPmeds$Fuc[order(PLIPmeds$Fuc, decreasing = T)], pch = 19, col = featColors[1:11][order(PLIPmeds$Fuc, decreasing = T)])


# copied @ 1400 x 900
par(mfrow = c(3,5))
for (i in c(1,3,2,4:15)){
  plot(medFeatPercentiles[order(medFeatPercentiles[, i], decreasing = T), i],
       pch = 19,
       col = featColors[order(medFeatPercentiles[,i], decreasing = T)],
       ylab = 'Median importance percentile', xlab = 'Ranked features', cex.lab = 1.5)
  title(main = ligNames[i], col.main = ligColors[i], cex.main = 1.5)
}

par(mfrow = c(3,5))
for (i in c(1,3,2,4:15)){
  plot(RESImeds[order(RESImeds[, i], decreasing = T), i],
       pch = 19,
       col = featColors[resiFeatTag][order(RESImeds[,i], decreasing = T)],
       ylab = 'Median importance percentile', xlab = 'Ranked features', cex.lab = 1.5,
       ylim = c(0,1))
  abline(h = 0.75)
  title(main = ligNames[i], col.main = ligColors[i], cex.main = 1.5)
}

par(mfrow = c(3,5))
for (i in c(1,3,2,4:15)){
  plot(POCKmeds[order(POCKmeds[, i], decreasing = T), i],
       pch = 19,
       col = featColors[pocketFeatTag][order(POCKmeds[,i], decreasing = T)],
       ylab = 'Median importance percentile', xlab = 'Ranked features', cex.lab = 1.5,
       ylim = c(0,1))
  abline(h = 0.75)
  title(main = ligNames[i], col.main = ligColors[i], cex.main = 1.5)
}

par(mfrow = c(3,5))
for (i in c(1,3,2,4:15)){
  plot(PLIPmeds[order(PLIPmeds[, i], decreasing = T), i],
       pch = 19,
       col = featColors[1:11][order(PLIPmeds[,i], decreasing = T)],
       ylab = 'Median importance percentile', xlab = 'Ranked features', cex.lab = 1.5,
       ylim = c(0,1))
  abline(h = 0.75)
  title(main = ligNames[i], col.main = ligColors[i], cex.main = 1.5)
}

#######################
# Shared features by sig & feat importance
#######################
# Re-do and filter by median feature importance percentile (>= nth percentile)

stats = read.delim(file = './analysis/training/weightedWMW_stats.tsv', sep = '\t', stringsAsFactors = F)
# all(round(stats,4) == round(stats_weighted,4))

mean(apply(stats[,grepl('_adj$', colnames(stats))] < 0.01, 2, sum))
sd(apply(stats[,grepl('_adj$', colnames(stats))] < 0.01, 2, sum))


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
dev.off()
par(mfrow = c(3,5))
for (i in c(1,3,2,4:15)){
  lig = colnames(medFeatPercentiles)[i]
  
  plot(abs(stats[,grepl('_effectSize$', colnames(stats))][,i]), stratAllFeatsImp[,i],
       xlab = '|Effect size|', ylab = 'Stratified feat imp percentile', main = ligNames[i], col.main = ligColors[i],
       cex.main = 1.5,
       cex.lab = 1.5,
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



# Global correlation
cor.test(unlist(abs(stats[,grepl('_effectSize$', colnames(stats))])), unlist(stratAllFeatsImp))


mean(diag(cor(abs(stats[,grepl('_effectSize$', colnames(stats))][,c(1,6,14)]), stratAllFeatsImp[,c(1,6,14)])))

##############################################
### Figure 3 panels B-D
##############################################

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



Feature_Type <- unique(featColors)
names(Feature_Type) <- levels(annot$Feature_Type)
anno_colors <- list(Feature_Type = Feature_Type)

pheatmap(t(stats[,grepl('_effectSize$', colnames(stats))]),
         color = colorRampPalette(c("royalblue1", "grey90", "gold1"))(length(breakLst)),
         clustering_distance_rows = 'correlation',
         # clustering_distance_cols = 'correlation',
         display_numbers = ifelse(t(topImpfeats & sigFeats), "+", ""), fontsize_number = 18,
         border_color = 'white',
         labels_row = ligNames,
         annotation_col = annot, annotation_colors = anno_colors,
         main = '',
         breaks = breakLst,
         show_colnames = F)

# pdf(file = paste('./manuscript/figures/subplots/', 
#                  'featLegend',
#                  '.pdf', sep = ''),
#     width = 24,
#     height = 7)
# pheatmap(t(stats[1:10,grepl('_effectSize$', colnames(stats))]),
#          color = colorRampPalette(c("royalblue1", "grey90", "gold1"))(length(breakLst)),
#          clustering_distance_rows = 'correlation',
#          # clustering_distance_cols = 'correlation',
#          #display_numbers = ifelse(t(stats[,grepl('_adj$', colnames(stats))]) < 0.01, "*", ""), fontsize_number = 18,
#          cellwidth = 3,
#          border_color = 'white',
#          labels_row = ligNames,
#          annotation_col = annot, annotation_colors = anno_colors,
#          main = '',
#          breaks = breakLst,
#          show_colnames = F)
# dev.off()

# Set-up for panel A from descriptive.R

allFeatCorrs = cor(stats[,grepl('_effectSize$', colnames(stats))], stats[,grepl('_effectSize$', colnames(stats))], method = 'pearson')
row.names(allFeatCorrs)  = gsub('_effectSize$', '', row.names(allFeatCorrs))

r_annot = data.frame(Terminal_Sugar = rep("", nrow(allFeatCorrs)))
row.names(r_annot) = row.names(allFeatCorrs)

r_annot$Terminal_Sugar[ligColors %in% c('purple2', 'darkviolet')] = 'NeuAc'
r_annot$Terminal_Sugar[ligColors %in% c('forestgreen', 'darkgreen')] = 'Man'
r_annot$Terminal_Sugar[ligColors %in% c('red1', 'firebrick3')] = 'Fuc'
r_annot$Terminal_Sugar[ligColors %in% c('goldenrod2', 'darkgoldenrod3')] = 'Gal'
r_annot$Terminal_Sugar[ligColors %in% c('mediumblue', 'royalblue2')] = 'Glc'

r_annot$Terminal_Sugar = factor(r_annot$Terminal_Sugar, levels = unique(r_annot$Terminal_Sugar))

Terminal_Sugar = c('purple2', 'forestgreen', 'red1', 'goldenrod2', 'mediumblue')
names(Terminal_Sugar) = levels(r_annot$Terminal_Sugar)


r_annot$Sugar_Cnt = rep("", nrow(allFeatCorrs))
r_annot$Sugar_Cnt[c(1:3, 9)] = '3+'
r_annot$Sugar_Cnt[c(4,12,13,15)] = '2'
r_annot$Sugar_Cnt[c(5,6,7,8,10,11,14)] = '1'

r_annot$Sugar_Cnt <- factor(r_annot$Sugar_Cnt, levels = c('1', '2', '3+'))

Sugar_Cnt = colorspace::sequential_hcl(3)[3:1]
names(Sugar_Cnt) = levels(r_annot$Sugar_Cnt)


annot_cols = list(Terminal_Sugar = Terminal_Sugar, Sugar_Cnt = Sugar_Cnt)

dev.off()
pdf(file = paste('./manuscript/figures/subplots/', 
                 'allFeats_MWM_corplots',
                 '.pdf', sep = ''),
    width = 11.5,
    height = 8.25)
breakLst = seq(-1,1,0.05)
pheatmap(allFeatCorrs,
         color = colorRampPalette(c("firebrick2", "ivory", "dodgerblue3"))(length(breakLst)),
         border_color = 'white',
         cellwidth = cHeight,
         cellheight = cHeight,
         legend_breaks = c(-1, -0.5, 0, 0.5, 1, 0.97),
         legend_labels = c("-1", "-0.5", "0", "0.5", "1", "Pearson Corr\n\n"),
         labels_row = ligNames,
         main = ' ',
         breaks = breakLst,
         show_colnames = F,
         cutree_rows = 5,
         # cutree_cols = 4,
         treeheight_col = 0,
         fontsize_row = 13,
         
         annotation_row = r_annot,
         annotation_colors = annot_cols)
# grid.text(label = 'Pearson correlations between feature-specific effect sizes across ligands',x = .415, y=0.985, gp=gpar(col="black", cex = 1.5))
dev.off() 

####
# PANEL B - adjusted from descriptive.R
## Residue features only
####


cWidth = 10
cHeight = 20 


breakLst = seq(-0.5,0.5,0.01)


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

resiAnnot_cols <- list(Feature_Type = Feature_Type, Bin = Bin, Terminal_Sugar = Terminal_Sugar, Sugar_Cnt = Sugar_Cnt)

resiFeat_stats = t(resiFeat_stats[resiFeatTag,grepl('_effectSize$', colnames(resiFeat_stats))])
row.names(resiFeat_stats) = gsub('_effectSize$', '', row.names(resiFeat_stats))

dev.off()
# pdf(file = paste('./manuscript/figures/subplots/', 
#                  'resiFeats_MWM_heatmap',
#                  '.pdf', sep = ''),
#     width = 28,
#     height = 7)
CairoPDF(file = paste('./manuscript/figures/subplots/', 
                      'resiFeats_MWM_heatmap',
                      '.pdf', sep = ''),
         width = 28,
         height = 7)
pheatmap(resiFeat_stats,
         color = colorRampPalette(c("royalblue1", "ivory", "gold1"))(length(breakLst)),
         border_color = 'white',
         cellwidth = cWidth,
         cellheight = cHeight,
         clustering_distance_rows = 'correlation',
         display_numbers = ifelse(t(topImpfeats & sigFeats)[,resiFeatTag] , "\u2022", ""),
         fontsize_number = 20, number_color = 'black',
         labels_row = ligNames,
         annotation_col = resiAnnot,
         annotation_row = r_annot,
         annotation_colors = resiAnnot_cols,
         # legend_breaks = c(1,5),
         cutree_rows = 4,
         main = '',
         breaks = breakLst,
         show_colnames = T,
         fontsize_col = 7,
         angle_col = 45,
         fontsize_row = 14)
dev.off()


####
# PANEL C - adjusted from descriptive.R
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

pockAnnot_cols <- list(Feature_Type = Feature_Type, Threshold = Threshold, Terminal_Sugar = Terminal_Sugar, Sugar_Cnt = Sugar_Cnt)


pocketFeat_stats = t(stats[pocketFeatTag,grepl('_effectSize$', colnames(stats))])
row.names(pocketFeat_stats) = gsub('_effectSize$', '', row.names(pocketFeat_stats))

dev.off()

CairoPDF(file = paste('./manuscript/figures/subplots/', 
                 'pocketFeats_MWM_heatmap',
                 '.pdf', sep = ''),
    width = 24,
    height = 7)
pheatmap(pocketFeat_stats,
         color = colorRampPalette(c("royalblue1", "ivory", "gold1"))(length(breakLst)),
         border_color = 'ivory',
         cellwidth = cWidth+3,
         cellheight = cHeight,
         clustering_distance_rows = 'correlation',
         display_numbers = ifelse(t(topImpfeats & sigFeats)[,pocketFeatTag] , "\u2022", ""), fontsize_number = 20, number_color = 'black',
         labels_row = ligNames,
         annotation_col = pockAnnot,
         annotation_row = r_annot,
         annotation_colors = pockAnnot_cols,
         main = '',
         cutree_rows = 5,
         breaks = breakLst,
         show_colnames = T,
         fontsize_col = 9,
         angle_col = 45,
         fontsize_row = 12)
dev.off()


####
# PANEL D - adjusted from descriptive.R
## PLIP interaction counts
####

annot_cols = list(Terminal_Sugar = Terminal_Sugar, Sugar_Cnt = Sugar_Cnt)

plipFeat_stats = t(stats[1:11,grepl('_effectSize$', colnames(stats))])
row.names(plipFeat_stats) = gsub('_effectSize$', '', row.names(plipFeat_stats))

dev.off()

CairoPDF(file = paste('./manuscript/figures/subplots/', 
                 'plipFeats_MWM_heatmap',
                 '.pdf', sep = ''),
    width = 24,
    height = 7)
pheatmap(plipFeat_stats,
         color = colorRampPalette(c("royalblue1", "ivory", "gold1"))(length(breakLst)),
         border_color = 'ivory',
         cellwidth = cWidth+3,
         cellheight = cHeight,
         clustering_distance_rows = 'correlation',
         display_numbers = ifelse(t(topImpfeats & sigFeats)[,1:11] , "\u2022", ""), fontsize_number = 20, number_color = 'black',
         labels_row = ligNames,
         annotation_row = r_annot,
         annotation_colors = annot_cols,
         main = '',
         cutree_rows = 4,
         breaks = breakLst,
         show_colnames = T,
         fontsize_col = 9,
         angle_col = 45,
         fontsize_row = 12)
dev.off()



########################
# Venn Diagrams
########################

# # Sialic acid features
# siaBindingFeats_pos = list(row.names(stats)[sigFeats$Sialic_Acid & stats$Sialic_Acid_effectSize > 0 & topImpfeats$Sialic_Acid],
#                            row.names(stats)[sigFeats$NeuAc & stats$NeuAc_effectSize > 0 & topImpfeats$NeuAc],
#                            row.names(stats)[sigFeats$NeuAc.a2.3.Gal.b1.4.Glc & stats$NeuAc.a2.3.Gal.b1.4.Glc_effectSize > 0 & topImpfeats$NeuAc.a2.3.Gal.b1.4.Glc])
# 
# siaBindingFeats_neg = list(row.names(stats)[sigFeats$Sialic_Acid & stats$Sialic_Acid_effectSize < 0 & topImpfeats$Sialic_Acid],
#                            row.names(stats)[sigFeats$NeuAc & stats$NeuAc_effectSize < 0 & topImpfeats$NeuAc],
#                            row.names(stats)[sigFeats$NeuAc.a2.3.Gal.b1.4.Glc & stats$NeuAc.a2.3.Gal.b1.4.Glc_effectSize < 0 & topImpfeats$NeuAc.a2.3.Gal.b1.4.Glc])
# 
# intersect(intersect(siaBindingFeats_pos[[1]], siaBindingFeats_pos[[2]]), siaBindingFeats_pos[[3]])
# intersect(intersect(siaBindingFeats_neg[[1]], siaBindingFeats_neg[[2]]), siaBindingFeats_neg[[3]])
# 
# venn.diagram(
#   x = siaBindingFeats_pos,
#   category.names = c("Sialic acid" , "NeuAc" , "NeuAc(a2-3)Gal(b1-4)Glc"),
#   filename = './manuscript/figures/subplots/neuac_feats_up.png',
#   output=F,
#   imagetype="png" ,
#   height = 480 , 
#   width = 480 , 
#   resolution = 300,
#   compression = "lzw",
#   lwd = 1,
#   col=c("#440154ff", '#21908dff', '#DAA520FF'),
#   fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#DAA520FF',0.3)),
#   cex = 0.5,
#   fontfamily = "sans",
#   cat.cex = 0.3,
#   cat.default.pos = "outer",
#   cat.pos = c(-27, 27, 135),
#   cat.dist = c(0.055, 0.055, 0.085),
#   cat.fontfamily = "sans",
#   cat.col = c("#440154ff", '#21908dff', '#DAA520FF'),
#   rotation = 1
# )
# venn.diagram(
#   x = siaBindingFeats_neg,
#   category.names = c("Sialic acid" , "NeuAc" , "NeuAc(a2-3)Gal(b1-4)Glc"),
#   filename = './manuscript/figures/subplots/neuac_feats_down.png',
#   output=F,
#   imagetype="png" ,
#   height = 480 , 
#   width = 480 , 
#   resolution = 300,
#   compression = "lzw",
#   lwd = 1,
#   col=c("#440154ff", '#21908dff', '#DAA520FF'),
#   fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#DAA520FF',0.3)),
#   cex = 0.5,
#   fontfamily = "sans",
#   cat.cex = 0.3,
#   cat.default.pos = "outer",
#   cat.pos = c(-27, 27, 135),
#   cat.dist = c(0.055, 0.055, 0.085),
#   cat.fontfamily = "sans",
#   cat.col = c("#440154ff", '#21908dff', '#DAA520FF'),
#   rotation = 1
# )
# 
# 
# # Fucose binding
# fucBindingFeats_pos = list(row.names(stats)[stats$Fuc_effectSize > 0 & sigFeats$Fuc & topImpfeats$Fuc],
#                            row.names(stats)[stats$Terminal_Fucose_effectSize > 0 & sigFeats$Terminal_Fucose & topImpfeats$Terminal_Fucose])
# 
# fucBindingFeats_neg = list(row.names(stats)[stats$Fuc_effectSize < 0 & sigFeats$Fuc & topImpfeats$Fuc],
#                            row.names(stats)[stats$Terminal_Fucose_effectSize < 0 & sigFeats$Terminal_Fucose & topImpfeats$Terminal_Fucose])
# 
# intersect(fucBindingFeats_pos[[1]], fucBindingFeats_pos[[2]])
# intersect(fucBindingFeats_neg[[1]], fucBindingFeats_neg[[2]])
# 
# dev.off()
# draw.pairwise.venn(area1 = length(fucBindingFeats_pos[[1]]),
#                    area2 = length(fucBindingFeats_pos[[2]]),
#                    cross.area = length(intersect(fucBindingFeats_pos[[1]], fucBindingFeats_pos[[2]])),
#                    euler.d = T, scaled = T, ext.text = F, cex = 2,
#                    category = c('Fucose (monosacc.)', 'Terminal fucose'), cat.cex = 2,
#                    cat.pos = c(340,20),
#                    cat.dist = c(0.04, 0.05),
#                    fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
#                    col = c("#440154ff", '#21908dff'),
#                    cat.col = c("#440154ff", '#21908dff'),
#                    cat.fontfamily = "sans",
#                    fontfamily = "sans"
# )
# dev.off()
# draw.pairwise.venn(area1 = length(fucBindingFeats_neg[[1]]),
#                    area2 = length(fucBindingFeats_neg[[2]]),
#                    cross.area = length(intersect(fucBindingFeats_neg[[1]], fucBindingFeats_neg[[2]])),
#                    euler.d = T, scaled = T, ext.text = F, cex = 2,
#                    category = c('Fucose (monosacc.)', 'Terminal fucose'), cat.cex = 2,
#                    cat.pos = c(340,30),
#                    fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
#                    col = c("#440154ff", '#21908dff'),
#                    cat.col = c("#440154ff", '#21908dff'),
#                    cat.fontfamily = "sans",
#                    fontfamily = "sans"
# )
# 
# # Mannose binding
# manBindingFeats_pos = list(row.names(stats)[stats$Man_effectSize > 0 & sigFeats$Man & topImpfeats$Man],
#                            row.names(stats)[stats$High_Mannose_effectSize > 0 & sigFeats$High_Mannose & topImpfeats$High_Mannose],
#                            row.names(stats)[stats$Man.a1.2.Man_effectSize > 0 & sigFeats$Man.a1.2.Man & topImpfeats$Man.a1.2.Man])
# 
# manBindingFeats_neg = list(row.names(stats)[stats$Man_effectSize < 0 & sigFeats$Man & topImpfeats$Man],
#                            row.names(stats)[stats$High_Mannose_effectSize < 0 & sigFeats$High_Mannose & topImpfeats$High_Mannose],
#                            row.names(stats)[stats$Man.a1.2.Man_effectSize < 0 & sigFeats$Man.a1.2.Man & topImpfeats$Man.a1.2.Man])
# 
# intersect(intersect(manBindingFeats_pos[[1]], manBindingFeats_pos[[2]]), manBindingFeats_pos[[3]])
# intersect(intersect(manBindingFeats_neg[[1]], manBindingFeats_neg[[2]]), manBindingFeats_neg[[3]])
# 
# dev.off()
# venn.diagram(
#   x = manBindingFeats_pos,
#   category.names = c("Man" , "High_mannose" , "Man(a1-2)Man"),
#   filename = './manuscript/figures/subplots/man_feats_up.png',
#   output=F,
#   imagetype="png" ,
#   height = 480 , 
#   width = 480 , 
#   resolution = 300,
#   compression = "lzw",
#   lwd = 1,
#   col=c("#440154ff", '#21908dff', '#DAA520FF'),
#   fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#DAA520FF',0.3)),
#   cex = 0.5,
#   fontfamily = "sans",
#   cat.cex = 0.3,
#   cat.default.pos = "outer",
#   cat.pos = c(-27, 27, 135),
#   cat.dist = c(0.055, 0.055, 0.085),
#   cat.fontfamily = "sans",
#   cat.col = c("#440154ff", '#21908dff', '#DAA520FF'),
#   rotation = 1
# )
# 
# 
# venn.diagram(
#   x = manBindingFeats_neg,
#   category.names = c("Man" , "High_mannose" , "Man(a1-2)Man"),
#   filename = './manuscript/figures/subplots/man_feats_down.png',
#   output=F,
#   imagetype="png" ,
#   height = 480 , 
#   width = 480 , 
#   resolution = 300,
#   compression = "lzw",
#   lwd = 1,
#   col=c("#440154ff", '#21908dff', '#DAA520FF'),
#   fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#DAA520FF',0.3)),
#   cex = 0.5,
#   fontfamily = "sans",
#   cat.cex = 0.3,
#   cat.default.pos = "outer",
#   cat.pos = c(-27, 27, 135),
#   cat.dist = c(0.055, 0.055, 0.085),
#   cat.fontfamily = "sans",
#   cat.col = c("#440154ff", '#21908dff', '#DAA520FF'),
#   rotation = 1
# )
# 
# # Galactose binding
# galBindingFeats_pos = list(row.names(stats)[stats$Gal.b1.4.Glc_effectSize > 0 & sigFeats$Gal.b1.4.Glc & topImpfeats$Gal.b1.4.Glc],
#                            row.names(stats)[stats$Gal_effectSize > 0 & sigFeats$Gal & topImpfeats$Gal],
#                            row.names(stats)[stats$GalNAc_effectSize > 0 & sigFeats$GalNAc & topImpfeats$GalNAc],
#                            row.names(stats)[stats$Gal.b1.4.GlcNAc_effectSize > 0 & sigFeats$Gal.b1.4.GlcNAc & topImpfeats$Gal.b1.4.GlcNAc])
# 
# galBindingFeats_neg = list(row.names(stats)[stats$Gal.b1.4.Glc_effectSize < 0 & sigFeats$Gal.b1.4.Glc & topImpfeats$Gal.b1.4.Glc],
#                            row.names(stats)[stats$Gal_effectSize < 0 & sigFeats$Gal & topImpfeats$Gal],
#                            row.names(stats)[stats$GalNAc_effectSize < 0 & sigFeats$Gal.b1.4.GlcNAc & topImpfeats$Gal.b1.4.GlcNAc],
#                            row.names(stats)[stats$Gal.b1.4.GlcNAc_effectSize < 0 & sigFeats$Gal.b1.4.GlcNAc & topImpfeats$Gal.b1.4.GlcNAc])
# 
# 
# intersect(intersect(intersect(galBindingFeats_pos[[1]], galBindingFeats_pos[[2]]), galBindingFeats_pos[[3]]), galBindingFeats_pos[[4]])
# intersect(intersect(intersect(galBindingFeats_neg[[1]], galBindingFeats_neg[[2]]), galBindingFeats_neg[[3]]), galBindingFeats_neg[[4]])
# 
# dev.off()
# draw.quad.venn(area1 = length(galBindingFeats_pos[[1]]),
#                area2 = length(galBindingFeats_pos[[2]]),
#                area3 = length(galBindingFeats_pos[[3]]),
#                area4 = length(galBindingFeats_pos[[4]]),
#                
#                n12 = length(intersect(galBindingFeats_pos[[1]], galBindingFeats_pos[[2]])),
#                n13 = length(intersect(galBindingFeats_pos[[1]], galBindingFeats_pos[[3]])),
#                n14 = length(intersect(galBindingFeats_pos[[1]], galBindingFeats_pos[[4]])),
#                n23 = length(intersect(galBindingFeats_pos[[2]], galBindingFeats_pos[[3]])),
#                n24 = length(intersect(galBindingFeats_pos[[2]], galBindingFeats_pos[[4]])),
#                n34 = length(intersect(galBindingFeats_pos[[3]], galBindingFeats_pos[[4]])),
#                
#                n123 = length(intersect(intersect(galBindingFeats_pos[[1]], galBindingFeats_pos[[2]]), galBindingFeats_pos[[3]])),
#                n124 = length(intersect(intersect(galBindingFeats_pos[[1]], galBindingFeats_pos[[2]]), galBindingFeats_pos[[4]])),
#                n134 = length(intersect(intersect(galBindingFeats_pos[[1]], galBindingFeats_pos[[3]]), galBindingFeats_pos[[4]])),
#                n234 = length(intersect(intersect(galBindingFeats_pos[[2]], galBindingFeats_pos[[3]]), galBindingFeats_pos[[4]])),
#                
#                n1234 = length(intersect(intersect(intersect(galBindingFeats_pos[[1]], galBindingFeats_pos[[2]]), galBindingFeats_pos[[3]]), galBindingFeats_pos[[4]])),
#                
#                category = c('Lac', 'Gal', 'GalNAc', 'LacNAc'), cat.cex = 2,
#                cat.fontfamily = "sans",
#                fontfamily = "sans",
#                euler.d = T, scaled = T, ext.text = F, cex = 2,
#                
#                col=c("#440154ff", '#21908dff', '#DAA520FF', 'red3'),
#                fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#DAA520FF',0.3), alpha('red3',0.3)),
#                cat.col = c("#440154ff", '#21908dff', '#DAA520FF', 'red3')
# )
# 
# 
# dev.off()
# draw.quad.venn(area1 = length(galBindingFeats_neg[[1]]),
#                area2 = length(galBindingFeats_neg[[2]]),
#                area3 = length(galBindingFeats_neg[[3]]),
#                area4 = length(galBindingFeats_neg[[4]]),
#                
#                n12 = length(intersect(galBindingFeats_neg[[1]], galBindingFeats_neg[[2]])),
#                n13 = length(intersect(galBindingFeats_neg[[1]], galBindingFeats_neg[[3]])),
#                n14 = length(intersect(galBindingFeats_neg[[1]], galBindingFeats_neg[[4]])),
#                n23 = length(intersect(galBindingFeats_neg[[2]], galBindingFeats_neg[[3]])),
#                n24 = length(intersect(galBindingFeats_neg[[2]], galBindingFeats_neg[[4]])),
#                n34 = length(intersect(galBindingFeats_neg[[3]], galBindingFeats_neg[[4]])),
#                
#                n123 = length(intersect(intersect(galBindingFeats_neg[[1]], galBindingFeats_neg[[2]]), galBindingFeats_neg[[3]])),
#                n124 = length(intersect(intersect(galBindingFeats_neg[[1]], galBindingFeats_neg[[2]]), galBindingFeats_neg[[4]])),
#                n134 = length(intersect(intersect(galBindingFeats_neg[[1]], galBindingFeats_neg[[3]]), galBindingFeats_neg[[4]])),
#                n234 = length(intersect(intersect(galBindingFeats_neg[[2]], galBindingFeats_neg[[3]]), galBindingFeats_neg[[4]])),
#                
#                n1234 = length(intersect(intersect(intersect(galBindingFeats_neg[[1]], galBindingFeats_neg[[2]]), galBindingFeats_neg[[3]]), galBindingFeats_neg[[4]])),
#                
#                category = c('Lac', 'Gal', 'GalNAc', 'LacNAc'), cat.cex = 2,
#                cat.fontfamily = "sans",
#                fontfamily = "sans",
#                euler.d = T, scaled = T, ext.text = F, cex = 2,
#                
#                col=c("#440154ff", '#21908dff', '#DAA520FF', 'red3'),
#                fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#DAA520FF',0.3), alpha('red3',0.3)),
#                cat.col = c("#440154ff", '#21908dff', '#DAA520FF', 'red3')
# )
# 
# 
# # Glucose binding
# glcBindingFeats_pos = list(row.names(stats)[stats$Glc_effectSize > 0 & sigFeats$Glc & topImpfeats$Glc],
#                            row.names(stats)[stats$GlcNAc_effectSize > 0 & sigFeats$GlcNAc & topImpfeats$GlcNAc])
# 
# glcBindingFeats_neg = list(row.names(stats)[stats$Glc_effectSize < 0 & sigFeats$Glc & topImpfeats$Glc],
#                            row.names(stats)[stats$GlcNAc_effectSize < 0 & sigFeats$GlcNAc & topImpfeats$GlcNAc])
# 
# intersect(glcBindingFeats_pos[[1]], glcBindingFeats_pos[[2]])
# intersect(glcBindingFeats_neg[[1]], glcBindingFeats_neg[[2]])
# 
# dev.off()
# draw.pairwise.venn(area1 = length(glcBindingFeats_pos[[1]]),
#                    area2 = length(glcBindingFeats_pos[[2]]),
#                    cross.area = length(intersect(glcBindingFeats_pos[[1]], glcBindingFeats_pos[[2]])),
#                    euler.d = T, scaled = T, ext.text = F, cex = 2,
#                    category = c('Glucose', 'GlcNAc'), cat.cex = 2,
#                    cat.pos = c(340,20),
#                    cat.dist = c(0.04, 0.05),
#                    fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
#                    col = c("#440154ff", '#21908dff'),
#                    cat.col = c("#440154ff", '#21908dff'),
#                    cat.fontfamily = "sans",
#                    fontfamily = "sans"
# )
# dev.off()
# draw.pairwise.venn(area1 = length(glcBindingFeats_neg[[1]]),
#                    area2 = length(glcBindingFeats_neg[[2]]),
#                    cross.area = length(intersect(glcBindingFeats_neg[[1]], glcBindingFeats_neg[[2]])),
#                    euler.d = T, scaled = T, ext.text = F, cex = 2,
#                    category = c('Glucose', 'GlcNAc'), cat.cex = 2,
#                    cat.pos = c(340,30),
#                    fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
#                    col = c("#440154ff", '#21908dff'),
#                    cat.col = c("#440154ff", '#21908dff'),
#                    cat.fontfamily = "sans",
#                    fontfamily = "sans"
# )

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
# previous FIGURE 4
###
pdf(file = paste('./manuscript/figures/subplots/', 
                 'important_feats',
                 '.pdf', sep = ''),
    width = 18,
    height = 7)
pheatmap(allSigFeats, cluster_rows = F, cluster_cols = F, # Drop empty columns after pruning rows
         color = cols,
         breaks = breakLst,
         border_color = 'white',
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

##################################
## Find prototypical binding sites for each ligand
row.names(stats) = gsub('^H_bin', 'a-helix_bin', row.names(stats))
row.names(stats) = gsub('^B_bin', 'b-bridge_bin', row.names(stats))
row.names(stats) = gsub('^E_bin', 'b-strand_bin', row.names(stats))
row.names(stats) = gsub('^G_bin', '3/10-helix_bin', row.names(stats))
row.names(stats) = gsub('^T_bin', 'turn_bin', row.names(stats))
row.names(stats) = gsub('^S_bin', 'bend_bin', row.names(stats))
row.names(stats) = gsub('^X._bin', 'loop_bin', row.names(stats))

ligSpefic_feat_means = stats[,grepl('_weightedFeatureMean', colnames(stats))] # Weighted means of scaled features
colnames(ligSpefic_feat_means) = gsub('_weightedFeatureMean$', '', colnames(ligSpefic_feat_means)) 

# ligSpefic_feat_means = ligSpefic_feat_means[colnames(allSigFeats),] # Limit to the set of 60 features that are significant and important in at least one of the 10 best performing ligands with RF

top60_scaledFeats = scaledFeats # Use all features
colnames(top60_scaledFeats) = row.names(stats) # carry over reformatted feature names
# top60_scaledFeats = top60_scaledFeats[,colnames(allSigFeats)] # Subset of 60 top features scaled from 0,1

for (i in 1:ncol(ligTags)){
  cat(colnames(ligTags[i]), '\n')
  boundDat = rbind(ligSpefic_feat_means[,i], top60_scaledFeats[ligTags[,i],])
  dists = distance(boundDat, method = 'euclidean', use.row.names = T)
  
  cat(row.names(boundDat)[which.min(dists[1,2:nrow(boundDat)]) +1],'\t')
  cat(min(dists[1,2:nrow(boundDat)]),'\n\n')
}





vols = predFeats[,grepl('^vol_',colnames(predFeats))]

vols$meanChng = (vols[,4]-vols[,3]) + (vols[,3]-vols[,2]) + (vols[,2]-vols[,1])

plot(density(vols$meanChng))


cor.test(vols$meanChng, predFeats$skew_4Ang)
cor.test(vols$meanChng, predFeats$skew_6Ang)
cor.test(vols$meanChng, predFeats$skew_8Ang)
cor.test(vols$meanChng, predFeats$skew_10Ang)


plot(vols$meanChng, predFeats$zern_PC4)
cor.test(vols$meanChng, predFeats$zern_PC4) # Moderate correlation

vols$pcntChng = (((vols[,4]-vols[,3])/vols[,3]) + ((vols[,3]-vols[,2])/vols[,2]) + ((vols[,2]-vols[,1])/vols[,1]))/3
vols$pcntChng[vols$vol_4Ang == 0] = 0 # Set NAs (where volume = 0) to 0
sum(is.na(vols$pcntChng))

plot(density(vols$pcntChng))



plot(vols$pcntChng, predFeats$zern_PC4)
cor.test(vols$pcntChng, predFeats$zern_PC4) # moderate correlation

vols$adjPcnt = vols$pcntChng / apply(vols[,1:4], 1, mean)
vols$adjPcnt[vols$pcntChng == 0] = 0 # Set NAs to 0
sum(is.na(vols$adjPcnt))

plot(density(vols$adjPcnt))


cor.test(vols$adjPcnt, predFeats$zern_PC4) # No correlation
plot(vols$adjPcnt, predFeats$zern_PC4)




# Binned D2 PC2
par(mfrow=c(1,3))
cor.test(vols$meanChng, predFeats$binnedD2_PC2) # No correlation
summary(lm(predFeats$binnedD2_PC2 ~ vols$meanChng))
plot(vols$meanChng, predFeats$binnedD2_PC2, xlab = 'Mean volume change across thresholds', ylab = 'Binned D2 PC2', main = 'Mean vol. change')
abline(a = 0.6696200, b = -0.0015655)
mtext('R = -0.15, p < 0.001', line = -2)

cor.test(vols$pcntChng, predFeats$binnedD2_PC2) # Very weak correlation
summary(lm(predFeats$binnedD2_PC2 ~ vols$pcntChng))
plot(vols$pcntChng, predFeats$binnedD2_PC2, xlab = 'Mean percent change across thresholds', ylab = 'Binned D2 PC2', main = 'Mean percent change')
abline(a = -1.0384, b = 4.3277)
mtext('R = 0.10, p < 0.001', line = -2)

cor.test(vols$adjPcnt, predFeats$binnedD2_PC2) # Moderate correlation
summary(lm(predFeats$binnedD2_PC2 ~ vols$adjPcnt))
plot(vols$adjPcnt, predFeats$binnedD2_PC2, xlab = 'Adjusted percent change across thresholds', ylab = 'Binned D2 PC2', main = 'Volume-adjusted percent change')
abline(a = -1.90697, b = 3927.39694)
mtext('R = 0.39, p < 0.001', line = -2)

# copied @ 900 x 550


# Zern PC4
par(mfrow=c(1,3))
cor.test(vols$meanChng, predFeats$zern_PC4) # Moderate correlation
summary(lm(predFeats$zern_PC4 ~ vols$meanChng))
plot(vols$meanChng, predFeats$zern_PC4, xlab = 'Mean volume change across thresholds', ylab = 'Zern PC4', main = 'Mean vol. change')
abline(a = -0.6966667, b = 0.0016288)
mtext('R = 0.24, p < 0.001', line = -2)

cor.test(vols$pcntChng, predFeats$zern_PC4) # Moderate correlation
summary(lm(predFeats$zern_PC4 ~ vols$pcntChng))
plot(vols$pcntChng, predFeats$zern_PC4, xlab = 'Mean percent change across thresholds', ylab = 'Zern PC4', main = 'Mean percent change')
abline(a = -1.5127, b = 6.3040)
mtext('R = 0.22, p < 0.001', line = -2)

cor.test(vols$adjPcnt, predFeats$zern_PC4) # No correlation
summary(lm(predFeats$zern_PC4 ~ vols$adjPcnt))
plot(vols$adjPcnt, predFeats$zern_PC4, xlab = 'Adjusted percent change across thresholds', ylab = 'Zern PC4', main = 'Volume-adjusted percent change')
abline(a = 0.16525, b = -340.33114)
mtext('R = -0.05, p < 0.01', line = -2)
# 900 x 550


ligSpefic_feat_means['vol_4Ang', 1]
ligSpefic_feat_means['vol_6Ang', 1]
ligSpefic_feat_means['vol_8Ang', 1]
ligSpefic_feat_means['vol_10Ang', 1]

ligSpefic_feat_means['vol_4Ang', 2]
ligSpefic_feat_means['vol_6Ang', 2]
ligSpefic_feat_means['vol_8Ang', 2]
ligSpefic_feat_means['vol_10Ang', 2]

ligSpefic_feat_means['vol_4Ang', 8]
ligSpefic_feat_means['vol_6Ang', 8]
ligSpefic_feat_means['vol_8Ang', 8]
ligSpefic_feat_means['vol_10Ang', 8]


vols$six8 = vols$vol_8Ang - vols$vol_6Ang
vols$six8_pcnt = vols$six8 / vols$vol_6Ang
vols$six8_pcnt[vols$vol_6Ang == 0] = 0
sum(is.na(vols$six8_pcnt))

plot(density(vols$six8_pcnt))



########################
## FIG 4 distributions
########################

d2Dists = read.delim('./analysis/d2_binnedMeasures.csv', header = T, sep =',', stringsAsFactors = F)
colnames(d2Dists) = gsub('^X','',colnames(d2Dists))

d2Dists_fig4 = d2Dists[c('1SID_BGC:I:1',
                         '1HGH_MNA:A:349',
                         '1CVN_MAN:E:1'),]




###################
dev.off()
pdf(file = paste('./manuscript/figures/subplots/', 
                 'sia_neuac_hiMan_d2plots',
                 '.pdf', sep = ''),
    #width = 14.25, height = 12)
    width = 9.5, height = 8)
par(mfrow = c(3,1),
    cex.lab = 1.5,
    cex.axis = 1.8,
    cex.main = 2.5,
    mai = c(0.51, 0.11, 0.11, .51),
    mar = c(3.1, 3.1, 4.1, 1.1))
xmax = 54 # plot values from first 54 bins (xlim of 27)
ymax = 70000
scale = 10000
threshStrings = c('_4Ang$', '_6Ang$', '_8Ang$', '_10Ang$')
mainTitles = c('Terminal NeuAc Group', 'NeuAc', 'High mannose group')
tCols = c('magenta', "firebrick3", 'darkorange', 'gold2')

for (i in 1:nrow(d2Dists_fig4)){
  for (t in 4:1){
    x = colnames(d2Dists_fig4)[grepl(threshStrings[t], colnames(d2Dists_fig4))]
    x = gsub(threshStrings[t], '', x)
    x = as.numeric(x)
    x = c(0, x[1:xmax])
    if (t ==4){
      barplot(height = c(0,as.numeric(d2Dists_fig4[i,grepl(threshStrings[t], colnames(d2Dists_fig4))])[1:xmax])/scale,
              width = 1, space = 0,
              ylab = '', xlab = '', main = mainTitles[i], # 'Distances between pairs of surface points (\uc5);
              ylim = c(0, ymax/scale),
              axes = F,
              names.arg = x,
              col = alpha(tCols[t],0.9))
      axis(side = 1, at = c(0:xmax), labels = FALSE)
      axis(side = 2, at = c(0:7), labels = c(0:7), line = -0.5)
      #mtext('Counts (10k)', side = 2, line = 1.5, cex = 1.5, col = alpha('black',1))  # Add y axis label
    } else{
      barplot(height = c(0,as.numeric(d2Dists_fig4[i,grepl(threshStrings[t], colnames(d2Dists_fig4))])[1:xmax])/scale, 
              width = 1, space = 0,
              ylim = c(0, ymax/scale),
              axes = F, ylab = '', xlab ='', main = '',
              col = alpha(tCols[t],0.9))
    }
    if (i == 2 & t == 4){
      legend(x = 46, y = 7.5,
             legend = c('4 \uc5', '6 \uc5', '8 \uc5', '10 \uc5'),
             pch = 15,
             col = tCols,
             cex = 3, pt.cex = 5,
             bty = 'n')
    }
    if (t != 1){par(new=T)}
  }
  volText = ''
  for (t in 4:1){
    med = predFeats[row.names(d2Dists_fig4)[i], grepl(paste0('^med',threshStrings[t]), colnames(predFeats))]
    abline(v = med*2, col = 'black', lwd = 5)
    abline(v = med*2, col = tCols[t], lwd = 2.5)
    volText = paste(volText, as.character())
  }
}



dev.off()
#####################






#####################

topLigs = row.names(allSigFeats)
topLigTag = apply(ligTags[,topLigs], 1, any)

top60.scaledPC = prcomp(x = top60_scaledFeats, scale. = T)

allLigColors = rep("grey40", nrow(predFeats))
allLigColors[apply(ligTags[,c(1,8,9)], 1, any)] = 'purple2'
allLigColors[apply(ligTags[,c(3,14)], 1, any)] = 'red1'
allLigColors[apply(ligTags[,c(2,6,12)], 1, any)] = 'forestgreen'
allLigColors[apply(ligTags[,c(4,5,7,11,15)], 1, any)] = 'goldenrod2'
allLigColors[apply(ligTags[,c(10,13)], 1, any)] = 'mediumblue'


plot(top60.scaledPC$x, col = allLigColors)

topLigShapes = rep(0, sum(topLigTag))
topLigShapes[apply(ligTags[topLigTag, c(14, 1, 4, 6, 10)], 1 , any)] = 15
topLigShapes[apply(ligTags[topLigTag, c(9, 5, 2, 13)], 1 , any)] = 16
topLigShapes[ligTags[topLigTag, 12]] = 17

plot(top60.scaledPC$x[topLigTag,], col = alpha(allLigColors[topLigTag], 0.6), pch = topLigShapes)
legend(x = 'center',
       legend = topLigs,
       col = c('red1', 
         'purple2', 'purple2', 
         'goldenrod2', 'goldenrod2', 
         'forestgreen', 'forestgreen', 'forestgreen', 
         'mediumblue', 'mediumblue'),
       pch = c(15,
               15, 16,
               15, 16,
               15, 16, 17,
               15, 16)
       )

topLig_top60.scaledPC = prcomp(x = top60_scaledFeats[topLigTag,], scale. = T)

plot(topLig_top60.scaledPC$x, col = alpha(allLigColors[topLigTag], 0.6), pch = topLigShapes)



all.umap = umap(top60_scaledFeats)
plot(all.umap$layout, #xlim = c(-5.5,7), ylim = c(-6.5,6.5),
     col = allLigColors,
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = 'UMAP vis for 60 top predictive features - Scaled')


plot(all.umap$layout[topLigTag,], #xlim = c(-5.5,7), ylim = c(-6.5,6.5),
     col = allLigColors[topLigTag], pch = topLigShapes,
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = 'UMAP vis for 60 top predictive features - Scaled - Top 10 gly subset')

##
top.umap = umap(top60_scaledFeats[topLigTag,])

plot(top.umap$layout,
     xlim = c(-11,11), ylim = c(-8,13),
     col = alpha(allLigColors[topLigTag], 0.6), pch = topLigShapes,
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = 'UMAP vis for 60 top predictive features - Scaled top 10 gly')

featMeans = predict(top.umap, t(ligSpefic_feat_means[,topLigs]))
par(new = T)
plot(featMeans, cex = 3.15,
     xlab = '', ylab = '', axes = F,
     xlim = c(-11,11), ylim = c(-8,13),
     col = 'black',
     pch = c(15,
             15, 16,
             15, 16,
             15, 16, 17,
             15, 16))
par(new = T)
plot(featMeans, cex = 3,
     xlab = '', ylab = '', axes = F,
     xlim = c(-11,11), ylim = c(-8,13),
     col = alpha(c('red1', 
                   'purple2', 'purple2', 
                   'goldenrod2', 'goldenrod2', 
                   'forestgreen', 'forestgreen', 'forestgreen', 
                   'mediumblue', 'mediumblue'), 1),
     pch = c(15,
             15, 16,
             15, 16,
             15, 16, 17,
             15, 16))


# Plot weighted distributions of all 221 features for all 15 ligands


# allTopFeatCols = rep('', 60)
# 
# for (i in 1:length(names(annot_cols$Feature_Type))){
#   feat = names(annot_cols$Feature_Type)[i]
#   allTopFeatCols[annot$Feature_Type == feat] = annot_cols$Feature_Type[i]
# }

allTopFeatCols = featColors

for (k in 1:nrow(stats)){
  featName = row.names(stats)[k]
  
  feat = featName
  feat = gsub('^b-strand_', 'E_', feat)
  feat = gsub('^loop_', 'X._', feat)
  feat = gsub('^3/10-helix_', 'G_', feat)
  feat = gsub('^bend_', 'S_', feat)
  feat = gsub('^turn_', 'T_', feat)
  feat = gsub('^a-helix_', 'H_', feat)
  feat = gsub('^b-bridge_', 'B_', feat)
  
  
  pdf(file = paste('./manuscript/figures/subplots/allSigFeatDists/', 
                   gsub('/','-',featName),
                   '.pdf', sep = ''),
      width = 11,
      height = 8.5)
  par(mfrow=c(3,5))
  for (i in 1:ncol(ligTags)){
    weightFlag = F
    lig = colnames(ligTags)[i]
    ligPattern = paste0('^', lig)
    
    des_w = subset(des, subset = ligTags[,lig])
    
    if (IQR(predFeats[ligTags[,lig], feat]) == 0){
      x = density(predFeats[ligTags[,lig], feat])
      yLim = c(0, max(max(density(predFeats[,feat])$y), max(x$y)))
    } else{
      x = svysmooth(formula = asOneSidedFormula(feat), design = des_w)
      yLim = c(0, max(max(density(predFeats[,feat])$y),max(x[[1]][[2]])))
      weightFlag = T
    }
    
    val = stats[featName,grepl(paste0(ligPattern, '_effectSize$'), colnames(stats))]
    
    if (val > 0){
      dir = 'enriched'
    } else{
      dir = 'depleted'
    }
    
    if (sigFeats[feat,i] & topImpfeats[feat,i]){ #ith lig
      meth = ', RF & WMW'
    } else if (sigFeats[feat,i] | topImpfeats[feat,i]){
      if (sigFeats[feat,i]){
        meth = ', WMW only'
      } else{
        meth = ', RF only' 
      }
    } else {
      meth = ', non-significant'
    }
    
    svyhist(formula = asOneSidedFormula(feat), design = des,
            xlim = range(predFeats[,feat]), ylim  = yLim,
            xlab = featName,
            main = paste(lig, '\n', dir, meth, sep = ''),
            col = allTopFeatCols[k])
    lines(x,
          xlim = range(predFeats[,feat]),
          ylim  = yLim,
          lwd = 3)
    if (! weightFlag){
      note = 'Unweighted\n(IQR == 0)'
      text(x = mean(range(predFeats[,feat])),
           y = yLim[2]*.9,
           note)
    }
  }
  dev.off()
}


###############################
## Fig 0 barplots
###############################

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


