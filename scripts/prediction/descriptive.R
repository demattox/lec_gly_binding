
library(pheatmap)
library(VennDiagram)
library(ggplot2)
library(corrplot)
library(reshape2)
library(RColorBrewer)
library(philentropy)
library(seqinr)
library(stringr)

library(survey)

colfunc = colorRampPalette(c("red","goldenrod","forestgreen","royalblue","darkviolet"))
threshColors = c('firebrick3', 'darkorange2', 'darkgoldenrod2', 'gold2')

###########################
# Read in data
###########################
homeDir = '/Users/dmattox/cbk/lec_gly_binding/'
setwd(homeDir)

# Load McCaldon aa frequencies
mccaldon = read.delim('~/episweep/data/mccaldon.csv', header = F, sep = ",", stringsAsFactors = F)
colnames(mccaldon) = c("aa", 'mcFreq')

ligTags <- read.delim(file = './analysis/training/data_in/ligTags.tsv', sep = '\t', stringsAsFactors = F)
predFeats <- read.delim(file = './analysis/training/data_in/predFeats.csv', sep = ',', stringsAsFactors = F)
bsResiDat <- read.delim(file = './analysis/training/data_in/bsResiDat.tsv', sep = '\t', stringsAsFactors = F)

###########################
# McCaldon analysis
###########################
mccaldon$aa = lapply(mccaldon$aa, FUN = aaa) # convert from 1-letter aminoacid code to 3, capitalized
mccaldon$aa = lapply(mccaldon$aa, FUN = str_to_upper)

mccaldon$freq_bin1 = 0
mccaldon$freq_bin2 = 0
mccaldon$freq_bin3 = 0
mccaldon$freq_bin4 = 0

clusLst50 = unique(bsResiDat$seqClust50)

for (i in 1:nrow(mccaldon)){
  b1 = b2 = b3 = b4 = 0
  for (j in 1:length(clusLst50)){
    tag = bsResiDat$seqClust50 == clusLst50[j]
    b1 = b1 + mean( bsResiDat[,grepl( paste('^',mccaldon$aa[i][[1]],'_bin1$', sep = ''), colnames(bsResiDat))] ) # Running sum of the average frequency for each cluster
    b2 = b2 + mean( bsResiDat[,grepl( paste('^',mccaldon$aa[i][[1]],'_bin2$', sep = ''), colnames(bsResiDat))] )
    b3 = b3 + mean( bsResiDat[,grepl( paste('^',mccaldon$aa[i][[1]],'_bin3$', sep = ''), colnames(bsResiDat))] )
    b4 = b4 + mean( bsResiDat[,grepl( paste('^',mccaldon$aa[i][[1]],'_bin4$', sep = ''), colnames(bsResiDat))] )
  }
  mccaldon$freq_bin1[i] = b1/length(clusLst50) # Take average across the clusters s.t. each cluster gets equal weight
  mccaldon$freq_bin2[i] = b2/length(clusLst50)
  mccaldon$freq_bin3[i] = b3/length(clusLst50)
  mccaldon$freq_bin4[i] = b4/length(clusLst50)
}

mcFreqs = mccaldon[,grep('[Ff]req', colnames(mccaldon))]
colnames(mcFreqs) = c('McCaldon', 'bin1', 'bin2', 'bin3', 'bin4')
mcFreqs$aa = unlist(mccaldon$aa)
mcFreqs = mcFreqs[,c(6,2:5,1)]

meltMc = melt(data = mcFreqs, id.vars = "aa", measure.vars = c("bin1", "bin2", "bin3", "bin4"))
colnames(meltMc) = c("Amino_acid", "Bin_Number", "Frequency")

# pdf(file = './analysis/prelimPlots/McCaldon_Frequencies.pdf', width = 16, height = 7)
# ggplot(meltMc, aes(fill = Bin_Number, x = Amino_acid, y = Frequency)) + geom_bar(stat="identity", color="black", position=position_dodge())+
#   scale_fill_manual(values=c(threshColors)) + 
#   theme_grey(base_size = 22) + 
#   geom_segment(aes(x = 0.6, y = mccaldon$mcFreq[1], xend = 1.4, yend = mccaldon$mcFreq[1]), cex=2) + 
#   geom_segment(aes(x = 1.6, y = mccaldon$mcFreq[2], xend = 2.4, yend = mccaldon$mcFreq[2]), cex=2) + 
#   geom_segment(aes(x = 2.6, y = mccaldon$mcFreq[3], xend = 3.4, yend = mccaldon$mcFreq[3]), cex=2) + 
#   geom_segment(aes(x = 3.6, y = mccaldon$mcFreq[4], xend = 4.4, yend = mccaldon$mcFreq[4]), cex=2) + 
#   geom_segment(aes(x = 4.6, y = mccaldon$mcFreq[5], xend = 5.4, yend = mccaldon$mcFreq[5]), cex=2) + 
#   geom_segment(aes(x = 5.6, y = mccaldon$mcFreq[6], xend = 6.4, yend = mccaldon$mcFreq[6]), cex=2) + 
#   geom_segment(aes(x = 6.6, y = mccaldon$mcFreq[7], xend = 7.4, yend = mccaldon$mcFreq[7]), cex=2) + 
#   geom_segment(aes(x = 7.6, y = mccaldon$mcFreq[8], xend = 8.4, yend = mccaldon$mcFreq[8]), cex=2) + 
#   geom_segment(aes(x = 8.6, y = mccaldon$mcFreq[9], xend = 9.4, yend = mccaldon$mcFreq[9]), cex=2) + 
#   geom_segment(aes(x = 9.6, y = mccaldon$mcFreq[10], xend = 10.4, yend = mccaldon$mcFreq[10]), cex=2) + 
#   geom_segment(aes(x = 10.6, y = mccaldon$mcFreq[11], xend = 11.4, yend = mccaldon$mcFreq[11]), cex=2) + 
#   geom_segment(aes(x = 11.6, y = mccaldon$mcFreq[12], xend = 12.4, yend = mccaldon$mcFreq[12]), cex=2) + 
#   geom_segment(aes(x = 12.6, y = mccaldon$mcFreq[13], xend = 13.4, yend = mccaldon$mcFreq[13]), cex=2) + 
#   geom_segment(aes(x = 13.6, y = mccaldon$mcFreq[14], xend = 14.4, yend = mccaldon$mcFreq[14]), cex=2) + 
#   geom_segment(aes(x = 14.6, y = mccaldon$mcFreq[15], xend = 15.4, yend = mccaldon$mcFreq[15]), cex=2) + 
#   geom_segment(aes(x = 15.6, y = mccaldon$mcFreq[16], xend = 16.4, yend = mccaldon$mcFreq[16]), cex=2) + 
#   geom_segment(aes(x = 16.6, y = mccaldon$mcFreq[17], xend = 17.4, yend = mccaldon$mcFreq[17]), cex=2) + 
#   geom_segment(aes(x = 17.6, y = mccaldon$mcFreq[18], xend = 18.4, yend = mccaldon$mcFreq[18]), cex=2) + 
#   geom_segment(aes(x = 18.6, y = mccaldon$mcFreq[19], xend = 19.4, yend = mccaldon$mcFreq[19]), cex=2) + 
#   geom_segment(aes(x = 19.6, y = mccaldon$mcFreq[20], xend = 20.4, yend = mccaldon$mcFreq[20]), cex=2)
# dev.off()

mcFreqs = mccaldon[,grep('[Ff]req', colnames(mccaldon))]
colnames(mcFreqs) = c('McCaldon', 'bin1', 'bin2', 'bin3', 'bin4')
mcFreqs$aa = unlist(mccaldon$aa)
mcFreqs = mcFreqs[,c(6,2:5,1)]

tag = (mcFreqs$bin1 > mcFreqs$McCaldon)*100
tag[!tag] = -100
mcFreqs$bin1 = (mcFreqs$bin1 / mcFreqs$McCaldon) * tag

tag = (mcFreqs$bin2 > mcFreqs$McCaldon)*100
tag[!tag] = -100
mcFreqs$bin2 = (mcFreqs$bin2 / mcFreqs$McCaldon) * tag

tag = (mcFreqs$bin3 > mcFreqs$McCaldon)*100
tag[!tag] = -100
mcFreqs$bin3 = (mcFreqs$bin3 / mcFreqs$McCaldon) * tag

tag = (mcFreqs$bin4 > mcFreqs$McCaldon)*100
tag[!tag] = -100
mcFreqs$bin4 = (mcFreqs$bin4 / mcFreqs$McCaldon) * tag

meltMc = melt(data = mcFreqs, id.vars = "aa", measure.vars = c("bin1", "bin2", "bin3", "bin4"))
colnames(meltMc) = c("Amino_acid", "Bin_Number", "PercentDiff_McCaldon")

ggplot(meltMc, aes(fill = Bin_Number, x = Amino_acid, y = PercentDiff_McCaldon)) + geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_hline(yintercept = c(0)) +
  scale_fill_manual(values=c(threshColors)) +
  theme_dark(base_size = 22) + 
  labs(title = 'Lectin binding site deviations from McCaldon frequency', x = "Amino Acids", y = "Percent change from McCaldon") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))

###########################
# Descriptive analysis
###########################
uniLigs = unique(bsResiDat$iupac)

scaledFeats = predFeats
for(i in 1:ncol(scaledFeats)){
  scaledFeats[,i] = (scaledFeats[,i] - min(scaledFeats[,i])) / (max(scaledFeats[,i]) - min(scaledFeats[,i]))
}

# # Ligand tags
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
# neuTag = neuCnt >= 1 # Has sialic acid
# fucTag = grepl('^Fuc',uniLigs) # Has a terminal fucose
# 
# 
# # get tags to indicate binding sites containing one of the 15 ligands of interest
# ligTags = as.data.frame(matrix(F, nrow = nrow(bsResiDat), ncol = 3+length(top50)))
# row.names(ligTags) =  row.names(bsResiDat)
# colnames(ligTags) = colnames(topLigOccurences)[2:ncol(topLigOccurences)] 
# 
# ligTags$High_Mannose = bsResiDat$iupac %in% uniLigs[manTag]
# ligTags$Sialic_Acid = bsResiDat$iupac %in% uniLigs[neuTag]
# ligTags$Terminal_Fucose = bsResiDat$iupac %in% uniLigs[fucTag]
# 
# 
# for (i in 1:length(top50)){
#   lig = top50[i]
#   lig = gsub('\\(', '\\\\(', lig)
#   lig = gsub('\\)', '\\\\)', lig)
#   lig = gsub('\\[', '\\\\]', lig)
#   lig = gsub('\\[', '\\\\]', lig)
#   lig = gsub('^', '\\^', lig)
#   lig = gsub('$', '\\$', lig)
#   
#   colInd = grep(lig, colnames(ligTags))
#   ligTags[,colInd] = grepl(lig, bsResiDat$iupac)
# }

# Add tags to bsResiDat
# bsResiDat = cbind(bsResiDat, ligTags)



# Statistical test for cluster independent feature enrichment by ligand
stats = as.data.frame(matrix(0, nrow = ncol(predFeats), ncol = (ncol(ligTags)*4)))
row.names(stats) = colnames(predFeats)
colnames(stats) = c(paste(colnames(ligTags), 'p', sep = '_'), paste(colnames(ligTags), '_FC', sep = '_'), paste(colnames(ligTags), 'effectSize', sep = '_'), paste(colnames(ligTags), 'adj', sep = '_'))

scaled_stats = stats # In parallel, try this with scaled features for more robust change metric, shouldn't impact WMW test

for (i in 1:ncol(ligTags)){ # for ligand i
  # meanClust50_data = as.data.frame(matrix(NA, nrow = length(clusLst50), ncol = (ncol(predFeats)*2))) # Rows for unique clusters, columns for predictive features from bsites WITH or WITHOUT current ligand
  # colnames(meanClust50_data) = c(paste(colnames(predFeats), 'with', sep = '_'),  paste(colnames(predFeats), 'without', sep = '_'))
  for(k in 1: ncol(predFeats)){ # for each feature k
    meanClust50_data = as.data.frame(matrix(NA, nrow = length(clusLst50), ncol = 2)) # Rows for unique clusters, columns for predictive feature from bsites WITH or WITHOUT current ligand
    colnames(meanClust50_data) = c('with', 'without')
    scaledClus50_data = meanClust50_data
    for (j in 1:length(clusLst50)){ # for each cluster j
      clustTag = bsResiDat$seqClust50 == clusLst50[j]
      withTag = clustTag & ligTags[,i]
      withoutTag = clustTag & !ligTags[,i]
      if ( sum(withTag) >= 1){
        meanClust50_data[j, 1] = mean( predFeats[withTag,k] )
        scaledClus50_data[j, 1] = mean( scaledFeats[withTag,k] )
      }
      if ( sum(withoutTag) >= 1){
        meanClust50_data[j, 2] = mean( predFeats[withoutTag,k] )
        scaledClus50_data[j, 2] = mean( scaledFeats[withoutTag,k] )
      }
    }
    # Raw Feats
    wDat = meanClust50_data[,1]
    wDat = wDat[!is.na(wDat)]
    woDat = meanClust50_data[,2]
    woDat = woDat[!is.na(woDat)]
    # Scaled feats
    sW = scaledClus50_data[!is.na(scaledClus50_data[,1]),1] # Non-NA feats from bSites w/ ligand
    sWO = scaledClus50_data[!is.na(scaledClus50_data[,2]),2] # Non-NA feats from bSites w/ ligand
    
    # ligTest = t.test(wDat[!is.na(wDat)], woDat[!is.na(woDat)])
    ligTest = wilcox.test(x= wDat, y = woDat) # Wilcoxon–Mann–Whitney test, sample sizes can be small (~5% of 230 clusters ~= 10), no reason to assume distribution is normal as it likely isn't, as well as sample sizes often being too small to test for normality
    sTest = wilcox.test(x = sW, y = sWO)
    
    stats[k, grepl('_p$', colnames(stats))][i] = ligTest$p.value # Raw p-value
    stats[k, grepl('_effectSize$', colnames(stats))][i] = ligTest$statistic / (length(wDat) * length(woDat)) # common language effect size
    scaled_stats[k, grepl('_p$', colnames(scaled_stats))][i] = sTest$p.value # Raw p-value
    scaled_stats[k, grepl('_effectSize$', colnames(scaled_stats))][i] = sTest$statistic / (length(sW) * length(sWO)) # common language effect size
    
    if (median(woDat) != 0){
      stats[k, grepl('_FC$', colnames(stats))][i] = median(wDat) / median(woDat) # Fold change in mean values  
    } else{
      stats[k, grepl('_FC$', colnames(stats))][i] = median(wDat)
    }
    if (median(sWO) != 0){
      scaled_stats[k, grepl('_FC$', colnames(scaled_stats))][i] = median(sW) / median(sWO) # Fold change in mean values  
    } else{
      scaled_stats[k, grepl('_FC$', colnames(scaled_stats))][i] = median(sW)
    }
  }
  stats[,grepl('_adj$', colnames(stats))][,i] = p.adjust(stats[grepl('_p$', colnames(stats))][,i], method = "BH") # Benjamini-Hochberg MHT correction (FDR)
  scaled_stats[,grepl('_adj$', colnames(scaled_stats))][,i] = p.adjust(scaled_stats[grepl('_p$', colnames(scaled_stats))][,i], method = "BH") # Benjamini-Hochberg MHT correction (FDR)
}


for (i in 1:ncol(ligTags)){
  lig = colnames(ligTags)[i]
  cat(lig)
  cat('\n----\n')
  cat('DOWN w/ ligand\n')
  cat(row.names(stats)[stats[,grepl('_adj$', colnames(stats))][,i] <= 0.1  &  stats[,grepl('_effectSize$', colnames(stats))][,i] < 0.5])
  cat('\nUP w/ ligand\n')
  cat(row.names(stats)[stats[,grepl('_adj$', colnames(stats))][,i] <= 0.1  &  stats[,grepl('_effectSize$', colnames(stats))][,i] > 0.5])
  cat('\n----------------\n')
}


# Colors to features for plotting
featColors = rep('', nrow(stats))
resiFeats = colorRampPalette(c("plum1","tomato", "firebrick4"))(4)
pocketFeats = colorRampPalette(c('turquoise', 'dodgerblue1', 'blue2'))(3)

featColors[1:11] = 'forestgreen'

featColors[grep('^vol_4Ang$', row.names(stats)) : grep('^leftskew_10Ang$', row.names(stats))] = pocketFeats[1] # features within the d2Feats range
featColors[grepl('^binnedD2', row.names(stats))] = pocketFeats[2] # PCs from the binned D2 measures
featColors[grepl('^zern', row.names(stats))] = pocketFeats[3] # PCs from the 3DZDs

featColors[grepl('^numBSresis', row.names(stats))] = resiFeats[1] # number of residues in binding site features
featColors[gsub('_bin\\d{1}', '', row.names(stats)) %in% c('H', 'B', 'E', 'G', 'T', 'S', 'X.')] = resiFeats[2] # secondary structure features
featColors[gsub('_bin\\d{1}', '', row.names(stats)) %in% c('nonpolar', 'polar', 'posCharge', 'negCharge', 'aromatic')] = resiFeats[3] # amino acid properties
featColors[grepl('^[[:upper:]]{3}_', row.names(stats)) | grepl('^CA$', row.names(stats))] = resiFeats[4] # amino acid identities

resiFeatTag = featColors %in% resiFeats
pocketFeatTag = featColors %in% pocketFeats

# Plot volcano plots from raw features

dev.off()
pdf(file = paste('./analysis/sfgPlots/', 
                 'commonEffect_volcanoes',
                 '.pdf', sep = ''),
    width = 21,
    height = 13)
par(mfrow=c(3,5))
xLim = c(0,1)
# yLim = c(0,10)
for(i in 1:ncol(ligTags)){
  
  yLim = c(0,max(-log10(0.1), -log10(min(stats[,grepl('_adj$', colnames(stats))][,i]))) + 1)
  
  tag = stats[,grepl('_adj$', colnames(stats))][,i] < 0.1
  
  # dev.off()
  plot(0,0,axes = F, main = '', xlab = '', ylab = '', pch = NA)
  bg = "seashell2"
  fg = "ivory"
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = bg)
  # abline(v = c(-1,-.5,0,.5,1), lwd = 6, col = fg)
  # abline(v = c(-1.25,-.75,-.25,.25,.75,1.25), lwd = 3, col = fg)
  # abline(h = c(0,1,2,3,4,5,6), lwd = 6, col = fg)
  # abline(h = c(0.5,1.5,2.5,3.5,4.5,5.5,6.5), lwd = 3, col = fg)
  par(new=T)
  
  plot(stats[,grepl('_effectSize$', colnames(stats))][,i], -log10(stats[,grepl('_adj$', colnames(stats))][,i]), # Plot all points w/ color @ alpha 0.5
       xlab = "Effect size", ylab = "-log10(FDR)", main = colnames(ligTags)[i],
       pch = 19, cex = 2, col = alpha(featColors, 0.5),
       cex.axis = 1.5, cex.main = 2, cex.lab = 1.5,
       xlim = xLim, ylim = yLim)
  par(new=T)
  plot(stats[,grepl('_effectSize$', colnames(stats))][tag,i], -log10(stats[,grepl('_adj$', colnames(stats))][tag,i]), # Plot stat sig points again with alpha 1
       pch = 19, col = featColors[tag], cex = 2,
       axes = F, xlab = "", ylab = "", main = "",
       xlim = xLim, ylim = yLim)
  par(new=T)
  plot(stats[,grepl('_effectSize$', colnames(stats))][tag,i], -log10(stats[,grepl('_adj$', colnames(stats))][tag,i]), # Outline stat sig points in black to highlight them
       col = 'black', cex = 2.05,
       xlab = "", ylab = "", axes = F, main = "",
       xlim = xLim, ylim = yLim)
  abline(h= -log10(0.1), lwd = 2)
  abline(v = 0.5, lty=2, lwd = 2)
}

dev.off()
plot(0,0,axes = F, main = '', xlab = '', ylab = '', pch = NA)
legend(x = 'center',
       col = c('forestgreen', resiFeats, pocketFeats),
       legend = c('PLIP interactions',
                  'bs resi cnt', 'sec struct', 'aa type', 'aa ident',
                  'D2 feats', 'D2 PCs', 'Zern PCs'),
       pch = 19)

dev.off()
par(mfrow=c(3,5))
xLim = c(0,1)
yLim = c(0,1)
# Compare raw feats to scaled feats
for (i in 1:ncol(ligTags)){
  plot(stats[,grepl('_adj$', colnames(stats))][,i], scaled_stats[,grepl('_adj$', colnames(scaled_stats))][,i],
       xlab = 'Raw feat adjusted p-values', ylab = 'Scaled feat adjusted p-values',
       main = colnames(ligTags)[i],
       xlim = xLim, ylim = yLim)
  abline(a = 0, b = 1)
  text(x=0.6, y=0.05, labels = paste( 'Pear. corr = ', cor(stats[,grepl('_adj$', colnames(stats))][,i], scaled_stats[,grepl('_adj$', colnames(scaled_stats))][,i]), sep = ''), cex = 2)
}

################################
# Calculate weighted MWM statistics instead of using means from each cluster
################################
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
des_alt <- svydesign(ids = ~1, data = wmwFeats, weights = 1/cWeights)
des_unweighted <- svydesign(ids = ~1, data = wmwFeats, weights = NULL)


# tst = svyranktest(formula = negCharge_bin1 ~ as.factor(ligTags$Sialic_Acid), design = des, test = 'wilcoxon')
# tst_alt = svyranktest(formula = negCharge_bin1 ~ as.factor(ligTags$Sialic_Acid), design = des_alt, test = 'wilcoxon')
# tst_unweighted = svyranktest(formula = negCharge_bin1 ~ ligTags$Sialic_Acid, design = des_unweighted, test = 'wilcoxon')

# wilcox.test(predFeats$negCharge_bin1[ligTags$Sialic_Acid], predFeats$negCharge_bin1[!ligTags$Sialic_Acid])
# wilcox.test(predFeats$negCharge_bin1[ligTags$Sialic_Acid], predFeats$negCharge_bin1[!ligTags$Sialic_Acid])$statistic / (sum(ligTags$Sialic_Acid) * sum(!ligTags$Sialic_Acid))

### Check formula construction for iterative tests across all features
# svyranktest(formula = negCharge_bin1 ~ Sialic_Acid, design = des_unweighted, test = 'wilcoxon')
# svyranktest(formula = as.formula(paste(colnames(predFeats)[26], ' ~ ', colnames(ligTags)[2], sep = '')),
#             design = des_unweighted, test = 'wilcoxon')
# 
# wilcox.test(predFeats$negCharge_bin1[ligTags$Sialic_Acid], predFeats$negCharge_bin1[!ligTags$Sialic_Acid])
# wilcox.test(predFeats$negCharge_bin1[ligTags$Sialic_Acid], predFeats$negCharge_bin1[!ligTags$Sialic_Acid])$statistic / (sum(ligTags$Sialic_Acid) * sum(!ligTags$Sialic_Acid))


# tmp = svyranktest(formula = as.formula(paste(colnames(predFeats)[1], ' ~ ', colnames(ligTags)[1], sep = '')),
#             design = des_unweighted, test = 'wilcoxon')
# 
# tmp2 = wilcox.test(predFeats$hydrophobics[ligTags$High_Mannose], predFeats$hydrophobics[!ligTags$High_Mannose])

### Check survey means without weights to validate syntax for weighted means
# i=2
# k=26
# des_w = subset(des_unweighted, subset = ligTags[,i]) # temporary design object holding all interactions with the ligand of interest
# des_wo = subset(des_unweighted, subset = !ligTags[,i]) # same as above but WITHOUT the ligand of interest
# 
# svymean(predFeats[ligTags[,i],k], des_w)
# mean(predFeats[ligTags[,i],k])
# 
# svymean(predFeats[!ligTags[,i],k], des_wo)
# mean(predFeats[!ligTags[,i],k])
# 
# des_w = subset(des, subset = ligTags[,i])
# des_wo = subset(des, subset = !ligTags[,i])
# svymean(predFeats[ligTags[,i],k], des_w)
# svymean(predFeats[!ligTags[,i],k], des_wo)




stats_weighted = as.data.frame(matrix(0, nrow = ncol(predFeats), ncol = (ncol(ligTags)*4)))
row.names(stats_weighted) = colnames(predFeats)
colnames(stats_weighted) = c(paste(colnames(ligTags), 'p', sep = '_'), paste(colnames(ligTags), '_FC', sep = '_'), paste(colnames(ligTags), 'effectSize', sep = '_'), paste(colnames(ligTags), 'adj', sep = '_'))

for (i in 1:ncol(ligTags)){ # for ligand i
  des_w = subset(des, subset = ligTags[,i]) # temporary design object holding all interactions with the ligand of interest
  des_wo = subset(des, subset = !ligTags[,i]) # same as above but WITHOUT the ligand of interest
  for(k in 1:ncol(predFeats)){ # for each feature k
    ligTest = svyranktest(formula = as.formula(paste(colnames(predFeats)[k], ' ~ ', colnames(ligTags)[i], sep = '')), # Wilcoxon–Mann–Whitney test, sample sizes can be small (~5% of 230 clusters ~= 10), no reason to assume distribution is normal as it likely isn't
                           design = des_unweighted, 
                           test = 'wilcoxon') 

    if (ligTest$p.value != 0){
      stats_weighted[k, grepl('_p$', colnames(stats_weighted))][i] = ligTest$p.value # Raw p-value
    } else{
      stats_weighted[k, grepl('_p$', colnames(stats_weighted))][i] = 5e-324 # If p-value is too small and R rounds to 0, reset as the smallest positive double as referenced by "?.Machine"
    }

    stats_weighted[k, grepl('_effectSize$', colnames(stats_weighted))][i] = ligTest$estimate # common language effect size (0.5 transformed to 0)
    
    stats_weighted[k, grepl('_FC$', colnames(stats_weighted))][i] = svymean(predFeats[ligTags[,i],k], des_w)[1] / svymean(predFeats[!ligTags[,i],k], des_wo)[1] # Fold change in weighted means
    
    
  }

  stats_weighted[,grepl('_adj$', colnames(stats_weighted))][,i] = p.adjust(stats_weighted[grepl('_p$', colnames(stats_weighted))][,i], method = "BH") # Benjamini-Hochberg MHT correction (FDR)
}

stats_weighted[,grepl('_adj$', colnames(stats_weighted))][stats_weighted[,grepl('_adj$', colnames(stats_weighted))] < 1e-16] <- 1e-16


# Plot volcanoes from weighted WMW test
dev.off()
pdf(file = paste('./manuscript/figures/', 
                 'weighted_commonEffect_volcanoes',
                 '.pdf', sep = ''),
    width = 21,
    height = 13)
par(mfrow=c(3,5))
xLim = c(-0.5,0.5)
# yLim = c(0,10)
for(i in 1:ncol(ligTags)){
  
  yLim = c(0,max(-log10(0.1), -log10(min(stats_weighted[,grepl('_adj$', colnames(stats_weighted))][,i]))) + 1)
  
  tag = stats_weighted[,grepl('_adj$', colnames(stats_weighted))][,i] < 0.01
  
  # dev.off()
  plot(0,0,axes = F, main = '', xlab = '', ylab = '', pch = NA)
  bg = "seashell2"
  fg = "ivory"
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = bg)
  # abline(v = c(-1,-.5,0,.5,1), lwd = 6, col = fg)
  # abline(v = c(-1.25,-.75,-.25,.25,.75,1.25), lwd = 3, col = fg)
  # abline(h = c(0,1,2,3,4,5,6), lwd = 6, col = fg)
  # abline(h = c(0.5,1.5,2.5,3.5,4.5,5.5,6.5), lwd = 3, col = fg)
  par(new=T)
  
  plot(stats_weighted[,grepl('_effectSize$', colnames(stats_weighted))][,i], -log10(stats_weighted[,grepl('_adj$', colnames(stats_weighted))][,i]), # Plot all points w/ color @ alpha 0.5
       xlab = "Effect size", ylab = "-log10(FDR)", main = colnames(ligTags)[i],
       pch = 19, cex = 2, col = alpha(featColors, 0.5),
       cex.axis = 1.5, cex.main = 2, cex.lab = 1.5,
       xlim = xLim, ylim = yLim)
  par(new=T)
  plot(stats_weighted[,grepl('_effectSize$', colnames(stats_weighted))][tag,i], -log10(stats_weighted[,grepl('_adj$', colnames(stats_weighted))][tag,i]), # Plot stat sig points again with alpha 1
       pch = 19, col = featColors[tag], cex = 2,
       axes = F, xlab = "", ylab = "", main = "",
       xlim = xLim, ylim = yLim)
  par(new=T)
  plot(stats_weighted[,grepl('_effectSize$', colnames(stats_weighted))][tag,i], -log10(stats_weighted[,grepl('_adj$', colnames(stats_weighted))][tag,i]), # Outline stat sig points in black to highlight them
       col = 'black', cex = 2.05,
       xlab = "", ylab = "", axes = F, main = "",
       xlim = xLim, ylim = yLim)
  abline(h= -log10(0.01), lwd = 2)
  abline(v = 0, lty=2, lwd = 2)
}

dev.off()
plot(0,0,axes = F, main = '', xlab = '', ylab = '', pch = NA)
legend(x = 'center',
       col = c('forestgreen', resiFeats, pocketFeats),
       legend = c('PLIP interactions',
                  'bs resi cnt', 'sec struct', 'aa type', 'aa ident',
                  'D2 feats', 'D2 PCs', 'Zern PCs'),
       pch = 19)

dev.off()



# Volcano plots with weighted mean fold change
dev.off()
pdf(file = paste('./manuscript/figures/', 
                 'weighted_FC_volcanoes',
                 '.pdf', sep = ''),
    width = 21,
    height = 13)
par(mfrow=c(3,5))
xLim = c(-10,10)
# yLim = c(0,10)
for(i in 1:ncol(ligTags)){
  
  yLim = c(0,max(-log10(0.1), -log10(min(stats_weighted[,grepl('_adj$', colnames(stats_weighted))][,i]))) + 1)
  
  tag = stats_weighted[,grepl('_adj$', colnames(stats_weighted))][,i] < 0.01
  
  # dev.off()
  plot(0,0,axes = F, main = '', xlab = '', ylab = '', pch = NA)
  bg = "seashell2"
  fg = "ivory"
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = bg)
  # abline(v = c(-1,-.5,0,.5,1), lwd = 6, col = fg)
  # abline(v = c(-1.25,-.75,-.25,.25,.75,1.25), lwd = 3, col = fg)
  # abline(h = c(0,1,2,3,4,5,6), lwd = 6, col = fg)
  # abline(h = c(0.5,1.5,2.5,3.5,4.5,5.5,6.5), lwd = 3, col = fg)
  par(new=T)
  
  plot(stats_weighted[,grepl('_FC$', colnames(stats_weighted))][,i], -log10(stats_weighted[,grepl('_adj$', colnames(stats_weighted))][,i]), # Plot all points w/ color @ alpha 0.5
       xlab = "log2(FC in weighted means)", ylab = "-log10(FDR)", main = colnames(ligTags)[i],
       pch = 19, cex = 2, col = alpha(featColors, 0.5),
       cex.axis = 1.5, cex.main = 2, cex.lab = 1.5,
       xlim = xLim, ylim = yLim)
  par(new=T)
  plot(stats_weighted[,grepl('_FC$', colnames(stats_weighted))][tag,i], -log10(stats_weighted[,grepl('_adj$', colnames(stats_weighted))][tag,i]), # Plot stat sig points again with alpha 1
       pch = 19, col = featColors[tag], cex = 2,
       axes = F, xlab = "", ylab = "", main = "",
       xlim = xLim, ylim = yLim)
  par(new=T)
  plot(stats_weighted[,grepl('_FC$', colnames(stats_weighted))][tag,i], -log10(stats_weighted[,grepl('_adj$', colnames(stats_weighted))][tag,i]), # Outline stat sig points in black to highlight them
       col = 'black', cex = 2.05,
       xlab = "", ylab = "", axes = F, main = "",
       xlim = xLim, ylim = yLim)
  abline(h= -log10(0.01), lwd = 2)
  abline(v = 0, lty=2, lwd = 2)
}

dev.off()





# Compare weighted test to mean-based approach
par(mfrow=c(3,5))
xLim = c(0,1)
yLim = c(0,1)
for (i in 1:ncol(ligTags)){
  plot(stats_weighted[,grepl('_adj$', colnames(stats_weighted))][,i], stats[,grepl('_adj$', colnames(stats))][,i],
       xlab = 'Weighted WMW test adjusted p-values', ylab = 'Cluster-means WMW test adjusted p-values',
       main = colnames(ligTags)[i],
       xlim = xLim, ylim = yLim)
  abline(a = 0, b = 1)
  text(x=0.6, y=0.05, labels = paste( 'Spearman corr = ', round(cor(stats_weighted[,grepl('_adj$', colnames(stats_weighted))][,i], stats[,grepl('_adj$', colnames(stats))][,i], method = 'spearman'), 4), sep = ''), cex = 2)
}

par(mfrow=c(3,5))
xLim = c(0,1)
yLim = c(0,1)
for (i in 1:ncol(ligTags)){
  plot(0.5 + stats_weighted[,grepl('_effectSize$', colnames(stats_weighted))][,i], stats[,grepl('_effectSize$', colnames(stats))][,i],
       xlab = 'Weighted WMW test effect size', ylab = 'Cluster-means WMW test effect size',
       main = colnames(ligTags)[i],
       xlim = xLim, ylim = yLim)
  abline(a = 0, b = 1)
  abline(h=0.5, lty =2)
  abline(v=0.5, lty =2)
  text(x=0.6, y=0.05, labels = paste( 'Pearson corr = ', round(cor(0.5 + stats_weighted[,grepl('_adj$', colnames(stats_weighted))][,i], stats[,grepl('_adj$', colnames(stats))][,i], method = 'spearman'), 4), sep = ''), cex = 2)
}


# # Plot volcano plots from scaled features
# dev.off()
# pdf(file = paste('./analysis/sfgPlots/', 
#                  'medFC_volcanoes',
#                  '.pdf', sep = ''),
#     width = 21,
#     height = 13)
# par(mfrow=c(3,5))
# xLim = c(-10,10)
# # yLim = c(0,10)
# for(i in 1:ncol(ligTags)){
#   
#   yLim = c(0,max(-log10(0.1), -log10(min(scaled_stats[,grepl('_adj$', colnames(scaled_stats))][,i]))) + 1)
#   
#   tag = scaled_stats[,grepl('_adj$', colnames(scaled_stats))][,i] < 0.1
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
#   par(new=T)
#   
#   plot(log2(scaled_stats[,grepl('_FC$', colnames(scaled_stats))][,i]), -log10(scaled_stats[,grepl('_adj$', colnames(scaled_stats))][,i]), # Plot all points w/ color @ alpha 0.5
#        xlab = "log2(median FC)", ylab = "-log10(FDR)", main = colnames(ligTags)[i],
#        pch = 19, cex = 2, col = alpha(featColors, 0.5),
#        cex.axis = 1.5, cex.main = 2, cex.lab = 1.5,
#        xlim = xLim, ylim = yLim)
#   par(new=T)
#   plot(log2(scaled_stats[,grepl('_FC$', colnames(scaled_stats))][tag,i]), -log10(scaled_stats[,grepl('_adj$', colnames(scaled_stats))][tag,i]), # Plot stat sig points again with alpha 1
#        pch = 19, col = featColors[tag], cex = 2,
#        axes = F, xlab = "", ylab = "", main = "",
#        xlim = xLim, ylim = yLim)
#   par(new=T)
#   plot(log2(scaled_stats[,grepl('_FC$', colnames(scaled_stats))][tag,i]), -log10(scaled_stats[,grepl('_adj$', colnames(scaled_stats))][tag,i]), # Outline stat sig points in black to highlight them
#        col = 'black', cex = 2.05,
#        xlab = "", ylab = "", axes = F, main = "",
#        xlim = xLim, ylim = yLim)
#   abline(h= -log10(0.1), lwd = 2)
#   abline(v = 0, lty=2, lwd = 2)
# }
# dev.off()

# Correlation between features for each ligand
pdf(file = paste('./analysis/sfgPlots/', 
                 'allFeat_MWM_corplot',
                 '.pdf', sep = ''),
    width = 11.5,
    height = 8.25)
breakLst = seq(-1,1,0.05)
pheatmap(cor(stats[,grepl('_effectSize$', colnames(stats))], stats[,grepl('_effectSize$', colnames(stats))], method = 'spearman'),
         color = colorRampPalette(c("royalblue1", "grey90", "gold1"))(length(breakLst)),
         labels_row = gsub('_effectSize$', '', colnames(stats[,grepl('_effectSize$', colnames(stats))])),
         main = 'Spearman correlation between feature-specific effect sizes across ligand classes',
         breaks = breakLst,
         show_colnames = F,
         cutree_rows = 1,
         treeheight_col = 0)
dev.off()

# Correlations by feature type
pdf(file = paste('./analysis/sfgPlots/', 
                 'resiFeat_MWM_corplot',
                 '.pdf', sep = ''),
    width = 11.5,
    height = 8.25)
breakLst = seq(-1,1,0.05)
pheatmap(cor(stats[resiFeatTag,grepl('_effectSize$', colnames(stats))], scaled_stats[resiFeatTag,grepl('_effectSize$', colnames(scaled_stats))], method = 'spearman'),
         color = colorRampPalette(c("royalblue1", "grey90", "gold1"))(length(breakLst)),
         labels_row = gsub('_effectSize$', '', colnames(stats[,grepl('_effectSize$', colnames(stats))])),
         main = 'Spearman correlation between RESIDUE feature effect sizes across ligand classes',
         breaks = breakLst,
         show_colnames = F,
         treeheight_col = 0, cutree_rows = 4)
dev.off()
pdf(file = paste('./analysis/sfgPlots/', 
                 'pocketFeat_MWM_corplot',
                 '.pdf', sep = ''),
    width = 11.5,
    height = 8.25)
pheatmap(cor(stats[pocketFeatTag,grepl('_effectSize$', colnames(stats))], scaled_stats[pocketFeatTag,grepl('_effectSize$', colnames(scaled_stats))], method = 'spearman'),
         color = colorRampPalette(c("royalblue1", "grey90", "gold1"))(length(breakLst)),
         labels_row = gsub('_effectSize$', '', colnames(stats[,grepl('_effectSize$', colnames(stats))])),
         main = 'Spearman correlation between POCKET feature effect sizes across ligand classes',
         breaks = breakLst,
         show_colnames = F,
         treeheight_col = 0)
dev.off()

pdf(file = paste('./analysis/sfgPlots/', 
                 'intFeat_MWM_corplot',
                 '.pdf', sep = ''),
    width = 11.5,
    height = 8.25)
pheatmap(cor(stats[1:11,grepl('_effectSize$', colnames(stats))], scaled_stats[1:11,grepl('_effectSize$', colnames(scaled_stats))], method = 'spearman'),
         color = colorRampPalette(c("royalblue1", "grey90", "gold1"))(length(breakLst)),
         labels_row = gsub('_effectSize$', '', colnames(stats[,grepl('_effectSize$', colnames(stats))])),
         main = 'Spearman correlation between INTERACTION feature effect sizes across ligand classes',
         breaks = breakLst,
         show_colnames = F,
         treeheight_col = 0)
dev.off()

# Shared significant features
# Sialic acid features
siaBindingFeats = list(row.names(stats)[stats$Sialic_Acid_adj < 0.1], row.names(stats)[stats$NeuAc_adj < 0.1], row.names(stats)[stats$`NeuAc(a2-3)Gal(b1-4)Glc_adj` < 0.1])
venn.diagram(
  x = siaBindingFeats,
  category.names = c("Sialic acid" , "NeuAc" , "NeuAc(a2-3)Gal(b1-4)Glc"),
  filename = './analysis/sfgPlots/neuac_feats.png',
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
intersect(intersect(siaBindingFeats[[1]], siaBindingFeats[[2]]), siaBindingFeats[[3]])

# Fucose binding
fucBindingFeats = list(row.names(stats)[stats$Fuc_adj < 0.1], row.names(stats)[stats$Terminal_Fucose_adj < 0.1])
intersect(fucBindingFeats[[1]], fucBindingFeats[[2]])
dev.off()
draw.pairwise.venn(area1 = length(fucBindingFeats[[1]]),
                   area2 = length(fucBindingFeats[[2]]),
                   cross.area = length(intersect(fucBindingFeats[[1]], fucBindingFeats[[2]])),
                   euler.d = T, scaled = T, ext.text = F, cex = 2,
                   category = c('Fucose (monosacc.)', 'Terminal fucose'), cat.cex = 2,
                   cat.pos = c(0,0),
                   fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
                   col = c("#440154ff", '#21908dff'),
                   cat.col = c("#440154ff", '#21908dff'),
                   cat.fontfamily = "sans",
                   fontfamily = "sans"
)

# Mannose binding
manBindingFeats = list(row.names(stats)[stats$Man_adj < 0.1], row.names(stats)[stats$High_Mannose_adj < 0.1], row.names(stats)[stats$`Man(a1-2)Man_adj` < 0.1])
intersect(intersect(manBindingFeats[[1]], manBindingFeats[[2]]), manBindingFeats[[3]])
dev.off()
venn.diagram(
  x = manBindingFeats,
  category.names = c("Man" , "High_mannose" , "Man(a1-2)Man"),
  filename = './analysis/sfgPlots/man_feats.png',
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




d2FeatCor = cor(d2Feats)
corrplot(d2FeatCor)
plot(density(d2FeatCor[upper.tri(d2FeatCor,diag = F)]),
     main = 'Density distribution of pairwise correlations b/w d2Feats',
     xlab = 'Pearson correlation')
abline(v = c(0.7, -0.7))

sum(abs(d2FeatCor[upper.tri(d2FeatCor,diag = F)]) > 0.7) / length(d2FeatCor[upper.tri(d2FeatCor,diag = F)])

subD2Feats = d2Feats[,grepl('_6Ang', colnames(d2Feats))]

corrplot(cor(subD2Feats))
sum(abs(cor(subD2Feats)[upper.tri(cor(subD2Feats),diag = F)]) > 0.7) / length(cor(subD2Feats)[upper.tri(cor(subD2Feats),diag = F)])