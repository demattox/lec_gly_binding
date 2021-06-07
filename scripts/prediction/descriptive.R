
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

# # Load McCaldon aa frequencies
# mccaldon = read.delim('~/episweep/data/mccaldon.csv', header = F, sep = ",", stringsAsFactors = F)
# colnames(mccaldon) = c("aa", 'mcFreq')

ligTags <- read.delim(file = './analysis/training/data_in/ligTags.tsv', sep = '\t', stringsAsFactors = F)
predFeats <- read.delim(file = './analysis/training/data_in/predFeats.csv', sep = ',', stringsAsFactors = F)
bsResiDat <- read.delim(file = './analysis/training/data_in/bsResiDat.tsv', sep = '\t', stringsAsFactors = F)

all(row.names(bsResiDat) == row.names(predFeats))
allOut = cbind(bsResiDat[,1:13], predFeats)
write.csv(allOut, 'interactionData/suppFile1.csv', quote = T)


uniLigs = unique(bsResiDat$iupac)

scaledFeats = predFeats
for(i in 1:ncol(scaledFeats)){
  scaledFeats[,i] = (scaledFeats[,i] - min(scaledFeats[,i])) / (max(scaledFeats[,i]) - min(scaledFeats[,i]))
}

# unique(bsResiDat$iupac[apply(ligTags, 1, sum) == 2])

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

###########################
## McCaldon analysis
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
## Mean-based WMW test
###########################

# Add tags to bsResiDat
# bsResiDat = cbind(bsResiDat, ligTags)


# 
# Statistical test for cluster independent feature enrichment by ligand
# stats = as.data.frame(matrix(0, nrow = ncol(predFeats), ncol = (ncol(ligTags)*4)))
# row.names(stats) = colnames(predFeats)
# colnames(stats) = c(paste(colnames(ligTags), 'p', sep = '_'), paste(colnames(ligTags), '_FC', sep = '_'), paste(colnames(ligTags), 'effectSize', sep = '_'), paste(colnames(ligTags), 'adj', sep = '_'))
# 
# scaled_stats = stats # In parallel, try this with scaled features for more robust change metric, shouldn't impact WMW test
# 
# for (i in 1:ncol(ligTags)){ # for ligand i
#   # meanClust50_data = as.data.frame(matrix(NA, nrow = length(clusLst50), ncol = (ncol(predFeats)*2))) # Rows for unique clusters, columns for predictive features from bsites WITH or WITHOUT current ligand
#   # colnames(meanClust50_data) = c(paste(colnames(predFeats), 'with', sep = '_'),  paste(colnames(predFeats), 'without', sep = '_'))
#   for(k in 1: ncol(predFeats)){ # for each feature k
#     meanClust50_data = as.data.frame(matrix(NA, nrow = length(clusLst50), ncol = 2)) # Rows for unique clusters, columns for predictive feature from bsites WITH or WITHOUT current ligand
#     colnames(meanClust50_data) = c('with', 'without')
#     scaledClus50_data = meanClust50_data
#     for (j in 1:length(clusLst50)){ # for each cluster j
#       clustTag = bsResiDat$seqClust50 == clusLst50[j]
#       withTag = clustTag & ligTags[,i]
#       withoutTag = clustTag & !ligTags[,i]
#       if ( sum(withTag) >= 1){
#         meanClust50_data[j, 1] = mean( predFeats[withTag,k] )
#         scaledClus50_data[j, 1] = mean( scaledFeats[withTag,k] )
#       }
#       if ( sum(withoutTag) >= 1){
#         meanClust50_data[j, 2] = mean( predFeats[withoutTag,k] )
#         scaledClus50_data[j, 2] = mean( scaledFeats[withoutTag,k] )
#       }
#     }
#     # Raw Feats
#     wDat = meanClust50_data[,1]
#     wDat = wDat[!is.na(wDat)]
#     woDat = meanClust50_data[,2]
#     woDat = woDat[!is.na(woDat)]
#     # Scaled feats
#     sW = scaledClus50_data[!is.na(scaledClus50_data[,1]),1] # Non-NA feats from bSites w/ ligand
#     sWO = scaledClus50_data[!is.na(scaledClus50_data[,2]),2] # Non-NA feats from bSites w/ ligand
# 
#     # ligTest = t.test(wDat[!is.na(wDat)], woDat[!is.na(woDat)])
#     ligTest = wilcox.test(x= wDat, y = woDat) # Wilcoxon–Mann–Whitney test, sample sizes can be small (~5% of 230 clusters ~= 10), no reason to assume distribution is normal as it likely isn't, as well as sample sizes often being too small to test for normality
#     sTest = wilcox.test(x = sW, y = sWO)
# 
#     stats[k, grepl('_p$', colnames(stats))][i] = ligTest$p.value # Raw p-value
#     stats[k, grepl('_effectSize$', colnames(stats))][i] = ligTest$statistic / (length(wDat) * length(woDat)) # common language effect size
#     scaled_stats[k, grepl('_p$', colnames(scaled_stats))][i] = sTest$p.value # Raw p-value
#     scaled_stats[k, grepl('_effectSize$', colnames(scaled_stats))][i] = sTest$statistic / (length(sW) * length(sWO)) # common language effect size
# 
#     if (median(woDat) != 0){
#       stats[k, grepl('_FC$', colnames(stats))][i] = median(wDat) / median(woDat) # Fold change in mean values
#     } else{
#       stats[k, grepl('_FC$', colnames(stats))][i] = median(wDat)
#     }
#     if (median(sWO) != 0){
#       scaled_stats[k, grepl('_FC$', colnames(scaled_stats))][i] = median(sW) / median(sWO) # Fold change in mean values
#     } else{
#       scaled_stats[k, grepl('_FC$', colnames(scaled_stats))][i] = median(sW)
#     }
#   }
#   stats[,grepl('_adj$', colnames(stats))][,i] = p.adjust(stats[grepl('_p$', colnames(stats))][,i], method = "BH") # Benjamini-Hochberg MHT correction (FDR)
#   scaled_stats[,grepl('_adj$', colnames(scaled_stats))][,i] = p.adjust(scaled_stats[grepl('_p$', colnames(scaled_stats))][,i], method = "BH") # Benjamini-Hochberg MHT correction (FDR)
# }
# 
# 
# for (i in 1:ncol(ligTags)){
#   lig = colnames(ligTags)[i]
#   cat(lig)
#   cat('\n----\n')
#   cat('DOWN w/ ligand\n')
#   cat(row.names(stats)[stats[,grepl('_adj$', colnames(stats))][,i] <= 0.1  &  stats[,grepl('_effectSize$', colnames(stats))][,i] < 0.5])
#   cat('\nUP w/ ligand\n')
#   cat(row.names(stats)[stats[,grepl('_adj$', colnames(stats))][,i] <= 0.1  &  stats[,grepl('_effectSize$', colnames(stats))][,i] > 0.5])
#   cat('\n----------------\n')
# }




# Plot volcano plots from raw features
# 
# dev.off()
# pdf(file = paste('./analysis/sfgPlots/', 
#                  'commonEffect_volcanoes',
#                  '.pdf', sep = ''),
#     width = 21,
#     height = 13)
# par(mfrow=c(3,5))
# xLim = c(0,1)
# # yLim = c(0,10)
# for(i in 1:ncol(ligTags)){
#   
#   yLim = c(0,max(-log10(0.1), -log10(min(stats[,grepl('_adj$', colnames(stats))][,i]))) + 1)
#   
#   tag = stats[,grepl('_adj$', colnames(stats))][,i] < 0.1
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
#   plot(stats[,grepl('_effectSize$', colnames(stats))][,i], -log10(stats[,grepl('_adj$', colnames(stats))][,i]), # Plot all points w/ color @ alpha 0.5
#        xlab = "Effect size", ylab = "-log10(FDR)", main = colnames(ligTags)[i],
#        pch = 19, cex = 2, col = alpha(featColors, 0.5),
#        cex.axis = 1.5, cex.main = 2, cex.lab = 1.5,
#        xlim = xLim, ylim = yLim)
#   par(new=T)
#   plot(stats[,grepl('_effectSize$', colnames(stats))][tag,i], -log10(stats[,grepl('_adj$', colnames(stats))][tag,i]), # Plot stat sig points again with alpha 1
#        pch = 19, col = featColors[tag], cex = 2,
#        axes = F, xlab = "", ylab = "", main = "",
#        xlim = xLim, ylim = yLim)
#   par(new=T)
#   plot(stats[,grepl('_effectSize$', colnames(stats))][tag,i], -log10(stats[,grepl('_adj$', colnames(stats))][tag,i]), # Outline stat sig points in black to highlight them
#        col = 'black', cex = 2.05,
#        xlab = "", ylab = "", axes = F, main = "",
#        xlim = xLim, ylim = yLim)
#   abline(h= -log10(0.1), lwd = 2)
#   abline(v = 0.5, lty=2, lwd = 2)
# }
# 
# dev.off()
# plot(0,0,axes = F, main = '', xlab = '', ylab = '', pch = NA)
# legend(x = 'center',
#        col = c('forestgreen', resiFeats, pocketFeats),
#        legend = c('PLIP interactions',
#                   'bs resi cnt', 'sec struct', 'aa type', 'aa ident',
#                   'D2 feats', 'D2 PCs', 'Zern PCs'),
#        pch = 19)
# 
# dev.off()
# par(mfrow=c(3,5))
# xLim = c(0,1)
# yLim = c(0,1)
# # Compare raw feats to scaled feats
# for (i in 1:ncol(ligTags)){
#   plot(stats[,grepl('_adj$', colnames(stats))][,i], scaled_stats[,grepl('_adj$', colnames(scaled_stats))][,i],
#        xlab = 'Raw feat adjusted p-values', ylab = 'Scaled feat adjusted p-values',
#        main = colnames(ligTags)[i],
#        xlim = xLim, ylim = yLim)
#   abline(a = 0, b = 1)
#   text(x=0.6, y=0.05, labels = paste( 'Pear. corr = ', cor(stats[,grepl('_adj$', colnames(stats))][,i], scaled_stats[,grepl('_adj$', colnames(scaled_stats))][,i]), sep = ''), cex = 2)
# }

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


for (i in 1:ncol(ligTags)){
  cat(colnames(ligTags)[i], '\n')
  cat('\t', sum(cuWeights[ligTags[,i]]), '\n')
}


wmwFeats = cbind(predFeats,ligTags)

des <- svydesign(ids = ~1, data = wmwFeats, weights = 1/cuWeights)
# des_alt <- svydesign(ids = ~1, data = wmwFeats, weights = 1/cWeights)
# des_unweighted <- svydesign(ids = ~1, data = wmwFeats, weights = NULL)


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

######################
scaledWMW_feats = cbind(scaledFeats,ligTags)

scaled_des <- svydesign(ids = ~1, data = scaledWMW_feats, weights = 1/cuWeights)

ligSpefic_feat_means = as.data.frame(matrix(0,nrow = ncol(predFeats), ncol = ncol(ligTags)))
colnames(ligSpefic_feat_means) = colnames(ligTags)
row.names(ligSpefic_feat_means) = colnames(predFeats)

# fuc = rep(0,ncol(scaledFeats))
# names(fuc) = colnames(scaledFeats)



######################


stats_weighted = as.data.frame(matrix(0, nrow = ncol(predFeats), ncol = (ncol(ligTags)*4)))
row.names(stats_weighted) = colnames(predFeats)
colnames(stats_weighted) = c(paste(colnames(ligTags), 'p', sep = '_'), paste(colnames(ligTags), '_FC', sep = '_'), paste(colnames(ligTags), 'effectSize', sep = '_'), paste(colnames(ligTags), 'adj', sep = '_'))

for (i in 1:ncol(ligTags)){ # for ligand i
  des_w = subset(des, subset = ligTags[,i]) # temporary design object holding all interactions with the ligand of interest
  des_wo = subset(des, subset = !ligTags[,i]) # same as above but WITHOUT the ligand of interest
  

  scaled_des_w = subset(scaled_des, subset = ligTags[,i]) # 

  for(k in 1:ncol(predFeats)){ # for each feature k
    ligTest = svyranktest(formula = as.formula(paste(colnames(predFeats)[k], ' ~ ', colnames(ligTags)[i], sep = '')), # Wilcoxon–Mann–Whitney test, sample sizes can be small (~5% of 230 clusters ~= 10), no reason to assume distribution is normal as it likely isn't
                           design = des, 
                           test = 'wilcoxon') 

    if (ligTest$p.value != 0){
      stats_weighted[k, grepl('_p$', colnames(stats_weighted))][i] = ligTest$p.value # Raw p-value
    } else{
      stats_weighted[k, grepl('_p$', colnames(stats_weighted))][i] = 5e-324 # If p-value is too small and R rounds to 0, reset as the smallest positive double as referenced by "?.Machine"
    }

    stats_weighted[k, grepl('_effectSize$', colnames(stats_weighted))][i] = ligTest$estimate # common language effect size (0.5 transformed to 0)
    
    stats_weighted[k, grepl('_FC$', colnames(stats_weighted))][i] = svymean(predFeats[ligTags[,i],k], des_w)[1] / svymean(predFeats[!ligTags[,i],k], des_wo)[1] # Fold change in weighted means
    
    ligSpefic_feat_means[k,i] = svymean(scaledFeats[ligTags[,i],k], scaled_des_w)[1] # Get the weighted mean for each feature in interactions with each glycan of interest
  }

  stats_weighted[,grepl('_adj$', colnames(stats_weighted))][,i] = p.adjust(stats_weighted[grepl('_p$', colnames(stats_weighted))][,i], method = "BH") # Benjamini-Hochberg MHT correction (FDR)
}

superSigTag = stats_weighted[,grepl('_adj$', colnames(stats_weighted))] < 1e-16
stats_weighted[,grepl('_adj$', colnames(stats_weighted))][superSigTag] <- 10**(-1*runif(sum(superSigTag), max = -log10(3e-19), min = -log10(1e-16))) # Sample from a log-uniform distribution

####################

colnames(ligSpefic_feat_means) = gsub('$', '_weightedFeatureMean', colnames(ligSpefic_feat_means))

# all(row.names(ligSpefic_feat_means) == row.names(stats_weighted))
stats_weighted = cbind(stats_weighted, ligSpefic_feat_means)

####################



# Colors to features for plotting
featColors = rep('', nrow(stats_weighted))
resiFeats = colorRampPalette(c("plum1","tomato", "firebrick4"))(4)
pocketFeats = colorRampPalette(c('paleturquoise3', 'deepskyblue', 'mediumblue'))(4)

featColors[1:11] = 'forestgreen'

featColors[grep('^vol_4Ang$', row.names(stats_weighted)) : grep('^leftskew_10Ang$', row.names(stats_weighted))] = pocketFeats[2] # features within the d2Feats range
featColors[grepl('^vol_', row.names(stats_weighted)) | grepl('^pcntSurf_', row.names(stats_weighted))] = pocketFeats[1] # General pocket descriptors
featColors[grepl('^binnedD2', row.names(stats_weighted))] = pocketFeats[3] # PCs from the binned D2 measures
featColors[grepl('^zern', row.names(stats_weighted))] = pocketFeats[4] # PCs from the 3DZDs

featColors[grepl('^numBSresis', row.names(stats_weighted))] = resiFeats[1] # number of residues in binding site features
featColors[gsub('_bin\\d{1}', '', row.names(stats_weighted)) %in% c('H', 'B', 'E', 'G', 'T', 'S', 'X.')] = resiFeats[2] # secondary structure features
featColors[gsub('_bin\\d{1}', '', row.names(stats_weighted)) %in% c('nonpolar', 'polar', 'posCharge', 'negCharge', 'aromatic')] = resiFeats[3] # amino acid properties
featColors[grepl('^[[:upper:]]{3}_', row.names(stats_weighted)) | grepl('^CA$', row.names(stats_weighted))] = resiFeats[4] # amino acid identities

resiFeatTag = featColors %in% resiFeats
pocketFeatTag = featColors %in% pocketFeats


# Plot volcanoes from weighted WMW test
dev.off()  
pdf(file = paste('./manuscript/figures/subplots/', 
                 'weighted_commonEffect_volcanoes',
                 '.pdf', sep = ''),
    width = 21,
    height = 13)
# par(mfrow=c(3,5))

par(mfrow=c(3,5), 
    mai = c(0.51, 0.11, 0.11, .51),
    mar = c(3.1, 3.1, 4.1, 2.1))

xLim = c(-0.5,0.5)
yLim = c(0,19)
for(i in 1:ncol(ligTags)){
  
  # yLim = c(0,max(-log10(0.1), -log10(min(stats_weighted[,grepl('_adj$', colnames(stats_weighted))][,i]))) + 1)
  
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
  
  abline(v = 0, lty=2, lwd = 6, col = 'white')
  
  par(new=T)
  
  plot(stats_weighted[,grepl('_effectSize$', colnames(stats_weighted))][,i], -log10(stats_weighted[,grepl('_adj$', colnames(stats_weighted))][,i]), # Plot all points w/ color @ alpha 0.5
       xlab = "", ylab = "", main = '',
       pch = 19, cex = 2, col = alpha(featColors, 0.33),
       cex.axis = 1.5,
       xlim = xLim, ylim = yLim)
  
  abline(h= -log10(0.01), lwd = 6, col = 'white')
  abline(h = -log10(1e-16), lwd = 3, lty = 2, col = 'white')
  
  title(main = ligNames[i], col.main = ligColors[i], cex.main = 2.5, font.main  = 2)
  # mtext(ligNames[i], col = ligColors[i],
  #       side=3, adj=0, outer = T,
  #       line=1.2, cex=1, font=2)

  par(new=T)
  plot(stats_weighted[,grepl('_effectSize$', colnames(stats_weighted))][tag,i], -log10(stats_weighted[,grepl('_adj$', colnames(stats_weighted))][tag,i]), # Plot stat sig points again with alpha 1
       pch = 19, col = featColors[tag], cex = 2,
       axes = F, xlab = "", ylab = "", main = "",
       xlim = xLim, ylim = yLim)
  par(new=T)
  plot(stats_weighted[,grepl('_effectSize$', colnames(stats_weighted))][tag,i], -log10(stats_weighted[,grepl('_adj$', colnames(stats_weighted))][tag,i]), # Outline stat sig points in black to highlight them
       col = alpha('black',0.2), cex = 2.05,
       xlab = "", ylab = "", axes = F, main = "",
       xlim = xLim, ylim = yLim)
  

}

dev.off()
# par(mar = c(10,5,10,5))
plot(0,0,axes = F, main = '', xlab = '', ylab = '', pch = NA)
legend(x = 'center',
       col = c('forestgreen', 'white','white',
               c(rbind(resiFeats,rep('white', length(resiFeats)))), 'white',
               c(rbind(pocketFeats,rep('white', length(pocketFeats)))), 'white'),
       legend = c('PLIP interaction\ncounts', '','',
                  'Residue counts', '', 'Residue sec.\nstructure', '', 'Amino acid\nproperties', '', 'Amino acid\nidentities', '','',
                  'Pocket descriptors', '', 'D2 distribution\nstatistics', '','D2 Principal\nComponents', '', '3DZD Principal\nComponents', ''),
       pch = 19,
       pt.cex = 2.5)

dev.off()




# Volcano plots with weighted mean fold change

# dev.off()
# pdf(file = paste('./manuscript/figures/subplots/', 
#                  'weighted_FC_volcanoes',
#                  '.pdf', sep = ''),
#     width = 21,
#     height = 13)
# par(mfrow=c(3,5))
# xLim = c(-10,10)
# # yLim = c(0,10)
# for(i in 1:ncol(ligTags)){
#   
#   yLim = c(0,max(-log10(0.1), -log10(min(stats_weighted[,grepl('_adj$', colnames(stats_weighted))][,i]))) + 1)
#   
#   tag = stats_weighted[,grepl('_adj$', colnames(stats_weighted))][,i] < 0.01
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
#   plot(stats_weighted[,grepl('_FC$', colnames(stats_weighted))][,i], -log10(stats_weighted[,grepl('_adj$', colnames(stats_weighted))][,i]), # Plot all points w/ color @ alpha 0.5
#        xlab = "log2(FC in weighted means)", ylab = "-log10(FDR)", main = colnames(ligTags)[i],
#        pch = 19, cex = 2, col = alpha(featColors, 0.5),
#        cex.axis = 1.5, cex.main = 2, cex.lab = 1.5,
#        xlim = xLim, ylim = yLim)
#   par(new=T)
#   plot(stats_weighted[,grepl('_FC$', colnames(stats_weighted))][tag,i], -log10(stats_weighted[,grepl('_adj$', colnames(stats_weighted))][tag,i]), # Plot stat sig points again with alpha 1
#        pch = 19, col = featColors[tag], cex = 2,
#        axes = F, xlab = "", ylab = "", main = "",
#        xlim = xLim, ylim = yLim)
#   par(new=T)
#   plot(stats_weighted[,grepl('_FC$', colnames(stats_weighted))][tag,i], -log10(stats_weighted[,grepl('_adj$', colnames(stats_weighted))][tag,i]), # Outline stat sig points in black to highlight them
#        col = 'black', cex = 2.05,
#        xlab = "", ylab = "", axes = F, main = "",
#        xlim = xLim, ylim = yLim)
#   abline(h= -log10(0.01), lwd = 2)
#   abline(v = 0, lty=2, lwd = 2)
# }
# 
# dev.off()




#######################
## Compare weighted test to mean-based approach
#######################
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
  text(x=0.6, y=0.05, labels = paste( 'Pearson corr = ', round(cor(0.5 + stats_weighted[,grepl('_adj$', colnames(stats_weighted))][,i], stats[,grepl('_adj$', colnames(stats))][,i], method = 'pearson'), 4), sep = ''), cex = 2)
}

#######################
## Correlation between features for each ligand
#######################
# Correlations by all features
pdf(file = paste('./manuscript/figures/subplots/', 
                 'allFeats_MWM_corplots',
                 '.pdf', sep = ''),
    width = 11.5,
    height = 8.25)
breakLst = seq(-1,1,0.05)
pheatmap(cor(stats_weighted[,grepl('_effectSize$', colnames(stats_weighted))], stats_weighted[,grepl('_effectSize$', colnames(stats_weighted))], method = 'pearson'),
         color = colorRampPalette(c("royalblue1", "grey90", "gold1"))(length(breakLst)),
         labels_row = gsub('_effectSize$', '', colnames(stats_weighted[,grepl('_effectSize$', colnames(stats_weighted))])),
         main = ' ',
         breaks = breakLst,
         show_colnames = F,
         cutree_rows = 1,
         treeheight_col = 0)
grid.text(label = 'Pearson correlation between feature-specific effect sizes across ligand classes',x = 0.5, y=0.985, gp=gpar(col="black", cex = 1.5))
dev.off()

# Correlations by Residue features
pdf(file = paste('./manuscript/figures/subplots/', 
                 'resiFeat_MWM_corplot',
                 '.pdf', sep = ''),
    width = 11.5,
    height = 8.25)
breakLst = seq(-1,1,0.05)
pheatmap(cor(stats_weighted[resiFeatTag,grepl('_effectSize$', colnames(stats_weighted))], stats_weighted[resiFeatTag,grepl('_effectSize$', colnames(stats_weighted))], method = 'pearson'),
         color = colorRampPalette(c("royalblue1", "grey90", "gold1"))(length(breakLst)),
         labels_row = gsub('_effectSize$', '', colnames(stats_weighted[,grepl('_effectSize$', colnames(stats_weighted))])),
         main = '',
         breaks = breakLst,
         show_colnames = F,
         treeheight_col = 0, 
         cutree_rows = 1)
grid.text(label = 'Pearson correlation between RESIDUE feature effect sizes across ligand classes',x = 0.5, y=0.985, gp=gpar(col="firebrick4", cex = 1.5))
dev.off()

# Correlations by Pocket features
pdf(file = paste('./manuscript/figures/subplots/', 
                 'pocketFeat_MWM_corplot',
                 '.pdf', sep = ''),
    width = 11.5,
    height = 8.25)
pheatmap(cor(stats_weighted[pocketFeatTag,grepl('_effectSize$', colnames(stats_weighted))], stats_weighted[pocketFeatTag,grepl('_effectSize$', colnames(stats_weighted))], method = 'pearson'),
         color = colorRampPalette(c("royalblue1", "grey90", "gold1"))(length(breakLst)),
         labels_row = gsub('_effectSize$', '', colnames(stats_weighted[,grepl('_effectSize$', colnames(stats_weighted))])),
         main = '',
         breaks = breakLst,
         show_colnames = F,
         treeheight_col = 0)
grid.text(label = 'Pearson correlation between POCKET feature effect sizes across ligand classes',x = 0.5, y=0.985, gp=gpar(col="blue2", cex = 1.5))
dev.off()

# Correlations by Interaction features
pdf(file = paste('./manuscript/figures/subplots/', 
                 'intFeat_MWM_corplot',
                 '.pdf', sep = ''),
    width = 11.5,
    height = 8.25)
pheatmap(cor(stats_weighted[1:11,grepl('_effectSize$', colnames(stats_weighted))], stats_weighted[1:11,grepl('_effectSize$', colnames(stats_weighted))], method = 'pearson'),
         color = colorRampPalette(c("royalblue1", "grey90", "gold1"))(length(breakLst)),
         labels_row = gsub('_effectSize$', '', colnames(stats_weighted[,grepl('_effectSize$', colnames(stats_weighted))])),
         main = '',
         breaks = breakLst,
         show_colnames = F,
         treeheight_col = 0)
grid.text(label = 'Pearson correlation between INTERACTION feature effect sizes across ligand classes',x = 0.5, y=0.985, gp=gpar(col="forestgreen", cex = 1.5))

dev.off()

#######################
# Heatmaps with features in columns
#######################

cWidth = 8
cHeight = 20 


## Ligand annotation


# All-feature based correlogram

allFeatCorrs = cor(stats_weighted[,grepl('_effectSize$', colnames(stats_weighted))], stats_weighted[,grepl('_effectSize$', colnames(stats_weighted))], method = 'pearson')
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
         
         annotation_row = r_annot,
         annotation_colors = annot_cols)
# grid.text(label = 'Pearson correlations between feature-specific effect sizes across ligands',x = .415, y=0.985, gp=gpar(col="black", cex = 1.5))
dev.off() 

cor.test(stats_weighted$Gal.b1.4.Glc_effectSize, stats_weighted$Gal.b1.4.GlcNAc_effectSize)


cor.test(stats_weighted$Man.a1.2.Man_effectSize, stats_weighted$Glc_effectSize)
cor.test(stats_weighted$Man_effectSize, stats_weighted$Glc_effectSize)

cor.test(stats_weighted$Man.a1.2.Man_effectSize, stats_weighted$High_Mannose_effectSize)


cor.test(stats_weighted$Sialic_Acid_effectSize, stats_weighted$NeuAc.a2.3.Gal.b1.4.Glc_effectSize)

cor.test(stats_weighted$Sialic_Acid_effectSize, stats_weighted$NeuAc_effectSize)
cor.test(stats_weighted$NeuAc_effectSize, stats_weighted$NeuAc.a2.3.Gal.b1.4.Glc_effectSize)


cor.test(stats_weighted$Man_effectSize, stats_weighted$Fuc_effectSize)
cor.test(stats_weighted$Man_effectSize, stats_weighted$Man.a1.2.Man_effectSize)
cor.test(stats_weighted$Man_effectSize, stats_weighted$Glc_effectSize)

cor.test(stats_weighted$Glc_effectSize, stats_weighted$GalNAc_effectSize)
cor.test(stats_weighted$Gal_effectSize, stats_weighted$GlcNAc_effectSize)

cor.test(stats_weighted$Gal_effectSize, stats_weighted$GalNAc_effectSize)
cor.test(stats_weighted$Glc_effectSize, stats_weighted$GlcNAc_effectSize)

# corrmat = cor(stats_weighted[,grepl('_effectSize$', colnames(stats_weighted))], stats_weighted[,grepl('_effectSize$', colnames(stats_weighted))], method = 'pearson')
# row.names(corrmat) = ligNames
# corrplot(corrmat, order = "hclust", method = "ellipse")




breakLst = seq(-0.5,0.5,0.01)

annot <- data.frame(Feature_Type = rep("", nrow(stats_weighted)))
row.names(annot) = row.names(stats_weighted)

annot$Feature_Type[featColors == 'forestgreen'] <- 'PLIP interaction counts'

annot$Feature_Type[grep('^vol_4Ang$', row.names(stats_weighted)) : grep('^leftskew_10Ang$', row.names(stats_weighted))] <- 'D2 distribution features'
annot$Feature_Type[grepl('^vol_', row.names(stats_weighted)) | grepl('^pcntSurf_', row.names(stats_weighted))] = 'Pocket descriptors' # General pocket descriptors
annot$Feature_Type[grepl('^binnedD2', row.names(stats_weighted))] <- 'D2 Principal Components'
annot$Feature_Type[grepl('^zern', row.names(stats_weighted))] <- '3DZD Principal Components'

annot$Feature_Type[grepl('^numBSresis', row.names(stats_weighted))] <- 'Residue counts/bin'
annot$Feature_Type[gsub('_bin\\d{1}', '', row.names(stats_weighted)) %in% c('H', 'B', 'E', 'G', 'T', 'S', 'X.')] <- 'Residue sec struct.'
annot$Feature_Type[gsub('_bin\\d{1}', '', row.names(stats_weighted)) %in% c('nonpolar', 'polar', 'posCharge', 'negCharge', 'aromatic')] <- 'Amino acid property counts'
annot$Feature_Type[grepl('^[[:upper:]]{3}_', row.names(stats_weighted)) | grepl('^CA$', row.names(stats_weighted))] <- 'Residue identities'

annot$Feature_Type <- factor(annot$Feature_Type, levels = unique(annot$Feature_Type))



Feature_Type <- unique(featColors)
names(Feature_Type) <- levels(annot$Feature_Type)
anno_colors <- list(Feature_Type = Feature_Type)

pheatmap(t(stats_weighted[,grepl('_effectSize$', colnames(stats_weighted))]),
         color = colorRampPalette(c("royalblue1", "grey90", "gold1"))(length(breakLst)),
         clustering_distance_rows = 'correlation',
         # clustering_distance_cols = 'correlation',
         display_numbers = ifelse(t(stats_weighted[,grepl('_adj$', colnames(stats_weighted))]) < 0.01, "*", ""), fontsize_number = 18,
         border_color = 'white',
         labels_row = ligNames,
         annotation_col = annot, annotation_colors = anno_colors,
         main = '',
         breaks = breakLst,
         show_colnames = F)


# Residue features only

resiAnnot = data.frame(Feature_Type = rep("", sum(resiFeatTag)))
row.names(resiAnnot) = row.names(stats_weighted)[resiFeatTag]

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

resiFeat_stats = stats_weighted
row.names(resiFeat_stats)[resiFeatTag] = row.names(resiAnnot)

Feature_Type <- unique(featColors[resiFeatTag])
names(Feature_Type) <- levels(resiAnnot$Feature_Type)

Bin <- c('firebrick3', 'darkorange2', 'darkgoldenrod2', 'gold2')
names(Bin) <- levels(resiAnnot$Bin)

resiAnnot_cols <- list(Feature_Type = Feature_Type, Bin = Bin, Terminal_Sugar = Terminal_Sugar, Sugar_Cnt = Sugar_Cnt)

resiFeat_stats = t(resiFeat_stats[resiFeatTag,grepl('_effectSize$', colnames(resiFeat_stats))])
row.names(resiFeat_stats) = gsub('_effectSize$', '', row.names(resiFeat_stats))

pdf(file = paste('./manuscript/figures/subplots/', 
                 'resiFeats_MWM_heatmap',
                 '.pdf', sep = ''),
    width = 24,
    height = 7)
pheatmap(resiFeat_stats,
         color = colorRampPalette(c("royalblue1", "ivory", "gold1"))(length(breakLst)),
         border_color = 'ivory',
         cellwidth = cWidth,
         cellheight = cHeight,
         clustering_distance_rows = 'correlation',
         # display_numbers = ifelse(t(stats_weighted[,grepl('_adj$', colnames(stats_weighted))]) < 0.01, "*", ""), fontsize_number = 18,
         labels_row = ligNames,
         annotation_col = resiAnnot,
         annotation_row = r_annot,
         annotation_colors = resiAnnot_cols,
         # legend_breaks = c(1,5),
         cutree_rows = 4,
         main = '',
         breaks = breakLst,
         show_colnames = T,
         fontsize_col = 6,
         angle_col = 45)
dev.off()

# Pocket features only

pockAnnot = data.frame(Feature_Type = rep("", sum(pocketFeatTag)))
row.names(pockAnnot) = row.names(stats_weighted)[pocketFeatTag]

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


pocketFeat_stats = t(stats_weighted[pocketFeatTag,grepl('_effectSize$', colnames(stats_weighted))])
row.names(pocketFeat_stats) = gsub('_effectSize$', '', row.names(pocketFeat_stats))

pdf(file = paste('./manuscript/figures/subplots/', 
                 'pocketFeats_MWM_heatmap',
                 '.pdf', sep = ''),
    width = 24,
    height = 7)
pheatmap(pocketFeat_stats,
         color = colorRampPalette(c("royalblue1", "ivory", "gold1"))(length(breakLst)),
         border_color = 'ivory',
         cellwidth = cWidth,
         cellheight = cHeight,
         clustering_distance_rows = 'correlation',
         # display_numbers = ifelse(t(stats_weighted[,grepl('_adj$', colnames(stats_weighted))]) < 0.01, "*", ""), fontsize_number = 18,
         labels_row = ligNames,
         annotation_col = pockAnnot,
         annotation_row = r_annot,
         annotation_colors = pockAnnot_cols,
         main = '',
         cutree_rows = 4,
         breaks = breakLst,
         show_colnames = T,
         fontsize_col = 6,
         angle_col = 45)
dev.off()

# PLIP features only

annot_cols = list(Terminal_Sugar = Terminal_Sugar, Sugar_Cnt = Sugar_Cnt)

plipFeat_stats = t(stats_weighted[1:11,grepl('_effectSize$', colnames(stats_weighted))])
row.names(plipFeat_stats) = gsub('_effectSize$', '', row.names(plipFeat_stats))


pdf(file = paste('./manuscript/figures/subplots/', 
                 'plipFeats_MWM_heatmap',
                 '.pdf', sep = ''),
    width = 24,
    height = 7)
pheatmap(plipFeat_stats,
         color = colorRampPalette(c("royalblue1", "ivory", "gold1"))(length(breakLst)),
         border_color = 'ivory',
         cellwidth = cWidth,
         cellheight = cHeight,
         clustering_distance_rows = 'correlation',
         # display_numbers = ifelse(t(stats_weighted[,grepl('_adj$', colnames(stats_weighted))]) < 0.01, "*", ""), fontsize_number = 18,
         labels_row = ligNames,
         annotation_row = r_annot,
         annotation_colors = annot_cols,
         main = '',
         cutree_rows = 4,
         breaks = breakLst,
         show_colnames = T,
         fontsize_col = 6,
         angle_col = 45)
dev.off()

#######################
# Shared significant features for similar ligands
#######################

write.table(stats_weighted, file = './analysis/training/weightedWMW_stats.tsv', quote = F, sep = '\t')

save(des, file = './analysis/training/surveyObject.RData')

# Sialic acid features
siaBindingFeats_pos = list(row.names(stats_weighted)[stats_weighted$Sialic_Acid_adj < 0.01 & stats_weighted$Sialic_Acid_effectSize > 0], row.names(stats_weighted)[stats_weighted$NeuAc_adj < 0.01 & stats_weighted$NeuAc_effectSize > 0], row.names(stats_weighted)[stats_weighted$NeuAc.a2.3.Gal.b1.4.Glc_adj < 0.01 & stats_weighted$NeuAc.a2.3.Gal.b1.4.Glc_effectSize > 0])
siaBindingFeats_neg = list(row.names(stats_weighted)[stats_weighted$Sialic_Acid_adj < 0.01 & stats_weighted$Sialic_Acid_effectSize < 0], row.names(stats_weighted)[stats_weighted$NeuAc_adj < 0.01 & stats_weighted$NeuAc_effectSize < 0], row.names(stats_weighted)[stats_weighted$NeuAc.a2.3.Gal.b1.4.Glc_adj < 0.01 & stats_weighted$NeuAc.a2.3.Gal.b1.4.Glc_effectSize < 0])
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
intersect(intersect(siaBindingFeats_pos[[1]], siaBindingFeats_pos[[2]]), siaBindingFeats_pos[[3]])
intersect(intersect(siaBindingFeats_neg[[1]], siaBindingFeats_neg[[2]]), siaBindingFeats_neg[[3]])


# Fucose binding
fucBindingFeats_pos = list(row.names(stats_weighted)[stats_weighted$Fuc_adj < 0.01 & stats_weighted$Fuc_effectSize > 0], row.names(stats_weighted)[stats_weighted$Terminal_Fucose_adj < 0.01 & stats_weighted$Terminal_Fucose_effectSize > 0])
fucBindingFeats_neg = list(row.names(stats_weighted)[stats_weighted$Fuc_adj < 0.01 & stats_weighted$Fuc_effectSize < 0], row.names(stats_weighted)[stats_weighted$Terminal_Fucose_adj < 0.01 & stats_weighted$Terminal_Fucose_effectSize < 0])

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
manBindingFeats_pos = list(row.names(stats_weighted)[stats_weighted$Man_adj < 0.1 & stats_weighted$Man_effectSize > 0], row.names(stats_weighted)[stats_weighted$High_Mannose_adj < 0.1 & stats_weighted$High_Mannose_effectSize > 0], row.names(stats_weighted)[stats_weighted$Man.a1.2.Man_adj < 0.1 & stats_weighted$Man.a1.2.Man_effectSize > 0])
manBindingFeats_neg = list(row.names(stats_weighted)[stats_weighted$Man_adj < 0.1 & stats_weighted$Man_effectSize < 0], row.names(stats_weighted)[stats_weighted$High_Mannose_adj < 0.1 & stats_weighted$High_Mannose_effectSize < 0], row.names(stats_weighted)[stats_weighted$Man.a1.2.Man_adj < 0.1 & stats_weighted$Man.a1.2.Man_effectSize < 0])

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










# d2FeatCor = cor(d2Feats)
# corrplot(d2FeatCor)
# plot(density(d2FeatCor[upper.tri(d2FeatCor,diag = F)]),
#      main = 'Density distribution of pairwise correlations b/w d2Feats',
#      xlab = 'Pearson correlation')
# abline(v = c(0.7, -0.7))
# 
# sum(abs(d2FeatCor[upper.tri(d2FeatCor,diag = F)]) > 0.7) / length(d2FeatCor[upper.tri(d2FeatCor,diag = F)])
# 
# subD2Feats = d2Feats[,grepl('_6Ang', colnames(d2Feats))]
# 
# corrplot(cor(subD2Feats))
# sum(abs(cor(subD2Feats)[upper.tri(cor(subD2Feats),diag = F)]) > 0.7) / length(cor(subD2Feats)[upper.tri(cor(subD2Feats),diag = F)])