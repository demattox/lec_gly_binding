rm (list = ls())

library(ggplot2)
library(corrplot)
library(reshape2)
#library(factoextra)
library(philentropy)
library(protr)
library(ggbiplot)
library(devtools)
library(RColorBrewer)


source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

setwd("/Users/dmattox/cbk/glycan_binding/analysis/prelim/prelim2/")

dat = read.delim(file = '/Users/dmattox/cbk/glycan_binding/analysis/prelim/prelim2/prelimData2.tsv', sep = '\t', header = T, stringsAsFactors = F)
length(unique(dat$pdb))
length(unique(dat$uniparc))

missingInteractions = rep(F, nrow(dat))

for (i in 1:nrow(dat)) {
  if ((dat$total[i] == 0)
    || (dat$numBSresis_bin1 == 0)) {
      missingInteractions[i] = T
    }
}

cnt = 0
for (i in 1:length(unique(dat$pdb))){
  pdbID = unique(dat$pdb)[i]
  tag = dat$pdb == pdbID
  if (length(unique(dat$longname[tag])) > 1){
    cnt = cnt + 1
  }
}
cnt

dat = dat[!missingInteractions,]
length(unique(dat$pdb))
length(unique(dat$uniparc))

for (i in 1:ncol(dat)){
  if (all(dat[,i] == 0) == TRUE){
    print(colnames(dat)[i])
    print(i)
  }
}
print(grep('I_bin4',colnames(dat))) # pi helix in bin 4 has a single binding site with any pi helix resiudes

bsDat = dat[,c(11, 13:31, 33:44, 46:57, 59:79)] # filter out non-binding site columns and columns 32, 45, & 58 for having 0 recorded pi helices in bins 1-3

pca = prcomp(x = bsDat, scale. = T)
plot(pca$x)

colfunc = colorRampPalette(c("red","goldenrod","forestgreen","royalblue","darkviolet"))


# PCA by lectin origin
xLim = c(min(pca$x[,1])+0.2, max(pca$x[,1])+0.2)
yLim = c(min(pca$x[,2])+0.2, max(pca$x[,2])+0.2)

colors = colfunc(length(unique(dat$origine)))
origins = sort(unique(dat$origine))
plot(0,0, col = 'white', xlab = 'PC1', ylab = 'PC2', xlim = xLim, ylim = yLim, main = 'PCA of all ligand sites')
par(new=T)
for (i in 1:length(origins)){
  orig = origins[i]
  plot(pca$x[dat$origine == orig, 1:2], pch = 19, col = alpha(colors[i],0.4), xlab = '', ylab = '', axes = F, xlim = xLim, ylim = yLim)
  par(new=T)
}
par(new=F)
legend(x = xLim[2]-4, y=yLim[2], col = colors, legend = origins, pch=15, bty = 'n')


# PCA by lectin fold

# colors = colfunc(length(unique(dat$new_fold.Not.Published)))
# plot(0,0, col = 'white', xlab = 'PC1', ylab = 'PC2', xlim = c(-6.2,6), ylim = c(-12,5.5), main = 'PCA of all ligand sites')
# par(new=T)
# for (i in 1:length(unique(dat$new_fold.Not.Published))){
#   fold = unique(dat$new_fold.Not.Published)[i]
#   plot(pca$x[dat$new_fold.Not.Published == fold, 1:2], col = colors[i], xlab = '', ylab = '', axes = F, xlim = c(-6.2,6), ylim = c(-12,5.5))
#   par(new=T)
# }
# par(new=F)

# PCA by monosaccharide binding
plot(pca$x, pch=16, col = alpha('black',alpha=0.4),xlim = xLim, ylim = yLim, main = "PCA of all binding sites")
par(new=T)
plot(pca$x[dat$iupac == 'Man',1:2], pch = 19, col = alpha('forestgreen',0.6), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F)
par(new=T)
plot(pca$x[dat$iupac == 'NeuAc',1:2], pch = 23, col = alpha('darkorchid',0.6), bg = 'darkorchid', xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F)
par(new=T)
plot(pca$x[dat$iupac == 'Gal',1:2], pch = 19, col = alpha('gold',0.6), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F)
par(new=T)
plot(pca$x[dat$iupac == 'Glc',1:2], pch = 19, col = alpha('dodgerblue',0.6), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F)
par(new=T)
plot(pca$x[dat$iupac == 'Fuc',1:2], pch = 17, col = alpha('red',0.6), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F)
legend(x = xLim[1]-0.6, y=-4.5, col = c('forestgreen', 'dodgerblue', 'gold', 'red', 'darkorchid'), pt.bg	 = c(NA, NA, NA, NA, 'darkorchid'), legend = c('Man','Glc','Gal', 'Fuc', 'NeuAc'), pch=c(19,19,19,17,23), bty = 'n')


tag = grepl('-', dat$iupac) # tag ligands with anotated with multiple sugar groups

monoDat = dat[!tag,] # All monosaccharide binding sites
monoBSDat = monoDat[,c(11, 13:31, 33:44, 46:57, 59:79)] # Filter out known all-zero columns

for (i in 1:ncol(monoBSDat)){
  if (all(monoBSDat[,i] == 0) == TRUE){
    print(colnames(monoBSDat)[i])
    print(i)
  }
}
monoBSDat = monoBSDat[,-57] # Drop I_bin4 that is all non-zero with single non-mono binding site


#######################
### FIGURE 4 in glycoT poster
######################

m.pca = prcomp(monoBSDat, scale. = T)

pc1var = round(m.pca$sdev[1]**2/sum(m.pca$sdev**2), 4) * 100
pc2var = round(m.pca$sdev[2]**2/sum(m.pca$sdev**2), 4) * 100


xLim = c(min(m.pca$x[,1])+0.2, max(m.pca$x[,1])+0.2)
yLim = c(min(m.pca$x[,2])+0.2, max(m.pca$x[,2])+0.2)
a = 0.5
pSize = 1.4


tag = monoDat$iupac %in% c('Man','Glc','Gal', 'Fuc', 'NeuAc', 'GalNAc', 'GlcNAc')
plot(m.pca$x[!tag, 1:2], pch=16, col = alpha('black',alpha=0.4),xlim = xLim, ylim = yLim, main = "PCA of all binding sites containing a monosaccharide",
     xlab = paste('PC1 (',as.character(pc1var),'% var)',sep = ""),
     ylab = paste('PC2 (',as.character(pc2var),'% var)',sep = ""), 
     cex = pSize)
par(new=T)
plot(m.pca$x[monoDat$iupac == 'Man',1:2], pch = 19, col = alpha('forestgreen',a), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F, cex = pSize)
par(new=T)
plot(m.pca$x[monoDat$iupac == 'NeuAc',1:2], pch = 23, col = alpha('darkorchid',a), bg = 'darkorchid', xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F, cex = pSize)
par(new=T)
plot(m.pca$x[monoDat$iupac == 'Gal',1:2], pch = 19, col = alpha('gold',a), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F, cex = pSize)
par(new=T)
plot(m.pca$x[monoDat$iupac == 'GalNAc',1:2], pch = 15, col = alpha('gold',a), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F, cex = pSize)
par(new=T)
plot(m.pca$x[monoDat$iupac == 'Glc',1:2], pch = 19, col = alpha('dodgerblue',a), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F, cex = pSize)
par(new=T)
plot(m.pca$x[monoDat$iupac == 'GlcNAc',1:2], pch = 15, col = alpha('dodgerblue',a), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F, cex = pSize)
par(new=T)
plot(m.pca$x[monoDat$iupac == 'Fuc',1:2], pch = 17, col = alpha('red',a), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F, cex = pSize)
legend(x = xLim[1], y=yLim[1]+5, col = c('forestgreen', 'dodgerblue', 'dodgerblue', 'gold', 'gold','red', 'darkorchid'), pt.bg	 = c(NA, NA, NA, NA, NA, NA, 'darkorchid'), legend = c('Man','Glc','GlcNAc','Gal','GalNAc', 'Fuc', 'NeuAc'), pch=c(19,19,15,19,15,17,23), bty = 'n', cex = 1.4)

pc1Cols = m.pca$rotation[,1]
pc1Cols = sort(abs(pc1Cols))

pc2Cols = m.pca$rotation[,2]
    pc2Cols = sort(abs(pc2Cols))



# xLim = c(min(m.pca$x[,1])+0.2, max(m.pca$x[,1])+0.2)
# yLim = c(min(m.pca$x[,2])+0.2, max(m.pca$x[,2])+0.2)
# 
# plot(m.pca$x, col = alpha('black',0.4),pch=16, xlim = xLim, ylim = yLim, main = "PCA of subset of binding sites with a monosaccharide")
# par(new=T)
# plot(m.pca$x[monoDat$iupac == 'Man',1:2], pch = 19, col = alpha('forestgreen',0.6), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F)
# par(new=T)
# plot(m.pca$x[monoDat$iupac == 'NeuAc',1:2], pch = 23, col = alpha('darkorchid',0.6), bg = 'darkorchid', xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F)
# par(new=T)
# plot(m.pca$x[monoDat$iupac == 'Gal',1:2], pch = 19, col = alpha('gold',0.6), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F)
# par(new=T)
# plot(m.pca$x[monoDat$iupac == 'Glc',1:2], pch = 19, col = alpha('dodgerblue',0.6), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F)
# par(new=T)
# plot(m.pca$x[monoDat$iupac == 'Fuc',1:2], pch = 17, col = alpha('red',0.6), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F)
# legend(x = xLim[2]-2.5, y=-2, col = c('forestgreen', 'dodgerblue', 'gold', 'red', 'darkorchid'), pt.bg	 = c(NA, NA, NA, NA, 'darkorchid'), legend = c('Man','Glc','Gal', 'Fuc', 'NeuAc'), pch=c(19,19,19,17,23), bty = 'n')


# Save out first 50 PCs for tSNE
out = pca$x[,1:50]
labs = dat[,1:12]
write.csv(file = 'prelim_pc50.csv', out, quote = F, row.names = F, col.names = F)
write.table(file = 'prelim_pc50_labels.tsv', labs, sep = '\t', quote = F, row.names = F)

# out = m.pca$x[,1:50]
# labs = monoDat[,1:12]
# write.csv(file = '/Users/dmattox/cbk/glycan_binding/analysis/prelimMonos_pc50.csv', out, quote = F, row.names = F, col.names = F)
# write.table(file = '/Users/dmattox/cbk/glycan_binding/analysis/prelimMonos_pc50_labels.tsv', labs, sep = '\t', quote = F, row.names = F)

################
# Try limiting to one binding site/pdb file
# Pick the one with the most interactions, if tie, one with most resis in bin1, tie --> bin2, ... bin3/bin4

row.names(dat) = 1:nrow(dat) # fix row numbering after deleting rows

uniDat = as.data.frame(matrix(nrow = length(unique(dat$pdb)), ncol = ncol(dat)))
colnames(uniDat) = colnames(dat)
for (i in 1:nrow(uniDat)){
  pdb = unique(dat$pdb)[i]
  rInds = row.names(dat)[dat$pdb == pdb] # Row indices for a given PDB ID
  hiScores = rep(0,5)
  keep = NA
  for (j in 1:length(rInds)){
    ind = as.integer(rInds[j])
    datScores = as.numeric(c(dat$total[ind], dat$numBSresis_bin1[ind], dat$numBSresis_bin2[ind], dat$numBSresis_bin3[ind], dat$numBSresis_bin4[ind]))
    for (k in 1:5){
      if (datScores[k] > hiScores[k]){
        keep = ind
        hiScores = datScores
        break
      } else if (datScores[k] < hiScores[k]){ #else (==), continue with loop
        break
      }
    }
  }
  if (is.na(keep)){
    print('No pick!')
    print(pdb)
  }else{
    uniDat[i,] = dat[keep,]
  }
}

for (i in 1:ncol(uniDat)){
  if (all(uniDat[,i] == 0)){
    print(i)
    print(colnames(uniDat)[i])
  }
}

uniBsDat = uniDat[,c(11, 13:31, 33:44, 46:57, 59:79)] # filter out non-binding site columns and columns 32, 45, & 58 for having 0 record pi helices

pca = prcomp(x = uniBsDat, scale. = T)
plot(pca$x)


##
xLim = c(min(pca$x[,1])+0.2, max(pca$x[,1])+0.2)
yLim = c(min(pca$x[,2])+0.2, max(pca$x[,2])+0.2)

colors = colfunc(length(unique(uniDat$origine)))
origins = sort(unique(uniDat$origine))
plot(0,0, col = 'white', xlab = 'PC1', ylab = 'PC2', xlim = xLim, ylim = yLim, main = 'PCA of 1 BS/PDB')
par(new=T)
for (i in 1:length(origins)){
  orig = origins[i]
  plot(pca$x[uniDat$origine == orig, 1:2], pch = 19, col = alpha(colors[i],0.4), xlab = '', ylab = '', axes = F, xlim = xLim, ylim = yLim)
  par(new=T)
}
par(new=F)
legend(x = xLim[1], y=-6, col = colors, legend = origins, pch=15, bty = 'n')



# PCA by monosaccharide binding
tag = uniDat$iupac %in% c('Man','Glc','Gal', 'Fuc', 'NeuAc', 'GalNAc', 'GlcNAc')
plot(pca$x[!tag, 1:2], pch=16, col = alpha('black',alpha=0.4),xlim = xLim, ylim = yLim, main = "PCA of 1 BS/PDB")
par(new=T)
plot(pca$x[uniDat$iupac == 'Man',1:2], pch = 19, col = alpha('forestgreen',0.6), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F)
par(new=T)
plot(pca$x[uniDat$iupac == 'NeuAc',1:2], pch = 23, col = alpha('darkorchid',0.6), bg = 'darkorchid', xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F)
par(new=T)
plot(pca$x[uniDat$iupac == 'Gal',1:2], pch = 19, col = alpha('gold',0.6), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F)
par(new=T)
plot(pca$x[uniDat$iupac == 'GalNAc',1:2], pch = 15, col = alpha('gold',0.6), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F)
par(new=T)
plot(pca$x[uniDat$iupac == 'Glc',1:2], pch = 19, col = alpha('dodgerblue',0.6), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F)
par(new=T)
plot(pca$x[uniDat$iupac == 'GlcNAc',1:2], pch = 15, col = alpha('dodgerblue',0.6), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F)
par(new=T)
plot(pca$x[uniDat$iupac == 'Fuc',1:2], pch = 17, col = alpha('red',0.6), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F)
legend(x = xLim[1], y=-6, col = c('forestgreen', 'dodgerblue', 'dodgerblue', 'gold', 'gold','red', 'darkorchid'), pt.bg	 = c(NA, NA, NA, NA, NA, NA, 'darkorchid'), legend = c('Man','Glc','GlcNAc','Gal','GalNAc', 'Fuc', 'NeuAc'), pch=c(19,19,15,19,15,17,23), bty = 'n')

# subset of monosaccs
tag = grepl('-', uniDat$iupac) # tag ligands with anotated with multiple sugar groups

uniMonoDat = uniDat[!tag,] # All monosaccharide binding sites
uniMonoBSDat = uniMonoDat[,c(11, 13:31, 33:44, 46:57, 59:79)]

for (i in 1:ncol(uniMonoBSDat)){
  if (all(uniMonoBSDat[,i] == 0) == TRUE){
    print(colnames(uniMonoBSDat)[i])
    print(i)
  }
}
uniMonoBSDat = uniMonoBSDat[,-57]


m.pca = prcomp(uniMonoBSDat, scale. = T)
pc1var = round(m.pca$sdev[1]**2/sum(m.pca$sdev**2), 4) * 100
pc2var = round(m.pca$sdev[2]**2/sum(m.pca$sdev**2), 4) * 100



xLim = c(min(m.pca$x[,1])+0.2, max(m.pca$x[,1])+0.2)
yLim = c(min(m.pca$x[,2])+0.2, max(m.pca$x[,2])+0.2)
a = 0.75
pSize = 1.4


tag = uniMonoDat$iupac %in% c('Man','Glc','Gal', 'Fuc', 'NeuAc', 'GalNAc', 'GlcNAc')
plot(m.pca$x[!tag, 1:2], pch=16, col = alpha('black',alpha=0.6),xlim = xLim, ylim = yLim, main = "PCA of 1 BS/PDB",
     xlab = paste('PC1 (',as.character(pc1var),'% var)',sep = ""),
     ylab = paste('PC2 (',as.character(pc2var),'% var)',sep = ""), 
     cex = pSize)
par(new=T)
plot(m.pca$x[uniMonoDat$iupac == 'Man',1:2], pch = 19, col = alpha('forestgreen',a), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F, cex = pSize)
par(new=T)
plot(m.pca$x[uniMonoDat$iupac == 'NeuAc',1:2], pch = 23, col = alpha('darkorchid',a), bg = 'darkorchid', xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F, cex = pSize)
par(new=T)
plot(m.pca$x[uniMonoDat$iupac == 'Gal',1:2], pch = 19, col = alpha('gold',a), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F, cex = pSize)
par(new=T)
plot(m.pca$x[uniMonoDat$iupac == 'GalNAc',1:2], pch = 15, col = alpha('gold',a), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F, cex = pSize)
par(new=T)
plot(m.pca$x[uniMonoDat$iupac == 'Glc',1:2], pch = 19, col = alpha('dodgerblue',a), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F, cex = pSize)
par(new=T)
plot(m.pca$x[uniMonoDat$iupac == 'GlcNAc',1:2], pch = 15, col = alpha('dodgerblue',a), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F, cex = pSize)
par(new=T)
plot(m.pca$x[uniMonoDat$iupac == 'Fuc',1:2], pch = 17, col = alpha('red',a), xlim = xLim, ylim = yLim, xlab = '', ylab = '', axes = F, cex = pSize)
legend(x = xLim[2]-2.5, y=yLim[2], col = c('forestgreen', 'dodgerblue', 'dodgerblue', 'gold', 'gold','red', 'darkorchid'), pt.bg	 = c(NA, NA, NA, NA, NA, NA, 'darkorchid'), legend = c('Man','Glc','GlcNAc','Gal','GalNAc', 'Fuc', 'NeuAc'), pch=c(19,19,15,19,15,17,23), bty = 'n', cex = 1.4)






out = pca$x[,1:50]
labs = uniDat[,1:12]
write.csv(file = '/Users/dmattox/cbk/glycan_binding/analysis/uniprelim_pc50.csv', out, quote = F, row.names = F)
write.table(file = '/Users/dmattox/cbk/glycan_binding/analysis/uniprelim_pc50_labels.tsv', labs, sep = '\t', quote = F, row.names = F)




















############################
# Corr plot - bs similarity vs seq similarity
############################s

# bsSim = distance(bsDat, method = "cosine", use.row.names = TRUE) # Cosine similarity matrix
# plot(density(bsSim)) # majority above 0.5, median 0.86

bsSim = cor(t(bsDat), method = 'spearman') # Use Spearman Correlation to cluster together profiles with similar shapes or show similar general trends (e.g. increasing expression with time), but whose expression levels may be very different.
plot(density(bsSim)) # fairly normal dist with a small bump around 1, median = 0.58

# bSimPears = cor(t(bsDat), method = 'pearson')
# plot(density(bSimPears)) # very large peak towards 1, uniform with gradual decline below 0.5 (VERY SIMILAR to cosine similarity), median = 0.84


uniBsSim = cor(t(uniBsDat), method = 'spearman') # Use Spearman Correlation to cluster together profiles with similar shapes or show similar general trends (e.g. increasing expression with time), but whose expression levels may be very different.
plot(density(uniBsSim))

# cat(unique(uniDat$uniparc)) # Print out uniprot IDs to get seqeunces from uniprot api
# for (i in 1:length(unique(uniDat$uniparc))){
#   cat(unique(uniDat$uniparc)[i], " ")
# }

seqs = readFASTA('/Users/dmattox/cbk/glycan_binding/analysis/prelim/uniprot-uniDat_seqs.fasta')

for (i in 1:length(seqs)){
  names(seqs)[i] = strsplit(names(seqs)[i], '\\|')[[1]][2]
}
seqs

tag1 = names(seqs) %in% uniDat$uniparc
tag2 = uniDat$uniparc %in% names(seqs)

plotSeq = seqs[tag1] # all uniprot IDs maded back (good!)
pseudoplotDat = uniDat[tag2,c(11, 13:31, 33:44, 46:57, 59:79)] # Limit to informative columns, lost 4 pdbs/lectins where provided uniprot ID did not map to uniprot API service
row.names(pseudoplotDat) = 1:nrow(pseudoplotDat)
datIDs = uniDat$uniparc[tag2]

folds = rep('',length(unique(datIDs)))
plotDat = as.data.frame(matrix(nrow = length(unique(datIDs)), ncol = ncol(pseudoplotDat)))
colnames(plotDat) = colnames(pseudoplotDat)

for (i in 1:length(unique(datIDs))){
  tag = datIDs == unique(datIDs)[i]
  folds[i] = unique(uniDat$new_fold.Not.Published[tag2][tag])
  if (sum(tag) == 1){
    row.names(plotDat)[i] = unique(datIDs)[i]
    plotDat[i,] = pseudoplotDat[tag,]
  }else{
    resCnts = pseudoplotDat$numBSresis_bin1[tag]
    tagInd = which.max(resCnts)
    row.names(plotDat)[i] = unique(datIDs)[i]
    plotDat[i,] = pseudoplotDat[tag,][tagInd,]
  }
}

all(names(plotSeq) == row.names(plotDat))

#psimmat = parSeqSim(plotSeq, cores = 4, type = "local")
#globPsimmat = parSeqSim(plotSeq, cores = 4, type = "global")

plot(density(psimmat))
plot(density(globPsimmat))

stGlobPsimmat = (globPsimmat - min(globPsimmat))/(max(globPsimmat) - min(globPsimmat)) # Feature scaling between 0 and 1

plot(density(stGlobPsimmat))

row.names(psimmat) = row.names(plotDat)
colnames(psimmat) = row.names(plotDat)

dists = cor(t(plotDat), method = 'spearman')
plot(density(dists), main = "Density distribution of binding site similarities (spearman cor)")

bsSim = distance(plotDat, method = "cosine", use.row.names = TRUE) # Cosine similarity matrix
plot(density(bsSim), main = "Density distribution of binding site similarities (cosine sim)") # majority above 0.5, median 0.86

bsSim = distance(plotDat, method = "euclidean", use.row.names = TRUE) # Euclidian similarity matrix
plot(density(bsSim), main = "Raw (euc dist)") # looks like cosine/pearson flipped
#  feature scale
bsSim = (bsSim - min(bsSim))/(max(bsSim) - min(bsSim))
bsSim = 1 - bsSim
dists = bsSim

normedPlotDat = as.data.frame(matrix(0, nrow = nrow(plotDat), ncol = ncol(plotDat)))
row.names(normedPlotDat) = row.names(plotDat); colnames(normedPlotDat) = colnames(plotDat)
for (i in 1:ncol(plotDat)){
  normedPlotDat[,i] = plotDat[,i]/max(plotDat[,i])
}

bsSim = distance(normedPlotDat, method = 'euclidean', use.row.names = T)
plot(density(bsSim), main = "Scaled (euc dist)")
bsSim = bsSim = (bsSim - min(bsSim))/(max(bsSim) - min(bsSim))
bsSim = 1 - bsSim
plot(density(bsSim[bsSim != 1]), main = "Scaled (euc sim)")
dists = bsSim


foldCols = list(
  "b-sandwich / virus globular domain" = "#DDACC9",
  "b-sandwich / cytolysin-like" = "#FF88AD",
  "a-helix triplets / Duffy-like" = "#FFB8CE",
  "b-sandwich / viral protein domain" = "#DD6091",
  "b-trefoil" = "#FF7290",
  "b-prism II" = "#FFA388",
  "small protein / Knottin" = "#C77963",
  "b-hairpin stack" = "#9440F3",
  "a/b mixed / LysM domain" = "#9900B3",
  "peptide" = "#C266D1",
  "a/b mixed with b-sheet / Fibrinogen C-ter like" = "#6C00BF",
  "a/b hairpin / non-globular proline-rich" = "#A700FF",
  "a/b barrel / TIM" = "#CA66FF",
  "b-barrel" = "#7779BF",
  "b-sandwich / PA14 adhesin" = "#8194CC",
  "a/b mixed / C-type lectin-like" = "#533691",
  "b-prism I" = "#9189FF",
  "b-helix" = "#B09FFF",
  "a/b mixed with b-sheet / MAR domain" = "#756FB3",
  "b-prism III" = "#9FAAFF",
  "b-sandwich / TNF-like" = "#FF00FF",
  "b-sandwich / CUB-like" = "#AF00E6",
  "b-sandwich / Ig-like" = "#FF00B3",
  "a/b OB-fold" = "#B3128A",
  "b-propeller" = "#FF4DC1",
  "b-sandwich / 2 calcium lectin" = "#BD3D9A",
  "b-sandwich / ConA-like" = "#882E81",
  "b-sandwich / cyanovirin-like" = "#AD589A",
  "b-sandwich / pili and adhesins" = "#AC3491",
  "b-sandwich / viral coat and capsid protein" = "#FFFF00",
  "a/b mixed with b-sheet / not classified" = "#FFBB33",
  "small protein / APPLE domain" = "#804811",
  "b-sandwich / Galactose-binding domain-like" = "#B06411"
)


seqPoints = rep(0,(nrow(dists)**2 - nrow(dists))/2)
bindPoints = rep(0,(nrow(dists)**2 - nrow(dists))/2)
foldX = foldY = rep(0, length(seqPoints))
fold = rep('', length(seqPoints))

prev = 0
for (i in 1:nrow(dists)){ # Flatten out dist and seqsim mats together
  for(j in i:nrow(dists)){
    if ( i != j){
      prev = prev + 1
      seqPoints[prev] = psimmat[i,j] 
      bindPoints[prev] = dists[i,j]
      if (folds[i] == folds[j]){
        foldX[prev] = psimmat[i,j]
        foldY[prev] = dists[i,j]
        fold[prev] = folds[i]
      }
    }
  }
}

tag = foldY != 0
foldX = foldX[tag]
foldY = foldY[tag]
fold = fold[tag]

foldcol = rep('',length(fold))
for (i in 1:length(fold)){
  foldcol[i] = alpha(foldCols[[fold[i]]], 0.6)
}



textSize = 0.65

plot(seqPoints[!tag], bindPoints[!tag],pch=16, col = alpha('black',0.75), xlab = "Sequence Similarity", ylab = "Binding Site Similarity (Euclidean)", main = "Pairwise comparisons of sequence similarity and binding site characterizations",xlim = c(0,1), ylim = c(0,1), cex = 1.2)
lines(c(0,1),c(0,1), lty=2, lwd = 3)
# cor(seqPoints, bindPoints)
# text(x=0.7, y =0.1, "Cor = 0.24", cex = 1.5)
par(new=T)
plot(foldX, foldY, col = foldcol, axes = F, xlab="",ylab="",main="",xlim = c(0,1), ylim = c(0,1), pch = 16, cex =1.2)
legend(x = 0,y=1.02,col = unlist(foldCols[1:10][names(foldCols[1:10])],use.names=F), legend = names(foldCols[1:10]),pch = 16,bty = 'n', cex = textSize,  pt.cex = 1.5)
legend(x = 0.23,y=1.02,col = unlist(foldCols[11:20][names(foldCols[11:20])],use.names=F), legend = names(foldCols[11:20]),pch = 16,bty = 'n', cex = textSize, pt.cex = 1.5)
legend(x = 0.35,y=0.2,col = unlist(foldCols[20:28][names(foldCols[20:28])],use.names=F), legend = names(foldCols[20:28]),pch = 16,bty = 'n', cex = textSize, pt.cex = 1.5)
legend(x = 0.55,y=0.1,col = unlist(foldCols[29:length(foldCols)][names(foldCols[29:length(foldCols)])],use.names=F), legend = names(foldCols[29:length(foldCols)]),pch = 16,bty = 'n', cex = textSize, pt.cex = 1.5)


tsne = read.csv("cbk/glycan_binding/analysis/tSNE_prelim_p30.csv", header = F, stringsAsFactors = F)
which.min(tsne[,2])








##########################################################################
# Heatmap of raw data
bsDatScaled = as.data.frame(matrix(0, nrow = nrow(bsDat), ncol = ncol(bsDat)))
row.names(bsDatScaled) = row.names(bsDat); colnames(bsDatScaled) = colnames(bsDat)
for (i in 1:ncol(bsDat)){
  bsDatScaled[,i] = bsDat[,i]/max(bsDat[,i])
}
colnames(bsDatScaled) = gsub("^X._", "X_", colnames(bsDatScaled))


zScores_bsDat = as.data.frame(matrix(0, nrow = nrow(bsDat), ncol = ncol(bsDat)))
row.names(zScores_bsDat) = row.names(bsDat); colnames(zScores_bsDat) = colnames(bsDat)
for (i in 1:ncol(bsDat)){
  m = mean(bsDat[,i])
  std = sd(bsDat[,i])
  zScores_bsDat[,i] = (bsDat[,i] - m) / std
}
colnames(zScores_bsDat) = gsub("^X._", "X_", colnames(zScores_bsDat))


bcols = brewer.pal(5,"BuPu")
binColors = matrix("grey95", nrow = ncol(bsDatScaled), ncol = 1)
binColors[grepl("bin1$", colnames(bsDatScaled))] = bcols[2]
binColors[grepl("bin2$", colnames(bsDatScaled))] = bcols[3]
binColors[grepl("bin3$", colnames(bsDatScaled))] = bcols[4]
binColors[grepl("bin4$", colnames(bsDatScaled))] = bcols[5]

fcols = brewer.pal(4, "Set3")
featColors = matrix("white", nrow = ncol(bsDatScaled), ncol = 1)
featColors[2:5] = fcols[1] # Residue counts
featColors[6:16] = fcols[2] # Interaction types
featColors[grepl("^[A-Z]_bin", colnames(bsDatScaled))] = fcols[3] # Secondary structure
tagRemaining = featColors == "white"
tagRemaining[1] = F # Number of binding sites
featColors[tagRemaining] = fcols[4] # Residue types

clab = cbind(binColors, featColors)
colnames(clab) = c("Bin identity", "Type of feature")


colfunc = colorRampPalette(c("red","goldenrod","forestgreen","royalblue","darkviolet"))
scols = colfunc(6)
origins = sort(unique(dat$origine))
speciesColors = matrix("", ncol = nrow(bsDatScaled), nrow = 1)
for (i in 1:length(origins)){
  speciesColors[dat$origine == origins[i]] = scols[i]
}

mcols = c('forestgreen', 'dodgerblue', 'dodgerblue3', 'darkgoldenrod1', 'gold3', 'firebrick1', 'darkorchid')
monos = c('Man','Glc','GlcNAc','Gal','GalNAc', 'Fuc', 'NeuAc')
monosacColors = matrix("grey90", ncol = nrow(bsDatScaled), nrow = 1)
monosacColors[!grepl('-', dat$iupac)] = "grey" # Other monosaccharides
for (i in 1:length(monos)){
  monosacColors[dat$iupac == monos[i]] = mcols[i]
}

foldColors = matrix("", ncol = nrow(bsDatScaled), nrow = 1)
for (i in 1:length(foldCols)){
  foldColors[dat$new_fold.Not.Published == names(foldCols)[i]] = foldCols[[i]]
}

caColors =  matrix("darkolivegreen1", ncol = nrow(bsDatScaled), nrow = 1)
caColors[dat$metal > 0] = "chartreuse4"

rlab = rbind(speciesColors, foldColors, caColors, monosacColors)
row.names(rlab) = c("Lectin Origin", "Fold Type", "Ca2+ Coord.", "Monosacc. Id.")

mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}
mycol <- colorRampPalette(c("ivory", "cornflowerblue", "navy"))(n = 299)


pdf(file = '/Users/dmattox/Documents/qbs/cbk/glycan_binding/analysis/prelim/prelim2/prelim_heatmap.pdf')

heatmap.3(bsDatScaled, hclustfun = myclust, distfun = mydist, ColSideColors = clab,
          main = "All Lectin Binding Sites", col = mycol,labCol = F, labRow = F
          , dendrogram = "both",  ColSideColorsSize = 2, KeyValueName="Scaled Feature Values"
          , scale = 'none', density.info = 'none',margins=c(6,12), Rowv = T, Colv = T,
          RowSideColors=rlab, RowSideColorsSize = 3)
legend("topright",legend=c("Global", "Bin 1","Bin 2","Bin 3","Bin 4",
                           "",
                           "Resi Cnts", "Int Types", "Sec Struct", "Resi Types",
                           "", "",
                           "Animal", "Bacterial", "Fungal/yeast", "Plant", "Protist/parasite", "Viral",
                           "",
                           "No Ca2+", "Ca2+ coord.",
                           "",
                           "No monosacc", "Other monosacc", monos[1:7])
       , fill=c("grey85", bcols[2:5],
                "white",
                fcols[1:4],
                "white", "white",
                scols[1:6],
                "white",
                "darkolivegreen1", "chartreuse4",
                "white",
                "grey85", "grey", mcols[1:7])
       , border=F, bty="n",y.intersp = 0.7, cex=0.7)

dev.off()




lowC = -5 # Cutoff values for z scores (otherwise spans -60 to 45)
hiC = 5
tag1 = zScores_bsDat < lowC
tag2 = zScores_bsDat > hiC
print(sum(tag2))
zScores_bsDat[tag1] = lowC
zScores_bsDat[tag2] = hiC




mycol <- colorRampPalette(c("blue", "cornflowerblue", "grey5", "lightgoldenrod2",  "yellow"))(n = 299)

pdf(file = '/Users/dmattox/Documents/qbs/cbk/glycan_binding/analysis/prelim/prelim2/prelim_heatmap_zscores.pdf')

heatmap.3(zScores_bsDat, hclustfun = myclust, distfun = mydist, ColSideColors = clab,
          main = "All Lectin Binding Sites", col = mycol,labCol = F, labRow = F
          , dendrogram = "both",  ColSideColorsSize = 2, KeyValueName="Normalized Feature Values"
          , scale = 'none', density.info = 'none',margins=c(6,12), Rowv = T, Colv = T,
          RowSideColors=rlab, RowSideColorsSize = 3)
legend("topright",legend=c("Global", "Bin 1","Bin 2","Bin 3","Bin 4",
                           "",
                           "Resi Cnts", "Int Types", "Sec Struct", "Resi Types",
                           "", "",
                           "Animal", "Bacterial", "Fungal/yeast", "Plant", "Protist/parasite", "Viral",
                           "",
                           "No Ca2+", "Ca2+ coord.",
                           "",
                           "No monosacc", "Other monosacc", monos[1:7])
       , fill=c("grey85", bcols[2:5],
                "white",
                fcols[1:4],
                "white", "white",
                scols[1:6],
                "white",
                "darkolivegreen1", "chartreuse4",
                "white",
                "grey85", "grey", mcols[1:7])
       , border=F, bty="n",y.intersp = 0.7, cex=0.7)

dev.off()

###########################################
# 80% similarity CD-HIT clustering
cdhit = read.csv("/Users/dmattox/cbk/glycan_binding/analysis/prelim/1592235610.result/1592235610.fas.1.clstr.sorted", sep = "\t", row.names = NULL, header = FALSE, stringsAsFactors = FALSE)
prev =""
for (i in 1:nrow(cdhit)){ # Clean up CD-HIT cluster input
  if (grepl("^>Cluster",cdhit$V1[i])){
    print(cdhit$V1[i])
    prev = gsub(">","", cdhit$V1[i])
  }else{
    cdhit$V1[i] = prev
  }
}
cdhit = cdhit[!cdhit$V2 == "",]

clustCnts = rep(0,length(unique(cdhit$V1)))
for (i in 1:length(unique(cdhit$V1))){
  clustCnts[i] = sum(cdhit$V1 == unique(cdhit$V1)[i])
}
hist(clustCnts, breaks = 16, xlim = c(0,15), col = "dodgerblue", xlab = "Number of sequences in a cluster", main = "Cluster sizes")
plot(density(clustCnts), col = "dodgerblue", lwd = 2, xlab = "Number of sequences in a cluster", main = "Cluster sizes")

cdhit$V2 = gsub("^.*, >..\\|", "", cdhit$V2)
cdhit$V2 = gsub("\\|.*$", "", cdhit$V2)

clustFams = rep("", length(unique(cdhit$V1)))
for (i in 1:length(unique(cdhit$V1))) {
  clustFams[i] = unique(dat$new_fold.Not.Published[dat$uniparc %in% cdhit$V2[cdhit$V1 == unique(cdhit$V1)[i]]])
}
tmp = cbind(unique(cdhit$V1), clustCnts, clustFams)


all(cdhit$V2 %in% dat$uniparc)

dat$cdhitClust = ""
for (i in 1:length(unique(cdhit$V1))){
  clus = unique(cdhit$V1)[i]
  uniIDs = cdhit$V2[cdhit$V1 == unique(cdhit$V1)[i]]
  dat$cdhitClust[dat$uniparc %in% uniIDs] = clus
}

all(unique(cdhit$V1) %in% dat$cdhitClust) # All clusters represented in dat

bsClustCnts = rep(0,length(unique(cdhit$V1)))
for (i in 1:length(unique(cdhit$V1))){
  bsClustCnts[i] = sum(dat$cdhitClust == unique(cdhit$V1)[i])
}
hist(bsClustCnts, breaks = 20, xlim = c(0,200), col = "dodgerblue", xlab = "Number of binding sites in a cluster", main = "Cluster sizes")
plot(density(bsClustCnts), col = "dodgerblue", lwd = 2, xlab = "Number of binding sites  in a cluster", main = "Cluster sizes")



tmp = strsplit(dat$iupac, split = "\\([a|b][[:digit:]]-[[:digit:]]\\)")
siaTag = rep(F, length(tmp))
manTag = rep(F, length(tmp))
for (i in 1:length(tmp)){
  if (grepl("NeuAc",tmp[[i]][1])){
    siaTag[i] = T
  }
  manContent = mean(grepl("Man", tmp[[i]]))
  if (manContent >= 0.66){
    manTag[i] = T
  }
}


manDat = as.data.frame(matrix(nrow = length(unique(cdhit$V1)), ncol = ncol(bsDat)))
row.names(manDat) = unique(cdhit$V1)
colnames(manDat) = colnames(bsDat)
nonManDat = siaDat = nonSiaDat = manDat

for (i in 1:nrow(manDat)){
  clus = row.names(manDat)[i]
  clusTag = dat$cdhitClust == clus
  if (any((manTag & clusTag))){
    manDat[clus,] = apply(X = dat[(manTag & clusTag), colnames(manDat)], FUN = mean, MARGIN = 2)
  }
  if (any((!manTag & clusTag))){
    nonManDat[clus,] = apply(X = dat[(!manTag & clusTag), colnames(manDat)], FUN = mean, MARGIN = 2)
  }
  if (any((siaTag & clusTag))){
    siaDat[clus,] = apply(X = dat[(siaTag & clusTag), colnames(manDat)], FUN = mean, MARGIN = 2)
  }
  if (any((!siaTag & clusTag))){
    nonSiaDat[clus,] = apply(X = dat[(!siaTag & clusTag), colnames(manDat)], FUN = mean, MARGIN = 2)
  }
}

manDat = manDat[!is.na(manDat$numBsites),]
nonManDat = nonManDat[!is.na(nonManDat$numBsites),]
siaDat = siaDat[!is.na(siaDat$numBsites),]
nonSiaDat = nonSiaDat[!is.na(nonSiaDat$numBsites),]

stats = as.data.frame(matrix(nrow = ncol(manDat), ncol = 4))
row.names(stats) = colnames(manDat)
row.names(stats)  = gsub("^X._", "X_", row.names(stats))
colnames(stats) = c("ManRichp","NeuTermp","ManRichDiff", "NeuTermDiff")
for (i in 1:ncol(manDat)){
  mTst = t.test(manDat[,i], nonManDat[,i])
  nTst = t.test(siaDat[,i], nonSiaDat[,i])
  stats[i,"ManRichp"] = mTst$p.value
  stats[i,"ManRichDiff"] = (mTst$estimate[1]-mTst$estimate[2]) / mTst$estimate[2]
  stats[i,"NeuTermp"] = nTst$p.value
  stats[i,"NeuTermDiff"] = (nTst$estimate[1]-nTst$estimate[2]) / nTst$estimate[2]
}
stats$ManRichAdj = p.adjust(stats$ManRichp, method = "BH")
stats$NeuTermAdj = p.adjust(stats$NeuTermp, method = "BH")

fcols = brewer.pal(4, "Set3")
stats$plotCols = "white"
stats$plotCols[1:5] = fcols[1] # Residue counts
stats$plotCols[6:16] = fcols[2] # Interaction types
stats$plotCols[grepl("^[A-Z]_bin", row.names(stats))] = fcols[3] # Secondary structure
tagRemaining = stats$plotCols == "white"
stats$plotCols[tagRemaining] = fcols[4] # Residue types

manFeats = row.names(stats)[stats$ManRichAdj < 0.05]
neuFeats = row.names(stats)[stats$NeuTermAdj < 0.05]

tmpCols = rep("",nrow(stats))
for (i in 1:length(stats$plotCols)){
  tmpCols[i] = alpha(stats$plotCols[i], 0.75)
  stats$plotCols[i] = alpha(stats$plotCols[i], 1)
}

tag = stats$ManRichAdj < 0.05
plot(stats$ManRichDiff, -log10(stats$ManRichp), xlab = "Fold difference from background", ylab = "-log10(p)", main = "Differences found in features from mannose rich ligands",
     col = tmpCols, cex = 2, pch = 16,
     cex.axis = 1.5, cex.main = 2, cex.lab = 1.5,
     xlim = c(-1.2,1.2), ylim = c(0,12))
par(new=T)
plot(stats$ManRichDiff[tag], -log10(stats$ManRichp[tag]), pch = 19, col = stats$plotCols[tag], cex = 2,
     axes = F, xlab = "", ylab = "", 
     xlim = c(-1.2,1.2), ylim = c(0,12))
par(new=T)
plot(stats$ManRichDiff[tag], -log10(stats$ManRichp[tag]), col = 'black', cex = 2.05,
     xlab = "", ylab = "", axes = F,
     xlim = c(-1.2,1.2), ylim = c(0,12))
abline(v= 0, lty =2)
abline(h = -log10(0.02))

# max(stats$ManRichp[tag])
# min(stats$ManRichp[!tag])
# max(stats$ManRichAdj[tag])
# min(stats$ManRichAdj[!tag])

tag = stats$NeuTermAdj < 0.05
xLim = c(-1.2,1.2)
yLim = c(0,5)

bg = "seashell2"
fg = "ivory"
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = bg)
abline(v = c(-1,-.5,0,.5,1), lwd = 6, col = fg)
abline(v = c(-1.25,-.75,-.25,.25,.75,1.25), lwd = 3, col = fg)
abline(h = c(0,1,2,3,4,5,6), lwd = 6, col = fg)
abline(h = c(0.5,1.5,2.5,3.5,4.5,5.5,6.5), lwd = 3, col = fg)
par(new=T)
plot(stats$NeuTermDiff, -log10(stats$NeuTermp), xlab = "Fold difference from background", ylab = "-log10(p)", main = "Differences enriched in NeuAc binding lectins",
     pch = 19, cex = 2, col = tmpCols,
     cex.axis = 1.5, cex.main = 2, cex.lab = 1.5,
     xlim = xLim, ylim = yLim)
par(new=T)
plot(stats$NeuTermDiff[tag], -log10(stats$NeuTermp[tag]), pch = 19, col = stats$plotCols[tag], cex = 2,
     axes = F, xlab = "", ylab = "", main = "",
     xlim = xLim, ylim = yLim)
par(new=T)
plot(stats$NeuTermDiff[tag], -log10(stats$NeuTermp[tag]), col = 'black', cex = 2.05,
     xlab = "", ylab = "", axes = F, main = "",
     xlim = xLim, ylim = yLim)
abline(h= -log10(0.008), lwd = 2)
abline(v = 0, lty=2, lwd = 2)

# max(stats$NeuTermp[tag])
# min(stats$NeuTermp[!tag])


#############
# Boxplots?
neuFeats

siadatFeats = siaDat[, neuFeats]
siadatFeats$Ligand = 'TerminalNeuAc'

nonsiadatFeats = nonSiaDat[, neuFeats]
nonsiadatFeats$Ligand = 'Background'

neuAll = rbind(siadatFeats, nonsiadatFeats)
for (i in 1:length(neuFeats)){
  neuAll[,neuFeats[i]] = neuAll[,neuFeats[i]]/max(neuAll[,neuFeats[i]])
}


neuMelt = melt(neuAll, id.vars = "Ligand")
colnames(neuMelt)[2:3] = c("Feature", "Scaled_Value")
neuMelt$Ligand = factor(neuMelt$Ligand)
neuMelt$Feature = factor(neuMelt$Feature)
neuMelt$Ligand = factor(neuMelt$Ligand, levels = rev(levels(neuMelt$Ligand)), ordered = T)
neuMelt$Feature = factor(neuMelt$Feature, levels = levels(neuMelt$Feature)[c(1,2,3,4,6,8, 7, 5, 9)], ordered = T)

fcols = brewer.pal(4, "Set3")
boxCols = rep("", length(neuFeats))
boxCols[1] = fcols[1] # Residue counts
boxCols[2:3] = fcols[2] # Interaction types
boxCols[4] = fcols[3] # Secondary structure
tagRemaining = boxCols == ""
boxCols[tagRemaining] = fcols[4] # Residue types

p1 = ggplot(data = neuMelt, aes(Feature, y = Scaled_Value, color  = Ligand, fill = Feature)) + labs(title="Significant features for terminal NeuAc (FDR < 0.05)", x="Binding site features", y="Scaled Value")
p1 + geom_boxplot(cex = 1.3) + scale_color_manual(values=c("darkorchid", "grey40")) + scale_fill_manual(values=boxCols) +
  theme(text = element_text(size = 18),
      axis.text.x = element_text(color = "grey20", size = 20, angle = 45, hjust = 1, vjust = 1, face = "plain"),
      axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
      axis.title.x = element_text(color = "grey20", size = 22, angle = 0, hjust = .5, vjust = 0, face = "bold"),
      axis.title.y = element_text(color = "grey20", size = 22, angle = 90, hjust = .5, vjust = .5, face = "bold"),
      legend.position='none')


#################################################
# Link to microarrays!!!

zBindingFilt = read.delim2("~/cbk/allGlycans/pipeline/zFilt_bound_glycans.txt", stringsAsFactors = F)

newInfo = read.delim2("~/cbk/allGlycans/prelim/lectin_info.txt", stringsAsFactors = F)
newInfo = newInfo[newInfo$uniprotID %in% dat$uniparc,] # Filter to lectins in our lectin set
newInfo = newInfo[newInfo$Lectin %in% row.names(zBindingFilt),]

all(newInfo$uniprotID %in% row.names(plotDat))
newInfo = newInfo[!newInfo$Lectin == "GSL.1",]
row.names(newInfo) = newInfo$uniprotID

plotBind = zBindingFilt[newInfo$Lectin,] # order by newInfo
row.names(plotBind) = newInfo$uniprotID

bpDists = distance(plotBind, method = 'euclidean', use.row.names = T)
plot(density(bpDists))
bpDists = (bpDists)/(max(bpDists))
bpDists = 1-bpDists
plot(density(bpDists))



plotDat = plotDat[row.names(plotDat) %in% newInfo$uniprotID,]
plotDat = plotDat[row.names(newInfo),] # order by order in newInfo
tag = apply(plotDat, sum, MARGIN = 2)
tag = tag != 0
plotDat = plotDat[,tag]
normedPlotDat = as.data.frame(matrix(0, nrow = nrow(plotDat), ncol = ncol(plotDat)))
row.names(normedPlotDat) = row.names(plotDat); colnames(normedPlotDat) = colnames(plotDat)
for (i in 1:ncol(plotDat)){
  normedPlotDat[,i] = plotDat[,i]/max(plotDat[,i])
}

bsSim = distance(normedPlotDat, method = 'euclidean', use.row.names = T)
plot(density(bsSim), main = "Scaled (euc dist)")
bsSim =(bsSim)/(max(bsSim))
bsSim = 1 - bsSim
plot(density(bsSim[bsSim != 1]), main = "Scaled (euc sim)")
dists = bsSim

all(row.names(dists) == row.names(bpDists))



bsPoints = rep(0,(nrow(dists)**2 - nrow(dists))/2)
bpPoints = rep(0,(nrow(dists)**2 - nrow(dists))/2)

prev = 0
for (i in 1:nrow(dists)){ # Flatten out dist and seqsim mats together
  for(j in i:nrow(dists)){
    if ( i != j){
      prev = prev + 1
      bsPoints[prev] = dists[i,j] 
      bpPoints[prev] = bpDists[i,j]
    }
  }
}

# tag = foldY != 0
# foldX = foldX[tag]
# foldY = foldY[tag]
# fold = fold[tag]



bg = "seashell2"
fg = "ivory"
par(mar=c(5,5,5,5)+.1)
plot(0,0,col = "white",axes = F, main = "", xlab = '', ylab = '', xlim = c(0,1), ylim = c(0,1))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = bg)
abline(v = seq(0,1,by = 0.2), lwd = 6, col = fg)
abline(v = seq(0,1,by = 0.1), lwd = 3, col = fg)
abline(h = seq(0,1,by = 0.2), lwd = 6, col = fg)
abline(h = seq(0,1,by = 0.1), lwd = 3, col = fg)
par(new=T)
plot(bsPoints, bpPoints,pch=16, col = alpha('dodgerblue',0.8), xlab = "Binding Site Similarity (Euclidean)", ylab = "Microarray Binding Profile Similarity", main = "Lectin specificity vs. binding site",
     xlim = c(0,1), ylim = c(0,1), cex = 1.5,
     cex.axis = 1.5, cex.main = 2, cex.lab = 1.8
     )
lines(c(0,1),c(0,1), lty=2, lwd = 3)



cor(bsPoints, bpPoints)






