
# rm (list = ls())

library(ggplot2)
library(corrplot)
library(reshape2)
#library(factoextra)
library(philentropy)
library(protr)
library(ggbiplot)
library(devtools)
library(RColorBrewer)
library(umap)
library(seqinr)
library(stringr)
library(randomForest)
library(caret)


source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

##########################

homeDir = '/Users/dmattox/cbk/lec_gly_binding/'

setwd(homeDir)

##########################
# functions
###########################
invarNorm <- function(invarVec){
  if (sum(invarVec) > 0){
    return(invarVec/sum(invarVec))
  }
  else{
    return(invarVec)
  }
}

pdbInRowNames <- function(df, pdbID){
  # Returns the rows of dataframe df that contain the provided PDB ID in their row names
  return(df[grepl(pdbID, row.names(df)),])
}

readCDHIT <- function(cdhitDir){
  FH = dir(cdhitDir)[grepl('*.clstr$', dir(cdhitDir))]
  cdhit = read.csv(paste(cdhitDir, FH, sep = ''), sep = "\t", row.names = NULL, header = FALSE, stringsAsFactors = FALSE)
  prev =""
  for (i in 1:nrow(cdhit)){ # Clean up CD-HIT cluster input
    if (grepl("^>Cluster",cdhit$V1[i])){
      # print(cdhit$V1[i])
      prev = gsub(">","", cdhit$V1[i])
    }else{
      cdhit$V1[i] = prev
    }
  }
  cdhit = cdhit[!cdhit$V2 == "",]
  cdhit$V1 = as.numeric(gsub('Cluster ','',cdhit$V1))
  cdhit$V2 = gsub("^.*, >", "", cdhit$V2)
  cdhit$V2 = gsub("\\.\\.\\..*$", "", cdhit$V2)
  colnames(cdhit) = c('clusterID','pdbID')
  return(cdhit)
}

countClustMembers <- function(cdhit){
  clustCnts = rep(0,length(unique(cdhit$V1)))
  for (i in 1:length(unique(cdhit[,1]))){
    clustCnts[i] = sum(cdhit[,1] == unique(cdhit[,1])[i])
  }
  return(clustCnts)
}

colfunc = colorRampPalette(c("red","goldenrod","forestgreen","royalblue","darkviolet"))
threshColors = c('firebrick3', 'darkorange2', 'darkgoldenrod2', 'gold2')


###########################
# Load dataframes
###########################
# Load McCaldon aa frequencies
mccaldon = read.delim('~/episweep/data/mccaldon.csv', header = F, sep = ",", stringsAsFactors = F)
colnames(mccaldon) = c("aa", 'mcFreq')

# Load residue-based features
bsResiDat = read.delim('./analysis/bsiteResiFeatures.tsv', header = T, sep ='\t', stringsAsFactors = F, row.names = 1)

# Load general voxelized pocket features & D2 measurement distribution features
d2Feats = read.delim('./analysis/d2_distributionFeats.csv', header = T, sep =',', stringsAsFactors = F, row.names = 1)

# load binned D2 measurements
if (dir.exists('./analysis/allD2binnedDists/')){
  d2Dists = read.delim(file = paste('./analysis/allD2binnedDists/', dir('./analysis/allD2binnedDists/')[1], sep = ''), sep = ',', header = T, stringsAsFactors = F)
  colnames(d2Dists) = gsub('^X','',colnames(d2Dists))
  for (i in 2:length(dir('./analysis/allD2binnedDists/'))){
    tmp = read.delim(file = paste('./analysis/allD2binnedDists/', dir('./analysis/allD2binnedDists/')[i], sep = ''), sep = ',', header = T, stringsAsFactors = F)
    colnames(tmp) = gsub('^X','',colnames(tmp))
    if (ncol(tmp) <= ncol(d2Dists)){
      d2Dists = dplyr::bind_rows(d2Dists, tmp)
    }else{
      d2Dists = dplyr::bind_rows(tmp, d2Dists)
    }
    d2Dists[is.na(d2Dists)] <- 0 # Fill columns that did not exist in other files as 0
  }
  rm(tmp)
  row.names(d2Dists) = d2Dists$bsite
  d2Dists$bsite <- NULL
  write.table(d2Dists, file = './analysis/d2_binnedMeasures.csv', quote = F, sep = ',') # save concatonated distances to a csv and remove the copied batch results (originals still in data dir)
  unlink('./analysis/allD2binnedDists/', recursive = T)
}else{
  d2Dists = read.delim('./analysis/d2_binnedMeasures.csv', header = T, sep =',', stringsAsFactors = F)
  colnames(d2Dists) = gsub('^X','',colnames(d2Dists))
}

all(row.names(bsResiDat) %in% row.names(d2Feats)) & nrow(d2Feats) == nrow(bsResiDat) # All three dataframes have the same rows names
all(row.names(bsResiDat) %in% row.names(d2Dists)) & nrow(d2Dists) == nrow(bsResiDat)

# Get dataframes in the same order
d2Feats = d2Feats[row.names(bsResiDat),]
d2Dists = d2Dists[row.names(bsResiDat),]

all(row.names(bsResiDat) == row.names(d2Feats)) & all(row.names(bsResiDat) == row.names(d2Dists)) # double check

# Read in 3DZD features and clean
zern3d = read.delim('./analysis/3DZDs_20thOrder.tsv', header = F, sep ='\t', stringsAsFactors = F, row.names = 1)
colnames(zern3d) = c("F00", "F11", "F20", "F22", "F31", "F33", "F40", "F42", "F44", "F51", "F53", "F55", "F60", "F62", "F64", "F66", "F71", "F73", "F75", "F77", "F80", "F82", "F84", "F86", "F88", "F91", "F93", "F95", "F97", "F99", "F100", "F102", "F104", "F106", "F108", "F1010", "F111", "F113", "F115", "F117", "F119", "F1111", "F120", "F122", "F124", "F126", "F128", "F1210", "F1212", "F131", "F133", "F135", "F137", "F139", "F1311", "F1313", "F140", "F142", "F144", "F146", "F148", "F1410", "F1412", "F1414", "F151", "F153", "F155", "F157", "F159", "F1511", "F1513", "F1515", "F160", "F162", "F164", "F166", "F168", "F1610", "F1612", "F1614", "F1616", "F171", "F173", "F175", "F177", "F179", "F1711", "F1713", "F1715", "F1717", "F180", "F182", "F184", "F186", "F188", "F1810", "F1812", "F1814", "F1816", "F1818", "F191", "F193", "F195", "F197", "F199", "F1911", "F1913", "F1915", "F1917", "F1919", "F200", "F202", "F204", "F206", "F208", "F2010", "F2012", "F2014", "F2016", "F2018", "F2020")
row.names(zern3d) = gsub('-',':',row.names(zern3d))
row.names(zern3d) = gsub('_20$','',row.names(zern3d))

zern4 = zern3d[grepl('_4$', row.names(zern3d)),]
zern6 = zern3d[grepl('_6$', row.names(zern3d)),]
zern8 = zern3d[grepl('_8$', row.names(zern3d)),]
zern10 = zern3d[grepl('_10$', row.names(zern3d)),]

row.names(zern4) = gsub('_4$', '', row.names(zern4))
row.names(zern6) = gsub('_6$', '', row.names(zern6))
row.names(zern8) = gsub('_8$', '', row.names(zern8))
row.names(zern10) = gsub('_10$', '', row.names(zern10))

###########################
# Read in clusters of sequences from CD-HIT at varied identity thresholds (sequences are unique sequences extracted from PDBs, ie only one copy of identical chains)
###########################

# write(unique(bsResiDat$pdb), file = paste("./data/structures/holo/seqs/pdbList.txt", sep = ""), ncol = length(unique(bsResiDat$pdb)), sep = " ") # Space delimited list of PDB IDs to get seqeunces for (cp & paste into ./scripts/misc/getSeqs.sh)

clusterDir = './analysis/seqClustering/'

id50 = readCDHIT(cdhitDir = "./analysis/seqClustering/id50/")
cntID50 = sort(countClustMembers(id50))
id60 = readCDHIT(cdhitDir = "./analysis/seqClustering/id60/")
cntID60 = sort(countClustMembers(id60))
id70 = readCDHIT(cdhitDir = "./analysis/seqClustering/id70/")
cntID70 = sort(countClustMembers(id70))
id80 = readCDHIT(cdhitDir = "./analysis/seqClustering/id80/")
cntID80 = sort(countClustMembers(id80))
id90 = readCDHIT(cdhitDir = "./analysis/seqClustering/id90/")
cntID90 = sort(countClustMembers(id90))

mycol =  colorRampPalette(c("navy","cyan3","darkgreen"))(5)
jitFact = 0.0001
alph = 0.9

# par(mfrow=c(2,2))
plot(0,0, xlim = c(0,400), ylim = c(0,max(c(cntID50, cntID60, cntID70, cntID80, cntID90))), pch = '', cex = 1.3, xlab = '', ylab = '')
title(xlab = 'Cluster (sorted by size)', ylab = 'Number of members in cluster', main = "Distributions of cluster sizes",
      cex.lab = 1.4, cex.main = 1.25)
lines(x = sort(rep(1:length(cntID50),3)), y = c(rbind(rep(0,length(cntID50)), cntID50, rep(0,length(cntID50)))), col = alpha(mycol[1],alph))
# points(x = 1:length(cntID50), y = jitter(cntID50,amount = jitFact), pch = 19, col = alpha(mycol[1],alph))
# plot(0,0, xlim = c(0,max(c(length(cntID50), length(cntID60), length(cntID70), length(cntID80)))), ylim = c(0,max(c(cntID50, cntID60, cntID70, cntID80))), pch = '', xlab = '', ylab = '', main = "60% id")
lines(x = sort(rep(1.2:(length(cntID60)+0.2),3)), y = c(rbind(rep(0,length(cntID60)), cntID60, rep(0,length(cntID60)))), col = alpha(mycol[2],alph))
# points(x = 1:length(cntID60), y = jitter(cntID60,amount = jitFact), pch = 19, col = alpha(mycol[2],alph))
# plot(0,0, xlim = c(0,max(c(length(cntID50), length(cntID60), length(cntID70), length(cntID80)))), ylim = c(0,max(c(cntID50, cntID60, cntID70, cntID80))), pch = '', xlab = 'Cluster ID (sorted)', ylab = 'Number of members in cluster', main = "70% id")
lines(x = sort(rep(1.4:(length(cntID70)+0.4),3)), y = c(rbind(rep(0,length(cntID70)), cntID70, rep(0,length(cntID70)))), col = alpha(mycol[3],alph))
# points(x = 1:length(cntID70), y = jitter(cntID70,amount = jitFact), pch = 19, col = alpha(mycol[3],alph))
# plot(0,0, xlim = c(0,max(c(length(cntID50), length(cntID60), length(cntID70), length(cntID80)))), ylim = c(0,max(c(cntID50, cntID60, cntID70, cntID80))), pch = '', xlab = 'Cluster ID (sorted)', ylab = '', main = "80% id")
lines(x = sort(rep(1.6:(length(cntID80)+0.6),3)), y = c(rbind(rep(0,length(cntID80)), cntID80, rep(0,length(cntID80)))), col = alpha(mycol[4],alph))
# points(x = 1:length(cntID80), y = jitter(cntID80,amount = jitFact), pch = 19, col = alpha(mycol[4],alph))
lines(x = sort(rep(1.8:(length(cntID90)+0.8),3)), y = c(rbind(rep(0,length(cntID90)), cntID90, rep(0,length(cntID90)))), col = alpha(mycol[5],alph))
legend(x = 'topleft', 
       legend = c('50% Identity', '60% Identity', '70% Identity', '80% Identity', '90% Identity'),
       cex = 1.3,
       col = mycol,
       pch = 19, pt.cex = 2)

dev.off()


# hist(clustCnts, breaks = 40, xlim = c(0,40), col = "dodgerblue", xlab = "Number of sequences in a cluster", main = "Cluster sizes")
# plot(density(clustCnts), col = "dodgerblue", lwd = 2, xlab = "Number of sequences in a cluster", main = "Cluster sizes")

clustDF = as.data.frame(matrix('', nrow = length(unique(bsResiDat$pdb)), ncol = 8), stringsAsFactors = F)
colnames(clustDF) = c('pdb','uniprot','fold','clus50','clus60','clus70','clus80','clus90')
clustDF$pdb = unique(bsResiDat$pdb) 
for (i in 1:nrow(clustDF)){
  pdb = clustDF$pdb[i]
  clustDF$uniprot[i] = unique(bsResiDat$uniparc[bsResiDat$pdb == pdb])
  clustDF$fold[i] = unique(bsResiDat$new_fold.Not.Published[bsResiDat$pdb == pdb])
  clustDF$clus50[i] = id50$clusterID[id50$pdbID == pdb]
  clustDF$clus60[i] = id60$clusterID[id60$pdbID == pdb]
  clustDF$clus70[i] = id70$clusterID[id70$pdbID == pdb]
  clustDF$clus80[i] = id80$clusterID[id80$pdbID == pdb]
  clustDF$clus90[i] = id90$clusterID[id90$pdbID == pdb]
}

#uniprot ids per cluster
upc50 = rep(0,length(unique(clustDF$clus50)))
for (i in 1:length(upc50)){upc50[i] = length(unique(clustDF$uniprot[clustDF$clus50 == unique(clustDF$clus50)[i]]))}
upc60 = rep(0,length(unique(clustDF$clus60)))
for (i in 1:length(upc60)){upc60[i] = length(unique(clustDF$uniprot[clustDF$clus60 == unique(clustDF$clus60)[i]]))}
upc70 = rep(0,length(unique(clustDF$clus70)))
for (i in 1:length(upc70)){upc70[i] = length(unique(clustDF$uniprot[clustDF$clus70 == unique(clustDF$clus70)[i]]))}
upc80 = rep(0,length(unique(clustDF$clus80)))
for (i in 1:length(upc80)){upc80[i] = length(unique(clustDF$uniprot[clustDF$clus80 == unique(clustDF$clus80)[i]]))}
upc90 = rep(0,length(unique(clustDF$clus90)))
for (i in 1:length(upc90)){upc90[i] = length(unique(clustDF$uniprot[clustDF$clus90 == unique(clustDF$clus90)[i]]))}

# xLim = c(0,28)
# yLim = c(0,1.1)
# plot(density(upc50), col = mycol[1], lwd = 1.5,  main = 'UniProt IDs per Cluster', xlab = 'Number of unique UniProt IDs',
#      xlim = xLim, ylim = yLim)
# par(new=T)
# plot(density(upc60), col = mycol[2], lwd = 1.5,  main = '', xlab = '', ylab = '', axes = F,
#      xlim = xLim, ylim = yLim)
# par(new=T)
# plot(density(upc70), col = mycol[3], lwd = 1.5,  main = '', xlab = '', ylab = '', axes = F,
#      xlim = xLim, ylim = yLim)
# par(new=T)
# plot(density(upc80), col = mycol[4], lwd = 1.5,  main = '', xlab = '', ylab = '', axes = F,
#      xlim = xLim, ylim = yLim)
# par(new=T)
# plot(density(upc90), col = mycol[5], lwd = 1.5,  main = '', xlab = '', ylab = '', axes = F,
#      xlim = xLim, ylim = yLim)

xLim = c(0,20)
yLim = c(0,1)

par(mfrow=c(2,3),oma = c(4, 4, 0.2, 0.2), mai = c(0.5, 0.1, 0.1, 0.5))
hist(upc50,probability = T, breaks = max(upc50), col = mycol[1], main = '', xlab = '',
     xlim = xLim, ylim = yLim)
hist(upc60,probability = T, breaks = max(upc60), col = mycol[2], main = '', xlab = '',
     xlim = xLim, ylim = yLim)
hist(upc70,probability = T, breaks = max(upc70), col = mycol[3], main = '', xlab = '',
     xlim = xLim, ylim = yLim)
hist(upc80,probability = T, breaks = max(upc80), col = mycol[4], main = '', xlab = '',
     xlim = xLim, ylim = yLim)
hist(upc90,probability = T, breaks = max(upc90), col = mycol[5], main = '', xlab = '',
     xlim = xLim, ylim = yLim)
plot(0,0, col = 'white',
     axes =F, xlab = '', ylab ='')
legend(x='center', col = mycol, legend = c('50% Identity', '60% Identity', '70% Identity', '80% Identity', '90% Identity'), cex = 3, pch = 19, pt.cex = 4.5, bty = 'n')
mtext('Histograms of counts of unique UniProt IDs in each cluster', side = 3, line = -1.5, outer = TRUE, cex = 1.5)
mtext('Number of unique UniProt IDs', side = 1, line = 0, outer = TRUE, cex = 1.25)
mtext('Density', side = 2, line = 2, outer = TRUE, cex = 1.25)

dev.off()

#clusters per uniprot ID???
cpu50 = rep(0,length(unique(clustDF$uniprot)))
for (i in 1:length(cpu50)){cpu50[i] = length(unique(clustDF$clus50[clustDF$uniprot == unique(clustDF$uniprot)[i]]))}
cpu60 = rep(0,length(unique(clustDF$uniprot)))
for (i in 1:length(cpu60)){cpu60[i] = length(unique(clustDF$clus60[clustDF$uniprot == unique(clustDF$uniprot)[i]]))}
cpu70 = rep(0,length(unique(clustDF$uniprot)))
for (i in 1:length(cpu70)){cpu70[i] = length(unique(clustDF$clus70[clustDF$uniprot == unique(clustDF$uniprot)[i]]))}
cpu80 = rep(0,length(unique(clustDF$uniprot)))
for (i in 1:length(cpu80)){cpu80[i] = length(unique(clustDF$clus80[clustDF$uniprot == unique(clustDF$uniprot)[i]]))}
cpu90 = rep(0,length(unique(clustDF$uniprot)))
for (i in 1:length(cpu90)){cpu90[i] = length(unique(clustDF$clus90[clustDF$uniprot == unique(clustDF$uniprot)[i]]))}

# par(mfrow=c(2,3),oma = c(4, 4, 0.2, 0.2), mai = c(0.5, 0.1, 0.1, 0.5))
# plot(density(cpu50), col = mycol[1], lwd = 3,  main = '')
# plot(density(cpu60), col = mycol[2], lwd = 3,  main = '')
# plot(density(cpu70), col = mycol[3], lwd = 3,  main = '')
# plot(density(cpu80), col = mycol[4], lwd = 3,  main = '')
# plot(density(cpu90), col = mycol[5], lwd = 3,  main = '')
# dev.off()

xLim = c(1,3)
yLim = c(0,2)

par(mfrow=c(2,3),oma = c(4, 4, 0.2, 0.2), mai = c(0.5, 0.1, 0.1, 0.5))
hist(cpu50,probability = T, breaks = 3, col = mycol[1], main = '', xlab = '',
     xlim = xLim, ylim = yLim)
hist(cpu60,probability = T, breaks = 3, col = mycol[2], main = '', xlab = '',
     xlim = xLim, ylim = yLim)
hist(cpu70,probability = T, breaks = 3, col = mycol[3], main = '', xlab = '',
     xlim = xLim, ylim = yLim)
hist(cpu80,probability = T, breaks = 3, col = mycol[4], main = '', xlab = '',
     xlim = xLim, ylim = yLim)
hist(cpu90,probability = T, breaks = 3, col = mycol[5], main = '', xlab = '',
     xlim = xLim, ylim = yLim)
plot(0,0, col = 'white',
     axes =F, xlab = '', ylab ='')
legend(x='center', col = mycol, legend = c('50% Identity', '60% Identity', '70% Identity', '80% Identity', '90% Identity'), cex = 3, pch = 19, pt.cex = 4.5, bty = 'n')
mtext('Histograms of numbers of clusters each UniprotID can be found in', side = 3, line = -2, outer = TRUE, cex = 1.5)
mtext('Number of clusters per UniProt ID', side = 1, line = 0, outer = TRUE, cex = 1.25)
mtext('Density', side = 2, line = 2, outer = TRUE, cex = 1.25)

dev.off()


# unique(clustDF$uniprot)[cpu80 == 3]
# unique(clustDF$uniprot)[cpu90 == 3]

# clustDF[clustDF$pdb %in% c('4MBZ', '4FMI', '4FMH','4FMJ'),c(1,2,4:8)]

bsResiDat$seqClust50 = bsResiDat$seqClust80 = 0
for (i in 1:nrow(bsResiDat)){
  bsResiDat$seqClust50[i] = id50$clusterID[id50$pdbID == bsResiDat$pdb[i]]
  bsResiDat$seqClust80[i] = id80$clusterID[id80$pdbID == bsResiDat$pdb[i]]
}

clusLst50 = unique(bsResiDat$seqClust50)
clusLst80 = unique(bsResiDat$seqClust80)

###########################
# Find ligands well represented in clusters
###########################

# For each cluster, how many different ligands are present?
lpc50 = lpc50Len = rep(0, length(clusLst50))
lpc80 = lpc80Len = rep(0, length(clusLst80))

for (i in 1:length(lpc50)){
  lpc50[i] = length(unique(bsResiDat$iupac[bsResiDat$seqClust50 == clusLst50[i]]))
  lpc50Len[i] = sum(id50$clusterID == clusLst50[i])
}
for (i in 1:length(lpc80)){
  lpc80[i] = length(unique(bsResiDat$iupac[bsResiDat$seqClust80 == clusLst80[i]]))
  lpc80Len[i] = sum(id80$clusterID == clusLst80[i])
}

par(mfrow=c(1,2))
hist(lpc50, col = mycol[1], breaks = 12, xlim = c(0,20), probability = T, ylim = c(0,0.8), xlab = 'ligands per cluster', main = '50% id cutoff')
hist(lpc80, col = mycol[4], breaks = 12, xlim = c(0,20), probability = T, ylim = c(0,0.8), xlab = 'ligands per cluster', main = '80% id cutoff')
mtext('Histograms of numbers of unique ligands in each cluster', side = 3, line = -1, outer = TRUE, cex = 1.5)

dev.off()


# tag = (lpc50 < 3 & lpc50Len >15)
# clusLst50[tag]
# bsResiDat[bsResiDat$seqClust50 == 14,c(1,2,3,4,5,7,8,9,176)]
# bsResiDat[bsResiDat$seqClust50 == 30,c(1,2,3,4,5,7,8,9,176)]

lpc50Len = lpc50Len[order(lpc50, decreasing = T)]
lpc50 = sort(lpc50, decreasing = T)
for (i in 1:length(unique(lpc50))){
  tag = lpc50 == unique(lpc50)[i]
  lpc50Len[tag] = sort(lpc50Len[tag], decreasing = T)
}
lpc80Len = lpc80Len[order(lpc80, decreasing = T)]
lpc80 = sort(lpc80, decreasing = T)
for (i in 1:length(unique(lpc80))){
  tag = lpc80 == unique(lpc80)[i]
  lpc80Len[tag] = sort(lpc80Len[tag], decreasing = T)
}


par(mfrow=c(2,1), mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(0,0,pch='', xlim = c(0,250), ylim = c(0,17), xlab = '', ylab = '')
title(main="Number of ligands per cluster (with cluster size)", col.main="black", cex.main = 2,
      xlab='Individual clusters (sorted)',
      col.lab='black', cex.lab=1.5)
title(ylab='Number of unique ligands in cluster',
      col.lab = mycol[1], cex.lab = 1.5)
lines(x = sort(rep(1:length(lpc50),3)), y = c(rbind(rep(0,length(lpc50)), lpc50, rep(0,length(lpc50)))), col = mycol[1], lwd = 2)
par(new=T)
plot(0,0,pch='', xlim = c(0,250), ylim = c(0,60),
     xlab = '', ylab = '', main = '', axes = F)
lines(x = sort(rep(1.5:(length(lpc50Len)+0.5),3)), y = c(rbind(rep(0,length(lpc50Len)), lpc50Len, rep(0,length(lpc50Len)))), col = 'firebrick', lwd = 2)
axis(side = 4, at = pretty(c(0,60)))      # Add second axis
mtext("Number of structures per cluster", side = 4, line = 3, cex = 1.5, col = 'firebrick')             # Add second axis label
legend(x= 'topright', legend = '50% Identity', pch = 15, col = mycol[1], bty = 'n', cex = 2, pt.cex = 3)

plot(0,0,pch='', xlim = c(0,350), ylim = c(0,17), xlab = '', ylab = '')
title(main="", col.main="black", cex.main = 2,
      xlab='Individual clusters (sorted)',
      col.lab='black', cex.lab=1.5)
title(ylab='Number of unique ligands in cluster',
      col.lab = mycol[4], cex.lab = 1.5)
lines(x = sort(rep(1:length(lpc80),3)), y = c(rbind(rep(0,length(lpc80)), lpc80, rep(0,length(lpc80)))), col = mycol[4], lwd = 2)
par(new=T)
plot(0,0,pch='', xlim = c(0,350), ylim = c(0,60),
     xlab = '', ylab = '', main = '', axes = F)
lines(x = sort(rep(1.5:(length(lpc80Len)+0.5),3)), y = c(rbind(rep(0,length(lpc80Len)), lpc80Len, rep(0,length(lpc80Len)))), col = 'firebrick', lwd = 2)
axis(side = 4, at = pretty(c(0,60)))      # Add second axis
mtext("Number of structures per cluster", side = 4, line = 3, cex = 1.5, col = 'firebrick')             # Add second axis label
legend(x= 'topright', legend = '80% Identity', pch = 15, col = mycol[4], bty = 'n', cex = 2, pt.cex = 3)

dev.off()




# For each ligand, how many different clusters can it be found in?
uniLigs = unique(bsResiDat$iupac)
# length(uniLigs)
cpl50  = cpl80 = rep(0, length(uniLigs))
for (i in 1:length(uniLigs)){
  cpl50[i] = length(unique(bsResiDat$seqClust50[bsResiDat$iupac == uniLigs[i]]))
  cpl80[i] = length(unique(bsResiDat$seqClust80[bsResiDat$iupac == uniLigs[i]]))
}

plot(density(cpl50))
plot(density(cpl80))

ligSort50 = uniLigs[order(cpl50, decreasing = T)]
cpl50 = sort(cpl50/length(clusLst50)*100, decreasing = T)
ligSort80 = uniLigs[order(cpl80, decreasing = T)]
cpl80 = sort(cpl80/length(clusLst80)*100, decreasing = T)

topLigOccurences = as.data.frame(matrix(0,nrow = 2, ncol = 16))
colnames(topLigOccurences) = c('id_cutoff', 'High_Mannose', 'Sialic_Acid', 'Terminal_Fucose', rep('', 12))
topLigOccurences$id_cutoff = c('50%', '80%')

## 50% id

parenCnt = bracCnt = manCnt = neuCnt = bracCnt = rep(0,length(ligSort50))
for (i in 1:length(ligSort50)){
  lig = ligSort50[i]
  parenCnt[i] = lengths(regmatches(lig, gregexpr("\\(", lig)))
  bracCnt[i] = lengths(regmatches(lig, gregexpr("\\[", lig)))
  manCnt[i] = lengths(regmatches(lig, gregexpr("Man", lig)))
  neuCnt[i] = lengths(regmatches(lig, gregexpr("Neu", lig)))
}

mTag = parenCnt == 0 & bracCnt == 0 # Monosaccharides
dTag = parenCnt == 1 & bracCnt == 0 # Disaccharides
tTag = (parenCnt == 2 & bracCnt == 0) | (parenCnt == 2 & bracCnt == 1) # Trisaccharides
qTag = (parenCnt == 3 & bracCnt == 0) | (parenCnt == 3 & bracCnt == 1) | (parenCnt == 1 & bracCnt == 2) # Tetrasaccharides
pTag = !(mTag | dTag | tTag | qTag) # 5+ sugars
bTag = bracCnt >= 1 # Branched glycans

hist(manCnt)

manTag = manCnt > 3 # High mannose
neuTag = neuCnt >= 1 # Has sialic acid
fucTag = grepl('^Fuc',ligSort50) # Has a terminal fucose

# ligSort50[mTag]
# ligSort50[dTag]
# ligSort50[tTag]
# ligSort50[qTag]
# ligSort50[pTag]
# 
# ligSort50[bTag]
# 
# ligSort50[manTag]
# ligSort50[neuTag]
# ligSort50[fucTag]

sacc_col =  colorRampPalette(c("plum1","tomato", "firebrick4"))(5)

plot(0,0,pch='', xlim = c(0.5,230), ylim = c(-1,20), xlab = '', ylab = '')
title(main="Number of clusters with any binding for each unique ligand (50% id)", col.main="black", cex.main = 2,
      xlab='Individual ligands (sorted)',
      col.lab='black', cex.lab=1.5)
title(ylab='Clusters containing ligand (%)',
      col.lab = 'black', cex.lab = 1.5)
lines(x = sort(rep(1:length(cpl50),3)), y = c(rbind(rep(0,length(cpl50)), cpl50, rep(0,length(cpl50)))), col = mycol[1], lwd = 5)
points(x=(1:length(mTag))[mTag], y = -0.25* mTag[mTag], pch = 15, col = sacc_col[1])
points(x=(1:length(dTag))[dTag], y = -0.25* dTag[dTag], pch = 15, col = sacc_col[2])
points(x=(1:length(tTag))[tTag], y = -0.25* tTag[tTag], pch = 15, col = sacc_col[3])
points(x=(1:length(qTag))[qTag], y = -0.25* qTag[qTag], pch = 15, col = sacc_col[4])
points(x=(1:length(pTag))[pTag], y = -0.25* pTag[pTag], pch = 15, col = sacc_col[5])
# text(x = 232, y =0, labels = '# of sugars', cex =1.1)
points(x=(1:length(bTag))[bTag], y = -0.5* bTag[bTag], pch = 15, col = 'navy')
points(x=(1:length(bTag))[!bTag], y = -0.5+ bTag[!bTag], pch = 15, col = 'gray85')
# text(x = 232, y =-0.75, labels = 'Branching?', cex =1.1)
points(x=(1:length(manTag))[manTag], y = -0.75* manTag[manTag], pch = 21, bg = 'forestgreen', col ='forestgreen')
points(x=(1:length(neuTag))[neuTag], y = -0.75* neuTag[neuTag], pch = 23, bg = 'darkorchid', col ='darkorchid')
points(x=(1:length(fucTag))[fucTag], y = -0.75* fucTag[fucTag], pch = 24, bg = 'firebrick1', col ='firebrick1')
# text(x = 232, y =-1.5, labels = 'Composition', cex =1.1)
legend(x = 'topright', pch = c(rep(15,9),21,23, 24),
       legend = c('Monosacc.','Disacc.','Trisacc.','Tetasacc.','5+ sugars',
                  '',
                  'Branching','No Branching',
                  '',
                  'High Mannose', 'NeuAc', 'Terminal Fucose'),
       col = c(sacc_col,
               'white',
               'navy', 'grey85',
               'white',
               'white', 'white', 'white'),
       pt.bg = c(rep(NA,9),
              'forestgreen', 'darkorchid', 'firebrick1'),
       pt.cex =1.8)

# ligSort50[manTag]
topLigOccurences$High_Mannose[1] = sum(cpl50[manTag])
# ligSort50[neuTag]
topLigOccurences$Sialic_Acid[1] = sum(cpl50[neuTag])
# ligSort50[fucTag]
topLigOccurences$Terminal_Fucose[1] = sum(cpl50[fucTag])

## 80% id

parenCnt = bracCnt = manCnt = neuCnt = bracCnt = rep(0,length(ligSort80))
for (i in 1:length(ligSort80)){
  lig = ligSort80[i]
  parenCnt[i] = lengths(regmatches(lig, gregexpr("\\(", lig)))
  bracCnt[i] = lengths(regmatches(lig, gregexpr("\\[", lig)))
  manCnt[i] = lengths(regmatches(lig, gregexpr("Man", lig)))
  neuCnt[i] = lengths(regmatches(lig, gregexpr("Neu", lig)))
}

mTag = parenCnt == 0 & bracCnt == 0 # Monosaccharides
dTag = parenCnt == 1 & bracCnt == 0 # Disaccharides
tTag = (parenCnt == 2 & bracCnt == 0) | (parenCnt == 2 & bracCnt == 1) # Trisaccharides
qTag = (parenCnt == 3 & bracCnt == 0) | (parenCnt == 3 & bracCnt == 1) | (parenCnt == 1 & bracCnt == 2) # Tetrasaccharides
pTag = !(mTag | dTag | tTag | qTag) # 5+ sugars
bTag = bracCnt >= 1 # Branched glycans

hist(manCnt)

manTag = manCnt > 3 # High mannose
neuTag = neuCnt >= 1 # Has sialic acid
fucTag = grepl('^Fuc',ligSort80) # Has a terminal fucose

plot(0,0,pch='', xlim = c(0.5,230), ylim = c(-1,20), xlab = '', ylab = '')
title(main="Number of clusters with any binding for each unique ligand (80% id)", col.main="black", cex.main = 2,
      xlab='Individual ligands (sorted)',
      col.lab='black', cex.lab=1.5)
title(ylab='Clusters containing ligand (%)',
      col.lab = 'black', cex.lab = 1.5)
lines(x = sort(rep(1:length(cpl80),3)), y = c(rbind(rep(0,length(cpl80)), cpl80, rep(0,length(cpl80)))), col = mycol[4], lwd = 5)
points(x=(1:length(mTag))[mTag], y = -0.25* mTag[mTag], pch = 15, col = sacc_col[1])
points(x=(1:length(dTag))[dTag], y = -0.25* dTag[dTag], pch = 15, col = sacc_col[2])
points(x=(1:length(tTag))[tTag], y = -0.25* tTag[tTag], pch = 15, col = sacc_col[3])
points(x=(1:length(qTag))[qTag], y = -0.25* qTag[qTag], pch = 15, col = sacc_col[4])
points(x=(1:length(pTag))[pTag], y = -0.25* pTag[pTag], pch = 15, col = sacc_col[5])
# text(x = 232, y =0, labels = '# of sugars', cex =1.1)
points(x=(1:length(bTag))[bTag], y = -0.5* bTag[bTag], pch = 15, col = 'navy')
points(x=(1:length(bTag))[!bTag], y = -0.5+ bTag[!bTag], pch = 15, col = 'gray85')
# text(x = 232, y =-0.75, labels = 'Branching?', cex =1.1)
points(x=(1:length(manTag))[manTag], y = -0.75* manTag[manTag], pch = 21, bg = 'forestgreen', col ='forestgreen')
points(x=(1:length(neuTag))[neuTag], y = -0.75* neuTag[neuTag], pch = 23, bg = 'darkorchid', col ='darkorchid')
points(x=(1:length(fucTag))[fucTag], y = -0.75* fucTag[fucTag], pch = 24, bg = 'firebrick1', col ='firebrick1')
# text(x = 232, y =-1.5, labels = 'Composition', cex =1.1)
legend(x = 'topright', pch = c(rep(15,9),21,23, 24),
       legend = c('Monosacc.','Disacc.','Trisacc.','Tetasacc.','5+ sugars',
                  '',
                  'Branching','No Branching',
                  '',
                  'High Mannose', 'NeuAc', 'Terminal Fucose'),
       col = c(sacc_col,
               'white',
               'navy', 'grey85',
               'white',
               'white', 'white', 'white'),
       pt.bg = c(rep(NA,9),
                 'forestgreen', 'darkorchid', 'firebrick1'),
       pt.cex =1.8)

# ligSort80[manTag]
topLigOccurences$High_Mannose[2] = sum(cpl80[manTag])
# ligSort80[neuTag]
topLigOccurences$Sialic_Acid[2] = sum(cpl80[neuTag])
# ligSort80[fucTag]
topLigOccurences$Terminal_Fucose[2] = sum(cpl80[fucTag])

top50 = ligSort50[cpl50 > 5]
top80 = ligSort80[cpl80 > 5]

# # Ligands present in >5% of clusters are consistent between ID thresholds
# all(top50 %in% top80)
# all(top80 %in% top50)

colnames(topLigOccurences)[5:16] = top50
topLigOccurences[1, 5:16] = cpl50[cpl50 > 5]

for (i in 5:16){
  lig = colnames(topLigOccurences)[i]
  lig = gsub('\\(', '\\\\(', lig)
  lig = gsub('\\)', '\\\\)', lig)
  lig = gsub('\\[', '\\\\]', lig)
  lig = gsub('\\[', '\\\\]', lig)
  lig = gsub('^', '\\^', lig)
  lig = gsub('$', '\\$', lig)
  topLigOccurences[2,i] = cpl80[grepl(lig, ligSort80)]
}

melt_ligOccur = melt(data =topLigOccurences, id.vars = 'id_cutoff')
colnames(melt_ligOccur) = c('Clust_ID', 'Ligand', 'Cluster_Percent_w_ligand')


ggplot(melt_ligOccur, aes(fill = Clust_ID, x = Ligand, y = Cluster_Percent_w_ligand)) + geom_bar(stat="identity", color="black", position=position_dodge())+
  scale_fill_manual(values=mycol[c(1,4)]) +
  geom_hline(yintercept = c(5)) +
  theme_linedraw(base_size = 14) +
  labs(title = "Well-represented ligands in lectin clusters", x = "IUPAC Ligands", y = "Percent of clusters with ligand") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))


###########################
# McCaldon analysis
###########################
mccaldon$aa = lapply(mccaldon$aa, FUN = aaa) # convert from 1-letter aminoacid code to 3, capitalized
mccaldon$aa = lapply(mccaldon$aa, FUN = str_to_upper)

mccaldon$freq_bin1 = 0
mccaldon$freq_bin2 = 0
mccaldon$freq_bin3 = 0
mccaldon$freq_bin4 = 0

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
# Extract features from d2Dists and d2Feats
###########################
corrplot(cor(d2Feats,d2Feats))

tag = rep(F,ncol(d2Dists)) # Filter out d2Dist columns with all zeroes
for (i in 1:ncol(d2Dists)){
  if (all(d2Dists[,i] == 0)){
    tag[i] = T
  }
}
if (sum(tag) > 0){ d2Dists = d2Dists[,!tag]}

pca = prcomp(x = d2Dists, scale. = T)
for( i in 1:10){ # How much variance is captured by the top principle components
  v = as.character(round(pca$sdev[i]**2/sum(pca$sdev**2), 4) * 100)
  print( paste('PC', as.character(i), ' accounts for ', v, '% of the total variance', sep = '') )
}
plot(pca$x, main = 'PCA of d2Dists (raw bin counts)', xlab = 'PC1 (61.52% of var.)', ylab = 'PC2 (20.67% of var.)')

corrplot(cor(pca$x[,1:10], d2Feats))

thresholds = c('_4Ang', '_6Ang', '_8Ang', '_10Ang')
for (i in 1:nrow(d2Dists)){
  for (t in 1:length(thresholds)){
    tag = grepl(thresholds[t], colnames(d2Dists)) # get columns for the threshold
    if (! all(d2Dists[i,tag] == 0)){ # avoid /0 error
      d2Dists[i,tag] = d2Dists[i,tag] / sum(d2Dists[i,tag])
    }
  }
}

pca.d2Dists = prcomp(x = d2Dists, scale. = T)
runSum = 0
for( i in 1:10){ # How much variance is captured by the top principle components
  v = round(pca.d2Dists$sdev[i]**2/sum(pca.d2Dists$sdev**2), 4) * 100
  runSum = runSum + v
  print( paste('PC', as.character(i), ' accounts for ', as.character(v), '% of the total variance', sep = '') )
}
print( paste('Top 10 PCs account for a total of ', as.character(runSum), '% of the total variance', sep = '') )

plot(pca.d2Dists$x, main = 'PCA of d2Dists (% comp per pocket)', xlab = 'PC1 (41.41% of var.)', ylab = 'PC2 (24.33% of var.)')

corrplot(cor(pca.d2Dists$x[,1:10], d2Feats))

# Map to roadmark pockets
smoothScatter(pca.d2Dists$x, main = 'PCA of d2Dists (% comp per pocket)', xlab = 'PC1 (41.41% of var.)', ylab = 'PC2 (24.33% of var.)')
pdb = '3LL2'
points(pca.d2Dists$x[grepl(pdb,row.names(d2Dists)), 1], pca.d2Dists$x[grepl(pdb,row.names(d2Dists)), 2], col = alpha('green',0.7), pch=19)
pdb = '2CL8'
points(pca.d2Dists$x[grepl(pdb,row.names(d2Dists)), 1], pca.d2Dists$x[grepl(pdb,row.names(d2Dists)), 2], col = alpha('orange2',0.7), pch=19) # 1 deep medium sized pockets only
pdb = '4URR'
points(pca.d2Dists$x[grepl(pdb,row.names(d2Dists)), 1], pca.d2Dists$x[grepl(pdb,row.names(d2Dists)), 2], col = alpha('red',0.7), pch=19) # deep long pocket
pdb = '1HJV'
points(pca.d2Dists$x[grepl(pdb,row.names(d2Dists)), 1], pca.d2Dists$x[grepl(pdb,row.names(d2Dists)), 2], col = alpha('purple',0.7), pch=19) # one 
pdb = '2WGC'
points(pca.d2Dists$x[grepl(pdb,row.names(d2Dists)), 1], pca.d2Dists$x[grepl(pdb,row.names(d2Dists)), 2], col = alpha('blue',0.7), pch=19) # two 0 pockets only
legend(x ="topright", legend = c('mGRFT', '2CL8', '2WGC', '4URR', '1HJV'), col = c('green','orange2', 'blue', 'red', 'purple'), pch = 19)

# volMod = lm(d2Feats$vol_4Ang ~ pca.d2Dists$x[,1] + pca.d2Dists$x[,2])
# summary(volMod)


# PCA of d2 distribution features
pca = prcomp(x = d2Feats, scale. = T)
runSum = 0
for( i in 1:10){ # How much variance is captured by the top principle components
  v = round(pca$sdev[i]**2/sum(pca$sdev**2), 4) * 100
  runSum = runSum + v
  print( paste('PC', as.character(i), ' accounts for ', as.character(v), '% of the total variance', sep = '') )
}
print( paste('Top 10 PCs account for a total of ', as.character(runSum), '% of the total variance', sep = '') )

plot(pca$x, main = 'PCA of d2Feats', xlab = 'PC1 (56.96% of var.)', ylab = 'PC2 (14.03% of var.)')


###########################
# Explore predictive feature set
###########################

row.names(zern8)[ ! (row.names(zern8) %in% row.names(zern4)) ]
d2Feats[grepl('2WBW',row.names((d2Feats))),] # No volume in smallest pocket threshold (probably below 15 Ang^3)

zern4 = rbind(zern4, rep(0,ncol(zern4))) # Add a row of zeros to Zern4 for the missing pocket at 4Ang threshold
row.names(zern4)[nrow(zern4)] = row.names(zern8)[ ! (row.names(zern8) %in% row.names(zern4)) ] # Name the row of zeros
zern4 = zern4[row.names(zern8),] # Re-order zern4 to match other zerns


# zern6 = zern6[!grepl('2WBW', row.names(zern6)),]
# zern8 = zern8[!grepl('2WBW', row.names(zern8)),]
# zern10 = zern10[!grepl('2WBW', row.names(zern10)),]

all(row.names(zern4) == row.names(zern6))
all(row.names(zern4) == row.names(zern8))
all(row.names(zern4) == row.names(zern10))

zern4Norm = t(apply(zern4, 1, invarNorm))
zern6Norm = t(apply(zern6, 1, invarNorm))
zern8Norm = t(apply(zern8, 1, invarNorm))
zern10Norm = t(apply(zern10, 1, invarNorm))

allZernNorm = as.data.frame(cbind(zern4Norm, zern6Norm, zern8Norm, zern10Norm))

all(row.names(allZernNorm) %in% row.names(bsResiDat))
colors = colfunc(length(unique(bsResiDat$seqClust50)))
labels = bsResiDat[row.names(allZernNorm), 'seqClust50']
labels.u = sort(unique(labels))
labels.int = rep(0,length(labels))
for (i in 1:length(labels.u)){
  labels.int[labels == labels.u[i]] = i
}


pca = prcomp(x = allZernNorm, scale. = T)
runSum = 0
for( i in 1:10){ # How much variance is captured by the top principle components
  v = round(pca$sdev[i]**2/sum(pca$sdev**2), 4) * 100
  runSum = runSum + v
  print( paste('PC', as.character(i), ' accounts for ', as.character(v), '% of the total variance', sep = '') )
}
print( paste('Top 10 PCs account for a total of ', as.character(runSum), '% of the total variance', sep = '') )

plot(pca$x, xlab = 'PC1 (13.88% of variance)', ylab = 'PC2 (7.75% of variance)', main = 'PCA of 3DZDs colored by cluster membership (50% id)', pch = 19, col = alpha(colors[labels.int], 0.6))

corrplot(cor(pca.d2Dists$x[row.names(d2Dists) %in% row.names(allZernNorm),1:10],pca$x[row.names(d2Dists)[row.names(d2Dists) %in% row.names(allZernNorm)],1:10]))

# zern.umap = umap(allZernNorm)

plot(x = zern.umap$layout[,1], y = zern.umap$layout[,2], xlim = c(-8,6), ylim = c(-8,6),
     pch = 19, col = alpha(colors[labels.int], 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = 'UMAP vis for 3DZDs (colored by 50% id clusters)')
pdb = '3LL2'
points(zern.umap$layout[grepl(pdb,row.names(allZernNorm)), 1], zern.umap$layout[grepl(pdb,row.names(allZernNorm)), 2], bg = alpha('green',0.9),col='white', pch=22, cex=3)
pdb = '2CL8'
points(zern.umap$layout[grepl(pdb,row.names(allZernNorm)), 1], zern.umap$layout[grepl(pdb,row.names(allZernNorm)), 2], bg = alpha('orange2',0.9),col='white', pch=22, cex=3) # 1 deep medium sized pockets only
pdb = '4URR'
points(zern.umap$layout[grepl(pdb,row.names(allZernNorm)), 1], zern.umap$layout[grepl(pdb,row.names(allZernNorm)), 2], bg = alpha('red',0.9),col='white', pch=22, cex=3) # deep long pocket
pdb = '1HJV'
points(zern.umap$layout[grepl(pdb,row.names(allZernNorm)), 1], zern.umap$layout[grepl(pdb,row.names(allZernNorm)), 2], bg = alpha('purple',0.9),col='white', pch=22, cex=3) # half large deep, half very shallow 

plot(x = zern.umap$layout[,1], y = zern.umap$layout[,2],# xlim = c(-8,6), ylim = c(-8,6),
     pch = 19, col = alpha(colors[labels.int], 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = 'UMAP vis for 3DZDs (colored by 50% id clusters)')

zern.umap$layout[(zern.umap$layout[,2] < -6),]

# tag = (zern.umap$layout[,1] < -20) & (zern.umap$layout[,2] < -20)
# zern.umap$layout[tag,]
# outlierLst = row.names(zern.umap$layout)[tag]
# tmp = apply(X = d2Feats, 2, mean)
# tmpO = apply(X = d2Feats[outlierLst,], 2, mean)
# barplot((tmpO-tmp)/tmp)
# abline(a = 0, b =1)
# summary(d2Feats)

missingRows = row.names(bsResiDat)[!(row.names(bsResiDat) %in% row.names(allZernNorm))]
zeroes2Zern = as.data.frame(matrix(0,nrow = length(missingRows), ncol = ncol(allZernNorm)))
row.names(zeroes2Zern) = missingRows
colnames(zeroes2Zern) = colnames(allZernNorm)

allZernNorm = rbind(allZernNorm, zeroes2Zern)

all(row.names(allZernNorm) %in% row.names(bsResiDat)) & all(row.names(bsResiDat) %in% row.names(allZernNorm))

allZernNorm = allZernNorm[row.names(bsResiDat),]

all(row.names(bsResiDat) == row.names(allZernNorm))
all(row.names(bsResiDat) == row.names(d2Feats))
all(row.names(bsResiDat) == row.names(pca.d2Dists$x))


# Clean up residue based features
dropCols = rep(F,ncol(bsResiDat))
for (i in 1:ncol(bsResiDat)){
  if (all(bsResiDat[,i] == 0)){
    dropCols[i] = T
    print(colnames(bsResiDat)[i])
  }
}

sum(! bsResiDat$ALT_bin1 == 0)
sum(! bsResiDat$ALT_bin2 == 0)
sum(! bsResiDat$I_bin4 == 0)

sum(dropCols)
dropCols[grep("ALT", colnames(bsResiDat))] = T
dropCols[grep("^I_bin", colnames(bsResiDat))] = T
sum(dropCols)

bsResiDat = bsResiDat[,!dropCols]

resiPredFeats = bsResiDat[,grepl('_bin', colnames(bsResiDat))]
resiPredFeats$numBSresis_bin1 = sqrt(resiPredFeats$numBSresis_bin1)/max(sqrt(resiPredFeats$numBSresis_bin1))


