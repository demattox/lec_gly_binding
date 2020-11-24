
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
library(pheatmap)
library(VennDiagram)

library(sjstats)

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

outlierPlots <- function(df, tag){
  for (i in (1:ncol(df))){
    df[,i] = (df[,i] - min(df[,i])) / (max(df[,i]) - min(df[,i]))
  }
  df$outlier = 'other'
  df$outlier[tag] = 'outlier'
  mdf = melt(df, id.vars = 'outlier')
  ggplot(data = mdf, aes(x = variable, y = value, color = outlier, fill = variable)) +
    geom_boxplot() + 
    scale_color_manual(values=c("black", "magenta")) + 
    labs(title = 'Outlier investigation', x = 'variables', y = 'Values (scaled from [0,1])') + 
    guides(fill=FALSE) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = 'bold.italic', color = 'black'))
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
# Separate skew into postive and negative
tag4 = d2Feats$skew_4Ang < 0
tag6 = d2Feats$skew_6Ang < 0
tag8 = d2Feats$skew_8Ang < 0
tag10 = d2Feats$skew_10Ang < 0

d2Feats$leftskew_4Ang = d2Feats$skew_4Ang # Copy skew column to negative skew
d2Feats$leftskew_4Ang[!tag4] = 0 # Zero out non-negative skew values in left skew column
d2Feats$leftskew_4Ang = d2Feats$leftskew_4Ang*-1 # # Set left skew values to postive
d2Feats$skew_4Ang[tag4] = 0 # zero out left skewed values in the original skew column

d2Feats$leftskew_6Ang = d2Feats$skew_6Ang
d2Feats$leftskew_6Ang[!tag6] = 0
d2Feats$leftskew_6Ang = d2Feats$leftskew_6Ang*-1
d2Feats$skew_6Ang[tag6] = 0

d2Feats$leftskew_8Ang = d2Feats$skew_8Ang
d2Feats$leftskew_8Ang[!tag8] = 0
d2Feats$leftskew_8Ang = d2Feats$leftskew_8Ang*-1
d2Feats$skew_8Ang[tag8] = 0

d2Feats$leftskew_10Ang = d2Feats$skew_10Ang
d2Feats$leftskew_10Ang[!tag10] = 0
d2Feats$leftskew_10Ang = d2Feats$leftskew_10Ang*-1
d2Feats$skew_10Ang[tag10] = 0


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

# Read in D2 distances binned into an equal number of bins for each shape
d2ScaledBins = read.delim('./analysis/d2_scaled_bins.csv', header = T, sep =',', stringsAsFactors = F, row.names = 1)
row.names(bsResiDat)[!(row.names(bsResiDat) %in% row.names(d2ScaledBins))] # Rows missing from d2ScaledBins
d2Feats[row.names(bsResiDat)[!(row.names(bsResiDat) %in% row.names(d2ScaledBins))],] # Zero volume binding sites

zeroRows = as.data.frame(matrix(0,nrow = sum(! row.names(bsResiDat) %in% row.names(d2ScaledBins) ), ncol = ncol(d2ScaledBins)))
row.names(zeroRows) = row.names(bsResiDat)[!(row.names(bsResiDat) %in% row.names(d2ScaledBins))]
colnames(zeroRows) = colnames(d2ScaledBins)
d2ScaledBins =  rbind(d2ScaledBins, zeroRows)

all(row.names(bsResiDat) %in% row.names(d2Feats)) & nrow(d2Feats) == nrow(bsResiDat) # All dataframes have the same rows names
all(row.names(bsResiDat) %in% row.names(d2Dists)) & nrow(d2Dists) == nrow(bsResiDat)
all(row.names(bsResiDat) %in% row.names(d2ScaledBins)) & nrow(d2ScaledBins) == nrow(bsResiDat)

# Get dataframes in the same order
d2Feats = d2Feats[row.names(bsResiDat),]
d2Dists = d2Dists[row.names(bsResiDat),]
d2ScaledBins = d2ScaledBins[row.names(bsResiDat),]


all(row.names(bsResiDat) == row.names(d2Feats)) & all(row.names(bsResiDat) == row.names(d2Dists)) & all(row.names(bsResiDat) == row.names(d2ScaledBins)) # double check

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

# mycol =  colorRampPalette(c("navy","cyan3","darkgreen"))(5)
mycol = colfunc(5)
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


pdf(file = paste('./analysis/sfgPlots/', 
                 'ligands_per_clust',
                 '.pdf', sep = ''),
    width = 12,
    height = 10)
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
lines(x = sort(rep(1.5:(length(lpc50Len)+0.5),3)), y = c(rbind(rep(0,length(lpc50Len)), lpc50Len, rep(0,length(lpc50Len)))), col = 'black', lwd = 2)
axis(side = 4, at = pretty(c(0,60)))      # Add second axis
mtext("Number of structures per cluster", side = 4, line = 3, cex = 1.5, col = 'black')             # Add second axis label
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
lines(x = sort(rep(1.5:(length(lpc80Len)+0.5),3)), y = c(rbind(rep(0,length(lpc80Len)), lpc80Len, rep(0,length(lpc80Len)))), col = 'black', lwd = 2)
axis(side = 4, at = pretty(c(0,60)))      # Add second axis
mtext("Number of structures per cluster", side = 4, line = 3, cex = 1.5, col = 'black')             # Add second axis label
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
  neuCnt[i] = lengths(regmatches(lig, gregexpr("NeuAc", lig)))
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

pdf(file = paste('./analysis/sfgPlots/', 
                 'clusts_with_ligands_id50',
                 '.pdf', sep = ''),
    width = 18,
    height = 10)
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
                  'High Mannose', 'NeuAc containing', 'Terminal Fucose'),
       col = c(sacc_col,
               'white',
               'navy', 'grey85',
               'white',
               'white', 'white', 'white'),
       pt.bg = c(rep(NA,9),
              'forestgreen', 'darkorchid', 'firebrick1'),
       pt.cex =1.8)
dev.off()

topLigOccurences$High_Mannose[1] = 0
topLigOccurences$Sialic_Acid[1] = 0
topLigOccurences$Terminal_Fucose[1] = 0
for (i in 1:length(clusLst50)){
  if(any(bsResiDat$iupac[bsResiDat$seqClust50 == clusLst50[i]] %in% ligSort50[manTag])){
    topLigOccurences$High_Mannose[1] = topLigOccurences$High_Mannose[1] + 1
  }
  if(any(bsResiDat$iupac[bsResiDat$seqClust50 == clusLst50[i]] %in% ligSort50[neuTag])){
    topLigOccurences$Sialic_Acid[1] = topLigOccurences$Sialic_Acid[1] + 1
  }
  if(any(bsResiDat$iupac[bsResiDat$seqClust50 == clusLst50[i]] %in% ligSort50[fucTag])){
    topLigOccurences$Terminal_Fucose[1] = topLigOccurences$Terminal_Fucose[1] + 1
  }
}

topLigOccurences[1,2:4] = (topLigOccurences[1,2:4] / length(clusLst50) )* 100

# ligSort50[manTag]
# ligSort50[neuTag]
# ligSort50[fucTag]

## 80% id

parenCnt = bracCnt = manCnt = neuCnt = bracCnt = rep(0,length(ligSort80))
for (i in 1:length(ligSort80)){
  lig = ligSort80[i]
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

hist(manCnt)

manTag = manCnt > 3 # High mannose
neuTag = neuCnt >= 1 # Has sialic acid
fucTag = grepl('^Fuc',ligSort80) # Has a terminal fucose

pdf(file = paste('./analysis/sfgPlots/', 
                 'clusts_with_ligands_id80',
                 '.pdf', sep = ''),
    width = 18,
    height = 10)
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
                  'High Mannose', 'NeuAc containing', 'Terminal Fucose'),
       col = c(sacc_col,
               'white',
               'navy', 'grey85',
               'white',
               'white', 'white', 'white'),
       pt.bg = c(rep(NA,9),
                 'forestgreen', 'darkorchid', 'firebrick1'),
       pt.cex =1.8)
dev.off()

topLigOccurences$High_Mannose[2] = 0
topLigOccurences$Sialic_Acid[2] = 0
topLigOccurences$Terminal_Fucose[2] = 0
for (i in 1:length(clusLst80)){
  if(any(bsResiDat$iupac[bsResiDat$seqClust80 == clusLst80[i]] %in% ligSort80[manTag])){
    topLigOccurences$High_Mannose[2] = topLigOccurences$High_Mannose[2] + 1
  }
  if(any(bsResiDat$iupac[bsResiDat$seqClust80 == clusLst80[i]] %in% ligSort80[neuTag])){
    topLigOccurences$Sialic_Acid[2] = topLigOccurences$Sialic_Acid[2] + 1
  }
  if(any(bsResiDat$iupac[bsResiDat$seqClust80 == clusLst80[i]] %in% ligSort80[fucTag])){
    topLigOccurences$Terminal_Fucose[2] = topLigOccurences$Terminal_Fucose[2] + 1
  }
}

topLigOccurences[2,2:4] = (topLigOccurences[2,2:4] / length(clusLst80) )* 100


# ligSort80[manTag]
# ligSort80[neuTag]
# ligSort80[fucTag]

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

pdf(file = paste('./analysis/sfgPlots/', 
                 'clust_pcts_w_ligands',
                 '.pdf', sep = ''),
    width = 15,
    height = 10.5)
ggplot(melt_ligOccur, aes(fill = Clust_ID, x = Ligand, alpha = Ligand, y = Cluster_Percent_w_ligand)) + geom_bar(stat="identity", color="black", position=position_dodge())+
  scale_fill_manual(values=mycol[c(1,4)]) +
  scale_alpha_manual(values=c(rep(0.6,3), rep(1,12)), guide = F) +
  ylim(c(0, 30)) +
  geom_hline(yintercept = c(5)) +
  geom_vline(xintercept = c(3.5), lty = 2) +
  theme_linedraw(base_size = 22) +
  labs(title = "Well-represented ligands in lectin clusters", x = "IUPAC Ligands", y = "Percent of clusters with ligand") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))
dev.off()

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

pc1 = pc2 = 0
runSum = 0
for( i in 1:10){ # How much variance is captured by the top principle components
  v = round(pca$sdev[i]**2/sum(pca$sdev**2), 4) * 100
  if ( i == 1){pc1 = as.character(v)}
  if ( i == 2){pc2 = as.character(v)}
  runSum = runSum + v
  print( paste('PC', as.character(i), ' accounts for ', as.character(v), '% of the total variance', sep = '') )
}
print( paste('Top 10 PCs account for a total of ', as.character(runSum), '% of the total variance', sep = '') )
plot(pca$x, xlab = paste('PC1 (', pc1, '% of variance)', sep =''), ylab = paste('PC2 (', pc2, '% of variance)', sep =''), main = 'PCA of d2Dists (raw bin counts)')



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
pc1 = pc2 = 0
runSum = 0
for( i in 1:10){ # How much variance is captured by the top principle components
  v = round(pca.d2Dists$sdev[i]**2/sum(pca.d2Dists$sdev**2), 4) * 100
  if ( i == 1){pc1 = as.character(v)}
  if ( i == 2){pc2 = as.character(v)}
  runSum = runSum + v
  print( paste('PC', as.character(i), ' accounts for ', as.character(v), '% of the total variance', sep = '') )
}
print( paste('Top 10 PCs account for a total of ', as.character(runSum), '% of the total variance', sep = '') )
plot(pca.d2Dists$x, xlab = paste('PC1 (', pc1, '% of variance)', sep =''), ylab = paste('PC2 (', pc2, '% of variance)', sep =''), main = 'PCA of d2Dists (% comp per pocket)')

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

# Scaled bins s.t. each shape has 40 equal bins of varied distance
pca = prcomp(x = d2ScaledBins, scale. = T)

pc1 = pc2 = 0
runSum = 0
for( i in 1:10){ # How much variance is captured by the top principle components
  v = round(pca$sdev[i]**2/sum(pca$sdev**2), 4) * 100
  if ( i == 1){pc1 = as.character(v)}
  if ( i == 2){pc2 = as.character(v)}
  runSum = runSum + v
  print( paste('PC', as.character(i), ' accounts for ', as.character(v), '% of the total variance', sep = '') )
}
print( paste('Top 10 PCs account for a total of ', as.character(runSum), '% of the total variance', sep = '') )
plot(pca$x, xlab = paste('PC1 (', pc1, '% of variance)', sep =''), ylab = paste('PC2 (', pc2, '% of variance)', sep =''), main = 'PCA of d2Dists in scaled bins')

corrplot(cor(pca$x[,1:10], d2Feats))

thresholds = c('_4Ang', '_6Ang', '_8Ang', '_10Ang')
for (i in 1:nrow(d2ScaledBins)){
  for (t in 1:length(thresholds)){
    tag = grepl(thresholds[t], colnames(d2ScaledBins)) # get columns for the threshold
    if (! all(d2ScaledBins[i,tag] == 0)){ # avoid /0 error
      d2ScaledBins[i,tag] = d2ScaledBins[i,tag] / sum(d2ScaledBins[i,tag])
    }
  }
}

pca.scaleBins = prcomp(x = d2ScaledBins, scale. = T)

pc1 = pc2 = 0
runSum = 0
for( i in 1:10){ # How much variance is captured by the top principle components
  v = round(pca.scaleBins$sdev[i]**2/sum(pca.scaleBins$sdev**2), 4) * 100
  if ( i == 1){pc1 = as.character(v)}
  if ( i == 2){pc2 = as.character(v)}
  runSum = runSum + v
  print( paste('PC', as.character(i), ' accounts for ', as.character(v), '% of the total variance', sep = '') )
}
print( paste('Top 10 PCs account for a total of ', as.character(runSum), '% of the total variance', sep = '') )
plot(pca.scaleBins$x, xlab = paste('PC1 (', pc1, '% of variance)', sep =''), ylab = paste('PC2 (', pc2, '% of variance)', sep =''), main = 'PCA of d2Dists in scaled bins (% per pocket)')

corrplot(cor(pca.scaleBins$x[,1:10], d2Feats))

smoothScatter(pca.scaleBins$x, main = 'PCA of d2ScaledBins (% comp per pocket)', xlab = paste('PC1 (', pc1, '% of variance)', sep =''), ylab = paste('PC2 (', pc2, '% of variance)', sep =''))
pdb = '3LL2'
points(pca.scaleBins$x[grepl(pdb,row.names(d2Dists)), 1], pca.scaleBins$x[grepl(pdb,row.names(d2Dists)), 2], col = alpha('green',0.7), pch=19)
pdb = '2CL8'
points(pca.scaleBins$x[grepl(pdb,row.names(d2Dists)), 1], pca.scaleBins$x[grepl(pdb,row.names(d2Dists)), 2], col = alpha('orange2',0.7), pch=19) # 1 deep medium sized pockets only
pdb = '4URR'
points(pca.scaleBins$x[grepl(pdb,row.names(d2Dists)), 1], pca.scaleBins$x[grepl(pdb,row.names(d2Dists)), 2], col = alpha('red',0.7), pch=19) # deep long pocket
pdb = '1HJV'
points(pca.scaleBins$x[grepl(pdb,row.names(d2Dists)), 1], pca.scaleBins$x[grepl(pdb,row.names(d2Dists)), 2], col = alpha('purple',0.7), pch=19) # one 
pdb = '2WGC'
points(pca.scaleBins$x[grepl(pdb,row.names(d2Dists)), 1], pca.scaleBins$x[grepl(pdb,row.names(d2Dists)), 2], col = alpha('blue',0.7), pch=19) # two 0 pockets only
legend(x ="topright", legend = c('mGRFT', '2CL8', '2WGC', '4URR', '1HJV'), col = c('green','orange2', 'blue', 'red', 'purple'), pch = 19)



# PCA of d2 distribution features
pca = prcomp(x = d2Feats, scale. = T)
runSum = 0
for( i in 1:10){ # How much variance is captured by the top principle components
  v = round(pca$sdev[i]**2/sum(pca$sdev**2), 4) * 100
  runSum = runSum + v
  print( paste('PC', as.character(i), ' accounts for ', as.character(v), '% of the total variance', sep = '') )
}
print( paste('Top 10 PCs account for a total of ', as.character(runSum), '% of the total variance', sep = '') )

plot(pca$x, main = 'PCA of d2Feats', xlab = 'PC1 (51.56% of var.)', ylab = 'PC2 (13.41% of var.)')


###########################
# Explore 3DZDs
###########################

row.names(zern8)[ ! (row.names(zern8) %in% row.names(zern4)) ]
d2Feats[grepl('2WBW',row.names((d2Feats))),] # No volume in smallest pocket threshold (probably below 15 Ang^3)

zern4 = rbind(zern4, rep(0,ncol(zern4))) # Add a row of zeros to Zern4 for the missing pocket at 4Ang threshold
row.names(zern4)[nrow(zern4)] = row.names(zern8)[ ! (row.names(zern8) %in% row.names(zern4)) ] # Name the row of zeros
zern4 = zern4[row.names(zern8),] # Re-order zern4 to match other zerns

all(row.names(zern4) == row.names(zern6)) & all(row.names(zern4) == row.names(zern8)) & all(row.names(zern4) == row.names(zern10))

# Normalize each shape s.t. each descriptor becomes the percentage of the sum of all descriptors for that shape
zern4Norm = t(apply(zern4, 1, invarNorm))
zern6Norm = t(apply(zern6, 1, invarNorm))
zern8Norm = t(apply(zern8, 1, invarNorm))
zern10Norm = t(apply(zern10, 1, invarNorm))

allZernNorm = as.data.frame(cbind(zern4Norm, zern6Norm, zern8Norm, zern10Norm))


# Add rows of zeros for the pockets with no volume
missingRows = row.names(bsResiDat)[!(row.names(bsResiDat) %in% row.names(allZernNorm))]
zeroes2Zern = as.data.frame(matrix(0,nrow = length(missingRows), ncol = ncol(allZernNorm)))
row.names(zeroes2Zern) = missingRows
colnames(zeroes2Zern) = colnames(allZernNorm)

allZernNorm = rbind(allZernNorm, zeroes2Zern)

# All dataframes have the same row names
all(row.names(allZernNorm) %in% row.names(bsResiDat)) & all(row.names(bsResiDat) %in% row.names(allZernNorm))

# order 3DZDs df to match others
allZernNorm = allZernNorm[row.names(bsResiDat),]

# All dataframes have the same row names and are in the same order
all(row.names(bsResiDat) == row.names(allZernNorm)) & all(row.names(bsResiDat) == row.names(d2Feats)) & all(row.names(bsResiDat) == row.names(pca.d2Dists$x))

# Plot the 3DZDs
colors = colfunc(length(unique(bsResiDat$seqClust50)))
labels = bsResiDat[row.names(allZernNorm), 'seqClust50']
labels.u = sort(unique(labels))
labels.int = rep(0,length(labels))
for (i in 1:length(labels.u)){
  labels.int[labels == labels.u[i]] = i
}


pca = prcomp(x = allZernNorm, scale. = T)
pc1 = pc2 = 0
runSum = 0
for( i in 1:10){ # How much variance is captured by the top principle components
  v = round(pca$sdev[i]**2/sum(pca$sdev**2), 4) * 100
  if ( i == 1){pc1 = as.character(v)}
  if ( i == 2){pc2 = as.character(v)}
  runSum = runSum + v
  print( paste('PC', as.character(i), ' accounts for ', as.character(v), '% of the total variance', sep = '') )
}
print( paste('Top 10 PCs account for a total of ', as.character(runSum), '% of the total variance', sep = '') )

plot(pca$x, xlab = paste('PC1 (', pc1, '% of variance)', sep =''), ylab = paste('PC2 (', pc2, '% of variance)', sep =''), main = 'PCA of 3DZDs colored by cluster membership (50% id)', pch = 19, col = alpha(colors[labels.int], 0.6))

corrplot(cor(pca.d2Dists$x[row.names(d2Dists) %in% row.names(allZernNorm),1:10],pca$x[row.names(d2Dists)[row.names(d2Dists) %in% row.names(allZernNorm)],1:10]))


##
zern.umap = umap(allZernNorm)
##


plot(x = zern.umap$layout[,1], y = zern.umap$layout[,2], xlim = c(-6,7), ylim = c(-6,7),
     pch = 19, col = alpha(colors[labels.int], 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = 'UMAP vis for 3DZDs (colored by 50% id clusters)')
pdb = '3LL2'
points(zern.umap$layout[grepl(pdb,row.names(allZernNorm)), 1], zern.umap$layout[grepl(pdb,row.names(allZernNorm)), 2], bg = alpha('green',0.75),col='white', pch=22, cex=3)
pdb = '2CL8'
points(zern.umap$layout[grepl(pdb,row.names(allZernNorm)), 1], zern.umap$layout[grepl(pdb,row.names(allZernNorm)), 2], bg = alpha('orange2',0.75),col='white', pch=22, cex=3) # 1 deep medium sized pockets only
pdb = '4URR'
points(zern.umap$layout[grepl(pdb,row.names(allZernNorm)), 1], zern.umap$layout[grepl(pdb,row.names(allZernNorm)), 2], bg = alpha('red',0.75),col='white', pch=22, cex=3) # deep long pocket
pdb = '1HJV'
points(zern.umap$layout[grepl(pdb,row.names(allZernNorm)), 1], zern.umap$layout[grepl(pdb,row.names(allZernNorm)), 2], bg = alpha('purple',0.75),col='white', pch=22, cex=3) # half large deep, half very shallow 
pdb = '2WGC'
points(zern.umap$layout[grepl(pdb,row.names(allZernNorm)), 1], zern.umap$layout[grepl(pdb,row.names(allZernNorm)), 2], bg = alpha('blue',0.75),col='white', pch=22, cex=3) # 2 zero volume pockets
pdb = '2WBW'
points(zern.umap$layout[grepl(pdb,row.names(allZernNorm)), 1], zern.umap$layout[grepl(pdb,row.names(allZernNorm)), 2], bg = alpha('pink',0.5),col='white', pch=22, cex=3) # 2 zero volume pockets 


# Checking outliers
plot(x = zern.umap$layout[,1], y = zern.umap$layout[,2],#xlim = c(-5.5,7), ylim = c(-6.5,6.5),
     pch = 19, col = alpha(colors[labels.int], 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = 'UMAP vis for 3DZDs (colored by 50% id clusters)')

out1tag = unname((zern.umap$layout[,2] < -7) & (zern.umap$layout[,2] > -15))
zern.umap$layout[out1tag,]
outlierPlots(d2Feats, out1tag)
outlierPlots(allZernNorm, out1tag)


out2tag = (zern.umap$layout[,1] > 15)
zern.umap$layout[out2tag ,]
outlierPlots(d2Feats, out2tag)
outlierPlots(allZernNorm, out2tag)

out3tag = (zern.umap$layout[,1] < -15)
zern.umap$layout[out3tag ,]
outlierPlots(d2Feats, out3tag)
outlierPlots(allZernNorm, out3tag)


## 10th order polynomials?
zern4Norm = t(apply(zern4[,1:36], 1, invarNorm))
zern6Norm = t(apply(zern6[,1:36], 1, invarNorm))
zern8Norm = t(apply(zern8[,1:36], 1, invarNorm))
zern10Norm = t(apply(zern10[,1:36], 1, invarNorm))

zern10thNorm = as.data.frame(cbind(zern4Norm, zern6Norm, zern8Norm, zern10Norm))


# Add rows of zeros for the pockets with no volume
missingRows = row.names(bsResiDat)[!(row.names(bsResiDat) %in% row.names(zern10thNorm))]
zeroes2Zern = as.data.frame(matrix(0,nrow = length(missingRows), ncol = ncol(zern10thNorm)))
row.names(zeroes2Zern) = missingRows
colnames(zeroes2Zern) = colnames(zern10thNorm)

zern10thNorm = rbind(zern10thNorm, zeroes2Zern)

all(row.names(zern10thNorm) %in% row.names(bsResiDat)) & all(row.names(bsResiDat) %in% row.names(zern10thNorm))# All dataframes have the smae row names
zern10thNorm = zern10thNorm[row.names(bsResiDat),]  # order 3DZDs df to match others
all(row.names(bsResiDat) == row.names(zern10thNorm)) & all(row.names(bsResiDat) == row.names(d2Feats)) & all(row.names(bsResiDat) == row.names(pca.d2Dists$x))# All dataframes have the smae row names and are in the same order

pca.z10 = prcomp(x = zern10thNorm, scale. = T)
pc1 = pc2 = 0
runSum = 0
for( i in 1:10){ # How much variance is captured by the top principle components
  v = round(pca.z10$sdev[i]**2/sum(pca.z10$sdev**2), 4) * 100
  if ( i == 1){pc1 = as.character(v)}
  if ( i == 2){pc2 = as.character(v)}
  runSum = runSum + v
  print( paste('PC', as.character(i), ' accounts for ', as.character(v), '% of the total variance', sep = '') )
}
print( paste('Top 10 PCs account for a total of ', as.character(runSum), '% of the total variance', sep = '') )

plot(pca.z10$x, xlab = paste('PC1 (', pc1, '% of variance)', sep =''), ylab = paste('PC2 (', pc2, '% of variance)', sep =''), main = 'PCA of 3DZDs colored by cluster membership (50% id)', pch = 19, col = alpha(colors[labels.int], 0.6))

###
zern.umap = umap(zern10thNorm)
###

# Checking outliers
plot(x = zern.umap$layout[,1], y = zern.umap$layout[,2],#xlim = c(-5.5,7), ylim = c(-6.5,6.5),
     pch = 19, col = alpha(colors[labels.int], 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = 'UMAP vis for 3DZDs (colored by 50% id clusters)')

tag = unname((zern.umap$layout[,1] < -20))
zern.umap$layout[tag,]
outlierPlots(d2Feats, tag)
outlierPlots(allZernNorm, tag)

plot(x = zern.umap$layout[,1], y = zern.umap$layout[,2], xlim = c(-6.5,6.5), ylim = c(-6.5,6.5),
     pch = 19, col = alpha(colors[labels.int], 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = 'UMAP vis for 3DZDs (colored by 50% id clusters)')



# Check reconstructions

# randPnts = row.names(zern10thNorm)[round(runif(n = 7, min = 1, max = nrow(zern10thNorm)))]
randPnts = c("4Z4Z_GLA:D:1","4LK7_04G:D:202","2JDH_TA5:B:1117","4FMH_GAL:X:1","1BOS_GAL:R:4370","4X07_GLA:D:1","4E52_GMH:E:1")

randCols = randClusts = rep(0, length(randPnts))
randLabs = rep('', length(randPnts))
for (i in 1: length(randPnts)){
  randClusts[i] = bsResiDat[randPnts[i], 'seqClust50']
  randLabs[i] = paste(randPnts[i], randClusts[i], sep = ', clust# ')
}
for (i in 1:length(labels.u)){
  randCols[randClusts == labels.u[i]] = i
}

plot(x = zern.umap$layout[,1], y = zern.umap$layout[,2], xlim = c(-6.5,6.5), ylim = c(-6.5,6.5),
     pch = 19, col = alpha(colors[labels.int], 0.3),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = 'UMAP vis for 10th order 3DZDs') # 1100x800
par(new = T)
plot(zern.umap$layout[randPnts,], xlim = c(-6.5,6.5), ylim = c(-6.5,6.5),
     cex = 1.5, lwd = 3, pch = 21, bg = colors[randCols], col = alpha('black', 0.8),
     axes = F, xlab = '', ylab = '', main = '')
text(x = zern.umap$layout[randPnts,1], y = zern.umap$layout[randPnts,2], labels = randLabs, cex = 1.3, font = 2, 
     pos = c(3, 3, 3, 3, 1, 3, 3))

##
z20.umap = umap(allZernNorm)
##
plot(z20.umap$layout, xlim = c(-7,6), ylim = c(-7,6),
     pch = 19, col = alpha(colors[labels.int], 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = 'UMAP vis for 3DZDs (colored by 50% id clusters)')


plot(z20.umap$layout, xlim = c(-7,6), ylim = c(-7,6),
     pch = 19, col = alpha(colors[labels.int], 0.3),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = 'UMAP vis for 20th order 3DZDs')
par(new = T)
plot(z20.umap$layout[randPnts,], xlim = c(-7,6), ylim = c(-7,6),
     cex = 1.5, lwd = 3, pch = 21, bg = colors[randCols], col = alpha('black', 0.8),
     axes = F, xlab = '', ylab = '', main = '')
text(x = z20.umap$layout[randPnts,1], y = z20.umap$layout[randPnts,2], labels = randLabs, cex = 1.3, font = 2, 
     pos = c(3, 3, 2, 3, 3, 3, 3))

###########################
# Assemble predictive feature set
###########################

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

# Limit to columns with number of residues, AA-specifc frequencies, secStruct frequencies, and Ca2+ ion count
resiPredFeats = bsResiDat[,grepl('(_bin)', colnames(bsResiDat)) | grepl('(^CA$)', colnames(bsResiDat))]

# Limit to columns with counts of the various interactions
intPredFeats = bsResiDat[,18:28]

# Combine resi features with pocket features for premlinary feature set
all(row.names(resiPredFeats) == row.names(pca.z10)) & all(row.names(resiPredFeats) == row.names(d2Feats)) & all(row.names(resiPredFeats) == row.names(pca.scaleBins$x)) & all(row.names(resiPredFeats) == row.names(intPredFeats))

runSum = rep(0, ncol(zern10thNorm))
sum = 0
for( i in 1:ncol(zern10thNorm)){ # How much variance is captured by the top principle components
  v = pca.z10$sdev[i]**2/sum(pca.z10$sdev**2) * 100
  sum = sum + v
  runSum[i] = sum
}

topZernPCs = pca.z10$x[,runSum <= 80]
colnames(topZernPCs) = gsub('^', 'zern_', colnames(topZernPCs))

runSum = rep(0, ncol(d2ScaledBins))
sum = 0
for( i in 1:ncol(d2ScaledBins)){ # How much variance is captured by the top principle components
  v = pca.scaleBins$sdev[i]**2/sum(pca.scaleBins$sdev**2) * 100
  sum = sum + v
  runSum[i] = sum
}

topBinPCs = pca.scaleBins$x[, runSum <= 80]
colnames(topBinPCs) = gsub('^', 'binnedD2_', colnames(topBinPCs))


predFeats = cbind(intPredFeats, resiPredFeats, d2Feats, topBinPCs, topZernPCs) # Preliminary predictive feature set

###########################
# Explore predictive feature set
###########################

colors = colfunc(length(unique(bsResiDat$seqClust50)))
labels = bsResiDat$seqClust50
labels.u = sort(unique(labels))
labels.int = rep(0,length(labels))
for (i in 1:length(labels.u)){
  labels.int[labels == labels.u[i]] = i
}

# initial UMAP
umap.rawfeats = umap(predFeats)


plot(umap.rawfeats$layout, #xlim = c(-5.5,7), ylim = c(-6.5,6.5),
     pch = 19, col = alpha(colors[labels.int], 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = 'UMAP vis for 229 predictive features (colored by 50% id clusters)')

# corrplot(cor(umap.rawfeats$layout, d2Feats))

# Try scaling features between [0,1]
scaledFeats = predFeats
for(i in 1:ncol(scaledFeats)){
  scaledFeats[,i] = (scaledFeats[,i] - min(scaledFeats[,i])) / (max(scaledFeats[,i]) - min(scaledFeats[,i]))
}

umap.scaledfeats = umap(scaledFeats)

pdf(file = paste('./analysis/sfgPlots/', 
                 'predFeats_UMAP_50',
                 '.pdf', sep = ''),
    width = 11,
    height = 11.25)
plot(umap.scaledfeats$layout, #xlim = c(-5.5,7), ylim = c(-6.5,6.5),
     pch = 19, col = alpha(colors[labels.int], 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = 'UMAP vis for 221 predictive features - Scaled (colored by 50% id clusters)')
dev.off()

colors = colfunc(length(unique(bsResiDat$seqClust80)))
labels = bsResiDat$seqClust80
labels.u = sort(unique(labels))
labels.int = rep(0,length(labels))
for (i in 1:length(labels.u)){
  labels.int[labels == labels.u[i]] = i
}

pdf(file = paste('./analysis/sfgPlots/', 
                 'predFeats_UMAP_80',
                 '.pdf', sep = ''),
    width = 11,
    height = 11.25)
plot(umap.scaledfeats$layout, #xlim = c(-5.5,7), ylim = c(-6.5,6.5),
     pch = 19, col = alpha(colors[labels.int], 0.6),
     xlab = 'UMAP 1', ylab = 'UMAP 2', main = 'UMAP vis for 221 predictive features - Scaled (colored by 80% id clusters)')
dev.off()

### Variance within clusters vs between clusters

# bwClus = cor(t(scaledFeats), method = 'spearman')
bwClus = distance(x = scaledFeats) # Precompute all pairwise distances between binding sites

inClus50Cor = rep(-1, length(clusLst50))
inClus50All = rep(-1, ((nrow(bsResiDat)**2)/2))
outClus50All = rep(-1, ((nrow(bsResiDat)**2)/2))
outClus50Mean = rep(-1, length(clusLst50))
ind1 = 1
ind2 = 1
for (i in 1:length(clusLst50)) {
  tag = bsResiDat$seqClust50 == clusLst50[i]
  if (sum(tag) > 1) { # Check correlations between cluster members if cluster has more than one member
    # corMat = cor(t(scaledFeats[bsResiDat$seqClust50 == clusLst50[i],]), method = 'spearman')
    # corMat = distance(x = scaledFeats[bsResiDat$seqClust50 == clusLst50[i], ])
    corMat = bwClus[tag,tag][upper.tri(bwClus[tag,tag])]
    inClus50Cor[i] = mean(corMat)
    inClus50All[ind1:(ind1+length(corMat)-1)] = corMat
    ind1 = ind1+length(corMat)
  }
  inds = (1:nrow(bsResiDat))[tag]
  outLst = matrix(0, ncol = (nrow(bsResiDat)-length(inds)), nrow = length(inds))
  for(j in 1:length(inds)){
    outLst[j,] = bwClus[inds[j],!tag]
  }
  outLst = as.vector(outLst)
  outClus50Mean[i] = mean(outLst)
  outClus50All[ind2:(ind2+length(outLst)-1)] =  outLst
  ind2 = ind2 + length(outLst)
}
inClus50Cor = inClus50Cor[inClus50Cor != -1] # drop values from clusters with only one member
inClus50All = inClus50All[inClus50All != -1] # drop unchanged values from initialized vector for all pairwise distances
outClus50All = outClus50All[outClus50All != -1]


inClus80Cor = rep(-1, length(clusLst80))
inClus80All = rep(-1, ((nrow(bsResiDat)**2)/2))
outClus80All = rep(-1, ((nrow(bsResiDat)**2)/2))
outClus80Mean = rep(-1, length(clusLst80))
ind1 = 1
ind2 = 1
for (i in 1:length(clusLst80)) {
  tag = bsResiDat$seqClust80 == clusLst80[i]
  if (sum(tag) > 1) { # Check correlations between cluster members if cluster has more than one member
    # corMat = cor(t(scaledFeats[bsResiDat$seqClust80 == clusLst80[i],]), method = 'spearman')
    # corMat = distance(x = scaledFeats[bsResiDat$seqClust80 == clusLst80[i], ])
    corMat = bwClus[tag,tag][upper.tri(bwClus[tag,tag])]
    inClus80Cor[i] = mean(corMat)
    inClus80All[ind1:(ind1+length(corMat)-1)] = corMat
    ind1 = ind1+length(corMat)
  }
  inds = (1:nrow(bsResiDat))[tag]
  outLst = matrix(0, ncol = (nrow(bsResiDat)-length(inds)), nrow = length(inds))
  for(j in 1:length(inds)){
    outLst[j,] = bwClus[inds[j],!tag]
  }
  outLst = as.vector(outLst)
  outClus80Mean[i] = mean(outLst)
  outClus80All[ind2:(ind2+length(outLst)-1)] =  outLst
  ind2 = ind2 + length(outLst)
}
inClus80Cor = inClus80Cor[inClus80Cor != -1] # drop values from clusters with only one member
inClus80All = inClus80All[inClus80All != -1] # drop unchanged values from initialized vector for all pairwise distances
outClus80All = outClus80All[outClus80All != -1]


# Plot inter/intra cluster distance distributions
bwClusPlot = bwClus[upper.tri(bwClus)]

plot(density(inClus50Cor))
plot(density(inClus50All))

plot(density(outClus50Mean))
plot(density(outClus50All))

plot(density(inClus80Cor))
plot(density(inClus80All))

plot(density(outClus80Mean))
plot(density(outClus80All))

plot(density(bwClusPlot), xlab = 'Euclidean distance', main = 'Distribution of all distances between all binding sites')
allMed = median(bwClusPlot)
abline(v = allMed, lty = 2, lwd = 3, col = 'red2')
text(x = allMed + 2, y = 1,
     labels = paste('Median = ', as.character(allMed), sep = ''),
     cex = 1.5,
     col = 'red1')


mclusdf = rbind(data.frame('type' = c('All binding sites'), 'seqID' = c('NA'),'value' = bwClusPlot),
                data.frame('type' = c('Cluster means (shared clusts)'), 'seqID' = c('80'), 'value' = inClus80Cor),
                data.frame('type' = c('All w/in shared clusts'), 'seqID' = c('80'), 'value' = inClus80All),
                data.frame('type' = c('Cluster means (shared clusts)'), 'seqID' = c('50'), 'value' = inClus50Cor),
                data.frame('type' = c('All w/in shared clusts'), 'seqID' = c('50'), 'value' = inClus50All),
                data.frame('type' = c('All (no shared clusts)'), 'seqID' = c('50'), 'value' = outClus50All),
                data.frame('type' = c('All (no shared clusts)'), 'seqID' = c('80'), 'value' = outClus80All),
                data.frame('type' = c('Cluster means (no shared clusts)'), 'seqID' = c('50'), 'value' = outClus50Mean),
                data.frame('type' = c('Cluster means (no shared clusts)'), 'seqID' = c('80'), 'value' = outClus80Mean))

pdf(file = paste('./analysis/sfgPlots/', 
                 'interclust_distances',
                 '.pdf', sep = ''),
    width = 18,
    height = 11.3)
ggplot(data = mclusdf, aes(x = type, y = value, col = seqID, fill = seqID)) + geom_boxplot(notch = T, outlier.alpha = 0.5) +
  scale_fill_manual(values = rep('grey80',3), guide = F) + 
  scale_color_manual(values = c('black', mycol[1], mycol[4])) +
  labs(title = 'Inter- & intracluster distances (from predictive features)', x = "Pairwise grouping", y = "Euclidean distance") +
  theme_dark(base_size = 22) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), title = element_text(face = "bold.italic", color = "black"))
dev.off()

# # Cosine similarity based on one-hot encoding of binned BS residues
# oneHotTag = grepl('^[[:upper:]]{3}_', colnames(predFeats))

  
# Sampling diversity w/in binding sites from the same sequence clusters
set.seed(27)

sampCnt = rep(0,length(clusLst50))
totCnt = rep(0, length(clusLst50))
for (i in 1:length(clusLst50)){
  tag = bsResiDat$seqClust50 == clusLst50[i]
  totCnt[i] = sum(tag)
  if (sum(tag) == 1){
    sampCnt[i] = 1
  } else if (sum(tag) ==2){
    if(distance(scaledFeats[tag,]) > allMed){
      sampCnt[i] = 2
    } else{
      sampCnt[i] = 1
    }
  } else {
    cnter = 1
    
    allInds = (1:nrow(bsResiDat))[tag] # Get the row indices for binding sites within a seq cluster
    sampInds = rep(0, length(allInds)) # Initialize a vector of zeros to hold the randomly sampled indices
    
    # tmp = prcomp(scaledFeats[allInds,], center = T)
    # tmpCol = colfunc(10)
    # plot(tmp$x[,1:2], main = "PCA of 113 binding sites in cluster 15", pch = 19, col = 'black', cex = 1.5,
    #      xlim = c(-1.6,1.5), ylim = c(-2,0.75))
    
    while( length(allInds) >= 2 ){
      sampInds[cnter] = sample(allInds, size = 1) # Sample an index randomly
      
      # points(x = tmp$x[row.names(scaledFeats)[sampInds[cnter]],1], y = tmp$x[row.names(scaledFeats)[sampInds[cnter]],2],
      #        pch = 21, bg = tmpCol[cnter], col = 'black' , cex = 2)
      
      cnter = cnter + 1 # increase the sample count
      allInds = allInds[! allInds %in% sampInds] # Drop sampled indices from vector of remaining indices
      
      distMat = distance(x = rbind(scaledFeats[sampInds[sampInds != 0],], scaledFeats[allInds,])) # Find all pairwise distances b/w binding sites in cluster
      distMat = distMat[1:sum(sampInds != 0),-1*(1:sum(sampInds != 0))] # Get the top n rows of the distance matrix dropping the first n columns, where n = the number of indices already sampled
      if (is.matrix(distMat)){
        distMat = apply(X = distMat, MARGIN = 2, FUN = min) # For each of the remaining binding sites, take the minimum pairwise distance to any sampled binding site
      }
      dropInds = allInds[distMat < allMed ]
      
      # if (! any(is.na(dropInds))){
      #   points(x = tmp$x[row.names(scaledFeats)[dropInds],1], y = tmp$x[row.names(scaledFeats)[dropInds],2],
      #          pch = 21, col = alpha(tmpCol[cnter-1], .8), bg = alpha(tmpCol[cnter-1], .4), cex = 1.8)
      # }
      
      allInds = allInds[! allInds %in% dropInds]
      
      if(length(allInds) == 1){
        sampInds[cnter] = allInds
        # points(x = tmp$x[row.names(scaledFeats)[sampInds[cnter]],1], y = tmp$x[row.names(scaledFeats)[sampInds[cnter]],2],
        #        pch = 21, bg = tmpCol[cnter], col = 'black' , cex = 2)
      }
      
    }
    # legend(x = 'bottomright', legend = c(1:length(tmpCol)), col = tmpCol, pch = 19, pt.cex = 2)
    
    sampInds = sampInds[sampInds != 0]
    sampCnt[i] = length(sampInds)
  }
}


plot(density(sampCnt))
plot(density(sampCnt/totCnt))

plot(jitter(totCnt, amount = 0.5), sampCnt, xlab = 'Number of binding sites in a seq cluster', ylab = 'Number of sampled binding sites', col = colfunc(3)[1],
     xlim = c(0,200), ylim = c(0,10))
par(new= T)
plot(jitter(totCnt, amount = 0.5), bk1, col = colfunc(3)[2],
     xlim = c(0,200), ylim = c(0,10),
     axes = F, xlab = '', ylab = '', main = '')
par(new= T)
plot(jitter(totCnt, amount = 0.5), bk2, col = colfunc(3)[3],
     xlim = c(0,200), ylim = c(0,10),
     axes = F, xlab = '', ylab = '', main = '')

###########################
# Descriptive analysis
###########################

# Ligand tags
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
neuTag = neuCnt >= 1 # Has sialic acid
fucTag = grepl('^Fuc',uniLigs) # Has a terminal fucose


# get tags to indicate binding sites containing one of the 15 ligands of interest
ligTags = as.data.frame(matrix(F, nrow = nrow(bsResiDat), ncol = 3+length(top50)))
row.names(ligTags) =  row.names(bsResiDat)
colnames(ligTags) = colnames(topLigOccurences)[2:ncol(topLigOccurences)] 

ligTags$High_Mannose = bsResiDat$iupac %in% uniLigs[manTag]
ligTags$Sialic_Acid = bsResiDat$iupac %in% uniLigs[neuTag]
ligTags$Terminal_Fucose = bsResiDat$iupac %in% uniLigs[fucTag]


for (i in 1:length(top50)){
  lig = top50[i]
  lig = gsub('\\(', '\\\\(', lig)
  lig = gsub('\\)', '\\\\)', lig)
  lig = gsub('\\[', '\\\\]', lig)
  lig = gsub('\\[', '\\\\]', lig)
  lig = gsub('^', '\\^', lig)
  lig = gsub('$', '\\$', lig)
  
  colInd = grep(lig, colnames(ligTags))
  ligTags[,colInd] = grepl(lig, bsResiDat$iupac)
}

# Add tags to bsResiDat
# bsResiDat = cbind(bsResiDat, ligTags)



# Statistical test for cluster independent feature enrichment by ligand
stats = as.data.frame(matrix(0, nrow = ncol(predFeats), ncol = (ncol(ligTags)*4)))
row.names(stats) = colnames(predFeats)
colnames(stats) = c(paste(colnames(ligTags), 'p', sep = '_'), paste(colnames(ligTags), '_FC', sep = '_'), paste(colnames(ligTags), 'effectSize', sep = '_'), paste(colnames(ligTags), 'adj', sep = '_'))

scaled_stats = stats # In parallel, try this with scaled features for more robust change metric, shouldn't impact WMW test much

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



# Plot volcano plots from scaled features

dev.off()
pdf(file = paste('./analysis/sfgPlots/', 
                 'medFC_volcanoes',
                 '.pdf', sep = ''),
    width = 21,
    height = 13)
par(mfrow=c(3,5))
xLim = c(-10,10)
# yLim = c(0,10)
for(i in 1:ncol(ligTags)){
  
  yLim = c(0,max(-log10(0.1), -log10(min(scaled_stats[,grepl('_adj$', colnames(scaled_stats))][,i]))) + 1)
  
  tag = scaled_stats[,grepl('_adj$', colnames(scaled_stats))][,i] < 0.1
  
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
  
  plot(log2(scaled_stats[,grepl('_FC$', colnames(scaled_stats))][,i]), -log10(scaled_stats[,grepl('_adj$', colnames(scaled_stats))][,i]), # Plot all points w/ color @ alpha 0.5
       xlab = "log2(median FC)", ylab = "-log10(FDR)", main = colnames(ligTags)[i],
       pch = 19, cex = 2, col = alpha(featColors, 0.5),
       cex.axis = 1.5, cex.main = 2, cex.lab = 1.5,
       xlim = xLim, ylim = yLim)
  par(new=T)
  plot(log2(scaled_stats[,grepl('_FC$', colnames(scaled_stats))][tag,i]), -log10(scaled_stats[,grepl('_adj$', colnames(scaled_stats))][tag,i]), # Plot stat sig points again with alpha 1
       pch = 19, col = featColors[tag], cex = 2,
       axes = F, xlab = "", ylab = "", main = "",
       xlim = xLim, ylim = yLim)
  par(new=T)
  plot(log2(scaled_stats[,grepl('_FC$', colnames(scaled_stats))][tag,i]), -log10(scaled_stats[,grepl('_adj$', colnames(scaled_stats))][tag,i]), # Outline stat sig points in black to highlight them
       col = 'black', cex = 2.05,
       xlab = "", ylab = "", axes = F, main = "",
       xlim = xLim, ylim = yLim)
  abline(h= -log10(0.1), lwd = 2)
  abline(v = 0, lty=2, lwd = 2)
}
dev.off()

# Correlation between features for each ligand
pdf(file = paste('./analysis/sfgPlots/', 
                 'allFeat_MWM_corplot',
                 '.pdf', sep = ''),
    width = 11.5,
    height = 8.25)
breakLst = seq(-1,1,0.05)
pheatmap(cor(stats[,grepl('_effectSize$', colnames(stats))], scaled_stats[,grepl('_effectSize$', colnames(scaled_stats))], method = 'spearman'),
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

###########################
# Save out files for prediction
###########################

write.table(bsResiDat, file = './analysis/training/data_in/bsResiDat.tsv', quote = F, sep = '\t')
write.table(predFeats, file = './analysis/training/data_in/predFeats.csv', quote = F, sep = ',')
write.table(ligTags, file = './analysis/training/data_in/ligTags.tsv', quote = F, sep = '\t')

# test = read.delim(file = './analysis/training/data_in/ligTags.tsv', sep = '\t', stringsAsFactors = F)
# all(test == ligTags)
# test = read.delim(file = './analysis/training/data_in/predFeats.csv', sep = ',', stringsAsFactors = F)
# all(round(test,4) == round(predFeats,4))
# test = read.delim(file = './analysis/training/data_in/bsResiDat.tsv', sep = '\t', stringsAsFactors = F)

