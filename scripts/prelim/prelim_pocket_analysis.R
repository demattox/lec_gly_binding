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

setwd("/Users/dmattox/cbk/glycan_binding/analysis/")

dat = read.delim(file = './prelim/prelim2/testPairwiseBins.csv', sep = ',', header = T, stringsAsFactors = F)

row.names(dat) = dat$bsite
dat$bsite <- NULL

dropColumns = rep(F, ncol(dat)) 

for (i in 1:ncol(dat)){
  if (sum(dat[,i]) == 0){
    dropColumns[i] = T
  }
}

dat = dat[,!dropColumns]

drop = c('1NWS_NAG:D:-5', '1NWS_NAG:C:-5', '1HJV_NAG:D:-5', '1NWS_NAG:B:-5')
dat = dat[!(row.names(dat) %in% drop), ]

distMat = distance(x = dat, method = 'cosine',use.row.names = TRUE)

heatmap(distMat, scale = 'none', col = colorRampPalette(brewer.pal(8, "Blues"))(25))

percentDat = dat
threshes = c('_4Ang', '_6Ang','_8Ang','_10Ang')
for (i in 1:length(threshes)){
  tag = grepl(pattern = threshes[i], colnames(percentDat))
  for (j in 1:nrow(percentDat)){
    tot = sum(percentDat[j,tag])
    percentDat[j,tag] = percentDat[j,tag]/tot
  }
}


######################################################
pca = prcomp(x = percentDat, scale. = T)
plot(pca$x)

pc1var = round(pca$sdev[1]**2/sum(pca$sdev**2), 4) * 100
pc2var = round(pca$sdev[2]**2/sum(pca$sdev**2), 4) * 100

out = pca$x[,1:50]
labs = row.names(dat)
write.csv(file = './prelim/prelim2/pocket_pc50.csv', out, quote = F, row.names = F, col.names = F)
write.table(file = './prelim/prelim2/pocket_pc50_labels.tsv', labs, sep = '\t', quote = F, row.names = F)

mydist=function(c) {dist(c,method="euclidean")}
myclust=function(c) {hclust(c,method="average")}
mycol <- colorRampPalette(c("ivory", "cornflowerblue", "navy"))(n = 299)

pdf(file = '/Users/dmattox/Documents/qbs/cbk/glycan_binding/analysis/prelim/prelim2/prelim_pocket_heatmap.pdf')

heatmap.3(percentDat, hclustfun = myclust, distfun = mydist,
          main = "All Pairwise Dists b/w surface pocket points", col = mycol,labCol = F, labRow = F
          , dendrogram = "both",  ColSideColorsSize = 2, KeyValueName="Scaled Feature Values"
          , scale = 'none', density.info = 'none',margins=c(6,12), Rowv = T, Colv = T)

dev.off()

#############################################
allDat = read.delim(file = './prelim/prelim2/testRandPairwiseBins.csv', sep = ',', header = T, stringsAsFactors = F)
clDat = read.delim(file = './prelim/prelim2/testRandCLbins.csv', sep = ',', header = T, stringsAsFactors = F)

row.names(allDat) = allDat$bsite
allDat$bsite <- NULL

row.names(clDat) = clDat$bsite
clDat$bsite <- NULL

dropColumns = rep(F, ncol(allDat)) 
for (i in 1:ncol(allDat)){
  if (sum(allDat[,i]) == 0){
    dropColumns[i] = T
  }
}
allDat = allDat[,!dropColumns]

dropColumns = rep(F, ncol(clDat)) 
for (i in 1:ncol(clDat)){
  if (sum(clDat[,i]) == 0){
    dropColumns[i] = T
  }
}
clDat = clDat[,!dropColumns]


# drop = c('1NWS_NAG:D:-5', '1NWS_NAG:C:-5', '1HJV_NAG:D:-5', '1NWS_NAG:B:-5')
# dat = dat[!(row.names(dat) %in% drop), ]

allDistMat = distance(x = allDat, method = 'cosine',use.row.names = TRUE)
heatmap(allDistMat, scale = 'none', col = colorRampPalette(brewer.pal(8, "Blues"))(25))

clDistMat = distance(x = clDat, method = 'cosine',use.row.names = TRUE)
heatmap(clDistMat, scale = 'none', col = colorRampPalette(brewer.pal(8, "Blues"))(25))


####################################################
binDat = read.delim(file = '/Users/dmattox/Documents/qbs/cbk/glycan_binding/analysis/prelim/prelim2/test2RandPairwiseBins.csv', sep = ',', header = T, stringsAsFactors = F)
featDat = read.delim(file = '/Users/dmattox/Documents/qbs/cbk/glycan_binding/analysis/prelim/prelim2/test2RandPairwiseFeatures.csv', sep = ',', header = T, stringsAsFactors = F)

row.names(binDat) = binDat$bsite
binDat$bsite <- NULL
row.names(featDat) = featDat$bsite
featDat$bsite <- NULL

dropColumns = rep(F, ncol(binDat)) 
for (i in 1:ncol(binDat)){
  if (sum(binDat[,i]) == 0){
    dropColumns[i] = T
  }
}
binDat = binDat[,!dropColumns]

dropColumns = rep(F, ncol(featDat)) 
for (i in 1:ncol(featDat)){
  if (sum(featDat[,i]) == 0){
    dropColumns[i] = T
  }
}
featDat = featDat[,!dropColumns]

binDistMat = distance(x = binDat, method = 'cosine',use.row.names = TRUE)
heatmap(binDistMat, scale = 'none', col = colorRampPalette(brewer.pal(8, "Blues"))(25))

featDistMat =  distance(x = featDat, method = 'cosine',use.row.names = TRUE)
heatmap(featDistMat, scale = 'none', col = colorRampPalette(brewer.pal(8, "Blues"))(25))


# Split skew features with negative values into two columns and normalize between 0 and 1 for euclidean distance


featDatScaled = featDat

tag4 = featDatScaled$skew_4Ang < 0
tag6 = featDatScaled$skew_6Ang < 0
tag8 = featDatScaled$skew_8Ang < 0
tag10 = featDatScaled$skew_10Ang < 0

featDatScaled$leftskew_4Ang = featDatScaled$skew_4Ang
featDatScaled$leftskew_4Ang[!tag4] = 0
featDatScaled$leftskew_4Ang = featDatScaled$leftskew_4Ang*-1
featDatScaled$skew_4Ang[tag4] = 0

featDatScaled$leftskew_6Ang = featDatScaled$skew_6Ang
featDatScaled$leftskew_6Ang[!tag6] = 0
featDatScaled$leftskew_6Ang = featDatScaled$leftskew_6Ang*-1
featDatScaled$skew_6Ang[tag6] = 0

featDatScaled$leftskew_8Ang = featDatScaled$skew_8Ang
featDatScaled$leftskew_8Ang[!tag8] = 0
featDatScaled$leftskew_8Ang = featDatScaled$leftskew_8Ang*-1
featDatScaled$skew_8Ang[tag8] = 0

featDatScaled$leftskew_10Ang = featDatScaled$skew_10Ang
featDatScaled$leftskew_10Ang[!tag10] = 0
featDatScaled$leftskew_10Ang = featDatScaled$leftskew_10Ang*-1
featDatScaled$skew_10Ang[tag10] = 0

for(i in 1:ncol(featDatScaled)){
  featDatScaled[,i] = featDatScaled[,i] / max(featDatScaled[,i])
}


featDistMat =  distance(x = featDat, method = 'cosine',use.row.names = TRUE)
heatmap(featDistMat, scale = 'none', col = colorRampPalette(brewer.pal(8, "Blues"))(25))
featDistMat =  distance(x = featDat, method = 'euclidean',use.row.names = TRUE)
heatmap(featDistMat, scale = 'none', col = rev(colorRampPalette(brewer.pal(8, "Blues"))(25)))

