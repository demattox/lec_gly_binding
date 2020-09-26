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

setwd("/Users/dmattox/cbk/glycan_binding/analysis/prelim/prelim3/")

################################
# Look at pockets only
featDat = read.delim(file = './pairwiseDistFeatures.csv', sep = ',', header = T, stringsAsFactors = F)

row.names(featDat) = featDat$bsite
featDat$bsite <- NULL

dropColumns = rep(F, ncol(featDat)) 
for (i in 1:ncol(featDat)){
  if (sum(featDat[,i]) == 0){
    dropColumns[i] = T
  }
}
featDat = featDat[,!dropColumns]

featDat = featDat[order(row.names(featDat)),]

dropRows = c('1NWS_NAG:C:-5', '1NWS_NAG:D:-5', '1HJV_NAG:C:-5')

featDatScaled = featDat[! (row.names(featDat) %in% dropRows), ]

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



mydist=function(c) {dist(c,method="euclidean")}
myclust=function(c) {hclust(c,method="average")}
mycol <- colorRampPalette(c("ivory", "cornflowerblue", "navy"))(n = 299)

pdf(file = '/Users/dmattox/Documents/qbs/cbk/glycan_binding/analysis/prelim/prelim3/prelim_heatmap_labels_dropBugs.pdf')

heatmap.3(featDatScaled, hclustfun = myclust, distfun = mydist,
          main = "All Lectin Binding Sites", col = mycol,labCol = gsub(pattern = "Ang", "", colnames(featDatScaled)), labRow = F
          , dendrogram = "both", KeyValueName="Scaled Feature Values"
          , scale = 'none', density.info = 'none',margins=c(6,12), Rowv = T, Colv = T)
dev.off()



pca = prcomp(x = featDat, scale. = T)
plot(pca$x)

out = pca$x[,1:36]
labs = row.names(featDat)
write.csv(file = 'distFeats_pc50.csv', out, quote = F, row.names = F)
write.table(file = 'distFeats_pc50_labels.tsv', labs, sep = '\t', quote = F, row.names = F)










################################

dat = read.delim(file = '/Users/dmattox/cbk/glycan_binding/analysis/prelim/prelim3/prelimData3.tsv', sep = '\t', header = T, stringsAsFactors = F)
length(unique(dat$pdb))
length(unique(dat$uniparc))

missingInteractions = rep(F, nrow(dat))

for (i in 1:nrow(dat)) {
  if ((dat$total[i] == 0)
      || (dat$numBSresis_bin1 == 0)) {
    missingInteractions[i] = T
  }
}

dat[missingInteractions,]

# cnt = 0
# for (i in 1:length(unique(dat$pdb))){
#   pdbID = unique(dat$pdb)[i]
#   tag = dat$pdb == pdbID
#   if (length(unique(dat$longname[tag])) > 1){
#     cnt = cnt + 1
#   }
# }
# cnt

# dat = dat[!missingInteractions,]
# length(unique(dat$pdb))
# length(unique(dat$uniparc))

row.names(dat) = dat$bSite
dat$bSite <- NULL

dat = dat[order(row.names(dat)),]

all(row.names(dat) == row.names(featDat))

for (i in 1:ncol(dat)){
  if(all(is.na(dat[,i]))){
    print(i)
  }
  else if (all(dat[,i] == 0) == TRUE){
    print(colnames(dat)[i])
    print(i)
  }
}
print(grep('I_bin4',colnames(dat))) # pi helix in bin 4 has a single binding site with any pi helix resiudes

bsDat = dat[,c(11, 13:31, 33:44, 46:57, 59:79 )] # filter out non-binding site columns and columns 32, 45, & 58 for having 0 recorded pi helices in bins 1-3

pca = prcomp(x = bsDat, scale. = T)
plot(pca$x)

out = pca$x[,1:50]
labs = cbind(row.names(bsDat),dat[,1:12])
write.csv(file = 'resiFeats_pc50.csv', out, quote = F, row.names = F)
write.table(file = 'resiFeats_pc50_labels.tsv', labs, sep = '\t', quote = F, row.names = F)




allDat = merge(x = bsDat, y= featDat, by =0, all=TRUE)
row.names(allDat) = allDat$Row.names
allDat$Row.names <- NULL 

pca = prcomp(x = allDat, scale. = T)
plot(pca$x)

out = pca$x[,1:50]
labs = dat[,1:12]
write.csv(file = 'allFeats_pc50.csv', out, quote = F, row.names = F)
write.table(file = 'allFeats_pc50_labels.tsv', labs, sep = '\t', quote = F, row.names = F)

########################################
zernFeats = read.delim(file = '/Users/dmattox/cbk/glycan_binding/analysis/prelim/prelim3/3DZD_test.csv', sep = ',', header = T, stringsAsFactors = F)

row.names(zernFeats) = zernFeats$bsite
zernFeats$bsite <- NULL

distMat = distance(x = zernFeats, method = 'cosine',use.row.names = TRUE)
heatmap(distMat, scale = 'none', col = colorRampPalette(brewer.pal(8, "Blues"))(25))

# distMat = distance(x = zernFeats, method = 'euclidean',use.row.names = TRUE)
# heatmap(distMat, scale = 'none', col = rev(colorRampPalette(brewer.pal(8, "Blues"))(25)))

zern20 = read.delim(file = '/Users/dmattox/cbk/glycan_binding/analysis/prelim/prelim3/3DZD_test20.csv', sep = ',', header = T, stringsAsFactors = F)

row.names(zern20) = zern20$bsite
zern20$bsite <- NULL

distMat = distance(x = zern20, method = 'cosine',use.row.names = TRUE)
heatmap(distMat, scale = 'none', col = colorRampPalette(brewer.pal(8, "Blues"))(25))


zern5 = read.delim(file = '/Users/dmattox/cbk/glycan_binding/analysis/prelim/prelim3/3DZD_test5.csv', sep = ',', header = T, stringsAsFactors = F)
row.names(zern5) = zern5$bsite
zern5$bsite <- NULL
distMat = distance(x = zern5, method = 'cosine',use.row.names = TRUE)
heatmap(distMat, scale = 'none', col = colorRampPalette(brewer.pal(8, "Blues"))(25))



moments = c('F00_','F11_') # 1st order
moments = gsub(pattern = '4Ang', replacement = '', x = colnames(zern20)[1:12]) # 5th order
moments = gsub(pattern = '4Ang', replacement = '', x = colnames(zern20)[1:30]) # 9th order
moments = gsub(pattern = '4Ang', replacement = '', x = colnames(zern20)[1:49]) # 12th order
moments = gsub(pattern = '4Ang', replacement = '', x = colnames(zern20)[1:81]) # 16th order
moments = gsub(pattern = '4Ang', replacement = '', x = colnames(zern20)[1:121]) # 20th order


slicedZern = zern20[,grepl(paste(moments,collapse="|"), colnames(zern20))]
distMat = distance(x = slicedZern, method = 'cosine',use.row.names = TRUE)
heatmap(distMat, scale = 'none', col = colorRampPalette(brewer.pal(8, "Blues"))(25))














