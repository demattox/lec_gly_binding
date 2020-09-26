rm(list = ls())


library(reshape2)
library(ggplot2)



dat = read.delim(file = '/Users/dmattox/cbk/glycan_binding/analysis/prelimData.tsv', sep = '\t', header = T, stringsAsFactors = F)

rawInts = read.delim(file = '/Users/dmattox/cbk/allGlycans/pipeline/Clean_Lectin_intensities_67uni.txt', sep = '\t', header = T, stringsAsFactors = F)

comps = as.data.frame(matrix(nrow= nrow(rawInts), ncol = 2))
colnames(comps) = c('gly', 'link')

pat = "(^.*)-(Sp[[:digit:]]+$)"
for (i in 1:nrow(rawInts)){
  combo = rawInts$glyName[i]
  comps[i,1] = gsub(pat, '\\1', combo) # glycan
  comps[i,2] = gsub(pat, '\\2', combo) # linker
}

# Unique glycans
length(unique(comps[,1])) # 553 glycans from 610 array components

n_occur = as.data.frame(table(comps$gly), stringsAsFactors = F)
sum(n_occur$Freq > 1) # 53 unique glycans that appear multiple times

n_occur$Var1[n_occur$Freq >1]


