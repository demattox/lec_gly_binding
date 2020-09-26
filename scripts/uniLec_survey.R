rm (list = ls())

library(ggplot2)
library(corrplot)
library(reshape2)
library(seqinr)
library(stringr)

gbDir = "/Users/dmattox/cbk/lec_gly_binding/"

dat = read.delim(paste(gbDir,"data/unilectin/UnilLectin_data.csv", sep = ""), sep=';', header = T, stringsAsFactors = F)

cfg = read.delim(paste(gbDir,"data/unilectin/CFG-array-data.csv", sep = ""), sep=';', header = T, stringsAsFactors = F)

dat$cfg = 0
for (i in 1:nrow(dat)){
  if (dat$pdb[i] %in% cfg$pdb){
    dat$cfg[i] = 1
  }
}

lecIDs = unique(dat$uniparc)
length(lecIDs)
sum(dat$cfg)

sum(dat$iupac != "")
sum(dat$ligand != "")
cnt = 0
for (i in 1:length(lecIDs)){
  lec = lecIDs[i]
  tag = dat$uniparc == lec
  if (any(dat$iupac[tag]!= "")){
    cnt = cnt + 1
  }
}
cnt # Number of uniProt IDs with at least one PDB file containing a glycan

# Cases with nothing in ligand column but entries in sucres/iupac columns

tag = dat$ligand == ""
sum(dat$iupac[tag] != "")
dat$pubmed[tag][dat$iupac[tag] != ""] # Seems to be the cases primarily where no pubmed ID associated


# Clean up PDB IDs

# PDBid was 43983
### From pubmed ID in dataframe, corrected PDB ID should be 1JUI
dat$pdb = gsub("43983", "1JUI", dat$pdb)

# Some PDB IDs autoformatted to exponential annotation
dat$pdb[grepl(pattern = "\\+", x = dat$pdb)]
dat$pdb = gsub("\\+", "", dat$pdb)

cfg$pdb[grepl(pattern = "\\+", x = cfg$pdb)] # Fix in CFG mapping dataframe too
cfg$pdb = gsub("\\+", "", cfg$pdb)

# No PDB IDs, structures achrived by PDB or replaced by more recent structures (publications provided by François in email)
dat$pdb[grepl(pattern = "^ND", x = dat$pdb)]
dat$pdb[grepl(pattern = "Model", x = dat$pdb)]
indsOut = c(grep(pattern = "^ND", x = dat$pdb), grep(pattern = "Model", x = dat$pdb))
excluded = dat[indsOut, ]
dat = dat[-indsOut, ] # Remove hose rows from dataframe, save them in separate dataframe

cfg = cfg[!cfg$pdb %in% excluded$pdb,] # Remove 3 instances from cfg dataframe (if no pdb file, no reason to link to cfg)

# Other non-PDB id formatted values? nchar != 4, 1st char non-numeric
for (i in 1:nrow(dat)){
  if (nchar(dat$pdb[i]) != 4){
    print(dat$pdb[i])
  }
}
dat$pdb[!grepl("^[[:digit:]]",dat$pdb)]
###

# PDB ID 3IYK file only contains alpha carbons from cryoEM at 7 Ang resolution. Drop from analysis
dat = dat[!dat$pdb == '3IYK',]


all((dat$sucres != "") == (dat$iupac != ""))
tag = dat$iupac != ""
complexedPDBs = dat$pdb[tag]
complexedPDBs = sort(complexedPDBs)
write(complexedPDBs, file = paste(gbDir,"data/unilectin/boundPDBs.csv", sep = ""), ncol = length(complexedPDBs), sep = " ") # Space delimited list of PDB IDs to copy/paste into runPlip.pbs


### Data characterization/visualization
dat$resolution = as.numeric(dat$resolution)

#PDB ID 4XTP has resolution of 211.97? RSCB: 1.97 Ang
dat$resolution[dat$pdb == "4XTP"] = 1.97

dat$resolution[dat$resolution == 0] =  NA
# plot(density(dat$resolution, na.rm = T), xlim = c(0,8)) # All structures (with listed resolutions)
jpeg(filename = paste(gbDir,'analysis/QC/structuralResolution.jpg', sep =''))
plot(density(dat$resolution[dat$iupac != ""], na.rm = T), xlab = 'Resolution (Angstoms)', main = 'Resolutions of lectin structures bound to glycans')
dev.off()
# Save out reformatted, updated dataframe with all PDB IDs for structures bound to glycans
boundOut = dat[dat$iupac != "",]
# Update missing Uniprot IDs manually (from rcsb.org page for each structure, sequence tab)
sum(boundOut$uniparc == "")
boundOut$uniparc[boundOut$pdb == '6V1C'] = 'Q07654'

boundOut$uniparc[boundOut$pdb == '5Y97'] = 'U3KRF8'

boundOut$uniparc[boundOut$pdb == '6RYM'] = 'P23805'
boundOut$uniparc[boundOut$pdb == '6RYN'] = 'P23805'

boundOut$uniparc[boundOut$pdb == '6Q40'] = 'F9XHX3'

boundOut$uniparc[boundOut$pdb == '5MZD'] = 'P22747'
boundOut$uniparc[boundOut$pdb == '6R41'] = 'P22747'
boundOut = boundOut[! boundOut$pdb %in% c('5MZB', '5MZ9'),]

boundOut$uniparc[boundOut$pdb == '4MQ0'] = 'P83304'

boundOut$uniparc[boundOut$pdb == '6IS5'] = 'Q66296'
 
boundOut$uniparc[boundOut$pdb == '5XG5'] = 'B3EWR1'

boundOut$uniparc[boundOut$pdb == '6TID'] = 'B4EH86'
boundOut$uniparc[boundOut$pdb == '6TIG'] = 'B4EH86'

length(unique(boundOut$uniparc)) #413 unique lectin UniProt IDs

# for (i in 1:length(unique(boundOut$uniparc))){ # make sure lectin info ~~generally~~  agrees for each uniprot tag
#   uniprot = unique(boundOut$uniparc)[i]
#   cat(uniprot, '\n')
#   on = length(unique(boundOut$origine[boundOut$uniparc == uniprot]))
#   en = length(unique(boundOut$espece[boundOut$uniparc == uniprot]))
#   nn = length(unique(boundOut$new_fold.Not.Published[boundOut$uniparc == uniprot]))
#   if ( (on + en + nn) > 3 ){
#     cat('\tMismatch!\n')
#     print(boundOut[boundOut$uniparc == uniprot, 1:7])
#   }
# }
# for (i in 1:sum(boundOut$uniparc == 'P18670')){
#   cat(boundOut$pdb[boundOut$uniparc == 'P18670'][i], ' ')
# }
# cat('\n')
  
  
write.table(boundOut,file = paste(gbDir,"data/unilectin/boundUnilectinData.tsv", sep = ""), sep = "\t", quote = F, row.names = F)

# 
# bs = read.delim(file = '/Users/dmattox/cbk/glycan_binding/data/prelimBS.tsv', sep = "\t", stringsAsFactors = F)
# 
# ints = bs[,9:19]
# plot(density(ints$hydrophobics), main = "Distribution of num of hydrophobic ints")
# plot(density(ints$hbonds), main = "Distribution of num of H bond ints")     
# plot(density(ints$total), main = "Distribution of num of total ints")
# plot(density(ints$metal), main = "Distribution of num of metal ints")
# plot(density(ints$sbridges), main = "Distribution of num of salt bridges")
# 
# # Distances between centroids in missing coavelnt bond linkages
# dists = c(5.338466648743167, 5.101439544210736, 5.173914946065158, 5.1792303875632015, 5.293956647855823, 5.139388489976182, 5.177267017354217, 0.1960114271764359, 5.248436341165829, 5.9450798449829145, 0.4136922426263439, 0.049512309676105495, 0.2221522308588242, 0.519004977431586, 6.2193579494518065, 6.8628272158775285, 7.0329599639995255, 0.8242511947648928, 5.800524885664845, 5.8785345257470025, 5.916874044381344, 0.13604711689164317, 0.19362870143137698, 4.975672240491402, 6.994698596460249, 5.402555433879537, 5.284046337739415, 5.418568668348148, 5.269287826942683, 5.376568504850945, 7.154716749182277)
# 
# hist(dists, breaks = 30, col = 'dodgerblue1', xlim = c(0,8), main = "Centroid distances for missing bonds", xlab = "Distance (Ang)")
# abline(v = 1.2, col = 'red')
# 
# allDists = read.delim(file = '/Users/dmattox/cbk/glycan_binding/data/supp/allCentDists.csv', sep = ',', header = F, nrows = 1)
# allDists = as.numeric(allDists)
# 
# plot(density(allDists), xlim = c(0,10), ylim = c(0,1), col = 'firebrick1', xlab = 'Distance (Ang)', main = 'Density of centroid distances')
# par(new = T)
# plot(density(dists), xlim = c(0,10), ylim = c(0,1), col ='dodgerblue1', axes = F, main = '', ylab = '', xlab = '')
# legend(x = 0, y = 1, fill = c('dodgerblue1', 'firebrick1'), legend = c('Missing links','All in thresh'))  
# 
# 
# allDists0_3 = allDists[allDists <= 3]
# hist(allDists0_3, breaks = 15, col = 'firebrick1', xlab = 'Distance (Ang)', main = 'Centroid dists <=3 for all within threshold')
# 
