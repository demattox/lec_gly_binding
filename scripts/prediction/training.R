
library(randomForest)
library(ggplot2)
library(caret)
# library(readxl)

# indicate home dir for projet
homeDir = '/Users/dmattox/cbk/lec_gly_binding/'
setwd(homeDir)

##########################
# functions
###########################




##########################
# Set up models
###########################
# Set a random seed
set.seed(27)

# Read in data
ligTags = read.delim(file = './analysis/training/data_in/ligTags.tsv', sep = '\t', stringsAsFactors = F)
predFeats = read.delim(file = './analysis/training/data_in/predFeats.csv', sep = ',', stringsAsFactors = F)
bsResiDat = read.delim(file = './analysis/training/data_in/bsResiDat.tsv', sep = '\t', stringsAsFactors = F)

# Prep RF model specifications
mtry = floor(sqrt(ncol(predFeats))) # 209 --> 14
train.control <- trainControl(method = "repeatedcv", number = 10, repeats = 3)



