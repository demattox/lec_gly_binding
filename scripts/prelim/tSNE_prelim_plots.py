#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 20:58:18 2020

@author: dmattox
"""

import numpy as np
import pylab as plt
import sys; sys.path.append('/Users/dmattox/FIt-SNE')
from fast_tsne import fast_tsne

datIn = '/Users/dmattox/cbk/glycan_binding/analysis/prelim_pc50.csv'
labsIn = '/Users/dmattox/cbk/glycan_binding/analysis/prelim_pc50_labels.tsv'
        
X = np.genfromtxt(datIn, delimiter = ',',  skip_header = 1)

origins = []
iupac = []
fold = []
uniprot = []
with open(labsIn, 'r') as inFH:
    for line in inFH:
        line = [l.strip() for l in line.split('\t')]
        if line[0] == 'pdb':
            oriInd = line.index('origine')
            iupacInd = line.index('iupac')
            foldInd = line.index('new_fold.Not.Published')
            uniInd = line.index('uniparc')
            continue
        uniprot.append(line[uniInd])
        origins.append(line[oriInd])
        iupac.append(line[iupacInd])
        fold.append(line[foldInd])

len(set(iupac))
len(set(fold))


colors = ['#DDACC9', '#FF88AD', '#FFB8CE', '#DD6091', '#FF7290', '#FFA388',
        '#C77963', '#9440F3', '#9900B3', '#C266D1', '#6C00BF', '#A700FF',
        '#CA66FF', '#7779BF', '#8194CC', '#533691', '#9189FF', '#B09FFF',
        '#756FB3', '#9FAAFF', '#FF00FF', '#AF00E6', '#FF00B3', '#B3128A',
        '#FF4DC1', '#BD3D9A', '#882E81', '#AD589A', '#AC3491', '#FFFF00',
        '#FFBB33', '#804811', '#B06411', '#BF480D', '#CC6D3D', '#FFDF11',
        '#D6C300', '#FF8011', '#FF9F2C', '#FFB307', '#D9C566', '#BF9F00',
        '#806B19', '#B95541', '#C77767', '#C11331', '#BF8219', '#994C00',
        '#802600', '#A81111', '#ED4C50', '#FF2F7E', '#FF4343', '#BC2B11',
        '#C94545', '#E62A5D', '#D6221D', '#E67A77', '#AF3F64', '#FF197F',
        '#D9F077', '#A6E6A9', '#7AE6AB', '#82AD7D', '#B8FFCA', '#ADE6A6',
        '#00979D', '#00DDC5', '#26FFF2', '#00A809', '#00FF00', '#26BF64',
        '#33A9CE', '#0094C2', '#00A79D', '#2F8C4D', '#00879D', '#669D6A',
        '#005C07', '#008F1F', '#53879D', '#53A39D', '#174F61', '#A19922',
        '#7F9922', '#C2E32C', '#96E32C', '#5100FF', '#0000FF', '#22737F',
        '#247740', '#00A863', '#29E043', '#53D385', '#1E806D', '#4A8044',
        '#008F39', '#3D9946', '#73CA95', '#47867A', '#00B8C3', '#006091',
        '#5C89CC', '#74A0FF', '#74CAFF', '#69A8E6', '#578EBF', '#1F6666',
        '#2B7880', '#388899', '#266180', '#494566', '#336D99', '#254566',
        '#335280', '#FF0000', '#00FF66', '#665C47', '#476655', '#475D4B',
        '#6B998D', '#6C8581', '#476662', '#664747', '#89997A', '#6B8059',
        '#74996B', '#665547', '#807059', '#997F7A', '#805F59', '#4C6647',
        '#598069']

colors = colors*4

print(X.shape)

foldDict = {}
for i,f in enumerate(set(fold)):
    foldDict[f] = colors[i]
foldCols = [foldDict[f] for f in fold]

origDict = {}
cols = ["#FF0000", "#E18419", "#6B9521", "#2E7D6E", "#5153DE", "#9400D3"]
origUni = list(set(origins))
origUni.sort()
for i,o in enumerate(origUni):
    origDict[o] = cols[i]
origCols = [origDict[o] for o in origins]

iupacDict = {'Man': 'green', 'Fuc': 'red', 'NeuAc': '#6C00BF', 'Gal': 'gold', 'Glc': 'blue'}
iupacCols = [iupacDict[c] if c in iupacDict.keys() else 'grey' for c in iupac]

uniDict = {}
for i,p in enumerate(set(uniprot)):
    uniDict[p] = colors[i]
uniCols = [uniDict[p] for p in uniprot]


# Check PCA
plt.figure(figsize=(4,4))
plt.scatter(X[:,0], X[:,1], s=1,
            color = foldCols)
plt.title('PCA colored by fold')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.tight_layout()

plt.figure(figsize=(4,4))
plt.scatter(X[:,0], X[:,1], s=1,
            color = origCols)
plt.title('PCA colored by origin')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.tight_layout()

plt.figure(figsize=(4,4))
plt.scatter(X[:,0], X[:,1], s=1,
            color = uniCols)
plt.title('PCA colored by uniprotID')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.tight_layout()

plt.figure(figsize=(4,4))
plt.scatter(X[:,0], X[:,1], s=1,
            color = iupacCols)
plt.title('PCA colored by iupac monosaccs')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.tight_layout()


# generic tSNE
%time tsne_default = fast_tsne(X, seed=42)

#save out tSNE
outFile = '/Users/dmattox/cbk/glycan_binding/analysis/tSNE_prelim_p30.csv'
with open(outFile, 'w') as outFH:
    for x,y in tsne_default:
        outFH.write(str(x) + ',' + str(y) + '\n')

plt.figure(figsize=(5,5))
plt.gca().set_aspect('equal', adjustable='datalim')
plt.scatter(tsne_default[:,0], tsne_default[:,1], s=1,
            color = foldCols)
plt.title('Lectin fold')
plt.tight_layout()

plt.figure(figsize=(5,5))
plt.gca().set_aspect('equal', adjustable='datalim')
plt.scatter(tsne_default[:,0], tsne_default[:,1], s=1,
            color = origCols)
plt.title('Lectin origin')
plt.tight_layout()

plt.figure(figsize=(5,5))
plt.gca().set_aspect('equal', adjustable='datalim')
plt.scatter(tsne_default[:,0], tsne_default[:,1], s=1,
           color = uniCols)
plt.title('Uniprot IDs')
plt.tight_layout()

plt.figure(figsize=(5,5))
plt.gca().set_aspect('equal', adjustable='datalim')
plt.scatter(tsne_default[:,0], tsne_default[:,1], s=1,
           color = iupacCols)
plt.title('Monosaccharide ligand identity')
plt.tight_layout()


tsne_10 = fast_tsne(X, seed=42, perplexity = 10)

plt.figure(figsize=(5,5))
plt.gca().set_aspect('equal', adjustable='datalim')
plt.scatter(tsne_10[:,0], tsne_10[:,1], s=1,
            color = foldCols)
plt.title('t-SNE, perplexity = 10 by fold')
plt.tight_layout()


tsne_50 = fast_tsne(X, seed=42, perplexity = 50)

plt.figure(figsize=(5,5))
plt.gca().set_aspect('equal', adjustable='datalim')
plt.scatter(tsne_50[:,0], tsne_50[:,1], s=1,
            color = foldCols)
plt.title('t-SNE, perplexity = 50 by fold')
plt.tight_layout()



# improved tSNE
%%time

pcaInit = X[:,:2] / np.std(X[:,0]) * 0.0001
tsne_ours = fast_tsne(X, initialization = pcaInit,
                      learning_rate = X.shape[0]/12,
                      perplexity_list = [30, int(X.shape[0]/100)])

plt.figure(figsize=(5,5))
plt.gca().set_aspect('equal', adjustable='datalim')
plt.scatter(tsne_ours[:,0], tsne_ours[:,1], s=1, color = foldCols)
plt.title('A better t-SNE by fold')
plt.tight_layout()

plt.figure(figsize=(5,5))
plt.gca().set_aspect('equal', adjustable='datalim')
plt.scatter(tsne_ours[:,0], tsne_ours[:,1], s=1, color = origCols)
plt.title('A better t-SNE by origin')
plt.tight_layout()

###########################################
#Repeat for monosaccharide subset

datIn = '/Users/dmattox/cbk/glycan_binding/analysis/prelimMonos_pc50.csv'
labsIn = '/Users/dmattox/cbk/glycan_binding/analysis/prelimMonos_pc50_labels.tsv'
        
X = np.genfromtxt(datIn, delimiter = ',',  skip_header = 1)

origins = []
iupac = []
fold = []
uniprot = []
with open(labsIn, 'r') as inFH:
    for line in inFH:
        line = [l.strip() for l in line.split('\t')]
        if line[0] == 'pdb':
            oriInd = line.index('origine')
            iupacInd = line.index('iupac')
            foldInd = line.index('new_fold.Not.Published')
            uniInd = line.index('uniparc')
            continue
        uniprot.append(line[uniInd])
        origins.append(line[oriInd])
        iupac.append(line[iupacInd])
        fold.append(line[foldInd])

len(set(iupac))
len(set(fold))


colors = ['#DDACC9', '#FF88AD', '#FFB8CE', '#DD6091', '#FF7290', '#FFA388',
        '#C77963', '#9440F3', '#9900B3', '#C266D1', '#6C00BF', '#A700FF',
        '#CA66FF', '#7779BF', '#8194CC', '#533691', '#9189FF', '#B09FFF',
        '#756FB3', '#9FAAFF', '#FF00FF', '#AF00E6', '#FF00B3', '#B3128A',
        '#FF4DC1', '#BD3D9A', '#882E81', '#AD589A', '#AC3491', '#FFFF00',
        '#FFBB33', '#804811', '#B06411', '#BF480D', '#CC6D3D', '#FFDF11',
        '#D6C300', '#FF8011', '#FF9F2C', '#FFB307', '#D9C566', '#BF9F00',
        '#806B19', '#B95541', '#C77767', '#C11331', '#BF8219', '#994C00',
        '#802600', '#A81111', '#ED4C50', '#FF2F7E', '#FF4343', '#BC2B11',
        '#C94545', '#E62A5D', '#D6221D', '#E67A77', '#AF3F64', '#FF197F',
        '#D9F077', '#A6E6A9', '#7AE6AB', '#82AD7D', '#B8FFCA', '#ADE6A6',
        '#00979D', '#00DDC5', '#26FFF2', '#00A809', '#00FF00', '#26BF64',
        '#33A9CE', '#0094C2', '#00A79D', '#2F8C4D', '#00879D', '#669D6A',
        '#005C07', '#008F1F', '#53879D', '#53A39D', '#174F61', '#A19922',
        '#7F9922', '#C2E32C', '#96E32C', '#5100FF', '#0000FF', '#22737F',
        '#247740', '#00A863', '#29E043', '#53D385', '#1E806D', '#4A8044',
        '#008F39', '#3D9946', '#73CA95', '#47867A', '#00B8C3', '#006091',
        '#5C89CC', '#74A0FF', '#74CAFF', '#69A8E6', '#578EBF', '#1F6666',
        '#2B7880', '#388899', '#266180', '#494566', '#336D99', '#254566',
        '#335280', '#FF0000', '#00FF66', '#665C47', '#476655', '#475D4B',
        '#6B998D', '#6C8581', '#476662', '#664747', '#89997A', '#6B8059',
        '#74996B', '#665547', '#807059', '#997F7A', '#805F59', '#4C6647',
        '#598069']

colors = colors*4

print(X.shape)

foldDict = {}
for i,f in enumerate(set(fold)):
    foldDict[f] = colors[i]
foldCols = [foldDict[f] for f in fold]

origDict = {}
cols = ["#FF0000", "#E86313", "#B59F20", "#469021", "#2E7D6E", "#4169E1"]
origUni = list(set(origins))
origUni.sort()
for i,o in enumerate(origUni):
    origDict[o] = cols[i]
origCols = [origDict[o] for o in origins]

iupacDict = {'Man': 'green', 'Fuc': 'red', 'NeuAc': '#6C00BF', 'Gal': 'gold', 'Glc': 'blue'}
iupacCols = [iupacDict[c] if c in iupacDict.keys() else 'grey' for c in iupac]

uniDict = {}
for i,p in enumerate(set(uniprot)):
    uniDict[p] = colors[i]
uniCols = [uniDict[p] for p in uniprot]


# Check PCA
plt.figure(figsize=(4,4))
plt.scatter(X[:,0], X[:,1], s=1,
            color = foldCols)
plt.title('PCA colored by fold')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.tight_layout()

plt.figure(figsize=(4,4))
plt.scatter(X[:,0], X[:,1], s=1,
            color = origCols)
plt.title('PCA colored by origin')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.tight_layout()

plt.figure(figsize=(4,4))
plt.scatter(X[:,0], X[:,1], s=1,
            color = uniCols)
plt.title('PCA colored by uniprotID')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.tight_layout()

plt.figure(figsize=(4,4))
plt.scatter(X[:,0], X[:,1], s=1,
            color = iupacCols)
plt.title('PCA colored by iupac monosaccs')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.tight_layout()


# generic tSNE
%time tsne_default = fast_tsne(X, seed=42)

plt.figure(figsize=(5,5))
plt.gca().set_aspect('equal', adjustable='datalim')
plt.scatter(tsne_default[:,0], tsne_default[:,1], s=1,
            color = foldCols)
plt.title('t-SNE, perplexity = 30 by fold')
plt.tight_layout()

plt.figure(figsize=(5,5))
plt.gca().set_aspect('equal', adjustable='datalim')
plt.scatter(tsne_default[:,0], tsne_default[:,1], s=1,
            color = origCols)
plt.title('t-SNE, perplexity = 30 by origin')
plt.tight_layout()

plt.figure(figsize=(5,5))
plt.gca().set_aspect('equal', adjustable='datalim')
plt.scatter(tsne_default[:,0], tsne_default[:,1], s=1,
           color = uniCols)
plt.title('t-SNE, perplexity = 30 by uniprot ID')
plt.tight_layout()

plt.figure(figsize=(5,5))
plt.gca().set_aspect('equal', adjustable='datalim')
plt.scatter(tsne_default[:,0], tsne_default[:,1], s=1,
           color = iupacCols)
plt.title('t-SNE, perplexity = 30 by iupac monosaccs')
plt.tight_layout()


#############################################################

# Repeat tSNE with single binding site from each PDB file 

datIn = '/Users/dmattox/cbk/glycan_binding/analysis/uniprelim_pc50.csv'
labsIn = '/Users/dmattox/cbk/glycan_binding/analysis/uniprelim_pc50_labels.tsv'
        
X = np.genfromtxt(datIn, delimiter = ',',  skip_header = 1)

origins = []
iupac = []
fold = []
uniprot = []
with open(labsIn, 'r') as inFH:
    for line in inFH:
        line = [l.strip() for l in line.split('\t')]
        if line[0] == 'pdb':
            oriInd = line.index('origine')
            iupacInd = line.index('iupac')
            foldInd = line.index('new_fold.Not.Published')
            uniInd = line.index('uniparc')
            continue
        uniprot.append(line[uniInd])
        origins.append(line[oriInd])
        iupac.append(line[iupacInd])
        fold.append(line[foldInd])

len(set(iupac))
len(set(fold))


colors = ['#DDACC9', '#FF88AD', '#FFB8CE', '#DD6091', '#FF7290', '#FFA388',
        '#C77963', '#9440F3', '#9900B3', '#C266D1', '#6C00BF', '#A700FF',
        '#CA66FF', '#7779BF', '#8194CC', '#533691', '#9189FF', '#B09FFF',
        '#756FB3', '#9FAAFF', '#FF00FF', '#AF00E6', '#FF00B3', '#B3128A',
        '#FF4DC1', '#BD3D9A', '#882E81', '#AD589A', '#AC3491', '#FFFF00',
        '#FFBB33', '#804811', '#B06411', '#BF480D', '#CC6D3D', '#FFDF11',
        '#D6C300', '#FF8011', '#FF9F2C', '#FFB307', '#D9C566', '#BF9F00',
        '#806B19', '#B95541', '#C77767', '#C11331', '#BF8219', '#994C00',
        '#802600', '#A81111', '#ED4C50', '#FF2F7E', '#FF4343', '#BC2B11',
        '#C94545', '#E62A5D', '#D6221D', '#E67A77', '#AF3F64', '#FF197F',
        '#D9F077', '#A6E6A9', '#7AE6AB', '#82AD7D', '#B8FFCA', '#ADE6A6',
        '#00979D', '#00DDC5', '#26FFF2', '#00A809', '#00FF00', '#26BF64',
        '#33A9CE', '#0094C2', '#00A79D', '#2F8C4D', '#00879D', '#669D6A',
        '#005C07', '#008F1F', '#53879D', '#53A39D', '#174F61', '#A19922',
        '#7F9922', '#C2E32C', '#96E32C', '#5100FF', '#0000FF', '#22737F',
        '#247740', '#00A863', '#29E043', '#53D385', '#1E806D', '#4A8044',
        '#008F39', '#3D9946', '#73CA95', '#47867A', '#00B8C3', '#006091',
        '#5C89CC', '#74A0FF', '#74CAFF', '#69A8E6', '#578EBF', '#1F6666',
        '#2B7880', '#388899', '#266180', '#494566', '#336D99', '#254566',
        '#335280', '#FF0000', '#00FF66', '#665C47', '#476655', '#475D4B',
        '#6B998D', '#6C8581', '#476662', '#664747', '#89997A', '#6B8059',
        '#74996B', '#665547', '#807059', '#997F7A', '#805F59', '#4C6647',
        '#598069']

colors = colors*4

print(X.shape)

foldDict = {}
for i,f in enumerate(set(fold)):
    foldDict[f] = colors[i]
foldCols = [foldDict[f] for f in fold]

origDict = {}
cols = ["#FF0000", "#E18419", "#6B9521", "#2E7D6E", "#5153DE", "#9400D3"]
origUni = list(set(origins))
origUni.sort()
for i,o in enumerate(origUni):
    origDict[o] = cols[i]
origCols = [origDict[o] for o in origins]

iupacDict = {'Man': 'green', 'Fuc': 'red', 'NeuAc': '#6C00BF', 'Gal': 'gold', 'Glc': 'blue'}
iupacCols = [iupacDict[c] if c in iupacDict.keys() else 'grey' for c in iupac]

uniDict = {}
for i,p in enumerate(set(uniprot)):
    uniDict[p] = colors[i]
uniCols = [uniDict[p] for p in uniprot]


# Check PCA
plt.figure(figsize=(4,4))
plt.scatter(X[:,0], X[:,1], s=1,
            color = foldCols)
plt.title('PCA colored by fold')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.tight_layout()

plt.figure(figsize=(4,4))
plt.scatter(X[:,0], X[:,1], s=1,
            color = origCols)
plt.title('PCA colored by origin')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.tight_layout()

plt.figure(figsize=(4,4))
plt.scatter(X[:,0], X[:,1], s=1,
            color = uniCols)
plt.title('PCA colored by uniprotID')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.tight_layout()

plt.figure(figsize=(4,4))
plt.scatter(X[:,0], X[:,1], s=1,
            color = iupacCols)
plt.title('PCA colored by iupac monosaccs')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.tight_layout()


# generic tSNE
%time tsne_default = fast_tsne(X, seed=42, perplexity = 30)

plt.figure(figsize=(5,5))
plt.gca().set_aspect('equal', adjustable='datalim')
plt.scatter(tsne_default[:,0], tsne_default[:,1], s=1,
            color = foldCols)
plt.title('t-SNE, perplexity = 30 by fold')
plt.tight_layout()

plt.figure(figsize=(5,5))
plt.gca().set_aspect('equal', adjustable='datalim')
plt.scatter(tsne_default[:,0], tsne_default[:,1], s=1,
            color = origCols)
plt.title('t-SNE, perplexity = 30 by origin')
plt.tight_layout()

plt.figure(figsize=(5,5))
plt.gca().set_aspect('equal', adjustable='datalim')
plt.scatter(tsne_default[:,0], tsne_default[:,1], s=1,
           color = uniCols)
plt.title('t-SNE, perplexity = 30 by uniprot ID')
plt.tight_layout()

plt.figure(figsize=(5,5))
plt.gca().set_aspect('equal', adjustable='datalim')
plt.scatter(tsne_default[:,0], tsne_default[:,1], s=1,
           color = iupacCols)
plt.title('t-SNE, perplexity = 30 by iupac monosaccs')
plt.tight_layout()














