#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 13:06:12 2020

@author: dmattox
"""

import os
import sys
import dill
import collections
import numpy as np

import lec_gly as LecGly

############################################

if __name__ == '__main__':
    
    os.chdir(LecGly.homeDir)
    
    unilecDataFH = './data/unilectin/boundUnilectinData.tsv'
    qcReport = './analysis/QC/ligandCheckReport.txt'
    
    uniprotMap = collections.defaultdict(list) # Dictionary to store PDB IDs associated with each UniProt ID
    unilecData = {} # dictionary for importing metadata for each lectin provided by UniLectin
    unilecCols = [] # for storing column headers of interest from
    
    with open(unilecDataFH, 'r') as dataIn:
        for line in dataIn:
            line = line.split('\t')
            if (line[0] == 'pdb'):
                    unilecCols = [line[1]] + line[2:5] + line[7:12] # Get column names of interest
            else:
                unilecData[line[0]] = [line[1]] + line[2:5] + line[7:12] # Get data from columns of interest
                uniprotMap[line[1]].append(line[0]) # Get uniprot-PDB associations
    
    knownInfo = {}
    for pdb in unilecData.keys():
        knownInfo[pdb] = {k:v for k,v in zip(unilecCols, unilecData[pdb])}
    
    with open("./data/unilectin/PLIPrecords.p", "rb") as pickleFH:
        allPLIP = dill.load(pickleFH)
    
    with open("./data/unilectin/icodeResi2PLIPmapping.p", "rb")as pickleFH:
        zeroResis = dill.load(pickleFH)
    
    ############################################
    
    # Binning distances
    bins = [3.5, 4.5, 6.5, 8]
    
    ############################################
    
    dropList = collections.defaultdict(lambda: collections.defaultdict(list)) # Nested distonary to store reasons for exclusion along with the pdb id and the list of bs ids to drop
    flgCnt = 0
    
    doManualFilter = input('Would you like to loop through all structures containing any potentially false binding sites (fewer than 3 residues within 4.5 Anstroms OR fewer than 3 PLIP interactions) ("y"/"n"):\n\t')
    if doManualFilter == 'n':
        useOldFilter = input('Would you like to refilter with a previously determined filter list? \n(stored as ./data/unilectin/manuallyFilteredBSites.csv)\n("y"/"n"):\n\t')
        if useOldFilter == 'y':
            try:
                with open('./data/unilectin/manuallyFilteredBSites.csv', 'r') as inFH:
                    for line in inFH:
                        line = line.split(',')
                        line = [l.strip() for l in line]
                        if len(line) == 1:
                            reason = line[0]
                        else:
                            pdb = line.pop(0)
                            dropList[reason][pdb].extend(line)
            except:
                print('The file "manuallyFilteredBSites.csv" was not found or was improperly formatted\n')
        else:
            print('No manual filter and no refiltering selecting. Saving input PLIPrecords.p out as cleanPLIPrecords.p ***without cleaning***\n')
            with open("./data/unilectin/cleanPLIPrecords.p", "wb") as pickleFH:
                dill.dump(allPLIP, pickleFH)
            sys.exit()
    else: # do manual filter!
    
        # pdb = '2J3U'
        
        for pdb in allPLIP.keys():
            # Get info on all binding sites for PDB
            bsiteInfo = collections.defaultdict(dict)
            for bs in allPLIP[pdb].bsites:
                bsiteInfo[bs]['longname'] = allPLIP[pdb].bsites[bs].longname
                bsiteInfo[bs]['binCnts'] = [0,0,0,0]
                bsiteInfo[bs]['totalInteractions'] = allPLIP[pdb].bsites[bs].counts['total']
                for r in allPLIP[pdb].bsites[bs].bs_res:
                    for bInd, bDist in enumerate(bins):
                        if r['min_dist'] <= bDist:
                            bsiteInfo[bs]['binCnts'][bInd] += 1
                            break
            ###
            #Check for issues such as a multiple longnames for the same PDB or no/few residues within first two bins
            allLongs = list(set([ bsiteInfo[bs]['longname'] for  bs in bsiteInfo]))
            bin1_2Cnts = np.array([ sum(bsiteInfo[bs]['binCnts'][:2]) for bs in bsiteInfo ])
            intCnts = np.array([ bsiteInfo[bs]['totalInteractions'] for bs in bsiteInfo ])
            
            if (any(bin1_2Cnts < 3) or any(intCnts < 3)):
                flgCnt += 1
                print('Flagged case number ' + str(flgCnt))
                print('\n####\nPotential false binding sites detected in PLIP record for', pdb + '\n####')
                print('\tUnilectin identified ligand (iupac):\n\t\t', knownInfo[pdb]['iupac'])
                print('\tAdditional Unilectin ligand info:\n\t\t', knownInfo[pdb]['ligand'], '\n\t\t', knownInfo[pdb]['sucres'],'\n')
                print('\tPLIP ligands by longname & PLIP bsID:')
                for bs in bsiteInfo:
                    print('\t\t',bsiteInfo[bs]['longname'], '\n\t\t[', bs,']')
                    print('\t\t\tTotal interaction count:', bsiteInfo[bs]['totalInteractions'])
                    print('\t\t\tResidues in bins 1-4 respectively:',', '.join([str(c) for c in bsiteInfo[bs]['binCnts']]))
                filt = input('Would you like to remove any of these binding sites? ("y"/"n"):\n')
                if filt == 'y':
                    dropCnt = 0
                    for bs in bsiteInfo:
                        print('Would you like to remove this binding site?')
                        print('\t\t',bsiteInfo[bs]['longname'], '\n\t\t[', bs,']')
                        print('\t\t\tTotal interaction count:', bsiteInfo[bs]['totalInteractions'])
                        print('\t\t\tResidues in bins 1-4 respectively:',', '.join([str(c) for c in bsiteInfo[bs]['binCnts']]))
                        drop = input('"y" or "n":\n')
                        if drop == 'y':
                            reason = input('Why is this binding site excluded? \n\t Missed glycosylation ("g"), Missed covalent bond ("c"), Other ("o") \n\t press any other key to cancel dropping this binding site\n\n\t')
                            if reason == 'g':
                                dropList['Missed_glycosylation'][pdb].append(bs)
                                dropCnt += 1
                            elif reason == 'c':
                                dropList['Missed_covalent_link'][pdb].append(bs)
                                dropCnt += 1
                            elif reason == 'o':
                                dropList['Other'][pdb].append(bs)
                                dropCnt += 1
                            else:
                                print('Input reason ' + reason + 'not a valid option, binding site kept')
                            
                    print('Dropped', dropCnt,'binding site(s) from', pdb,'\n\n')
                else:
                    print('\n\n')
        
        print('\n\nDONE!\n\n' + str(flgCnt) + ' PLIP records (out of ' + str(len(allPLIP.keys())) + ') flagged for manual inspection\n ')
        ###

    # Remove flagged binding sites
    zeroedOut = []
    totDropCnt = 0
    for reason in dropList.keys():
        for pdb in dropList[reason]:
            for bs in dropList[reason][pdb]:
                del allPLIP[pdb].bsites[bs]
                totDropCnt += 1
            allPLIP[pdb].num_bsites = len(allPLIP[pdb].bsites) # update binding site counter
            
            if allPLIP[pdb].num_bsites == 0: # Check if all ligands for a PDB have been removed
                zeroedOut.append(pdb)
                del allPLIP[pdb]
    
    pdbDropCnt = len(set( [ pdb for reason in dropList.keys() for pdb in dropList[reason].keys() ] ))
    
    
    print(totDropCnt, 'false binding sites removed from', pdbDropCnt, 'PLIP records')
    if zeroedOut:
        print('WARNING!\n\t All ligands removed for', len(zeroedOut),'PDB IDs')
        print('\t',zeroedOut)
    
    with open(qcReport, "a")as qcFH:
        qcFH.write('\n\n#################################\nAddendum, added on by manualBSfilter.py\n#################################\n\n')
        qcFH.write('Manual inspection of binding sites in lectin structures if any binding sites in the structure have:\n\tFewer than 3 interactions\n\t(OR)\n\tFewer than 3 proteins residues within 4.5 Angstoms\n')
        qcFH.write('Binding sites are dropped if the glycan in the binding site appears to be a glycosylation occurence, a separate binding site because of a missing covalent bond, or other (artifact/glycan not near protein/etc.)\n\n')
        qcFH.write('In total, ' + str(totDropCnt) + ' false binding sites removed from ' + str(pdbDropCnt) + ' PLIP records\n\n')
        if zeroedOut:
            qcFH.write('WARNING! All ligands removed for ' + str(len(zeroedOut)) + ' PDB IDs:\n')
            qcFH.write('\t' + '\n\t'.join(zeroedOut) + '\n\n')
        for reason in dropList.keys():
            qcFH.write(str(len(dropList[reason])) + ' PDB IDs with dropped binding sites for reason (' + reason.replace('_',' ').lower() + '):\n')
            for pdb in dropList[reason]:
                qcFH.write('\t' + pdb + ': ' + ', '.join(dropList[reason][pdb]) + '\n')
                
        pdbCnt = len(allPLIP)
        bsCnt = 0
        for pdb in allPLIP.keys():
            for bs in allPLIP[pdb].bsites:
                bsCnt += 1
        qcFH.write('#####\nThere remain ' + str(bsCnt) + ' binding sites from ' + str(pdbCnt) + ' structures\n#####\n\n')
    
    ###
    # Save out as cleaned
    with open("./data/unilectin/cleanPLIPrecords.p", "wb") as pickleFH:
        dill.dump(allPLIP, pickleFH)
    
    # Save out dropped binding sites
    with open("./data/unilectin/manuallyFilteredBSites.csv", "w")as outFH:
        for reason in dropList.keys():
            outFH.write(reason + '\n')
            for pdb in dropList[reason]:
                out = [pdb]
                out.extend(dropList[reason][pdb])
                outFH.write(', '.join(out) + '\n')
    
    print('DONE!\nCleaned PLIP records serialized to file cleanPLIPrecords.p\nDetails of filtering appended to QC log file in ./analysis/QC/\n\n')
    
        
        
        
        
        
        
    
    
    
    
    
    
    
    
    
    
