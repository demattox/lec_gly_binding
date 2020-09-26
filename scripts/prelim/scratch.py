#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 17:45:39 2020

@author: dmattox
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

lst = ['OC[C@H]1OC[C@@H]([C@H]([C@@H]1O)O)O', 'OC[C@H]1OC(O[C@@H]2[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O']

mols = [Chem.MolFromSmiles(smile, sanitize = False) for smile in lst]

mergeMol = Chem.CombineMols(mols[0], mols[1])
edMerge = Chem.EditableMol(mergeMol)
edMerge.AddBond(20, 4, order = Chem.rdchem.BondType.SINGLE)
mergeMol = edMerge.GetMol()

print(Chem.MolToSmiles(mergeMol))
print(Chem.inchi.MolToInchiKey(mergeMol))
print(Descriptors.MolLogP(mergeMol)) 


{'ACE:A:0': <pliptool.plip.modules.plipxml.BSite at 0x128700710>,
 'ACE:B:0': <pliptool.plip.modules.plipxml.BSite at 0x1281c0a10>,
 'BGC:A:122': <pliptool.plip.modules.plipxml.BSite at 0x1287cadd0>,
 'BGC:A:124': <pliptool.plip.modules.plipxml.BSite at 0x1287caed0>,
 'BGC:A:126': <pliptool.plip.modules.plipxml.BSite at 0x1287a8350>,
 'BGC:B:122': <pliptool.plip.modules.plipxml.BSite at 0x1287a8e10>,
 'BGC:B:124': <pliptool.plip.modules.plipxml.BSite at 0x1287a8410>,
 'BGC:B:126': <pliptool.plip.modules.plipxml.BSite at 0x1287a8ed0>,
 'EDO:A:601': <pliptool.plip.modules.plipxml.BSite at 0x1287a8210>,
 'EDO:A:603': <pliptool.plip.modules.plipxml.BSite at 0x1287a8890>,
 'EDO:B:602': <pliptool.plip.modules.plipxml.BSite at 0x12879ea90>,
 'GLC:A:123': <pliptool.plip.modules.plipxml.BSite at 0x12879ee10>,
 'GLC:A:125': <pliptool.plip.modules.plipxml.BSite at 0x1281b5950>,
 'GLC:A:127': <pliptool.plip.modules.plipxml.BSite at 0x12877aad0>,
 'GLC:B:123': <pliptool.plip.modules.plipxml.BSite at 0x12877af90>,
 'GLC:B:125': <pliptool.plip.modules.plipxml.BSite at 0x128785450>,
 'GLC:B:127': <pliptool.plip.modules.plipxml.BSite at 0x128785fd0>,
 'MG:A:501': <pliptool.plip.modules.plipxml.BSite at 0x128796090>,
 'MG:A:502': <pliptool.plip.modules.plipxml.BSite at 0x128796e10>,
 'SO4:A:401': <pliptool.plip.modules.plipxml.BSite at 0x1287969d0>,
 'SO4:A:402': <pliptool.plip.modules.plipxml.BSite at 0x128796a50>,
 'SO4:B:403': <pliptool.plip.modules.plipxml.BSite at 0x126dab410>}

{'ACE:A:0': <pliptool.plip.modules.plipxml.BSite at 0x1288274d0>,
 'ACE:B:0': <pliptool.plip.modules.plipxml.BSite at 0x128827f50>,
 'BGC:A:122': <pliptool.plip.modules.plipxml.BSite at 0x128827410>,
 'BGC:A:124': <pliptool.plip.modules.plipxml.BSite at 0x12885b450>,
 'BGC:A:126': <pliptool.plip.modules.plipxml.BSite at 0x12885b750>,
 'BGC:B:122': <pliptool.plip.modules.plipxml.BSite at 0x12885bb90>,
 'BGC:B:124': <pliptool.plip.modules.plipxml.BSite at 0x128832c10>,
 'BGC:B:126': <pliptool.plip.modules.plipxml.BSite at 0x128832bd0>,
 'EDO:A:601': <pliptool.plip.modules.plipxml.BSite at 0x128832290>,
 'EDO:A:603': <pliptool.plip.modules.plipxml.BSite at 0x128832cd0>,
 'EDO:B:602': <pliptool.plip.modules.plipxml.BSite at 0x126db2f10>,
 'MG:A:501': <pliptool.plip.modules.plipxml.BSite at 0x128820f10>,
 'MG:A:502': <pliptool.plip.modules.plipxml.BSite at 0x128820610>,
 'SO4:A:401': <pliptool.plip.modules.plipxml.BSite at 0x128839190>,
 'SO4:A:402': <pliptool.plip.modules.plipxml.BSite at 0x128839590>,
 'SO4:B:403': <pliptool.plip.modules.plipxml.BSite at 0x128839450>}



# tree = ET.parse(plipIN)
# root = tree.getroot()

# for child in root:
#     print(child.tag, child.attrib, child.text)
# [elem.tag for elem in root.iter()]
# print(ET.tostring(root, encoding='utf8').decode('utf8'))

#### Initial preview of binding sites
# with open('./data/unilectin/binding_sites.csv', "w") as outFH:
#     outFH.write(','.join(header) + "\n")
#     for name in tmp:
#         pdb = fname.rstrip(".xml")[-4:]
#         print('\n\n---------\n')
#         print(pdb,unilecLigs[pdb])
#         tree = ET.parse(fname)
#         for bs in tree.iter('bindingsite'):
#             if bs.attrib['has_interactions'] and (float(bs.find("lig_properties").find("molweight").text) > 115): # Erythrose (4C monosacc) has MW of 120, preliminary screen against non-glycans
#                 outFH.write(','.join([pdb, unilecSupp[pdb][0], unilecSupp[pdb][1], unilecLigs[pdb][2], unilecSupp[pdb][2]]) + ',')
#                 glyChain = bs.find('identifiers').find('chain').text
#                 longName = bs.find('identifiers').find('longname').text
#                 hetID = bs.find('identifiers').find('hetid').text
#                 print(glyChain, longName, hetID)

# #### Exploring residues deemed to be in the binding site by PLIP
                

# bs_resis = []
# tmp = {}
# intCnt = 0
# for bs in tree.iter('bindingsite'):
#     if bs.attrib['has_interactions'] and (float(bs.find("lig_properties").find("molweight").text) > 115): # Erythrose (4C monosacc) has MW of 120, preliminary screen against non-glycans
#         print(bs.find('identifiers').find('longname').text)
#         print(bs.find('identifiers').find('smiles').text)
#         for res in bs.find('bs_residues'):
#             if ( res.text[-1] == 'A' ):
#                 print(res.tag, res.attrib, res.text)
#                 if (res.attrib['contact'] == 'True'): intCnt += 1
#                 if float(res.attrib['min_dist']) <= 100:
#                     bs_resis.append((int(res.text[:-1]), d[res.attrib["aa"]]) )
#                     tmp[res.text] = res.attrib
#         # for interact in bs.find('interactions'):
#         #     print(interact.tag)
#     elif (bs.attrib['has_interactions'] and (float(bs.find("lig_properties").find("molweight").text) < 115)):
#         print('SMOL BOI!:')
#         print(bs.find('identifiers').find('longname').text)

# bs_resis.sort()
# bs_resis
# for r in bs_resis:
#     print(r[0],end = '+')
# print()
# print('All:',len(bs_resis))
# print('Interacting', intCnt)
# tmp2 = {}
# for k in tmp.keys():
#     tmp2[int(k[:-1])] = tmp[k]
# k2 = list(tmp2.keys())
# k2.sort()
# for k in k2:
#     print(k,tmp2[k])

########################
# Trying out plipxml parser
plipDat = px.PLIPXML(plipIN)

plipDat.num_bsites
plipDat.bsites

bs = plipDat.bsites['MAN:A:122']
bs.bsid
bs.hetid
bs.longname
bs.position
bs.smiles
bs.molweight

bs.bs_res
cnt = 0
for r in bs.bs_res:
    if r['contact']:
        print(r)
        cnt+=1
print(cnt)
bs.num_contacts # Number of interactions
bs.get_counts()

hit = False
for f in plipFiles:
    plipDat = px.PLIPXML(f)
    for bs in plipDat.bsites:
        if plipDat.bsites[bs].composite:
            print(bs, plipDat.bsites[bs].uniqueid)
            hit = True
    if hit == True:
        print(f)
        break
    
# Checking Ca2+ coordination
for site in plipDat.bsites:
    bs = plipDat.bsites[site]
    if bs.molweight > 115:
        print(site, bs.longname, 'composite:',bs.composite)
        print(len(bs.bs_res),'binding site residues')
        if bs.counts['total'] != 0:
           print(bs.counts['total'], 'interactions')
        print(bs.counts)
        print()
###############################
biolip_list = ['ACE', 'HEX', 'TMA', 'SOH', 'P25', 'CCN', 'PR', 'PTN', 'NO3', 'TCN', 'BU1', 'BCN', 'CB3', 'HCS', 'NBN',
               'SO2', 'MO6', 'MOH', 'CAC', 'MLT', 'KR', '6PH', 'MOS', 'UNL', 'MO3', 'SR', 'CD3', 'PB', 'ACM', 'LUT',
               'PMS', 'OF3', 'SCN', 'DHB', 'E4N', '13P', '3PG', 'CYC', 'NC', 'BEN', 'NAO', 'PHQ', 'EPE', 'BME', 'TB',
               'ETE', 'EU', 'OES', 'EAP', 'ETX', 'BEZ', '5AD', 'OC2', 'OLA', 'GD3', 'CIT', 'DVT', 'OC6', 'MW1', 'OC3',
               'SRT', 'LCO', 'BNZ', 'PPV', 'STE', 'PEG', 'RU', 'PGE', 'MPO', 'B3P', 'OGA', 'IPA', 'LU', 'EDO', 'MAC',
               '9PE', 'IPH', 'MBN', 'C1O', '1PE', 'YF3', 'PEF', 'GD', '8PE', 'DKA', 'RB', 'YB', 'GGD', 'SE4', 'LHG',
               'SMO', 'DGD', 'CMO', 'MLI', 'MW2', 'DTT', 'DOD', '7PH', 'PBM', 'AU', 'FOR', 'PSC', 'TG1', 'KAI', '1PG',
               'DGA', 'IR', 'PE4', 'VO4', 'ACN', 'AG', 'MO4', 'OCL', '6UL', 'CHT', 'RHD', 'CPS', 'IR3', 'OC4', 'MTE',
               'HGC', 'CR', 'PC1', 'HC4', 'TEA', 'BOG', 'PEO', 'PE5', '144', 'IUM', 'LMG', 'SQU', 'MMC', 'GOL', 'NVP',
               'AU3', '3PH', 'PT4', 'PGO', 'ICT', 'OCM', 'BCR', 'PG4', 'L4P', 'OPC', 'OXM', 'SQD', 'PQ9', 'BAM', 'PI',
               'PL9', 'P6G', 'IRI', '15P', 'MAE', 'MBO', 'FMT', 'L1P', 'DUD', 'PGV', 'CD1', 'P33', 'DTU', 'XAT', 'CD',
               'THE', 'U1', 'NA', 'MW3', 'BHG', 'Y1', 'OCT', 'BET', 'MPD', 'HTO', 'IBM', 'D01', 'HAI', 'HED', 'CAD',
               'CUZ', 'TLA', 'SO4', 'OC5', 'ETF', 'MRD', 'PT', 'PHB', 'URE', 'MLA', 'TGL', 'PLM', 'NET', 'LAC', 'AUC',
               'UNX', 'GA', 'DMS', 'MO2', 'LA', 'NI', 'TE', 'THJ', 'NHE', 'HAE', 'MO1', 'DAO', '3PE', 'LMU', 'DHJ',
               'FLC', 'SAL', 'GAI', 'ORO', 'HEZ', 'TAM', 'TRA', 'NEX', 'CXS', 'LCP', 'HOH', 'OCN', 'PER', 'ACY', 'MH2',
               'ARS', '12P', 'L3P', 'PUT', 'IN', 'CS', 'NAW', 'SB', 'GUN', 'SX', 'CON', 'C2O', 'EMC', 'BO4', 'BNG',
               'MN5', '__O', 'K', 'CYN', 'H2S', 'MH3', 'YT3', 'P22', 'KO4', '1AG', 'CE', 'IPL', 'PG6', 'MO5', 'F09',
               'HO', 'AL', 'TRS', 'EOH', 'GCP', 'MSE', 'AKR', 'NCO', 'PO4', 'L2P', 'LDA', 'SIN', 'DMI', 'SM', 'DTD',
               'SGM', 'DIO', 'PPI', 'DDQ', 'DPO', 'HCA', 'CO5', 'PD', 'OS', 'OH', 'NA6', 'NAG', 'W', 'ENC', 'NA5',
               'LI1', 'P4C', 'GLV', 'DMF', 'ACT', 'BTB', '6PL', 'BGL', 'OF1', 'N8E', 'LMT', 'THM', 'EU3', 'PGR', 'NA2',
               'FOL', '543', '_CP', 'PEK', 'NSP', 'PEE', 'OCO', 'CHD', 'CO2', 'TBU', 'UMQ', 'MES', 'NH4', 'CD5', 'HTG',
               'DEP', 'OC1', 'KDO', '2PE', 'PE3', 'IOD', 'NDG', 'CL', 'HG', 'F', 'XE', 'TL', 'BA', 'LI', 'BR', 'TAU',
               'TCA', 'SPD', 'SPM', 'SAR', 'SUC', 'PAM', 'SPH', 'BE7', 'P4G', 'OLC', 'OLB', 'LFA', 'D10', 'D12', 'DD9',
               'HP6', 'R16', 'PX4', 'TRD', 'UND', 'FTT', 'MYR', 'RG1', 'IMD', 'DMN', 'KEN', 'C14', 'UPL', 'CMJ', 'ULI',
               'MYS', 'TWT', 'M2M', 'P15', 'PG0', 'PEU', 'AE3', 'TOE', 'ME2', 'PE8', '6JZ', '7PE', 'P3G', '7PG', 'PG5',
               '16P', 'XPE', 'PGF', 'AE4', '7E8', '7E9', 'MVC', 'TAR', 'DMR', 'LMR', 'NER', '02U', 'NGZ', 'LXB', 'A2G',
               'BM3', 'NAA', 'NGA', 'LXZ', 'PX6', 'PA8', 'LPP', 'PX2', 'MYY', 'PX8', 'PD7', 'XP4', 'XPA', 'PEV', '6PE',
               'PEX', 'PEH', 'PTY', 'YB2', 'PGT', 'CN3', 'AGA', 'DGG', 'CD4', 'CN6', 'CDL', 'PG8', 'MGE', 'DTV', 'L44',
               'L2C', '4AG', 'B3H', '1EM', 'DDR', 'I42', 'CNS', 'PC7', 'HGP', 'PC8', 'HGX', 'LIO', 'PLD', 'PC2', 'PCF',
               'MC3', 'P1O', 'PLC', 'PC6', 'HSH', 'BXC', 'HSG', 'DPG', '2DP', 'POV', 'PCW', 'GVT', 'CE9', 'CXE', 'C10',
               'CE1', 'SPJ', 'SPZ', 'SPK', 'SPW', 'HT3', 'HTH', '2OP', '3NI', 'BO3', 'DET', 'D1D', 'SWE', 'SOG']

'NAG' in biolip_list
'NGA' in biolip_list

comps = '/Users/dmattox/cbk/allGlycans/structures/components_quotes.csv'
comp = {}
with open(comps, 'r') as inFH:
    for line in inFH:
        line = line.split('","')
        comp[line[0][1:]] = [line[1], line[2][:-2]]
        
print(len(biolip_list))
for c in biolip_list:
    if c in comp.keys():
        if 'SACCH' in comp[c][1] or 'sacch' in comp[c][1]:
            biolip_list.remove(c)
            print(c,comp[c])
        if 'PYRAN-3,4,5-TRIOL' in comp[c][0]:
            print('PYRAN-3,4,5-TRIOL','\n\t', c, comp[c])
print(len(biolip_list))

# with open('/Users/dmattox/cbk/glycan_binding/data/comps.txt', 'w') as outFH:
#     for c in biolip_list:
#         if c in comp.keys():
#             outFH.write(c+'\n')
#             outFH.write('\t' + comp[c][0]+'\n')
#             outFH.write('\t\t' + comp[c][1]+'\n')

excluded = []
for f in plipFiles:
    plipDat = px.PLIPXML(f)
    if len(plipDat.excluded) > 0:
        for lig in plipDat.excluded:
            if lig not in excluded:
                excluded.append(lig)
len(excluded)
for e in excluded:
    if e in comp.keys():
        print(e,comp[e])




######################
# Ligand identities
noLig = []
ligs = []
ligMap = {}
for pdb in unilecData.keys():
    
    plipIN = [f for f in plipFiles if f[-8:-4] == pdb]
    pdbIN = [f for f in pdbFiles if f[-8:-4] == pdb]
    
    if (len(plipIN) != 1) or (len(pdbIN) != 1):
        print("ERROR: ambigious or missing match to provided files for PDB ID in UniLectin data: " + pdb)
    else:
        plipIN = plipIN[0]
        pdbIN = pdbIN[0]
    
    plipDat = px.PLIPXML(plipIN)
    
    lig = ''
    bestmass = 115
    for site in plipDat.bsites:
        bs = plipDat.bsites[site]
        if bs.counts['total'] > 0 and bs.molweight > bestmass:
            lig = site
            bestmass = bs.molweight
    if lig == '':
        noLig.append((pdb, bestmass))
        continue
    else:
        ligMap[pdb] = lig[:3]
        if lig[:3] not in ligs: ligs.append(lig[:3])

print(noLig)
len(ligs)
print(ligs)

comps = '/Users/dmattox/cbk/allGlycans/structures/components_quotes.csv'
comp = {}
with open(comps, 'r') as inFH:
    for line in inFH:
        line = line.split('","')
        comp[line[0][1:]] = [line[1], line[2][:-2]]

ions = []
for c in comp:
    if ' ION' in comp[c][0]:
        ions.append(c)

for lig in ligs:
    if lig in comp.keys():
        if 'SACCH' in comp[lig][1] or 'sacch' in comp[lig][1]:
            print(lig,comp[lig])
        else:
            print('NOT A SACCH?')
            print('\t',lig,comp[lig])
    else:
        print('NOT IN COMP FILE')
        print('\t',lig)

check = '1PG'
for f in ligMap.keys():
    if check in ligMap[f]:
        print(f,ligMap[f])



for bs in bsIDs:
    print('\n\tLig: bs')
    bs = plipDat.bsites[bs]
    print(accBS.longname, bs.longname)
    print(accBS.ligtype, bs.ligtype)
    print(accBS.position, bs.position)
    print(accBS.chain, bs.chain)
    print(accBS.members, bs.members)
    print(accBS.composite, bs.composite)
    




# Print residues in binding site to select in pymol
print('+'.join([str(r['resnr']) for r in allPLIP[pdb].bsites[bs].bs_res if r['min_dist'] <= 15]))


missedLinks = collections.defaultdict(list)

for pdb in allPLIP.keys():
    plipDat = allPLIP[pdb]
    
    ligs = [] # List of lists of members of each PLIP binding site
    for site in plipDat.bsites:
        ligs.append(plipDat.bsites[site].members)
    allLigs = [lig for lst in ligs for lig in lst]

    for l1,l2 in plipDat.covlinks:
        for lig in plipDat.bsites:
            siteMembers = plipDat.bsites[lig].members
            if ((l1 in siteMembers) != (l2 in siteMembers)) and l1 in allLigs and l2 in allLigs:
                missedLinks[pdb].append((l1,l2))
            
    
with open("/Users/dmattox/pdbList.txt", "w") as ofh:
    ofh.write(" ")
    for pdb in allPLIP.keys():
        ofh.write(pdb + " ")
    

def recurObjPrint(obj):
    for attr, value in obj.__dict__.items():
        print("Attribute: " + str(attr or ""))
        print("\t\tAttrType: " + str(type(attr) or ""))
        print("\tValue: " + str(value or ""))
        print("\t\tValueType: " + str(type(value) or ""))
        print()
        #recurObjPrint(value)


####
[proSearcher.search(atm.coord, 1.6, level = 'R') for atm in model['L'][133].get_atoms()]


for pdb in dsspErrs.keys():
    pdbIN = [f for f in pdbFiles if f[-8:-4] == pdb][0]
    parser = Bio.PDB.PDBParser(QUIET=True)
    model = parser.get_structure(pdb, pdbIN)[0] # grab the first model in PDB file
    for r in dsspErrs[pdb]:
        n = int(r[:-1])
        c = r[-1]
        print(model[c][n])


###############################
# Example PLIP binding site
        
  </bindingsite>
  <bindingsite has_interactions="True" id="2">
    <identifiers>
      <longname>LBT</longname>
      <ligtype>SMALLMOLECULE</ligtype>
      <hetid>LBT</hetid>
      <chain>B</chain>
      <position>1300</position>
      <composite>False</composite>
      <members>
        <member id="1">LBT:B:1300</member>
      </members>
      <smiles>OC[C@H]1O[C@H](O)[C@@H]([C@H]([C@@H]1O[C@@H]1O[C@H](CO)[C@@H]([C@@H]([C@H]1O)O)O)O)O</smiles>
      <inchikey>GUBGYTABKSRVRQ-XLOQQCSPSA-N
</inchikey>
    </identifiers>
    <lig_properties>
      <num_heavy_atoms>23</num_heavy_atoms>
      <num_hbd>8</num_hbd>
      <num_unpaired_hbd>3</num_unpaired_hbd>
      <num_hba>11</num_hba>
      <num_unpaired_hba>4</num_unpaired_hba>
      <num_hal>0</num_hal>
      <num_unpaired_hal>0</num_unpaired_hal>
      <num_aromatic_rings>0</num_aromatic_rings>
      <num_rotatable_bonds>4</num_rotatable_bonds>
      <molweight>342.29648</molweight>
      <logp>-5.3972</logp>
    </lig_properties>
    <interacting_chains>
      <interacting_chain id="1">B</interacting_chain>
    </interacting_chains>
    <bs_residues>
      <bs_residue aa="ARG" contact="True" id="1" min_dist="2.8">64B</bs_residue>
      <bs_residue aa="VAL" contact="False" id="2" min_dist="4.1">72B</bs_residue>
      <bs_residue aa="THR" contact="False" id="3" min_dist="6.4">75B</bs_residue>
      <bs_residue aa="TRP" contact="True" id="4" min_dist="3.7">81B</bs_residue>
      <bs_residue aa="ARG" contact="True" id="5" min_dist="4.3">43B</bs_residue>
      <bs_residue aa="GLY" contact="False" id="6" min_dist="5.1">82B</bs_residue>
      <bs_residue aa="PRO" contact="False" id="7" min_dist="6.4">83B</bs_residue>
      <bs_residue aa="ASN" contact="True" id="8" min_dist="3.0">62B</bs_residue>
      <bs_residue aa="ASN" contact="True" id="9" min_dist="2.9">74B</bs_residue>
      <bs_residue aa="CYS" contact="False" id="10" min_dist="5.9">73B</bs_residue>
      <bs_residue aa="GLU" contact="True" id="11" min_dist="2.4">84B</bs_residue>
      <bs_residue aa="GLU" contact="False" id="12" min_dist="7.0">85B</bs_residue>
      <bs_residue aa="HIS" contact="False" id="13" min_dist="2.8">60B</bs_residue>
      <bs_residue aa="LYS" contact="False" id="14" min_dist="7.0">76B</bs_residue>
      <bs_residue aa="VAL" contact="False" id="15" min_dist="5.6">45B</bs_residue>
      <bs_residue aa="ARG" contact="True" id="16" min_dist="3.1">86B</bs_residue>
      <bs_residue aa="GLU" contact="False" id="17" min_dist="7.2">66B</bs_residue>
      <bs_residue aa="TYR" contact="False" id="18" min_dist="5.2">70B</bs_residue>
    </bs_residues>
    <interactions>
      <hydrophobic_interactions />
      <hydrogen_bonds>
        <hydrogen_bond id="1">
          <resnr>62</resnr>
          <restype>ASN</restype>
          <reschain>B</reschain>
          <resnr_lig>1300</resnr_lig>
          <restype_lig>LBT</restype_lig>
          <reschain_lig>B</reschain_lig>
          <sidechain>True</sidechain>
          <dist_h-a>2.54</dist_h-a>
          <dist_d-a>2.97</dist_d-a>
          <don_angle>106.52</don_angle>
          <protisdon>True</protisdon>
          <donoridx>1705</donoridx>
          <donortype>Nam</donortype>
          <acceptoridx>2457</acceptoridx>
          <acceptortype>O3</acceptortype>
          <ligcoo>
            <x>25.578</x>
            <y>3.129</y>
            <z>49.094</z>
          </ligcoo>
          <protcoo>
            <x>23.484</x>
            <y>5.236</y>
            <z>49.058</z>
          </protcoo>
        </hydrogen_bond>
        <hydrogen_bond id="2">
          <resnr>64</resnr>
          <restype>ARG</restype>
          <reschain>B</reschain>
          <resnr_lig>1300</resnr_lig>
          <restype_lig>LBT</restype_lig>
          <reschain_lig>B</reschain_lig>
          <sidechain>True</sidechain>
          <dist_h-a>1.89</dist_h-a>
          <dist_d-a>2.78</dist_d-a>
          <don_angle>149.23</don_angle>
          <protisdon>True</protisdon>
          <donoridx>1722</donoridx>
          <donortype>Ng+</donortype>
          <acceptoridx>2470</acceptoridx>
          <acceptortype>O3</acceptortype>
          <ligcoo>
            <x>22.262</x>
            <y>-0.280</y>
            <z>46.584</z>
          </ligcoo>
          <protcoo>
            <x>21.206</x>
            <y>2.250</y>
            <z>46.126</z>
          </protcoo>
        </hydrogen_bond>
        <hydrogen_bond id="3">
          <resnr>64</resnr>
          <restype>ARG</restype>
          <reschain>B</reschain>
          <resnr_lig>1300</resnr_lig>
          <restype_lig>LBT</restype_lig>
          <reschain_lig>B</reschain_lig>
          <sidechain>True</sidechain>
          <dist_h-a>2.23</dist_h-a>
          <dist_d-a>3.07</dist_d-a>
          <don_angle>142.29</don_angle>
          <protisdon>True</protisdon>
          <donoridx>1723</donoridx>
          <donortype>Ng+</donortype>
          <acceptoridx>2457</acceptoridx>
          <acceptortype>O3</acceptortype>
          <ligcoo>
            <x>25.578</x>
            <y>3.129</y>
            <z>49.094</z>
          </ligcoo>
          <protcoo>
            <x>23.046</x>
            <y>2.503</y>
            <z>47.481</z>
          </protcoo>
        </hydrogen_bond>
        <hydrogen_bond id="4">
          <resnr>74</resnr>
          <restype>ASN</restype>
          <reschain>B</reschain>
          <resnr_lig>1300</resnr_lig>
          <restype_lig>LBT</restype_lig>
          <reschain_lig>B</reschain_lig>
          <sidechain>True</sidechain>
          <dist_h-a>1.90</dist_h-a>
          <dist_d-a>2.86</dist_d-a>
          <don_angle>165.37</don_angle>
          <protisdon>True</protisdon>
          <donoridx>1800</donoridx>
          <donortype>Nam</donortype>
          <acceptoridx>2459</acceptoridx>
          <acceptortype>O3</acceptortype>
          <ligcoo>
            <x>24.769</x>
            <y>-0.826</y>
            <z>50.358</z>
          </ligcoo>
          <protcoo>
            <x>24.537</x>
            <y>-2.528</y>
            <z>52.646</z>
          </protcoo>
        </hydrogen_bond>
        <hydrogen_bond id="5">
          <resnr>84</resnr>
          <restype>GLU</restype>
          <reschain>B</reschain>
          <resnr_lig>1300</resnr_lig>
          <restype_lig>LBT</restype_lig>
          <reschain_lig>B</reschain_lig>
          <sidechain>True</sidechain>
          <dist_h-a>1.69</dist_h-a>
          <dist_d-a>2.62</dist_d-a>
          <don_angle>167.25</don_angle>
          <protisdon>True</protisdon>
          <donoridx>1880</donoridx>
          <donortype>O3</donortype>
          <acceptoridx>2459</acceptoridx>
          <acceptortype>O3</acceptortype>
          <ligcoo>
            <x>24.769</x>
            <y>-0.826</y>
            <z>50.358</z>
          </ligcoo>
          <protcoo>
            <x>22.665</x>
            <y>-1.122</y>
            <z>48.831</z>
          </protcoo>
        </hydrogen_bond>
        <hydrogen_bond id="6">
          <resnr>86</resnr>
          <restype>ARG</restype>
          <reschain>B</reschain>
          <resnr_lig>1300</resnr_lig>
          <restype_lig>LBT</restype_lig>
          <reschain_lig>B</reschain_lig>
          <sidechain>True</sidechain>
          <dist_h-a>2.35</dist_h-a>
          <dist_d-a>3.10</dist_d-a>
          <don_angle>132.78</don_angle>
          <protisdon>True</protisdon>
          <donoridx>1900</donoridx>
          <donortype>Ng+</donortype>
          <acceptoridx>2470</acceptoridx>
          <acceptortype>O3</acceptortype>
          <ligcoo>
            <x>22.262</x>
            <y>-0.280</y>
            <z>46.584</z>
          </ligcoo>
          <protcoo>
            <x>19.231</x>
            <y>-0.778</y>
            <z>46.157</z>
          </protcoo>
        </hydrogen_bond>
      </hydrogen_bonds>
      <water_bridges>
        <water_bridge id="1">
          <resnr>43</resnr>
          <restype>ARG</restype>
          <reschain>B</reschain>
          <resnr_lig>1300</resnr_lig>
          <restype_lig>LBT</restype_lig>
          <reschain_lig>B</reschain_lig>
          <dist_a-w>3.40</dist_a-w>
          <dist_d-w>2.92</dist_d-w>
          <don_angle>165.57</don_angle>
          <water_angle>93.07</water_angle>
          <protisdon>True</protisdon>
          <donor_idx>1545</donor_idx>
          <donortype>Ng+</donortype>
          <acceptor_idx>2455</acceptor_idx>
          <acceptortype>O3</acceptortype>
          <water_idx>2858</water_idx>
          <ligcoo>
            <x>27.543</x>
            <y>1.836</y>
            <z>45.512</z>
          </ligcoo>
          <protcoo>
            <x>24.232</x>
            <y>4.822</y>
            <z>44.651</z>
          </protcoo>
          <watercoo>
            <x>24.648</x>
            <y>2.073</y>
            <z>43.748</z>
          </watercoo>
        </water_bridge>
        <water_bridge id="2">
          <resnr>62</resnr>
          <restype>ASN</restype>
          <reschain>B</reschain>
          <resnr_lig>1300</resnr_lig>
          <restype_lig>LBT</restype_lig>
          <reschain_lig>B</reschain_lig>
          <dist_a-w>3.53</dist_a-w>
          <dist_d-w>2.93</dist_d-w>
          <don_angle>123.04</don_angle>
          <water_angle>90.47</water_angle>
          <protisdon>False</protisdon>
          <donor_idx>2457</donor_idx>
          <donortype>O3</donortype>
          <acceptor_idx>1704</acceptor_idx>
          <acceptortype>O2</acceptortype>
          <water_idx>2860</water_idx>
          <ligcoo>
            <x>25.578</x>
            <y>3.129</y>
            <z>49.094</z>
          </ligcoo>
          <protcoo>
            <x>23.581</x>
            <y>7.459</y>
            <z>49.389</z>
          </protcoo>
          <watercoo>
            <x>25.666</x>
            <y>5.482</y>
            <z>47.342</z>
          </watercoo>
        </water_bridge>
        <water_bridge id="3">
          <resnr>81</resnr>
          <restype>TRP</restype>
          <reschain>B</reschain>
          <resnr_lig>1300</resnr_lig>
          <restype_lig>LBT</restype_lig>
          <reschain_lig>B</reschain_lig>
          <dist_a-w>4.00</dist_a-w>
          <dist_d-w>3.16</dist_d-w>
          <don_angle>168.28</don_angle>
          <water_angle>87.44</water_angle>
          <protisdon>True</protisdon>
          <donor_idx>1855</donor_idx>
          <donortype>Nar</donortype>
          <acceptor_idx>2456</acceptor_idx>
          <acceptortype>O3</acceptortype>
          <water_idx>2924</water_idx>
          <ligcoo>
            <x>27.899</x>
            <y>3.688</y>
            <z>47.683</z>
          </ligcoo>
          <protcoo>
            <x>29.543</x>
            <y>-0.618</y>
            <z>49.948</z>
          </protcoo>
          <watercoo>
            <x>31.114</x>
            <y>1.334</y>
            <z>48.030</z>
          </watercoo>
        </water_bridge>
      </water_bridges>
      <salt_bridges>
        <salt_bridge id="1">
          <resnr>64</resnr>
          <restype>ARG</restype>
          <reschain>B</reschain>
          <resnr_lig>1300</resnr_lig>
          <restype_lig>LBT</restype_lig>
          <reschain_lig>B</reschain_lig>
          <dist>4.25</dist>
          <protispos>True</protispos>
          <lig_group>Carboxylate</lig_group>
          <lig_idx_list>
            <idx id="1">2458</idx>
            <idx id="2">2454</idx>
          </lig_idx_list>
          <ligcoo>
            <x>25.272</x>
            <y>0.347</y>
            <z>46.871</z>
          </ligcoo>
          <protcoo>
            <x>21.955</x>
            <y>3.011</y>
            <z>46.921</z>
          </protcoo>
        </salt_bridge>
      </salt_bridges>
      <pi_stacks />
      <pi_cation_interactions />
      <halogen_bonds />
      <metal_complexes />
    </interactions>
    <mappings>
      <smiles_to_pdb>1:2470,2:2465,3:2464,4:2469,5:2460,6:2466,7:2461,8:2462,9:2463,10:2454,11:2448,12:2458,13:2452,14:2453,15:2459,16:2451,17:2450,18:2449,19:2455,20:2456,21:2457,22:2468,23:2467</smiles_to_pdb>
    </mappings>
  </bindingsite>