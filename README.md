# Glycan binding

Project to identify and map features of lectin binding sites and associations with affinity or specificity towards particular glycans
---
## Overview
1. Package dependencies
2. Analysis notes

### Package dependencies
The scripts used to conduct this analysis were run with the following software packages:

**R version 3.3.2**
- ggplot2 2.2.1
- reshape2 1.4.3

**Python 3.7.6**
- rdkit 2019.09.3
- lxml 4.5.0
- biopython 1.76

[PLIP v1.4.5 was run on **Python 2.7.17**]
Using custom version of PLIP modified from v1.4.5 that can be found here: https://github.com/demattox/plip
- openbabel 2.4.1

(full list of dependencies to be finalized soon)

### Analysis notes
For the sake of transparency, all abandoned and previous versions of scripts are still included. However, the final analysis pipeline is (roughly) laid out below:
#### Pre-processing/general scripts
- lec_gly.py: helper functions used across python scripts
- ./bashScripts/runPlip.pbs + manualPlip.sh: bash scripts to run PLIP tool on all PDB IDs obtained from UniLectin3D
- ./bashScripts/extractPLIPreports.sh: extract PLIP reports from nest subdirectories, rename them, and place them in a common directory, along with "plipfixed" PDB files
- ./bashScripts/runDssp.pbs + manualDSSP.sh: run DSSP on all PDB/mmCIF files to get secondary structure, phi/psi angles, and solvent accesibility information
- checkPLIPligands.py: massive automated QC script to remove non-glycan ligands, repair/replace missing covalent bonds, exclude N-,O-, and C-linked glycans, recover residue numbering from plip"fixed" PDB files where inserted residues with letters appended to the numbers get set to zero, handle cases of duplicated saccharide residues in the same position (different anomericity), expand the set of residues in the binding site, and attached previously generated DSSP information
- manualBSfilter.py: semi-automated method to flag suspicious interactions (no/few protein residues close to the glycan means it probably isn't really bound by the protein) and provide the opportunity to manual inspect the flagged interactions and exclude unlikely interactions while providing the reason for exclusion to be compiled in the auto-generated QC report
- wholeSequence.py: returns the entire sequence of protein it the PDB file, ordered by chain identifier
- ./misc/getSeqs.sh: bash script to run wholeSequence.py script on every PDB file
- ./misc/runChainLvlCDHIT.sh: run CD-HIT at 90% seq id to collapse redundant chains from the same PDB file
- chainClusts2seqs.py: process output of CD-HIT from above script to write non-redundant protein sequences into the same file
- runProteinLvlCDHIT.sh: bash script to run CD-HIT at varied seq id thresholds on the non-redundant protein sequences, clustering lectins by various levels of homology


#### Generating features
- bSiteResiFeatures.py: generates the file containing residue-base features and features from PLIP report. Bins the binding site residues by distance to the closest heavy glycan atom and return statistics of residue identity, residue class, secondary structure for each bin and each interaction
- ./pocketAnalysis/bsResis2pdbs.py: copy all of the binding site residues to a separate .pdb file
- ./pocketAnalysis/pdb2wrl_surface.py: export van der Waals surface of binding site residues to a .wrl file
