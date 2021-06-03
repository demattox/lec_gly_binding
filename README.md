# Glycan binding

Project to identify and map features of lectin binding sites and associations with affinity or specificity towards particular glycans
---
## Overview
1. Analysis notes
2. Package dependencies

### Analysis notes
To re-run any portion of this analysis, the first step would be to edit scripts/lec_gly.py at line 10 and change the variable "homeDir" to the full path pointing to where this repository was cloned.

This analysis was primarily run on a Mac (iOS 10.15.6) as well as our local computing cluster running CentOS Linux 7 (release 7.9.2009). It should work with most unix/linux systems put would likely require modification to be compatible with Windows machines.

For the sake of transparency, all abandoned and previous versions of scripts are still included. However, the final analysis pipeline is (roughly) laid out below:
#### Pre-processing/general scripts
- uniLec_survey.R: calculate initial statistics from UniLectin3D data and manually address some formatting issues from Excel or missing labels
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
##### Residue-based features & PLIP features
- bSiteResiFeatures.py: generates the file containing residue-base features and features from PLIP report. Bins the binding site residues by distance to the closest heavy glycan atom and return statistics of residue identity, residue class, secondary structure for each bin and each interaction
##### Pocket features
- genSplitPDBlists.py: Split all PDB IDs into 28 lists and save each list to a separate file to facilitate parallel batch processing of binding sites
- ./pocketAnalysis/bsResis2pdbs.py: copy all of the binding site residues to a separate .pdb file
- ./pocketAnalysis/pdb2wrl_surface.py: export van der Waals surface of binding site residues to a .wrl file
- ./pocketAnalysis/generateVoxels.py: Script to voxelize the space in the binding pocket around the glycan, saving out the centroids of each voxel to mol2 files (working version of included legacy script ./pocketAnalysis/bs3Dcharacterization.py)
- ./misc/runVoxGen.pbs: Run ./pocketAnalysis/generateVoxels.py on the computing cluster
- ./pocketAnalysis/getPocketSurfacePoints.py: Pre-process voxelized pocket representations to split voxel centroids into "buried" and "surface" points, saving the surface points to a separate mol2 file
- ./misc/runSurfacePnts.pbs: Run ./pocketAnalysis/getPocketSurfacePoints.py on the computing cluster
- ./pocketAnalysis/featurizePockets.py: Process voxelized pocket representations, calculate D2 measures (all pairwise distances between surface points) for each threshold of each pocket, and return pocket descriptor features, D2 distribution statistics, and the scaled-bin-size D2 distributions
- ./misc/concatD2dists_scaled.sh & concatD2distFeats: Concatenate raw D2 files from each interaction together into the same file (each version of this script isused for its corresponding D2 output file)
- ./pocketAnalysis/getZernike.py: Translate pocket voxel centroid points to a defined grid centered at (0,0,0) and scale coordinates between [0,0.6] to fit entirely within the unit sphere of radius 1, processed voxel points are saved out as .vox files. This script will then execute **vox** (compiled binary from modifed source code of Daberdaku et al. 2019 previously implemented from https://github.com/sebastiandaberdaku/AntibodyInterfacePrediction) which returns the 3DZD invariants in .inv files. 20th order polynomials were found by default but only descriptors up to the 10th order were used for the final features.
- ./pocketAnalysis/vox: Not called directly but through getZernike.py above, returns a .inv file with the invariants describing the provided shape at the specified order. Accepts 4 positional arguments: (1) .vox file input (2) integer, which order of polynomial to calculate (3) output file name + path (no extension) (4) 0/1, whether or not to reconstruct the input voxels from the calculated polynomials, if 1 an additional .vox file is generated with the same name from (3) as the returned .inv file (warning, reconstruction adds significantly to the runtime)
- ./misc/concat3DZDs.sh: Concatenate all returned .inv files in the same directory into one .tsv file


#### Using features
##### Pre-processing & Feature generation above can be skipped by using the included data files and the below scripts (starting with descriptive.R) to repeat the analysis in this work
- ./prediction/preliminary.R: Read in residue + PLIP features, UniLectin3D-provided info, pocket descriptors, D2 statistics, raw D2 binned features, and 3DZDs, & returns the final feature set (also found **here**) and Booleans to index rows of the data correponding to interactions with the 15 glycans of interest. Along the way, this script was also used to investigate patterns and preliminary trends in the data, look at clustering sequences with identity from 50% to 90% (stayed with the most conservative seq ID level at 50%), & identify the ligands of interest.
- ./prediction/descriptive.R: With final set of features for each interaction, run univariate comparative WMW tests weighted by 50% seq ID homology clusters and investigate associations
- ./prediction/train_validate.R: Build RF models for a provided glycan of interest and run 10x repeated leave-one-(cluster-)out cross-validation, saving out feature importances, training (5x CV), and LO(C)O CV perfromances form. (Run on the cluster with `runTrainValidate.pbs` and jobs are sumbitted with `startSubJobs.sh` to generate nested subdirectories specifying the ligand of interest (1st directory) and the random seed to use for each repeat (2nd directory)
- ./prediction/train_validate.R: Build RF models for a provided glycan of interest and run 10x repeated leave-one-(cluster-)out cross-validation, saving out feature importances, training (5x CV), and LO(C)O CV perfromances form. (Run on the cluster with `./prediction/cluster/runTrainValidate.pbs` and jobs are sumbitted with `./prediction/cluster/startSubJobs.sh` to generate nested subdirectories specifying the ligand of interest (1st directory) and the random seed to use for each repeat (2nd directory)
- ./prediction/train_validate_random.R: Build RF models for a provided glycan of interest **(with those glycan labels shuffled between all interactions)** and run 10x repeated leave-one-(cluster-)out cross-validation, saving out feature importances, training (5x CV), and LO(C)O CV perfromances form. (Run on the cluster with `./prediction/cluster/runTrainValidateRand.pbs` and jobs are sumbitted with `./prediction/cluster/startSubJobs.sh` to generate nested subdirectories specifying the ligand of interest (1st directory) and the random seed to use for each repeat (2nd directory)
- ./prediction/results.R: Process and interpret RF results, combined with statisical associations from the weighted WMW
- fine_specificity-NeuAc23_26.R: Fine specificity analysis of 6' vs 3' NeuAc terminal glycans, including weighted WMW and clustering by binding site seqs or feature interactions.

### Package dependencies
The scripts used to conduct this analysis were run with the following software packages:

**R version 4.0.3**
- RColorBrewer_1.1-2
- devtools_2.3.2
- usethis_2.0.1
- stringr_1.4.0
- seqinr_4.2-5
- reshape2_1.4.4
- corrplot_0.84
- ggplot2_3.3.3
- umap_0.2.7.0
- doParallel_1.0.16
- iterators_1.0.13
- foreach_1.5.1
- protr_1.6-2
- Cairo_1.5-12.2
- pheatmap_1.0.12
- vioplot_0.3.5
- zoo_1.8-9
- sm_2.2-5.6
- philentropy_0.4.0
- colorspace_2.0-0
- scales_1.1.1
- survey_4.0
- survival_3.2-10
- Matrix_1.3-2

**Python 3.8.1**
- rdkit 2019.09.3
- lxml 4.5.0
- biopython 1.78
- numpy 1.19.2
- dill 0.3.2
- scikit-learn 0.23.2
- scipy 1.5.0

[PLIP v1.4.5 was run on **Python 2.7.18**]
Using custom version of PLIP modified from v1.4.5 that can be found here: https://github.com/demattox/plip
- openbabel 2.4.1

(full list of dependencies to be finalized soon)

