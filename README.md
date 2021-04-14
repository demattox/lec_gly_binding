# Glycan binding

Project to identify and map features of lectin binding sites and associations with affinity or specificity towards particular glycans
---
## Overview
1. Package dependencies
2. Analysis notes

### Package dependencies
The scripts used to conduct this analysis were run with the following software packages:

**R version 4.0.3**
R version 4.0.3 (2020-10-10)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Catalina 10.15.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
  [1] parallel  grid      stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
  [1] RColorBrewer_1.1-2 devtools_2.3.2     usethis_2.0.1      stringr_1.4.0      seqinr_4.2-5       reshape2_1.4.4     corrplot_0.84      ggplot2_3.3.3      umap_0.2.7.0       doParallel_1.0.16  iterators_1.0.13   foreach_1.5.1      protr_1.6-2
[14] Cairo_1.5-12.2     pheatmap_1.0.12    vioplot_0.3.5      zoo_1.8-9          sm_2.2-5.6         philentropy_0.4.0  colorspace_2.0-0   scales_1.1.1       survey_4.0         survival_3.2-10    Matrix_1.3-2

loaded via a namespace (and not attached):
  [1] pkgload_1.2.0      jsonlite_1.7.2     splines_4.0.3      assertthat_0.2.1   askpass_1.1        remotes_2.2.0      progress_1.2.2     sessioninfo_1.1.1  pillar_1.5.1       lattice_0.20-41    glue_1.4.2         reticulate_1.18    plyr_1.8.6
[14] pkgconfig_2.0.3    purrr_0.3.4        processx_3.5.0     RSpectra_0.16-0    tibble_3.1.0       openssl_1.4.3      generics_0.1.0     farver_2.1.0       ellipsis_0.3.1     cachem_1.0.4       withr_2.4.1        cli_2.3.1          magrittr_2.0.1
[27] crayon_1.4.1       memoise_2.0.0      ps_1.6.0           fs_1.5.0           fansi_0.4.2        MASS_7.3-53.1      pkgbuild_1.2.0     tools_4.0.3        prettyunits_1.1.1  hms_1.0.0          mitools_2.4        lifecycle_1.0.0    munsell_0.5.0
[40] callr_3.5.1        ade4_1.7-16        compiler_4.0.3     rlang_0.4.10       tcltk_4.0.3        testthat_3.0.2     gtable_0.3.0       codetools_0.2-18   DBI_1.1.1          R6_2.5.0           dplyr_1.0.5        fastmap_1.1.0      utf8_1.2.1
[53] rprojroot_2.0.2    KernSmooth_2.23-18 desc_1.3.0         stringi_1.5.3      Rcpp_1.0.6         vctrs_0.3.6        tidyselect_1.1.0

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

### Analysis notes
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
- bSiteResiFeatures.py: generates the file containing residue-base features and features from PLIP report. Bins the binding site residues by distance to the closest heavy glycan atom and return statistics of residue identity, residue class, secondary structure for each bin and each interaction
- ./pocketAnalysis/bsResis2pdbs.py: copy all of the binding site residues to a separate .pdb file
- ./pocketAnalysis/pdb2wrl_surface.py: export van der Waals surface of binding site residues to a .wrl file
