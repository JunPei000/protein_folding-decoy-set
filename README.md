# protein_folding-decoy-set
The decoy sets we used include the multiple decoy sets from the Decoys ‘R’ Us collection (http://compbio.buffalo.edu/dd/download.shtml), which include the 4state_reduced, fisa, fisa_casp3, hg_structal, ig_structal, ig_structal_hires, lattice_ssfit, lmds, and lmds_v2 decoy sets. The MOULDER decoy set1 was downloaded from https://salilab.org/decoys/; the I-TASSER decoy set-II was obtained from https://zhanglab.ccmb.med.umich.edu/decoys/decoy2.html; and the ROSETTA all-atom decoy set from https://zenodo.org/record/48780#.WvtCA63MzLF. 

All protein structures (including both native and decoys) were converted into their biological oligomerization state and prepared with the Protein Preparation Wizard, which adds missing atoms, optimizes the H-bond network, and performs energy minimization to clean up the structures for subsequent calculations. 

Please cite:
Random Forest Refinement of the KECSA2 Knowledge-Based Scoring Function for Protein Decoy Detection
Jun Pei, Zheng Zheng, and Kenneth M. Merz, Jr.
Journal of Chemical Information and Modeling Article ASAP
DOI: 10.1021/acs.jcim.8b00734
