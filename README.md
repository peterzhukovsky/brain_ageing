# Overview
A set of scripts used to analyse brain ageing using morphometric similarity networks (MSN).

# network_con_hphi.m
A script that constructs MSN from freesurfer-derived label statistics. The script outputs mean regional MS and also the full 360x360 matrix of pairwise similarity values.

# network_comp.m
A script that compares groups using independent samples t-tests with FDR and also uses permutation testing as alternative to multiple comparison correction. PFDR<0.01 Briefly, 60 regions were significant at Pperm<0.01 and 43 regions were significant at Pfdr<0.01, suggesting that FDR correction at the 0.01 level is more stringent. 

# dev_trajectories.m
A script constructing developmental "trajectories" from cross-sectional NKI-RS data. It uses these quadratic models to create residuals for "hold-out" data for individuals with stimulant use disorder (cases) and healthy controls. These residuals are then used to test for singificant deviation in cases and controls from the "normative" NKI-RS data. While the model provides good fit for controls, cases show significant deviations (PFDR<0.01); for many regions, the direction of the deviation is the same as the direction of the trajectory. Such MS deviations are interpreted as suggesting accelerated ageing of those brain regions. 

# network_PAL.m
Partial least squares analysis linking regional MS to paired associates learning (PAL) and age using *plsregress* and permutation testing following **Morgan et al 2019 (PNAS)**; see also  
https://github.com/SarahMorgan/Morphometric_Similarity_SZ/blob/master/Gene_analyses.md

# Results visualization:
**network_viz_fs.m**

**network_viz_fs2.m**

**network_viz_fs3_PLS.m**

These scripts visualize the resulting brain maps by creating a freesurfer color file similar to
https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT

Where the colors are representing the statistic of interest (average MS or case-control t-statistics for MS). Can be adapted to work with any statistics as long as they are in the Glasser parcellation. 

Note that a different order of the HCP 360 labels (Glasser et al 2016) was used in the analyses from the order of the regions in e.g. freesurfer aparc stats: HOA_VOL_SURF_GLASSER_labels.xlsx provides a mapping of these ids. It also includes a coarse mapping to Harvard-Oxford atlas regions based on manual annotation.
