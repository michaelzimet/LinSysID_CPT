<a href="https://doi.org/10.5281/zenodo.18048428"><img src="https://zenodo.org/badge/1122413425.svg" alt="DOI"></a>

# LinSysID_CPT
Linear System Identification of Heart Physiology and Pupilometry Time Series from the Cold Presser Task

This repository contains the Matlab code for performing linear system identification of time series data recorded during the Cold Presser Task.

# Experimental Data

From the experimental paper's data repository (https://ucsb.box.com/s/5i3zwm9jvry2eicxyddgdv7x1ouy8kji), the physiology data are in PHYSIO_MASTER_FINAL.mat and the pupilometry data are in CPT_EYE_Master.mat. The scripts PHYSIO_Plot_Main_Figs_Analyses.m and EYE_Plot_Main_Figs_Analyses.m generates plots of the experimental data with descriptive statistics. These scripts remove subjects with missing or corrupted data, and average over all subjects in the population. These require the resampling toolbox.

# System Identification Scripts

The scripts physio_sysID.m and eye_sysID.m perform linear system idenfication of these time series data. These generate Figs 1 through 6 of the paper, as well as Supplementary Figs SM1 through SM4. These also generate all the data for Tables 1, 2, and SM1 through SM3 as well as the sample model parameters in the Supplement.

# Analyses

The script lin_main.m with the accompanying function func_lin.m performs eigensystem analysis of the recovery trajectories in a planar pupilometry model. This generates Figs 7 and 8.

The script step_main.m with the accompanying function func_step.m applies the state-space models to generate predictions for Impulse and Ramp Response. This generates Fig 9.
