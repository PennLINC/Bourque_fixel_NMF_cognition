# NMF_fixel_psychopatho

This repository is divided into three (3) sets of scripts. Preprocessing takes the raw diffusion data to Fixel Density Cross-Section (FDC) metric using the QSIPrep preprocessing pipeline and the Fixel-Based Analysis pipeline from MRtrix3. NMF includes the steps required for the non-negative matrix factorization of the FDC metric. The figures and GAMs folder includes all scripts to generate the paper's figures and investigate the developmental and cognitive effects of the WM covariance networks delineated with NMF. Details below.

This project includes typically-developing PNC participants - LTN criteria was used. 

# Preprocessing

### 0.Run_QSIPrep
Raw DWI -> Denoising procedure; correction for motion, Gibb's ringing artifact, distortion with a B0 field map, and eddy current; coregistration to T1; quality assurance (mean framewise displacement, and number of bad slices); and voxel upsampling (1.25mm).

### 0.QSIPrep_download
Download the preprocessed DWI output data.


**The following steps to reconstruct FDC follow the recommended FBA pipeline (https://mrtrix.readthedocs.io/en/latest/fixel_based_analysis/mt_fibre_density_cross-section.html)**

### 1.RUN_nii2mif.sh
Convert nii.gz to .mif files (for MRtrix) and update the size of the anatomical image (T1) to match the size of the DWI images.

### 2.RUN_response_function.sh
Calculate the response function for the 3 tissue types (WM, GM, CSF) with single-shell 3-tissue MRtrix fork for our 30 template subjects. 

## 3.RUN_response_mean.sh
Average the response function of the template subjects.

## 4.RUN_FODs.sh
Estimate Fiber Orientation Distributions (FODs) for all individuals based on the tissue types' response functions.

## 5.RUN_normalise.sh
Bias correction and global intensity normalisation of all individual's FOD image.

## 6.RUN_template_creation.sh
Create the population FOD template based on the FOD image of the 30 template subjects.

## 7.RUN_FODreg.sh
Register all subjects FODs to the FOD template.

## 8.RUN_fixel_mask.sh
Create the fixel mask - the intersection of all subject masks in template space.

## 9.RUN_segment_fixel.sh
Segment fixel (calculate fiber density - FD) and reorient and assign fixels to the template fixels. Compute Fibre Cross-Section (FC), and Fiber Density Cross-Section (FDC) - our metric of interest.

## 10.RUN_tracks.sh
Perform whole brain fibre tractography of the FOD template.

## 11.RUN_smoothing.sh
Create the fixel-to-fixel connectivity to then smooth the fixel data accordingly.




