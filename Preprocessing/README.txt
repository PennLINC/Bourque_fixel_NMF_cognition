#Steps from QSIprep outputs to smoothed fixels for creating the NMF components
#Reconstruction with MRtrix3 according to this guideline with a few modification for single-shell 3 t-issue https://mrtrix.readthedocs.io/en/latest/fixel_based_analysis/mt_fibre_density_cross-section.html. This reconstruction workflow uses a lot of independent scripts to run sequentially because at the time of the analysis single-shell 3-tissue was an independent fork of MRtrix3, also some steps require to run on all sample, others on template subjects only. Some MRtrix3 functions used in these scripts may be outdated, a new release has made changes to only a few commands (functions).

#Unless specified all these steps run on the qsiprep-0.4.5.simg image which contains all code from the main MRtrix3 RC3

#1 Convert nii.gz to .mif files and update the size of the anatomical image (T1) to match the size of the DWI images
# 1.RUN_nii2mif.sh

#2 Calculate the response function for the 3 tissue types with ss3t fork for template subjects only
# 2.RUN_response_function.sh

#3 Average the response function
# 3.RUN_response_mean

#4 Estimate FODs based on the tissue types' response functions using the ss3t fork of MRtrix3
# 4.RUN_FODs.sh

#5 Bias correction and global intensity normalisation
# 5.RUN_normalise.sh

#6 Create the population FOD template
# 6.RUN_template_creation.sh

#7 Register all subjects FODs to the FOD template
# 7.RUN_FODreg.sh
 
#8 Create the fixel mask
# 8.RUN_fixel_mask.sh

#9 Segment fixel (calculate fiber density) and reorient and assign fixels to the template fixels. Compute Fibre Cross-Section
# 9.RUN_segment_fixel.sh
**might need to use instead RUN_segment_fixel_sequential.sh script because I ran into problems when running this script in parallel on cubic. 

#10 Perform whole brain fibre tractography of the FOD template
# 10.RUN_tracks.sh

#11 Create the fixel-to-fixel connectivity to then smooth the fixel data accordingly
# 11.RUN_smoothing.sh


