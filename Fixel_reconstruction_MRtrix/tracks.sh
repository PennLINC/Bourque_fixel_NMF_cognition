#!/bin/bash
export MRTRIX_NTHREADS=$1

tckgen \
	-angle 22.5 \
	-maxlen 250 \
	-minlen 10 \
	-power 1.0 \
	/work/data/MRtrix_recon/template/wmfod_template_new.mif \
	-seed_image /work/data/MRtrix_recon/template/template_mask_new.mif \
	-mask /work/data/MRtrix_recon/template/template_mask_new.mif \
	-select 100000 \
	-cutoff 0.06 \
	/work/data/MRtrix_recon/template/tracks_100k_new.tck \
	-force
 
tckgen \
        -angle 22.5 \
        -maxlen 250 \
        -minlen 10 \
        -power 1.0 \
        /work/data/MRtrix_recon/template/wmfod_template_new.mif \
        -seed_image /work/data/MRtrix_recon/template/template_mask_new.mif \
        -mask /work/data/MRtrix_recon/template/template_mask_new.mif \
        -select 20000000 \
        -cutoff 0.06 \
        /work/data/MRtrix_recon/template/tracks_20million_new.tck \
	-force

#Reduce biases in tractogram
tcksift \
	/work/data/MRtrix_recon/template/tracks_20million_new.tck \
	/work/data/MRtrix_recon/template/wmfod_template_new.mif \
	/work/data/MRtrix_recon/template/tracks_2million_new_sift.tck \
	-term_number 2000000 \
	-force 



