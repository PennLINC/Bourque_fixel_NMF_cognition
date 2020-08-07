#$ -pe threaded 8
#$ -l h_vmem=20G

#!/bin/bash
SIF=/cbica/projects/pnc_fixel_cs/scripts/qsiprep-0.4.5.simg

#Fixelcorrespondence 
singularity exec \
	-B /cbica/projects/pnc_fixel_cs:/work \
	$SIF bash /work/scripts/MRtrix_recon/segment_fixel_sequential.sh

