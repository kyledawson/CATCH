#!/bin/csh
#BSUB -W 12:00
#BSUB -n 4
#BSUB -q cos
source /usr/local/apps/MATLAB/matlab2016a.csh
#BSUB -o /gpfs_backup/meskhidze_data/DISCOVER_AQ/logfiles/output.%J
#BSUB -e /gpfs_backup/meskhidze_data/DISCOVER_AQ/logfiles/error.%J
matlab -nodisplay -r make_aerosol_id_masks

