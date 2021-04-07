#!/bin/bash

# Runs MorpholoPy using the following variables:
python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/03_SNE_4e51_L006N0188  \
-s 23 \
-n SNE4e51 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/03_SNE_4e51_L006N0188/output_023

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/03_SNE_4e51_L006N0188  \
-s 23 23 \
-n Reference SNE4e51 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_SNE_01_03

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/05_delay_40Myr_L006N0188 \
-s 23 \
-n SNTimeDelay \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/05_delay_40Myr_L006N0188/output_023

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/05_delay_40Myr_L006N0188 \
-s 23 23 \
-n Reference SNTimeDelay \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_SNE_TIME_DELAY_01_05