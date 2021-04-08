#!/bin/bash

# Runs MorpholoPy using the following variables:
python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/19_thermalSNE_2e51_L006N0188  \
-s 23 \
-n SNE2e51 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/19_thermalSNE_2e51_L006N0188/output_023

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/19_thermalSNE_2e51_L006N0188  \
-s 23 23 \
-n Reference SNE2e51 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_SNE_01_19

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/20_thermalSNE_3e51_L006N0188 \
-s 23 \
-n SNE3e51 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/20_thermalSNE_3e51_L006N0188/output_023

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/20_thermalSNE_3e51_L006N0188 \
-s 23 23 \
-n Reference SNE3e51 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_SNE_01_20
