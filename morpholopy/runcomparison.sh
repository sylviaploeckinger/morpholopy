#!/bin/bash

# Runs MorpholoPy using the following variables:
#python morpholopy.py \
#-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188_CHIMES \
#-s 23 23 \
#-n Reference Reference+CHIMES \
#-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/01_reference_L006N0188/comparison

python -i morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/39_only_UVB_L006N0188 \
-s 23 23 \
-n Reference Reference+OnlyUVB \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_only_UVB_01_39

