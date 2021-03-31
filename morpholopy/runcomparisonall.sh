#!/bin/bash

# Runs MorpholoPy using the following variables:
python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188_CHIMES \
-s 23 23 \
-n Reference Reference+CHIMES \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_TABLE_VS_CHIMES

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/09_steepEOS_1e3_pressurelaw_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/10_steepEOS_1e3_schmidtlaw_L006N0188 \
-s 23 23 \
-n 09_steepEOS_1e3_pressurelaw 10_steepEOS_1e3_schmidtlaw \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_EOS_VARY_STAR_FORMATION_LAW_09_10

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/29_SNE_densdep_2_4_L006N018 \
-s 23 23 \
-n Reference 29_SNE_densdep_2_4 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_E_SNII_DENS_DEP_VS_REF_01_29

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/31_SNE_metdep_2_4_L006N0188  \
-s 23 23 \
-n Reference 31_SNE_metdep_2_4 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_E_SNII_Z_DEP_VS_REF_01_31

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/15_alphavir_0p5_L006N0188 \
-s 23 23 \
-n Reference 15_alphavir_0p5 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_ALPHA_VIRIAL_01_15

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/40_alphavir_0p25_L006N0188 \
-s 23 23 \
-n Reference 40_alphavir_0p25 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_ALPHA_VIRIAL_01_40
