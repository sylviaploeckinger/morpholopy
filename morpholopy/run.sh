#!/bin/bash

# Runs MorpholoPy using the following variables:
python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/07_BHboost_m0p5dex_L006N0188 \
-s 23 \
-n BHboost \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/07_BHboost_m0p5dex_L006N0188/output_023

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/07_BHboost_m0p5dex_L006N0188  \
-s 23 23 \
-n Reference BHboost \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_BH_BOOST_01_07

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/33_BHboost_densdep_beta0p5_L006N0188 \
-s 23 \
-n BHboost \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/33_BHboost_densdep_beta0p5_L006N0188/output_023

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/33_BHboost_densdep_beta0p5_L006N0188 \
-s 23 23 \
-n Reference BHboostDensDep \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_BH_BOOST_INCL_DENS_DEP_01_33
