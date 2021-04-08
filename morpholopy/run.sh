#!/bin/bash

# Runs MorpholoPy using the following variables:
python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/35_SF_efficiency_0p003_L006N0188 \
-s 23 \
-n SFeff0p003 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/35_SF_efficiency_0p003_L006N0188/output_023

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/11_fkin_0p1_L006N0188/ \
-s 23 23 \
-n Reference fkin0p1 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_fkin_01_11

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/35_SF_efficiency_0p003_L006N0188 \
-s 23 23 \
-n SFeff0p003 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_SFeff_01_35

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/36_SF_efficiency_0p001_L006N0188 \
-s 23 23 \
-n Reference SFeff0p001 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_SFeff_01_36

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/36_SF_efficiency_0p001_L006N0188 \
-s 23 \
-n SFeff0p001 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/36_SF_efficiency_0p001_L006N0188/output_023

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/37_same_gr_softening_L006N0188 \
-s 23 23 \
-n Reference SameGrSoftening \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_DM_softening_01_37

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/38_AGN_dT_3e8_L006N0188 \
-s 23 23 \
-n Reference AGNdT3e8 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_AGNdT_01_38
