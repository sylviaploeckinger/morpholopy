#!/bin/bash

# Runs MorpholoPy using the following variables:
python galaxy_population.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/47_var_dT_agn_T_floor_1e9_L006N0188 \
-s 23 23 \
-n Reference AGNvardT9 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_AGNdT_01_47

python galaxy_population.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/46_var_dT_agn_T_floor_1e8_L006N0188 \
-s 23 23 \
-n Reference AGNvardT8 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_AGNdT_01_46

python galaxy_population.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/48_AGNdT_1e8_L006N0188 \
-s 23 23 \
-n Reference AGNdT8 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_AGNdT_01_48
