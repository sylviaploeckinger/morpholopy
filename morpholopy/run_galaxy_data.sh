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

python galaxy_population.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/38_AGNdT_3e8_L006N0188 \
-s 23 23 \
-n Reference AGNdT3e8 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_AGNdT_01_38

python galaxy_population.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/11_fkin_0p1_L006N0188 \
-s 23 23 \
-n Reference fkin0p1 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_fkin_01_11

python galaxy_population.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/27_SNdT_1e8_L006N0188 \
-s 23 23 \
-n Reference SNdT1e8 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_SNII_TEMPERATURE_01_27

python galaxy_population.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/54_SNdT_1e7_L006N0188 \
-s 23 23 \
-n Reference SNdT1e7 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_SNII_TEMPERATURE_01_54

python galaxy_population.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/05_delay_40Myr_L006N0188 \
-s 23 23 \
-n Reference SNTimeDelay \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_SNE_TIME_DELAY_01_05


python galaxy_population.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/51_noAGN_L006N0188 \
-s 23 23 \
-n Reference noAGN \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_noAGN_01_51

python galaxy_population.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/52_thermalSNE_2e51_no_AGN_L006N0188 \
-s 23 23 \
-n Reference thermalSNE2e51 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_thermalSNE_01_52

python galaxy_population.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/53_thermalSNE_3e51_no_AGN_L006N0188 \
-s 23 23 \
-n Reference thermalSNE3e51 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_thermalSNE_01_53

python galaxy_population.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/19_thermalSNE_2e51_L006N0188 \
-s 23 23 \
-n Reference SNE2e51 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_SNE_01_19

python galaxy_population.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/20_thermalSNE_3e51_L006N0188 \
-s 23 23 \
-n Reference SNE3e51 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_SNE_01_20
