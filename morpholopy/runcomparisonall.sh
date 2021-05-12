#!/bin/bash

# Runs MorpholoPy using the following variables:

#python morpholopy.py \
#-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/39_only_UVB_L006N0188 \
#-s 23 23 \
#-n Reference Reference+OnlyUVB \
#-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_only_UVB_01_39

#python morpholopy.py \
#-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/51_noAGN_L006N0188 \
#-s 23 23 \
#-n Reference noAGN \
#-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_noAGN_01_51
#
#python morpholopy.py \
#-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/47_var_dT_agn_T_floor_1e9_L006N0188 \
#-s 23 23 \
#-n Reference AGNvardT9 \
#-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_AGNdT_01_47
#
#python morpholopy.py \
#-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/46_var_dT_agn_T_floor_1e8_L006N0188 \
#-s 23 23 \
#-n Reference AGNvardT8 \
#-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_AGNdT_01_46
#
#python morpholopy.py \
#-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/48_AGNdT_1e8_L006N0188 \
#-s 23 23 \
#-n Reference AGNdT8 \
#-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_AGNdT_01_48
#
#python morpholopy.py \
#-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/38_AGN_dT_3e8_L006N0188 \
#-s 23 23 \
#-n Reference AGNdT3e8 \
#-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_AGNdT_01_38
#
#python morpholopy.py \
#-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/53_thermalSNE_3e51_no_AGN_L006N0188 \
#-s 23 23 \
#-n Reference thermalSNE3e51 \
#-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_thermalSNE_01_53
#
#python morpholopy.py \
#-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/52_thermalSNE_2e51_no_AGN_L006N0188 \
#-s 23 23 \
#-n Reference thermalSNE2e51 \
#-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_thermalSNE_01_52
#
#python morpholopy.py \
#-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/54_SNdT_1e7_L006N0188 \
#-s 23 23 \
#-n Reference SNdT1e7 \
#-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_SNII_TEMPERATURE_01_54
#
#python morpholopy.py \
#-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/27_SNdT_1e8_L006N0188 \
#-s 23 23 \
#-n Reference SNdT1e8 \
#-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_SNII_TEMPERATURE_01_27
#
#python morpholopy.p \
#-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/11_fkin_0p1_L006N0188 \
#-s 23 23 \
#-n Reference fkin0p1 \
#-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_fkin_01_11
#
#python morpholopy.py \
#-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/20_thermalSNE_3e51_L006N0188 \
#-s 23 23 \
#-n Reference SNE3e51 \
#-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_SNE_01_20

#python morpholopy.py \
#-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/15_alphavir_0p5_L006N0188 \
#-s 23 23 \
#-n Reference Alphavir0p5 \
#-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_ALPHA_VIRIAL_01_15
#
#python morpholopy.py \
#-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/40_alphavir_0p25_L006N0188 \
#-s 23 23 \
#-n Reference Alphavir0p25 \
#-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_ALPHA_VIRIAL_01_40
#
#python morpholopy.py \
#-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/59_alphavir_2p00_L006N0188 \
#-s 23 23 \
#-n Reference Alphavir2p0 \
#-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_ALPHA_VIRIAL_01_59
#
#python morpholopy.py \
#-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/35_SF_efficiency_0p003_L006N0188 \
#-s 23 23 \
#-n Reference SFeff0p003 \
#-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_SFeff_01_35
#
#python morpholopy.py \
#-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/36_SF_efficiency_0p001_L006N0188 \
#-s 23 23 \
#-n Reference SFeff0p001 \
#-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_SFeff_01_36
#
#python morpholopy.py \
#-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/49_SF_efficiency_0p03_L006N0188 \
#-s 23 23 \
#-n Reference SFeff0p03 \
#-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_VARY_SFeff_01_49

#python morpholopy.py \
#-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/39_only_UVB_L006N0188 \
#-s 19 19 \
#-n Reference OnlyUVB \
#-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z1.0_only_UVB_01_39

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/extra02_SNE_3e51_L012N0376 \
-s 23 23 \
-n Reference L12extra02SNE3e51 \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/comparisons/z0.0_extra02_SNE_3e51
