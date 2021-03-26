#!/bin/bash

# Runs MorpholoPy using the following variables:

folder="/snap7/scratch/dp004/dc-chai1/my_cosmological_box/XMAS2020_L012N376_FKIN03_NOEOS_VKICK050SLOPE00NORM00_FIXEDDELAY_BOOST01"
output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/XMAS2020_L012N376_FKIN03_NOEOS_VKICK050SLOPE00NORM00_FIXEDDELAY_BOOST01/output_0823"
snap="823"
python morpholopy.py -d=$folder -n=$snap -o=$output

#folder="/cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188"
#output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/01_reference_L006N0188/output_023"
#snap="23"
#python morpholopy.py -d=$folder -n=$snap -o=$output


#folder="/Users/camila/Dropbox/Science-projects/swift-COLIBRE/morphology_estimators/data"
#output="/Users/camila/Dropbox/Science-projects/swift-COLIBRE/morphology_estimators/data/plots"
#snap="34"
#python morpholopy.py -d=$folder -n=$snap -o=$output


