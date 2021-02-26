#!/bin/bash

# Runs MorpholoPy using the following variables:

#folder="/cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188"
#output="/cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/01_reference_L006N0188/output_023" 
#snap="23"

folder="/Users/Camila/Dropbox/Science-projects/swift-COLIBRE/morphology_estimators/data"
output="/Users/Camila/Dropbox/Science-projects/swift-COLIBRE/morphology_estimators/data/plots"
snap="34"

python -i morpholopy.py -d=$folder -n=$snap -o=$output


