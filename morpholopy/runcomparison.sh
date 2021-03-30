#!/bin/bash

# Runs MorpholoPy using the following variables:

#python -i morpholopy.py \
#-d /Users/camila/Dropbox/Science-projects/swift-COLIBRE/morphology_estimators/data \
#-s 34 \
#-n Reference \
#-o /Users/camila/Dropbox/Science-projects/swift-COLIBRE/morphology_estimators/data/plots

#python -i morpholopy.py \
#-d /Users/camila/Dropbox/Science-projects/swift-COLIBRE/morphology_estimators/data/Testing_Sims /Users/camila/Dropbox/Science-projects/swift-COLIBRE/morphology_estimators/data/HAWKReferenceL006N188 \
#-s 34 23 \
#-n RandomRun HAWKReference \
#-o /Users/camila/Dropbox/Science-projects/swift-COLIBRE/morphology_estimators/data/Testing_Sims/plots

python morpholopy.py \
-d /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188 /cosma7/data/dp004/dc-chai1/HAWK/01_reference_L006N0188_CHIMES \
-s 23 23 \
-n Reference Reference+CHIMES \
-o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/MorpholoPy/HAWK/01_reference_L006N0188/comparison

