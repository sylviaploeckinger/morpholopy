#!/bin/bash


#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDPowerLaw \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDPowerLaw_Fe05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDtau1p5nu1p6_Fe05O05Mg05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDtau2nu1p6 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L025NYSNIIlo8hi40DTDComparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDPowerLaw DTDtauPowerLawFe05 DTDtau1p5nu1p6Fe05O05Mg05 DTDtau2nu1p6 \
#                     -m 5e8


python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_StandardYields_SNIIlo8hi100_DTDtau2nu1p6 \
                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p6 \
                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDtau2nu1p6 \
                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L025CompareYields \
                     -c halo_0036.properties halo_0036.properties halo_0036.properties \
                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
                     -n StandardYields NewYieldsSNIIlo8hi100 NewYieldsSNIIlo8hi40 \
                     -m 5e8 

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDPowerLaw \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDPowerLaw_Fe05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLaw \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLaw_Fe05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L025NYSNIIMassPowerLawComparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n lo8hi40DTDPowerLaw lo8hi40DTDtauPowerLawFe05 lo10hi40DTDPowerLaw lo10hi40DTDPowerLawFe05 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLaw \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLaw_Fe05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLawbeta08 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLawbeta08_Fe05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L025NYSNIIlo10hi40DTDPowerLawComparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDPowerLaw DTDtauPowerLawFe05 DTDPowerLawBeta08 DTDPowerLawBeta08Fe05 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLaw \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLaw_Fe05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLawbeta08 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L025NYSNIIlo10hi40DTDPowerLawComparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDPowerLaw DTDtauPowerLawFe05 DTDPowerLawBeta08 \
#                     -m 5e8


