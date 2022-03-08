#!/bin/bash

python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Mass_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe2M8to10 \
                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Mass_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe2M8to20 \
                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Mass_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe2M8to30 \
                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe05 \
                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_comparisons/TestMassDepBoostFactorFe2 \
                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
                     -n DTDExptau2nu1p6Fe2M8to10 DTDExptau2nu1p6Fe2M8to20 DTDExptau2nu1p6Fe2M8to30 DTDExptau2nu1p6Fe05 \
                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDExptau2nu4 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDExptau4nu1_O05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDExptau4nu1p5 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDExptau4nu1p5_O05Mg15 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_comparisons/DTDExptau4NewComparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDExptau2nu4 DTDExptau4nu1O05 DTDExptau4nu1p5 DTDtau4nu1p5O05Mg15 \
#                     -m 5e8


#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_Fe05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_O15Mg15 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_Fe07 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_comparisons/DTDExptau2NewComparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDExptau2nu2 DTDExptau2nu2Fe05 DTDExptau2nu2O15Mg15 DTDtau2nu2Fe07 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu4 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDExptau2nu4 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu3 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu3_Fe05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_comparisons/DTDExpNewComparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n 8to40DTDExptau2nu4 10to40DTDExptau2nu4 8to40DTDExptau2nu3 8to40DTDExptau2nu3Fe05 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDtau1p5nu1p6_Fe05O05Mg05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDtau2nu1p6 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDtau2nu1_O05Fe05Si05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_comparisons/L025NYSNIIDTDExpTau2Comparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDtau1p5nu1p6Fe05O05 DTDtau2nu1p6 DTDtau2nu1Fe05O05 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDPowerLaw \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDPowerLaw_Fe05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDtau1p5nu1p6_Fe05O05Mg05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDtau2nu1p6 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L025NYSNIIlo8hi40DTDComparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDPowerLawBeta05nu1 DTDtauPowerLawBeta05nu1Fe05 DTDtau1p5nu1p6Fe05O05Mg05 DTDtau2nu1p6 \
#                     -m 5e8


#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_StandardYields_SNIIlo8hi100_DTDtau2nu1p6 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p6 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDtau2nu1p6 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L025CompareYields \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n StandardYields NewYieldsSNIIlo8hi100 NewYieldsSNIIlo8hi40 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDPowerLaw \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDPowerLaw_Fe05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLaw \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLaw_Fe05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L025NYSNIIMassPowerLawComparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n lo8hi40DTDPowerLawBeta05nu1 lo8hi40DTDtauPowerLawBeta05nu1Fe05 lo10hi40DTDPowerLawBeta05nu1 lo10hi40DTDPowerLawBeta05nu1Fe05 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLaw \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLaw_Fe05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLawbeta08 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLawbeta08_Fe05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L025NYSNIIlo10hi40DTDPowerLawComparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDPowerLawBeta05nu1 DTDtauPowerLawBeta05nu1Fe05 DTDPowerLawBeta08nu1 DTDPowerLawBeta08nu1Fe05 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLaw \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLaw_Fe05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLawbeta08 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L025NYSNIIlo10hi40DTDPowerLawComparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDPowerLaw DTDtauPowerLawFe05 DTDPowerLawBeta08 \
#                     -m 5e8


