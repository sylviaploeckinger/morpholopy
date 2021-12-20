#!/bin/bash

#python morpholopy.py -d /cosma7/data/dp004/dc-chai1/HAWK/prediction2_full_noZ \
#                     -s colibre_0023.hdf5 -c halo_0023.properties -n wave2ColibreDTDFe1Mg4 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d /cosma7/data/dp004/dc-chai1/HAWK/prediction1 \
#                     -s colibre_0023.hdf5 -c halo_0023.properties -n wave1ColibreDTDFe1Mg4 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Test_yield_tables/L12N188_SNIa2Gyr_Nomoto_boost_1/ \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Test_yield_tables/L12N188_SNIa2Gyr_Nomoto_boost_2/ \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Test_yield_tables/L12N188_newAGB_oldSNII_SNIa2Gyr/ \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 -c halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -n L12NomotoFe05Mg2 L12NomotoFe1Mg1 L12PortinariFe05Mg2 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L012N188_SNIa_model_1 \
#                     /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L012N188_SNIa_model_2 \
#                     /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L012N188_SNIa_model_3 \
#                     /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/L012N188_SNIa_model_4 \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -n L12SNIatau1Gyr L12SNIatau2Gyr L12SNIatau3Gyr L12SNIatau4Gyr \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Test_yield_tables/L12N188_newAGB_oldSNII_SNIa2Gyr/ \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/EAGLE_SNIa_Fe_05_Mg_2_L012N188/ \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 -c halo_0036.properties halo_0036.properties \
#                     -n L12newAGB L12oldyields \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d ../../ \
#                     -s colibre_0034.hdf5 -c halo_0034.properties -n test \
#                     -o ../../Morpholopy_Plots/


python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Test_yield_tables/L12N188_newAGB_oldSNII_SNIa2Gyr/ \
                     -s colibre_0036.hdf5 -c halo_0036.properties \
                     -n L12newAGB \
                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Test_yield_tables/L12N188_newAGB_oldSNII_SNIa2Gyr/ \
#                      /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/EAGLE_SNIa_Fe_05_Mg_2_L012N188/ \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 -c halo_0036.properties halo_0036.properties \
#                     -n L12newAGB L12oldyields \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d ../../ \
#                     -s colibre_0034.hdf5 -c halo_0034.properties -n test \
#                     -o ../../Morpholopy_Plots/


#python morpholopy.py -d /cosma7/data/dp004/dc-chai1/HAWK/123b2_121b2_SNII_Mg_2_Fe_1_L012N188/ \
#                     -s colibre_0023.hdf5 -c halo_0023.properties -n 123b2121b2L012N18 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/
#
#python morpholopy.py -d /cosma7/data/dp004/dc-chai1/HAWK/116b2_115b2_Fe1_new_cooling_L012N188/ \
#                     -s colibre_0023.hdf5 -c halo_0023.properties -n 116b2115b2L012N18 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/
#
#python morpholopy.py -d /cosma7/data/dp004/dc-chai1/HAWK/106b3_104b2_norm_0p3_new_cooling_L012N188/ \
#                     -s colibre_0023.hdf5 -c halo_0023.properties -n 106b3104b2L012N18 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/
#
#python morpholopy.py -d /cosma7/data/dp004/dc-chai1/HAWK/122b2_121b2_SNIa_best_fit_L012N188/ \
#                     -s colibre_0023.hdf5 -c halo_0023.properties -n 122b2121b2L012N18 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/EAGLE_SNIa_boost_1/ \
#                     -s colibre_0034.hdf5 -c halo_0034.properties -n SNIaFe05Mg2L06N94 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/EAGLE_SNIa_boost_2/ \
#                     -s colibre_0036.hdf5 -c halo_0036.properties -n SNIaFe025Mg2L06N94 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_11_COLIBRE/EAGLE_SNIa_diffusion_1/ \
#                     -s colibre_0036.hdf5 -c halo_0036.properties -n SNIaFe1Mg2L06N94 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/SNIa_variations/Model_2/ \
#                     -s colibre_0036.hdf5 -c halo_0036.properties -n SNIa4GyrFe1Mg2L06N94 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/SNIa_variations/Model_5/ \
#                     -s colibre_0036.hdf5 -c halo_0036.properties -n SNIa4GyrFe05Mg2L06N94 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/SNIa_variations/Model_6/ \
#                     -s colibre_0036.hdf5 -c halo_0036.properties -n SNIa4GyrFe025Mg2L06N94 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/EAGLE_SNIa_Fe_05_Mg_2_L012N188/ \
#                     -s colibre_0036.hdf5 -c halo_0036.properties -n SNIa2GyrFe05Mg2L12N188 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/
#
#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/EAGLE_SNIa_Fe_025_Mg_2_L012N188/ \
#                     -s colibre_0036.hdf5 -c halo_0036.properties -n SNIa2GyrFe025Mg2L12N188 \
#                     -o /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2021_12_COLIBRE/Morholopy_Plots/


