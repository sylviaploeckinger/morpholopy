MorpholoPy
=========

Python package that calculates the following morphological parameters :

kappa_co (see e.g. https://arxiv.org/pdf/1704.06283.pdf), specific angular momentum and axial ratios (see e.g https://arxiv.org/pdf/1805.03210.pdf ) for stars and HI+H2 gas. It also calculates HI+H2 surface densities and makes
particle projections of the ten most massive galaxies in the simulations.

Morphology outputs figures and tables in a html format.

Requirements
----------------

The morpholopy package requires:

+ `python3.6` or above
+ see requirements.txt

Usage
---------------

To run the script in the _single-run_ mode use
```bash
 python3 morpholopy.py -d run_directory \
                       -s snapshot_name \
                       -c catalogue_name \
                       -n name_of_the_run \
                       -g number_of_individual_galaxies_to_show \
                       -o path_to_output_directory \
                       -m minimal_stellar_mass_for_galaxies_to_be_analysed
```

To run the script in the _comparison_ mode use
```bash
 python3 morpholopy.py -d directory_of_run1 directory_of_run2 \
                       -s snapshot_name_run1 snapshot_name_run2 \
                       -c catalogue_name_run1 catalogue_name_run2 \
                       -n name_of_the_run1 name_of_the_run2 \
                       -g number_of_individual_galaxies_to_show \
                       -o path_to_output_directory \
                       -m minimal_stellar_mass_for_galaxies_to_be_analysed
```



