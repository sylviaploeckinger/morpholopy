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

How to use it
---------------

Specify in the file run.sh the folder where your simulation is stored (-d), the simulation
snapshot (e.g. -s) and the folder (-o) where plots will be stored.

Then type `bash run.sh` or `bash runcomparison.sh` or alternatively run it as

```bash
python morpholopy.py \
  -d /cosma7/data/Simulation1 /cosma7/data/Simulation2 ... \ # These are directories where sims outputs are
  -s 23 23 ... \ # Snapshot numbers from the different sims
  -n name_one name_two ... \ # Names for different sims (for legend)
  -o output_directory
```



