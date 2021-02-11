MorpholoPy
=========

Python package that calculates the following morphological parameters :

kappa_co (see e.g. )
specific angular momentum
axial ratios (see e.g )

for stars and HI+H2 gas.

Morphology outputs figures and tables in hdf5 format.

Requirements
----------------

The morpholopy package requires:

+ `python3.6` or above
+ see requirements.txt

Installation
------------

You can install this library from PyPI using:
```
pip3 install velociraptor
```

How to use it
---------------

Given the folder where your simulation is stored (simulation_folder) and the
specific snapshot of the simulation (e.g 36) type:

```python morpholopy simulation_folder snapshot
```


