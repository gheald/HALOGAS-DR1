# HALOGAS-DR1

This repository contains some scripts that are being used to prepare data for the first public data release from the HALOGAS survey.

Contents:

* `README.md`	This file
* `cube-locations.dat` Index of cube locations on local computer
* `get_masses.py` Determines total HI masses by integrating column density maps
* `get_masses_2.py` Determines total HI masses by integrating HI global profiles
* `get_vel_limits.py` Small script to determine velocity min/max of cube
* `make_cubes.py` Automates generation of HI cubes (fits format)
* `make_moments.py` Automates generation of moment maps from cubes using SoFiA
* `mask_moments.sh` Automates processing of moment maps to produce final products
* `plot_globprofs.py` Determines HI global profiles
* `plot_moments.py` Plots moment maps for atlas paper
* `sf_params.txt` SoFiA template file used by `make_moments.py`
