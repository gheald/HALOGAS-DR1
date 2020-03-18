# HALOGAS-DR1

This repository contains some scripts that are being used to prepare data for the first public data release from the HALOGAS survey.

Contents:

* `README.md`	This file
* `continuum-images.dat` Index of continuum image locations on local computer
* `cube-locations.dat` Index of cube locations on local computer
* `galaxy-values.dat` Data file containing key parameters of each galaxy, needed for plots
* `get_masses.py` Determines total HI masses by integrating column density maps
* `get_masses_2.py` Determines total HI masses by integrating HI global profiles
* `get_vel_limits.py` Small script to determine velocity min/max of cube
* `make_cubes.py` Automates generation of HI cubes (fits format)
* `make_mom2.sh` Generate moment-2 maps (may not be used, not released in DR1)
* `make_moments.py` Automates generation of moment maps from cubes using SoFiA
* `make_pv.py` Generates fits files containing major axis PV diagrams of each galaxy
* `mask_moments_n0891.sh` As `mask_moments.sh` but tuned for NGC 891
* `mask_moments_n2403.sh` As `mask_moments.sh` but tuned for NGC 2403
* `mask_moments.sh` Automates processing of moment maps to produce final products
* `optical-images.dat` Index of optical image locations on local computer
* `plot_globprofs.py` Determines HI global profiles
* `plot_grid2x2.py` Plots "atlas" images (2x2 grids with HI moment maps, optical and continuum images)
* `plot_moments.py` Plots moment maps for atlas paper
* `plot_pv.py` Plots the PV diagrams created with `make_pv.py` in a format suitable to fit in the "atlas" grids
* `sf_params.txt` SoFiA template file used by `make_moments.py`
