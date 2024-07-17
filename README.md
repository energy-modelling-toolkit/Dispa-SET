![dispaset logo](https://raw.githubusercontent.com/energy-modelling-toolkit/Dispa-SET/master/Docs/figures/logo.png)
===================
 ![Documentation](https://img.shields.io/badge/python-2.7,%203.7-blue.svg) [![License](https://img.shields.io/badge/License-EUPL--1.2-blue.svg)](https://opensource.org/licenses/EUPL-1.2) [![Documentation](https://readthedocs.org/projects/dispa-set/badge/?branch=master)](http://dispa-set.readthedocs.io/en/latest/) [![Build Status](https://travis-ci.org/energy-modelling-toolkit/Dispa-SET.svg?branch=master)](https://travis-ci.org/energy-modelling-toolkit/Dispa-SET)

### Description
The Dispa-SET model is a unit commitment and dispatch model developed within the “Joint Research Centre” and focused on the balancing and flexibility problems focusing on the European context. It is written in GAMS with advanced input/output data handling and visualization routines in Python.

Three different formulations are available offering a trade-off between accuracy and computational complexity ( Linear Programming (LP), Mixed-Integer Linear Programming (MILP)). This allows
 to model a power system at any level of detail e.g. micro-grid, region, country, continent. A Pan-European scenario is included with the model as of version 2.3.
 
### Features
The model is expressed as an optimization problem. 
Continuous variables include the individual unit dispatched power, the shedded load and the curtailed power generation. The binary variables are the commitment status of each unit. The main model features can be summarized as follows:

- Minimum and maximum power for each unit
- Power plant ramping limits
- Reserves up and down
- Minimum up/down times
- Load Shedding
- Curtailment
- Pumped-hydro storage
- Non-dispatchable units (e.g. wind turbines, run-of-river, etc.)
- Start-up, ramping and no-load costs
- Multi-nodes with capacity constraints on the lines (congestion)
- Constraints on the targets for renewables and/or CO2 emissions
- Yearly schedules for the outages (forced and planned) of each units
- CHP power plants and thermal storage

The demand is assumed to be inelastic to the price signal. The MILP objective function is therefore the total generation cost over the optimization period. 

### Quick start

If you want to download the latest version from github for use or development purposes, make sure that you have git and the [anaconda distribution](https://www.anaconda.com/distribution/) installed and type the following:

```bash
git clone https://github.com/energy-modelling-toolkit/Dispa-SET.git
cd Dispa-SET
conda env create  # Automatically creates environment based on environment.yml
conda activate dispaset # Activate the environment
pip install -e . # Install editable local version
```
The above commands create a dedicated environment so that your anaconda configuration remains clean from the required dependencies installed.

At this point, it is necessary to make sure that the GAMS API is properly installed in the newly created environment:
- Make sure to define an evironment variable GAMSDIR pointing to the gams installation folder (e.g. C:\GAMS\win62\47.7.0)
- For GAMS version 45 and higher, the api can be installed directly from pip, replacing the x values below with the installed version of gams:
```bash
pip install gamsapi==4x.x.x
```
- For older version of gams, the API can be compiled from the apifile/Python directory in the gams installation folder. The the path to the library must then be added to the PYTHONPATH environment variable. For example, in linux:
```bash
cd /path/to/gams/apifiles/Python/api_310
python gamssetup.py install
export PYTHONPATH=/path/to/gams/apifiles/Python/gams:/path/to/gams/apifiles/Python/api_310
```
- For very old versions of GAMS (e.g. 24.x, 25.x), it is possible to install the old low-level api from pip. Note that if this does not work, the api must be de-installed and compiled from the gams apifile folder.
```bash
pip install gdxcc, gamsxcc, optcc
```

To check that everything runs fine, you can build and run a test case by typing:
```bash
dispaset -c ConfigFiles/ConfigTest.xlsx build simulate
```

### Documentation
The documentation and the stable releases are available on the main Dispa-SET website: http://www.dispaset.eu
 
### Get involved
This project is an open-source project. Interested users are therefore invited to test, comment or [contribute](CONTRIBUTING.md) to the tool. Submitting issues is the best way to get in touch with the development team, which will address your comment, question, or development request in the best possible way. We are also looking for contributors to the main code, willing to contibute to its capabilities, computational-efficiency, formulation, etc. Finally, we are willing to collaborate with national agencies, reseach centers, or academic institutions on the use on the model for different data sets relative to EU countries.

### License
Dispa-SET is a free software licensed under the “European Union Public Licence" EUPL v1.2. It 
can be redistributed and/or modified under the terms of this license.

### Main developers
This software has been developed initially within the Directorate C Energy, Transport and Climate, which is one of the 7 scientific directorates of the Joint Research Centre (JRC) of the European Commission. Directorate C is based both in Petten, the Netherlands, and Ispra, Italy. 
Currently the main developers are the following:

- Sylvain Quoilin (KU Leuven, Belgium)
- Konstantinos Kavvadias (Joint Research Centre, European Commission)
- Matija Pavičević  (KU Leuven, Belgium)
- Matthias Zech (Deutsches Zentrum für Luft-und Raumfahrt, DLR)
- Matteo De Felice (Joint Research Centre, European Commission)

