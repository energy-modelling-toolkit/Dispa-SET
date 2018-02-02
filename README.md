Dispa-SET 
===================
 ![Documentation](https://img.shields.io/badge/python-2.7,%203.6-blue.svg) [![License](https://img.shields.io/badge/License-EUPL--1.1-blue.svg)](https://opensource.org/licenses/EUPL-1.1) [![Documentation](https://readthedocs.org/projects/dispa-set/badge/?branch=master)](http://dispa-set.readthedocs.io/en/latest/)

### Description
The Dispa-SET model is a unit commitment and dispatch model developed within the “Joint Research Centre” and focused on the balancing and flexibility problems in European grids. It is written in GAMS and coupled to Matlab and Excel for input/output data handling and visualization. The selected Mixed-Integer Linear Programming (MILP) solver is CPLEX.

It is written in GAMS an Python (Pyomo) and uses csv files for input data handling. The optimisation is defined as a Linear Programming (LP) or Mixed-Integer Linear Programming (MILP) problem, depending on the desired level of accuracy and complexity. 

 
### Features
The model is expressed as an optimization problem. Continuous variables include the individual unit dispatched power, the shedded load and the curtailed power generation. The binary variables are the commitment status of each unit. The main model features can be summarized as follows:

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

If you want to download the latest version from github for use or development purposes, make sure that you have git and the [anaconda distribution](https://www.continuum.io/downloads) installed and type the following:

```bash
git clone https://github.com/energy-modelling-toolkit/Dispa-SET.git
cd Dispa-SET
conda env create  # Automatically creates environment based on environment.yml
source activate dispaset # in Windows: activate dispaset
pip install -e . # Install editable local version
```

The above commands create a dedicated environment so that your anconda configuration remains clean from the required dependencies installed.
To check that everything runs fine, you can build and run a test case by typing:
```bash
dispaset -c ConfigFiles/ConfigTest.xlsx build simulate
```

### Documentation
The documentation and the stable releases are available on the main Dispa-SET website: http://www.dispaset.eu


### Ongoing developments
The Dispa-SET project is relatively recent, and a number of improvements will be brought to the project in a close future:

- Grid constraints (DC power-flow)
- Stochastic scenarios
- Extension of the modeled areas
- Modeling of the ancillary markets
- User interface

 
### Get involved
This project is an open-source project. Interested users are therefore invited to test, comment or contribute to the tool. Submitting issues is the best way to get in touch with the development team, which will address your comment, question, or development request in the best possible possible. We are also looking for contributors to the main code, willing to contibute to its capabilities, computational-efficiency, formulation, etc. Finally, we are willing to collaborate with national agencies, reseach centers, or academic institutions on the use on the model for different data sets relative to EU countries.

### Public administration reference
This software has been developed initially within the Directorate C Energy, Transport and Climate, which is one of the 7 scientific directorates of the Joint Research Centre (JRC) of the European Commission. Directorate C is based both in Petten, the Netherlands, and Ispra, Italy. 

### Licence
Dispa-SET is a free software licensed under the “European Union Public Licence" EUPL v1.1. It 
can be redistributed and/or modified under the terms of this license.

### Main developers
- Sylvain Quoilin (University of Liège)
- Konstantinos Kavvadias (Joint Research Centre, European Commission)
- Andreas Zucker (Joint Research Centre, European Commission)

