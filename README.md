Dispa-SET
===================

### Description
The Dispa-SET model is a unit commitment and dispatch model developed within the “Joint Research Centre” and focused on the balancing and flexibility problems in European grids. It is written in GAMS and coupled to Matlab and Excel for input/output data handling and visualization. The selected Mixed-Integer Linear Programming (MILP) solver is CPLEX.

 
### Features
The model is expressed as a MILP problem. Continuous variables include the individual unit dispatched power, the shedded load and the curtailed power generation. The binary variables are the commitment status of each unit. The main model features can be summarized as follows:

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

The demand is assumed to be inelastic to the price signal. The MILP objective function is therefore the total generation cost over the optimization period. 

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
This software is primarily developed and used within the Institute for Energy and Transport, which is one of the 7 scientific institutes of the Joint Research Centre (JRC) of the European Commission. The IET is based both in Petten, the Netherlands, and Ispra, Italy. The Dispa-SET model is developed in the framework of the "Energy Systems Modelling" (ESM) project.

### Licence
Dispa-SET is a free software licensed under the “European Union Public Licence" EUPL v1.1. It 
can be redistributed and/or modified under the terms of this license.

### Main developpers
- Sylvain Quoilin (European Commission, Institute for Energy and Transport)
- Ignacio Hidalgo (European Commission, Institute for Energy and Transport)
- Andreas Zucker (European Commission, Institute for Energy and Transport)

