.. _overview:

Overview
========
The model is expressed as a MILP problem. Continuous variables include the individual unit dispatched power, the shedded load and the curtailed power generation. The binary variables are the commitment status of each unit. The main model features can be summarized as follows:


Features
--------

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

Ongoing developments
--------------------
The Dispa-SET project is relatively recent, and a number of improvements will be brought to the project in a close future:

- Grid constraints (DC power-flow)
- Stochastic scenarios
- Extension of the modeled areas
- Modeling of the ancillary markets
- User interface

Public administration reference
-------------------------------
This software is primarily developed and used within the Institute for Energy and Transport, which is one of the 7 scientific institutes of the Joint Research Centre (JRC) of the European Commission. The IET is based both in Petten, the Netherlands, and Ispra, Italy. The Dispa-SET model is developed in the framework of the "Energy Systems Modelling" (ESM) project.


Licence
-------
Dispa-SET is a free software licensed under the “European Union Public Licence" EUPL v1.1. It 
can be redistributed and/or modified under the terms of this license.

Main Developers
---------------
- Sylvain Quoilin (European Commission, Institute for Energy and Transport)
- Ignacio Hidalgo (European Commission, Institute for Energy and Transport)
- Andreas Zucker (European Commission, Institute for Energy and Transport)
- Konstantinos Kavvadias (European Commission, Institute for Energy and Transport)

References
----------

.. [1] Hidalgo González, I., Quoilin, S., & Zucker, A. (2014). Dispa-SET 2.0: unit commitment and power dispatch model (EUR 27015 EN). Petten, Netherlands: European Commission. 
.. [2] Quoilin, S., Nijs, W., Hidalgo, I., & Thiel, C. (2015). Evaluation of simplified flexibility evaluation tools using a unit commitment model. IEEE Digital Library. 
.. [3] Quoilin, S., Gonzalez Vazquez, I., Zucker, A., & Thiel, C. (2014). Available technical flexibility for balancing variable renewable energy sources: case study in Belgium. Proceedings of the 9th Conference on Sustainable Development of Energy, Water and Environment Systems. 
