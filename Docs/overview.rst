.. _overview:

Overview
========

:Organization:  `Joint Research Centre`_, 
		`Institute of Energy and Transport`_,
                Energy Technology Policy Outlook Unit
:Version: |version|
:Date: |today|

The Dispa-SET model is mainly developed within the “Joint Research Centre” of the European Commission and focused on the balancing and flexibility problems in European grids [1]_. 

It is written in GAMS an Python (Pyomo) and uses csv files for input data handling. The optimisation is defined as a Linear Programming (LP) or Mixed-Integer Linear Programming (MILP) problem, depending on the desired level of accuracy and complexity. Continuous variables include the individual unit dispatched power, the shedded load and the curtailed power generation. The binary variables are the commitment status of each unit. The main model features can be summarized as follows:


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
- CHP power plants and thermal storage

The demand is assumed to be inelastic to the price signal. The MILP objective function is therefore the total generation cost over the optimization period. 


Libraries used
--------------

* `pyomo`_ Optimization object library, interface to LP solver (e.g. CPLEX)
* `pandas`_ for input and result data handling 
* `matplotlib`_ for plotting
* `GAMS_api`_ for the communication with GAMS


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
- Sylvain Quoilin (University of Liège)
- Andreas Zucker (European Commission, Institute for Energy and Transport)
- Konstantinos Kavvadias (European Commission, Institute for Energy and Transport)

References
----------
.. [1] Quoilin, S., Hidalgo Gonzalez, I., & Zucker, A. (2017). Modelling Future EU Power Systems Under High Shares of Renewables: The Dispa-SET 2.1 open-source models. Publications Office of the European Union.  
.. [2] Hidalgo González, I., Quoilin, S., & Zucker, A. (2014). Dispa-SET 2.0: unit commitment and power dispatch model (EUR 27015 EN). Petten, Netherlands: European Commission. 
.. [3] Quoilin, S., Nijs, W., Hidalgo, I., & Thiel, C. (2015). Evaluation of simplified flexibility evaluation tools using a unit commitment model. IEEE Digital Library. 
.. [4] Quoilin, S., Gonzalez Vazquez, I., Zucker, A., & Thiel, C. (2014). Available technical flexibility for balancing variable renewable energy sources: case study in Belgium. Proceedings of the 9th Conference on Sustainable Development of Energy, Water and Environment Systems. 

.. _matplotlib: http://matplotlib.org
.. _pandas: http://pandas.pydata.org
.. _pyomo: http://www.pyomo.org/
.. _GAMS_api: http://www.gams.com/help/index.jsp?topic=%2Fgams.doc%2Fapis%2Findex.html
.. _Institute of Energy and Transport: https://ec.europa.eu/jrc/en/institutes/iet
.. _Joint Research Centre: https://ec.europa.eu/jrc/en



