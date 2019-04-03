.. _caseeu28:

Dispa-SET for the EU28
======================

Description
-----------

Dispa-SET is provided with a ready-to-use dataset of the EU28 (+Norway +Switzerland) power system. A detailed description of the model and of the selected input data is available in [1]_ and [2]_.

The power plants in the model are represented as a cluster of units powered by the same fuel type and technology. They can be modelled together with a large number of RES units with separate hourly distribution curves.
 
Features
--------

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
- Integer clustering: a representative unit is considered for each technology and multiplied N times.

The demand is assumed to be inelastic to the price signal. The MILP objective function is therefore the total generation cost over the optimization period. 

Run the EU model
----------------

A specific config file is provided with the standard Dispa-SET installation (starting from v2.3). After installing Dispa-SET and checking that everything is fine, you can run the EU model in different ways:

From the command line (if Dispa-SET was properly installed)::

	dispaset -c ConfigFiles/ConfigEU.xlsx build simulate

Using the Dispa-SET Api and the provided example script:

	scripts\build_and_run_EU_model.py

  
Documentation
-------------
The general documentation of the Dispa-SET model and the stable releases are available on the main Dispa-SET website: http://www.dispaset.eu

Licence
-------
Dispa-SET is a free software licensed under the “European Union Public Licence" EUPL v1.2. It can be redistributed and/or modified under the terms of this license.

Important results
-----------------

.. .. image:: figures/Balkans_capacity.png

.. .. image:: figures/Balkans_generation.png

Main developpers
----------------
- Sylvain Quoilin (University of Liège, KU Leuven)
- Konstantinos Kavvadias (European Commission, Institute for Energy and Transport)
- Matija Pavičević (KU Leuven)

References
----------
More details regarding the model and its implementation are available in the following publications

.. [1] Kavvadias, K., Hidalgo Gonzalez, I., Zucker, A., & Quoilin, S. (2018). Integrated modelling of future EU power and heat systems - The Dispa-SET v2.2 open-source model (EUR 29085 EN). Luxembourg: European Commission.

.. [2] Matija Pavičević, Wouter Nijs, Konstantinos Kavvadias, Sylvain Quoilin,  Modelling flexible power demand and supply in the EU power system: soft-linking between JRC-EU-TIMES and the open-source Dispa-SET model, ECOS 2019 


