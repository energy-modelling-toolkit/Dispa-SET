.. _overview:

Overview
========

:Organization:  `Joint Research Centre`_,
		`European Commission`_
:Version: |version| (|release|)
:Date: |today|

The Dispa-SET model is mainly developed within the “Joint Research Centre” of the European Commission and focuses on the balancing and flexibility problems in European grids [1]_.

It is written in GAMS and uses csv files for input data handling. The optimisation is defined as a Linear Programming (LP) or Mixed-Integer Linear Programming (MILP) problem, depending on the desired level of accuracy and complexity. Continuous variables include the individual unit dispatched power, the shedded load and the curtailed power generation. The binary variables are the commitment status of each unit. The main model features can be summarized as follows:


Features
--------

- Minimum and maximum power for each unit
- Power plant constraints: minimum power, ramping limits, minimum up/down times, start-up and no-load costs
- Outages (forced and planned) for each unit
- Generalised reserve system with European reserve types (FFR, FCR, aFRR, mFRR) — both upward and downward
- System inertia requirements: synchronous inertia from conventional units and virtual inertia from inverter-based resources
- Under Frequency Load Shedding (UFLS) and Optional Downward Flexibility Management (OFDM) as emergency frequency services
- Load shedding
- Curtailment of variable renewable generation
- Pumped-hydro storage, batteries, and generic storage with alert/flood-control bounds
- Non-dispatchable units (e.g. wind turbines, run-of-river, etc.)
- Multi-node network with NTC or DC Power Flow (PTDF) constraints
- Emission caps per zone and pollutant
- CHP power plants and thermal storage
- Power-to-X and X-to-Power units for sector coupling (heat, hydrogen, Power-to-Liquid, etc.)
- Boundary sector model for representing non-electric energy sectors with storage, flexible demand and supply
- DSM-ready demand and flexible demand model
- Integrated mid-term scheduling (MTS) and short-term optimal dispatch
- Rolling horizon simulation with configurable horizon and look-ahead period
- Multiple model formulations (LP, MIP, integer clustering) generated from the same dataset
- Python-based solver using `linopy`_ as an alternative to the GAMS/Gurobi solver

The demand is assumed to be inelastic to the price signal. The MILP objective function is therefore the total generation cost over the optimization period. 


Libraries used and requirements
-------------------------------

* `Python 3.8+`_
* `pandas`_ for input and result data handling
* `matplotlib`_ for plotting
* `linopy`_ for the Python-native LP/MILP solver interface
* `GAMS_api`_ (optional) for the communication with GAMS/Gurobi

The above are auto-installed in a conda environment if you follow the Quick Start instructions.

Dispa-SET in the scientific literature
--------------------------------------

In the past years, Dispa-SET has been used in various scientific works covering different geographical areas and with different focus points. The works for which scientific articles have been published are summarized hereunder:


* Hydropower for flexibility services in the European power system [2]_.
* Generating stylized flexibility constraints for the JRC-EU-TIMES model [3]_ [4]_.
* Impact of Electric Vehicle deployment in The Netherlands [5]_.
* Open-source model of the Balkans area, with some simulations involving high shares of renewables [6]_ [7]_.
* Specific country studies for RES integration (Belgium, Greece) [3]_ [9]_
* Comparison between model formulations and levels of clustering [10]_ [11]_ [20]_
* Benders decomposition for capacity expansion [14]_
* The water-energy nexus in Greece and in Africa [8]_ [9]_ [15]_
* Soft-linking between JRC-EU-TIMES and Dispa-SET at the EU level [17]_
* Quantifying the flexibility provided by coupling the heating and power sectors [12]_ [13]_ [18]_ [19]_ [16]_ 
* Power systems adequacy and flexibility assessments in developing countrie (Africa, Bolivia) [15]_ [21]_ 



Ongoing developments
--------------------
The Dispa-SET project is under active development. Areas of ongoing work include:

- Stochastic scenarios
- Modelling of investment and capacity expansion
- Modelling of ancillary markets and frequency services
- Extended sector coupling (Power-to-Gas, Power-to-Liquid, district heating)
- Improved Python-native solver (linopy) feature parity with the GAMS formulation


Licence
-------
Dispa-SET is a free software licensed under the “European Union Public Licence" EUPL v1.2. It 
can be redistributed and/or modified under the terms of this license.

Main Developers
---------------
- `Sylvain Quoilin`_ (University of Liège, KU Leuven)
- `Konstantinos Kavvadias`_ (European Commission, Joint Research Centre)
- `Matija Pavičević`_ (KU Leuven, Belgium)


References
----------
.. [1] Quoilin, S., Hidalgo Gonzalez, I., & Zucker, A. (2017). Modelling Future EU Power Systems Under High Shares of Renewables: The Dispa-SET 2.1 open-source models. Publications Office of the European Union.
.. [2] Sánchez Pérez, A. (2017), Modelling Hydropower in detail to assess its contribution to flexibility services in the European power system. Master Thesis, University of Utrecht, Netherlands.
.. [3] Quoilin, S., Nijs, W., Gonzalez, I. H., Zucker, A. and Thiel, C. (2015), Evaluation of simplified flexibility evaluation tools using a unit commitment model, In 12th International Conference on the European Energy Market (EEM), pp. 1 5.
.. [4] Quoilin, S., Nijs, W. and Zucker, A. (2017), Evaluating flexibility and adequacy in future EU power systems: model coupling and long-term forecasting, In Proceedings of the 2017 ECOS Conference, San Diego.
.. [5] Beltramo, A., Julea, A., Refa, N., Drossinos, Y., Thiel, C. and Quoilin, S. (2017),`Using electric vehicles as flexible resource in power systems: A case study in the Netherlands, In 14th International Conference on the European Energy Market (EEM).
.. [6] Pavičević, M., Tomić, I., Quoilin, S., Zucker, A. and Pukšec, T. and Krajačić, G. (2017), Applying the Dispa-SET model on the Western Balkans power systems, In Proceedings of the 2017 SDEWES Conference
.. [7] Tomić, I., Pavičević, M., Quoilin, S., Zucker, A., Krajačić, G., Pukšec, T. and Duić, N. (2017), Applying the Dispa-SET model on the seven countries from the South East Europe, In 8th Energy Planning and Modeling of Energy Systems-Meeting, Belgrade
.. [8] Ricardo Fernandez Blanco Carramolino, Konstantinos Kavvadias, Ignacio Hidalgo Gonzalez (2017). Water-related modelling in electric power systems: WATERFLEX Exploratory Research Project.
.. [9] Ricardo Fernandez Blanco Carramolino, Konstantinos Kavvadias, I Hidalgo Gonzalez (2017). Quantifying the water-power linkage on hydrothermal power systems: A Greek case study. Applied Energy.
.. [10] Pavičević, M., Quoilin, S. and Pukšec, T., (2018). Comparison of Different Power Plant Clustering Approaches for Modeling Future Power Systems, Proceedings of the 3rd SEE SDEWES Conference, Novi Sad.
.. [11] Pavičević, M., Kavvadias, K. and Quoilin, S. (2018). Impact of model formulation on power system simulations - Example with the Dispa-SET Balkans model, EMP-E conference 2018: Modelling Clean Energy Pathways, Brussels.
.. [12] Juan Pablo Jiménez Navarro, Konstantinos Kavvadias, Sylvain Quoilin, Zucker Andreas (2018). The joint effect of centralised cogeneration plants and thermal storage on the efficiency and cost of the power system. Energy.
.. [13] Kavvadias, K., Jimenez Navarro, J.-P., Zucker, A., & Quoilin, S. (2018). Case study on the impact of cogeneration and thermal storage on the flexibility of the power system (KJ-NA-29082-EN-N). Netherlands: Publication Office of the European Commission.
.. [14] Matthias Zech, Acceleration strategies of the Generation Expansion Planning problem using Benders Decomposition, Master Thesis, Dresden University of Technology, 2018
.. [15] Matteo De Felice, Iratxe Gonzalez-Aparicio, Thomas Huld, Sebastian Busch, Ignacio Hidalgo-Gonzalez . Analysis of the water-power nexus in the West African power pool. JRC Technical Report, 2019.
.. [16] Matija Pavičević, Juan-Pablo Jimenez, Konstantinos Kavvadias, Sylvain Quoilin (2019). Modeling the flexibility offered by coupling the heating sector and the power sector: an assessment at the EU level. 5th International Conference On Smart Energy Systems.
.. [17] Matija Pavičević, Wouter Nijs, Konstantinos Kavvadias, Sylvain Quoilin (2019). Modelling flexible power demand and supply in the EU power system: soft-linking between JRC-EU-TIMES and the open-source Dispa-SET model. Proceedings of the 32nd International Conference on Efficiency, Cost, Optimization, Simulation and Environmental Impact of Energy Systems.
.. [18] Konstantinos Kavvadias, Georg Thomassen, Matija Pavičević, Sylvain Quoilin (2019). Electrifying the heating sector in Europe: The impact on the power sector. Proceedings of the 32nd International Conference on Efficiency, Cost, Optimization, Simulation and Environmental Impact of Energy Systems.
.. [19] Konstantinos Kavvadias, Juan Pablo Jimenez Navarro, Georg Thomassen (2019). Decarbonising the EU heating sector: Integration of the power and heating sector.
.. [20] Pavičević, M., Kavvadias, K., Pukšec, T., & Quoilin, S. (2019, June). Comparison of different model formulations for modelling future power systems with high shares of renewables – The Dispa-SET Balkans model. Applied Energy.
.. [21] Rojas Candia, R., Balderrama Subieta, S. L., Adhemar Araoz Ramos, J., Vicente Senosiain, M., Peña Balderrama, G., Jaldín Florero, H., & Quoilin, S. (2019). Techno-economic assessment of high variable renewable energy penetration in the Bolivian interconnected electric system. International Journal of Sustainable Energy Planning and Management, 22.



.. _Python 3.8+: https://www.anaconda.com/distribution/
.. _matplotlib: http://matplotlib.org
.. _pandas: http://pandas.pydata.org
.. _GAMS_api: https://github.com/kavvkon/gams-api
.. _European Commission: https://ec.europa.eu/
.. _Joint Research Centre: https://ec.europa.eu/jrc/en
.. _Sylvain Quoilin: http://squoilin.eu
.. _Konstantinos Kavvadias: http://kavvadias.eu
.. _Matija Pavičević: https://www.mpavicevic.com
.. _linopy: https://linopy.readthedocs.io


