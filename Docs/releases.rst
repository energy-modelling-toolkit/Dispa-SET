.. _releases:

Releases
========

Major stable releases:

*  Dispa-SET v3.0

*  Dispa-SET v2.5

* `Dispa-SET v2.4`_

* `Dispa-SET v2.3`_ 

* `Dispa-SET v2.2`_

* `Dispa-SET v2.1`_

* `Dispa-SET v2.0`_

Changelog
---------

Version 3.0
^^^^^^^^^^^
**New**

* Boundary Sectors (renamed to Sector X)
    - Full integration of multi-sector coupling through boundary sectors
    - Flexible demand for boundary sectors
    - Network Transfer Capacity (NTC) between boundary sectors
    - Storage integrated with boundary sectors
    - Sector X water management with flood control
    - Alert level and spillage for boundary sectors
    - Proper support for CHP and power-to-heat technologies

* DC Power Flow
    - Implemented DC power flow with PTDF matrix
    - Backward compatibility with NTC model
    - Power flow limits with maximum and minimum bidirectional NTC

* Inertia and PFR constraints
    - Added version 1.0, 1.1, and 1.2 of inertia constraints
    - Primary Frequency Response (PFR) constraints

* New Units and Technologies
    - Added electrolyser model
    - Added x2p (X-to-Power) units
    - Added boundary sector unit to the base case
    - Improved cycling storage conditions for medium-term storage
    - Fuel prices can now be defined per unit

* Model Engine
    - Basic structure for a Linopy implementation started

* Post-Processing
    - Added dispatch plot for boundary sectors
    - Added hide_storage feature in the plot_zone function

* User Interface
    - Replaced Excel configuration files with an HTML configuration editor
    - Improved user experience with a more intuitive interface
    - Better visualization of configuration parameters

**Improvements**

* Code Quality
    - Upgraded compatibility with newer versions of libraries
    - Fixed dependencies with the new GAMS API
    - Added CI/CD with GitHub Actions
    - Improved error handling

* Fixes and Enhancements
    - Better handling of demand balancing to prevent charging reservoirs with slack
    - Improved representation of heating technologies
    - Fixed various issues with boundary sector formulation
    - Updated startup and shutdown costs implementation to display properly
    - Improved water-energy nexus representation
    - Enhanced zone configuration
    - Optimized plot functions for better performance and visualization

Version 2.5
^^^^^^^^^^^
**New**

* Fuels/Technologies
	- AIR, AMO, WHT
	- WAVE, WHEN, GETH, HOBO, SOTH, ABHP, ASHP, GSHP, HYHP, REHE, WSHP

* GIS
	- Added ability to load geographic coordinates of indivdual zones inside config file (currently only applicable for electricity). An example file can be found inside the \test directory 

* Heating sector
	- Added Heating zones analogue to electricity zones
	- Demand in each zone can now be covered by multiple units
	- Added Heat only generation technologies to the mix
	- Thermal storage can be attached directly to a unit or used as a standalone technology
	- Heating sector dispatch plots
	- Heat sector installed capacity and anual generation plot
	
* Reserves
	- Added three new reserve forcasting methods: Exogenous, Percentage and Probablistic
	- CHP and power to heat units can now provide reserves

* Hydrogen sector
	- Added new hydrogen technologies
	- H2 zones designed analoguos to the heating sector (i.e. one demand per H2 node and multiple technologies patricipating to the same H2 market)

* Water-energy nexus
	- Water consumption and water withdrawal parameters used by the cooling systems inside the thermal generators are added to the database. These indicators are useful for the assessment of the water-energy nexus. 

* RES curtailment
	- RES curtailment is now an optimization variable
	- Curtailment cost can be specified in the cinfig file. If the cost is > 0 EUR/MWh and there is some curtailment, shadow prices will be negative (this alligns closer to the real world markets where exces wind generation results in negative market clearing prices) 

* New model outputs (i.e. Results)
	- OutputCurtailedHeat, OutputCurtailedPower - Hourly time series of curtailments
	- ShadowPrice_2U, ShadowPrice_2D, ShadowPrice_3U - Hourly timeseries of balancing prices
	- HeatShadowPrice, OutputH2ShadowPrice - Hourly timeseries of heat and H2 prices
	- StorageShadowPrice - Hourly timeseries of storage prices (i.e. intra temporal price of the next MWh dispatched by the storage units)
	- OutputPtLDemand, OutputH2Output - Hourly timeseries of the Power to liquid demand (fixed anual value distributed to hourly values through MTS) and Power output of H2 units
	- OutputReserve_2U, OutputReserve_2D, OutputReserve_3U - Hourly timeseries of the power reserved for the balancing
	- ShadowPrice_RampUp_TC, ShadowPrice_RampDown_TC - Hourly timeseries of ramping costs (i.e. if dispatchable unit needs to be switched on/off the price would reflect the start-up shut-down costs specified in the power plants database)
	- OutputRampRate - Hourly timeseries of ramping rates
	- OutputStartUp, OutputShutDown - Hourly timeseries of start-up and shut-down events (i.e. if one powerplant has several Nunits of which 2 are switched on or 7 are switched off this number would be visible here)
	- OutputCostStartUpH, OutputCostRampUpH - Hourly timseries of startup and tamping costs
	- OutputEmissions - Hourly timseries of zonal CO2 intesity (i.e. value will be >0 if fossil generators are on and 0 if demand is 100% satisfied by the RES)
	- OutputStorageSlack - Unsatisifed minimum storage constraint at the end of the optimization horizon (one value per horizon) 
	
* Post-processing
	- (plot_zone) New fuel/technology-types added to the plots (color codes are now slightly improved)
	- (plot_zone) Power consumption and storage charging (storage input) added to the energy generation plots.  This gives us more insigts into energy conversion/transmission  
	- (plot_zone) Added generation share plot (i.e. expressed in terms of % rather than absolute values)
	- (plot_zone_capacities) Peak power consumption line added to the installed capacities plot. This gives us more insights into sector-coupling links
	- (plot_tech_cap) Added installed storage capacity plot. It shows how much storage capacity is located in individual zones
	- (plot_co2) Added CO2 intensity violin plots. It provides insights into the distribution of CO2 intensity inisde individual zones
	- (plot_power_flow_tracing_matrix) Added power flow tracing matrix. Indication of local generation and exports to neighboring ones
	- (plot_net_flows_map) Added power feed plot. Indication of net importing/exporting zones (WRNING! this plots take a while to generate. Sometimes Stamen server is timed out and background pictures are not loaded properly resulting in "http.client.IncompleteRead: IncompleteRead(xxxxx bytes read, yyyy more expected)". Currently the only way to fix this is to try the plotting function at a later point in time.
	- (plot_line_congestion_map) Added line congestion plot. Indication of the congested lines and directions (i.e. line between Z1 and ZZ3 is congested 60% of the time, meaning that Flow/MaxFLow = 1) 
	- (get_result_analysis) Added a more detailed statistical system representation such as Total and peak/max values for individual variables of intrest, 
	- (get_result_analysis) Added a more detailed statistical zonal representation same as above but filtered by zone and summarized under ZoneData 
	- (get_result_analysis) Added a more detailed statistical unit representation same as above but filtered by units and summarized under UnitData
	- (get_result_analysis) Added a more detailed statistical fuel representation same as above but filtered by fuel and summarized under FuelData
	- (get_result_analysis) Added water consumption which can be filtered on per Zone/Unit level and is summarized under WaterConsumptionData 
	- (get_result_analysis) Added a more detailed storage analysis summarized under StorageData 

**Bugfixes**

* Variable time step
	- The pre-processing and the GAMS file have been updated to handle different time steps (not only one hour)
	- This is currently restricted to three time steps: 15min, 1h, 24h
	- The input data whose time step is lower than the desired one is averaged

* Miscellaneous
	- Improved error handling


Version 2.4
^^^^^^^^^^^
* Mid-term scheduling
	- The yearly storage level profiles can now be calculated internally (i.e. without providing exogenous profiles).
	- A first, simplified version of dispa-set is run over a whole year to generate these profiles during the pre-processing phase
	- This option is activated in the config file and is transparent for the user.

* Flexible Demand:
	- To model demand-side management, it is now possible to define a share of the demand curve as "flexible"
	- In this flexible demand, the load can be shifted from one hour to the other
	- The maximum flexibility is characterized by the equivalent number of storage hours for the shifted load, which is defined as parameter in the configuration file.

* Power-to-heat units
	- P2HT units (heat pumps, electrical heater) have now been added
	- They are coupled to a heat demand and possibly to a thermal storage capacity
	- COP can be defined as temperature-dependent. An additional input with temperature times for each zone has been defined.

* Transmission prices have been added to the pre-processing and can now be fully parametrized

* Fuel Prices can now be country-specific

* Input data in the csv files can now be defined with time stamps from any year or with a numerical index

* Post-processing:
	- Improved dispatch plot with shifted, shed loads and electricity consumption from P2HT units
	- Storage levels are now differentiated by technology

* Miscellaneous:
	- Multiple bug fixes, code improvement and usability improvement.
	- All config files and the example scripts have been checked and cleaned
	- New formulation of the clustering function with significant simulation time improvements
	- The Pyomo version of Dispa-SET has now been removed since it was no longer up-to-date
	- The end-of horizon reservoir level is no longer a firm constraint. A water value can be defined to impose a price on the unmet level requirements.
	- Excel configuration files are now subject to versioning, which ensures backward compatibility with older configuration files.
	- Countries are now renamed into "zones" in all API functions.
	- The option to cache csv file data when loading has been removed
	- Implemented a more robust versioning system

Version 2.3
^^^^^^^^^^^
* Input Data: 
	- A complete EU dataset has been included to the repository for the year 2016. 
	- More information: :ref:`caseeu28`.

* Reformulation of the reserve constraints:
	- Secondary reserves are now covered by spinning units only. 
	- Tertirary reserves can also be covered by quick start units. 
	- In total, three different reserve markets are now considered: Secondary up; Secondary down; and Tertiary up

* Implementation of a new formulation (integer clustering) for power plant related constraints. This formulation divides the simulation time by a factor higher than 10 and allows extending the geographical scope of Dispa-SET. There are now four standard model formulations, which can be run with the same input data:
	- Standard formulation: low capacity or highly flexible units are merged
	- No clustring: all units are considered individually
	- LP clustering: all units are aggregated by technology and binary constraints are removed
	- Integer clustering: a representative unit is considered for each technology and multiplied N times.

* Improved pre-processing:
	- Improved log message during input data checks
	- New config files to test the different clustering methods
	- Added functions to perform parametric studies
	- Example scripts for Monte Carlo analyses using lating hypercube samplings

* Improved post-processing:
	- Netting interconnections in dispatch plots
 	- New colour palette and polished dispatch plot
	- New fuels included
	- Improved representation of curtailment

* External dependencies:
	- Removed pre-compiled libraries for unix systems
	- Use of the low-level GAMS API (https://github.com/kavvkon/gams-api)

* Python 3.7: 
	- Dispa-SET now runs exclusively on Python 3.7. 
	- The compatibility with previous Python versions (2.7, 3.6) is not guaranteed anymore.

* Miscellaneous:
	- Unit tests on travis (https://travis-ci.org/energy-modelling-toolkit/Dispa-SET)
	- Bug fixes

Version 2.2
^^^^^^^^^^^

* Inclusion of CHP, power2heat and thermal storage (these new features can be tested by running the config file for Cyprus: 'ConfigFiles/ConfigCY.xlsx')

* Bug fixes

* Improved user interface


Version 2.1
^^^^^^^^^^^

* Major refactoring of the folder structure

* New data included in the database

* Inclusion of the LP formulation (in addition to the MILP)


Version 2.0
^^^^^^^^^^^

First public version of the Dispa-SET model.

.. _Dispa-SET v2.4: https://github.com/energy-modelling-toolkit/Dispa-SET/archive/v2.4.zip
.. _Dispa-SET v2.3: https://github.com/energy-modelling-toolkit/Dispa-SET/archive/v2.3.zip
.. _Dispa-SET v2.2: https://github.com/energy-modelling-toolkit/Dispa-SET/archive/v2.2.zip
.. _Dispa-SET v2.1: https://github.com/energy-modelling-toolkit/Dispa-SET/archive/v2.1.zip
.. _Dispa-SET v2.0: https://github.com/energy-modelling-toolkit/Dispa-SET/archive/v2.0.zip



