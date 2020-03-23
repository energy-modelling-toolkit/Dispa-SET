.. _releases:

Releases
========

Major stable releases:

* Dispa-SET v2.4

* `Dispa-SET v2.3`_ 

* `Dispa-SET v2.2`_

* `Dispa-SET v2.1`_

* `Dispa-SET v2.0`_

Changelog
---------

Version 2.x
^^^^^^^^^^^
* Variable time step
	- The pre-processing and the GAMS file have been update to handle different time steps (not only one hour)
	- This is currently restricted to three time steps: 15min, 1h, 24h
	- The input data whose time step is lower than the desired one is averaged

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


.. _Dispa-SET v2.3: https://github.com/energy-modelling-toolkit/Dispa-SET/archive/v2.3.zip
.. _Dispa-SET v2.2: https://github.com/energy-modelling-toolkit/Dispa-SET/archive/v2.2.zip
.. _Dispa-SET v2.1: https://github.com/energy-modelling-toolkit/Dispa-SET/archive/v2.1.zip
.. _Dispa-SET v2.0: https://github.com/energy-modelling-toolkit/Dispa-SET/archive/v2.0.zip



