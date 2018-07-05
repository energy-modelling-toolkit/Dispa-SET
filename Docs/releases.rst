.. _releases:

Releases
========

Major stable releases:

* `Dispa-SET v2.3`_

* `Dispa-SET v2.2`_

* `Dispa-SET v2.1`_

* `Dispa-SET v2.0`_

Changelog
---------

Version 2.3
^^^^^^^^^^^

* Reformulation of the reserve constraints. Secondary reserves are now covered by spinning units only. Tertirary reserves can also be covered by quick start units. In total, three different reserve markets are now considered: 
	- Secondary up 
	- Secondary down
	- Tertiary up

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
	- New fuels included
	- Improved representation of curtailment

* External dependencies:
	- Removed pre-compiled libraries for unix systems
	- Use of the low-level GAMS API (https://github.com/kavvkon/py-gdx)

* Miscellaneous:
	- Bug fixes
	- Compatibility with Python 3.6

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



