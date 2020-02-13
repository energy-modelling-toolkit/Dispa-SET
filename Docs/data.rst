.. _data:

Input Data
==========

In this section, "Input Data" refers to the data stored in the Dispa-SET database. The format of this data is pre-defined and imposed, in such a way that it can be read by the pre-processing tool.

Two important preliminary comments should be formulated:

* All the time series should be registered with their timestamps (e.g. '2013-02-20 02:00:00') or with a numerical index. Dispa-SET will issue an error if the day is located before the month. It is also advised to remove all time zone information from the time stamps. If the index is an integer, Dispa-SET will only recognize it if contains 8760 elements (one full year) or if it has exactly the same length as the simulation horizon.
* Although the optimisation model is designed to run with any technology or fuel name, the pre-processing and the post-processing tools of Dispa-SET use some hard-coded values. The Dispa-SET database should also comply with this convention (described in the next sections). Any non-recognized technology or fuel will be discarded in the pre-processing.

General simulation parameters
-----------------------------
A number of simulation options and parameters need to be defined in the configuration file. In order to obtain default values and a complete description of the options, it is commended to open the ConfigTest.xlsx configuration file, which is always kept up-to-date.

The options to be filled are summarized  hereunder.

.. table:: Dispa-SET Simulation Options

	=============================== ================================================== 
	General Options			Description				
	=============================== ================================================== 
	SimulationDirectory		Folder with all simulation files and input data
	WriteGDX			Write the inputs in a GDX file (required for gams)
	WritePickle			Write the inputs to a pickle file
	GAMS_folder			Path the GAMS installation folder
	cplex_path			Path to the cplex folder
	**Horizon Settings**	
	StartDate			Start date of the simulation 
	StopDate			End data of the simulation
	HorizonLength			Simulation horizon length in days
	Look ahead			Overlap period in days
	DataTimeStep			Time step of the date in the csv files
	SimulationTimeStep		Time step for the simulation
	**Simulation Options**
	SimulationType			Stanard/LP/LP clustered/Integer clustering
	ReserveCalculation		Generic (only available option for now)
	AllowCurtailment		True/False
	**Mid-term scheduling**
	HydroScheduling			Off/Zonal/Regional
	HydroSchedulingHorizon		"Annual"/"Stop-date driven"
	InitialFinalReservoirLevel	True/False (if False, use StorageProfile)
	ReservoirLevelInitial		Initial res. level if the above option is true
	ReservoirLevelFinal		Fainl reservoir level if the above option is true
	=============================== ================================================== 


Technologies
------------

The Dispa-SET input distinguishes between the technologies defined in the table below. The VRES column indicates the variable renewable technologies (set "tr" in the optimisation) and the Storage column indicates the technologies which can accumulate energy. 

.. table:: Dispa-SET technologies

	=============== ======================================= ======= ========
	Technology	Description				VRES	Storage
	=============== ======================================= ======= ========
	COMC		Combined cycle				N	N
	GTUR		Gas turbine				N	N
	HDAM		Conventional hydro dam			N	Y
	HROR		Hydro run-of-river			Y	N
	HPHS		Pumped hydro storage			N	Y
	ICEN 		Internal combustion engine		N	N
	PHOT		Solar photovoltaic			Y	N
	STUR		Steam turbine				N	N
	WTOF		Offshore wind turbine			Y	N
	WTON		Onshore wind turbine			Y	N
	CAES		Compressed air energy storage		N	Y
	BATS		Stationary batteries			N	Y
	BEVS		Battery-powered electric vehicles	N	Y
	THMS		Thermal storage				N	Y
	P2GS		Power-to-gas storage			N	Y
	P2HT		Power-to-heat				N	Y
	=============== ======================================= ======= ========

Fuels
-----

Dispa-SET only considers a limited number of different fuel types. They are summarised in the following table, together with some examples.

.. table:: Dispa-SET fuels

	======= =============
	Fuel	Examples
	======= =============
	BIO	Bagasse, Biodiesel, Gas From Biomass, Gasification, Biomass, Briquettes, Cattle Residues, Rice Hulls Or Padi Husk, Straw, Wood Gas (From Wood Gasification), Wood Waste Liquids Excl Blk Liq (Incl Red Liquor, Sludge, Wood,Spent Sulfite Liquor And Oth Liquids, Wood And Wood Waste
	GAS	Blast Furnace Gas, Boiler Natural Gas, Butane, Coal Bed Methane, Coke Oven Gas, Flare Gas, Gas (Generic), Methane, Mine Gas, Natural Gas, Propane, Refinery Gas, Sour Gas, Synthetic Natural Gas, Top Gas, Voc Gas & Vapor, Waste Gas, Wellhead Gas
	GEO	Geothermal steam
	HRD	Anthracite, Other Anthracite, Bituminous Coal, Coker By-Product, Coal Gas (From Coal Gasification), Coke, Coal (Generic), Coal-Oil Mixture, Other Coal, Coal And Pet Coke Mi, Coal Tar Oil, Anthracite Coal Waste, Coal-Water Mixture, Gob, Hard Coal / Anthracite, Imported Coal, Other Solids, Soft Coal, Anthracite Silt, Steam Coal, Subbituminous, Pelletized Synthetic Fuel From Coal, Bituminous Coal Waste)
	HYD	Hydrogen
	LIG	Lignite black, Lignite brown, lignite
	NUC	U, Pu
	OIL	Crude Oil, Distillate Oil, Diesel Fuel, No. 1 Fuel Oil, No. 2 Fuel Oil, No. 3 Fuel Oil, No. 4 Fuel Oil, No. 5 Fuel Oil, No. 6 Fuel Oil, Furnace Fuel, Gas Oil, Gasoline, Heavy Oil Mixture, Jet Fuel, Kerosene, Light Fuel Oil, Liquefied Propane Gas, Methanol, Naphtha, ,Gas From Fuel Oil Gasification, Fuel Oil, Other Liquid, Orimulsion, Petroleum Coke, Petroleum Coke Synthetic Gas, Black Liquor, Residual Oils, Re-Refined Motor Oil, Oil Shale, Tar, Topped Crude Oil, Waste Oil
	PEA	Peat Moss
	SUN	Solar energy
	WAT	Hydro energy
	WIN	Wind energy
	WST	Digester Gas (Sewage Sludge Gas), Gas From Refuse Gasification, Hazardous Waste, Industrial Waste, Landfill Gas, Poultry Litter, Manure, Medical Waste, Refused Derived Fuel, Refuse, Waste Paper And Waste Plastic, Refinery Waste, Tires, Agricultural Waste, Waste Coal, Waste Water Sludge, Waste
	======= =============

Different fuels may be used to power a given technology, e.g. steam turbines may be fired with almost any fuel type. In Dispa-SET, each unit must be defined with the pair of values (technology,fuel). The next tables is derived from a commercial power plant database and indicates the number of occurences of each combination. It appears clearly that, even through some combinations are irrelevant, both characteristics are needed to define a power plant type.

======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ==========
f/t	COMC	GTUR	HDAM	HPHS	HROR	ICEN	PHOT	STUR	WTOF	WTON	Total
======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ==========
BIO		2				10		79			91
GAS	485	188				28		97			798
GEO								10			10
HRD	4							389			393
HYD		1						1			2
LIG								249			249
NUC								138			138
OIL	7	94				27		146			274
PEA								17			17
SUN							20	7			27
UNK		2				1		1			4
WAT			33	23	21			1			78
WIN									9	27	36
WST		3				7		46			56
Total	496	290	33	23	21	73	20	1181	9	27	2173
======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ==========


Unit-specific or technology-specific inputs
-------------------------------------------

Some parameters, such as the availability factor, the outage factor or the inflows may be defined at the unit level or at the technology level. For that reason, the pre-processing tool first lookups the unit name in the database to assign it a value, and then lookups the technology or the fuel if no unit-specific information has been found.

Demand
------

Electricity demand is given per zone and the first row of each column with the time series should be the zone name.

Heat demand timeseries is needed where CHP or P2HT plants are used. In the current formulation, each CHP/P2HT unit is covering a heat load. In other words, one power plant is connected to a single district heating network. Therefore, in the heat demand input file, the first column has to be a time index and the following columns the heat demand in MW. The first row should contain the exact name of the power plant that will cover this demand.

It si possible to assume that a share of the demand is flexible (see model formulation for more information). In that case, this flexible share is provided as times series for each zone (see for example the tests/dummy_data/ShareFlexible.csv file), referencend in the "FlexibleDemand" field of the config file. It is also necessary to specify the number of hours of equivalent demand shifting capacity. This is achieved through the "DemandFlexibility" field of the config file and is expressed in hours (i.e. the number of hours during which the maximum flexible demand can be stored for shifting). An example of such configuration is proivded in the ConfigTest

Countries
---------
Although the nodes names can be freely user-defined in the database, for the Dispa-SET EU model, the ISO 3166-1 standard has been adopted to describe each country at the NUTS1 level (except for Greece and the United Kingdom, for which the abbreviations EL and UK are used according to `EU Interinstitutional style guide <http://publications.europa.eu/code/pdf/370000en.htm>`_ ). The list of countries is defined as:

======= =======
Code	Country
======= =======
AT	Austria
BE	Belgium
BG	Bulgaria
CH	Switzerland
CY	Cyprus
CZ	Czech Republic
DE	Germany
DK	Denmark
EE	Estonia
EL	Greece
ES	Spain
FI	Finland
FR	France
HR	Croatia
HU	Hungary
IE	Ireland
IT	Italy
LT	Lituania
LU	Luxembourg
LV	Latvia
MT	Malta
NL	Netherlands
NO	Norway
PL	Poland
PT	Portugal
RO	Romania
SE	Sweden	
SI	Slovenia
SK	Slovakia
UK      United Kingdom
======= =======


Power plant data
----------------
The power plant database may contain as many fields as desired, e.g. to ensure that the input data can be traced back, or to provide the id of this plant in another database. However, some fields are required by Dispa-SET and must therefore be defined in the database. 

Common fields
^^^^^^^^^^^^^

The following fields must be defined for all units:

.. table:: Common fields for all units

	=============================== =============== ===========
	Description			Field name	Units
	=============================== =============== ===========
	Unit name			Unit
	Power Capacity (for one unit) 	PowerCapacity	MW		
	Number of units			Nunits	
	Technology			Technology	
	Primary fuel			Fuel		
	Zone				Zone		
	Efficiency 			Efficiency	%
	Efficiency at minimum load 	MinEfficiency	%
	CO2 intensity 			CO2Intensity	TCO2/MWh
	Minimum load 			PartLoadMin	%
	Ramp up rate			RampUpRate	%/min
	Ramp down rate 			RampDownRate	%/min)
	Start-up time			StartUPTime	h
	Minimum up time 		MinUpTime	h
	Minimum down time		MinDownTime	h
	No load cost 			NoLoadCost	EUR/h
	Start-up cost 			StartUpCost	EUR
	Ramping cost			RampingCost	EUR/MW
	=============================== =============== ===========


NB: the fields indicated with % as unit must be entered in a non-dimensional way (i.e. 90% should be written 0.9).

Storage units
^^^^^^^^^^^^^

Some parameters must only be defined for the units equipped with storage. They can be left blank for all other units.

.. table:: Specific fields for storage units

	=============================== =======================	===========
	Description			Field name		Units
	=============================== =======================	===========
	Storage capacity 		STOCapacity		MWh
	Self-discharge rate		STOSelfDischarge	%/h
	Maximum charging power 		STOMaxChargingPower	MW
	Charging efficiency 		STOChargingEfficiency	%
	=============================== =======================	===========


In the case of a storage unit, the discharge efficiency should be assigned to the common field "Efficiency". Similarly, the common field "PowerCapacity" is the nominal power in discharge mode.

CHP units
^^^^^^^^^

Some parameters must only be defined for the units equipped with CHP. They can be left blank for all other units.

.. table:: Specific fields for CHP units

    ========================================= ================== ===========
    Description                               Field name         Units
    ========================================= ================== ===========
    CHP Type                                  CHPType            extraction/back-pressure/p2h
    Power-to-heat ratio                       CHPPowerToHeat     -
    Power Loss factor                         CHPPowerLossFactor -
    Maximum heat production                   CHPMaxHeat         MW(th)
    Capacity of heat Storage                  STOCapacity        MWh(th)
    % of storage heat losses per timestep     STOSelfDischarge   %
    ========================================= ================== ===========

In the current version of DispaSet three type of combined heat and power units are supported:

* Extraction/condensing units
* Backpressure units
* Power to heat 

For each of the above configurations the following fields must be filled:

.. table:: Mandatory fields per type of CHP unit (X: mandatory, o:optional)

    ================== =========== ============ =============
    Description        Extraction  Backpressure Power to heat
    ================== =========== ============ =============
    CHPType            X           X            X
    CHPPowerToHeat     X           X
    CHPPowerLossFactor X                        X
    CHPMaxHeat         o           o            X
    STOCapacity        o           o            o
    STOSelfDischarge   o           o            o
    ================== =========== ============ =============

There are numerous data checking routines to ensure that all data provided is consistent.

.. warning::
    For extraction/condensing CHP plants, the power plant capacity (*PowerCapacity*) must correspont to the nameplate capacity in the maximum heat and power mode. Internal Dispaset calculations will use the equivalent stand-alone plants capacity based on the parameters provided.


P2HT units
^^^^^^^^^^

Some parameters must only be defined for the power-to-heat units (heat pumps, electrical heaters). They can be left blank for all other units.

.. table:: Specific fields for P2HT units

    ========================================= ================== ===========
    Description                               Field name         Units
    ========================================= ================== ===========
    Nominal coefficient of performance	      COP                -
    Nominal temperature                       Tnominal           °C
    First coefficient                         coef_COP_a         -
    Second coefficient	                      coef_COP_b         - 
    Capacity of heat Storage                  STOCapacity        MWh(th)
    % of storage heat losses per timestep     STOSelfDischarge   %
    ========================================= ================== ===========

NB:

* Electrical heaters can be simulated by setting the nominal COP to 1 and the temperature coefficients to 0
* The two coefficients a and b aim at correcting the COP for the ambient temperatures. They are calculated as follows:

.. math::

	 \mathit{COP} = \mathit{COP}_{nom} + \mathit{coef}_{a} \cdot (T - T_{nom}) + \mathit{coef}_{b} \cdot (T - T_{nom})^2

where T is the atmospheric temperature which needs to be provided as a times sereis for each zone in a csv file. The first row of the csv file is the zone name and a proper time index is required. The csv file path must be provided in the "Temperatures" field of the configuration file (see ConfigTest.xlsx for an example)

.. warning::
    For power-to-heat units, the power plant capacity (*PowerCapacity*) must correspont to the nameplate nominal ELECTRICAL consumption, thus given by the thermal capacity divided by the nominal COP.


Renewable generation
--------------------
Variable renewable generation is defined as power generation from renewable source that cannot be stored: its is either fed to the grid or curtailed. The technologies falling under this definition are the ones described in the subset "tr" in the model definition. 

The time-dependent genration of for these technologies must be provided as an exogenous time series in the form of an "availability factor". The latter is defined as the proportion of the nominal power capacity that can be generated at each hour.

In the database, the time series are provided as column vectors with the technology name as header. After the pre-processing, an availability factor is attributed to each unit according to their technology. Non-renewable technologies are assigned an availability factor of 1. 



Storage and hydro data
----------------------

Storage units are an extension of the regular units, including additional constraints and parameters. In the power plant table, four additional parameters are required: storage capacity (in MWh), self-discharge (in %/h), discharge power (in MW) and discharge efficiency (in %). 

Some other parameters must be introduced in the form of time series in the "HydroData" section of the Dispa-SET database. There are described hereunder.

It should be noted that the nomenclature adopted for the modeling of storage units refers to the characteristics of hydro units with water reservoirs. However, these parameters (e.g. inflows, level) can easily be transposed to the case of alternative storage units such as batteries or CAES.

Inflows
^^^^^^^
The Inflows are defined as the contribution of exogenous sources to the level (or state of charge) or the reservoir. They are expressed in MWh of potential energy. If the inflows are provided as m³/h, they must be converted.

The input to dispaset is defined as "ScaledInflows". It is the normalized values of the inflow with respect to the nominal power of the storage unit (in discharge mode). As an example, if the inflow value at a certain time is 100MWh/h and if the turbining capacity of the hydro plant is 200 MW, the scaled inflow value must be defined as 0.5.

Scaled inflows should be provided in the form of time series with the unit name or the technology as columns header.


Storage level
^^^^^^^^^^^^^
Because emptying the storage has a zero marginal cost, a non-constrained optimization tends to leave the storage completely empty at the end of the optimisation horizon. For that reason, a minimum storage level is imposed at the last hour of each horizon. In Dispa-SET, a typical optimisation horizon is a few days. The model is therefore not capable of optimising the storage level e.g. for seasonal variations. The minimum storage level at the last hour is therefore an exogenous input. It can be selected from a historical level or obtained from a long-term hydro scheduling optimization.

The level input in the Dispa-SET database is normalized with respect to the storage capacity: its minimum value is zero and its maximum is one. 

Variable capacity storage
^^^^^^^^^^^^^^^^^^^^^^^^^
In special cases, it might be necessary to simulate a storage unit whose capacity varies in time. A typical example is the simulation of the storage capacity provided by electric vehicles: depending on the time of the day, the connected battery capacity varies. 

This special case can be simulated using the "AvailabilityFactor" input. In the case of a storage unit, reduces the avaiable capacity by a factor varying from 0 to 1. 


Power plant outages
-------------------
In the current version, Dispa-SET does not distinguish planned outages from unplanned outages. They are characterized for each unit by the "OutageFactor" parameter. This parameter varies from 0 (no outage) to 1 (full outage). The available unit power is thus given by its nominal capacity multiplied by (1-OutageFactor). 

The outages are provided in the dedicated section of the Database for each unit. They consist of a time series with the unit name as columns header.


Interconnections
----------------

Two case should be distinguished when considering interconnections:

* Interconnections occuring between the simulated zones
* Interconnections occuring between the simulated zones and the Rest of the World (RoW)

These two cases are addresses by two different datasets described here under.

Net transfer capacities
^^^^^^^^^^^^^^^^^^^^^^^
Dispa-SET indogenously models the internal exchanges between countries (or zones) using a commercial net transfer caapcity (NTC). It does not consider (yet) DC power flows or more complex grid simulations. 

Since the NTC values might vary in time, they must be supplied as time series, whose header include the origin country, the string ' -> ' and the destination country. As an example, the NTC from belgium to france must be provided with the header 'BE -> FR'. 

Because NTCs are not necessarily symetrical, they must be provided in both directions (i.e. 'BE -> FR' and 'FR -> BE'. Non-provided NTCs are considered to be zero (i.e. no interconnection).


Historical physical flows
^^^^^^^^^^^^^^^^^^^^^^^^^
In Dispa-SET, the flows between internal zones and the rest of the world cannot be modeled endogenously. They must be provided as exogenous inputs. These inputs are referred to as "Historical physical flows", although they can also be user-defined. 

In the input table of historical flows, the headers are similar to those of the NTCs (ie. 'XX -> YY'). All flows occuring an internal zone of the simulation and outside zones are considered as external flows and summed up. As an example, the historical flows 'FR -> XX', 'FR -> YY' and 'FR -> ZZ' will be aggregated in to a single interconnection flow 'FR -> RoW' if XX, YY and ZZ are not simulated zones. 

These aggregated historical flows are then imposed to the solver as exogenous inputs.

In Dispa-SET, the flows are defined as positive variables. For each zone, there will thus be a maximum of two vectors defining its exchanges with the rest of the world (e.g. 'FR -> RoW' and 'RoW -> FR').

As for the NTCs, undefined historical flows are considered to be zero, i.e. not provided any historical flows is equivalent to consider the system as islanded.


Fuel Prices
-----------
Fuel prices vary both geographically and in time. They must therefore be provided as a time series for each simulated zone. One table is provided per fuel type, with as column header the zone to which it applies. If no header is provided, the fuel price is applied to all the simulated zones.





