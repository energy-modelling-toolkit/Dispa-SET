.. _data:

Input Data
==========

In this section, "Input Data" refers to the data stored in the Dispa-SET database. The format of this data is pre-defined and imposed, in such a way that it can be read by the pre-processing tool.

Two important preliminary comments should be formulated:

* All the time series should be registered with their timestamps (e.g. '2013-02-20 02:00:00') or with a numerical index. Dispa-SET will issue an error if the day is located before the month. It is also advised to remove all time zone information from the time stamps. If the index is an integer, Dispa-SET will only recognize it if contains 8760 elements (one full year) or if it has exactly the same length as the simulation horizon.
* Although the optimisation model is designed to run with any technology or fuel name, the pre-processing and the post-processing tools of Dispa-SET use some hard-coded values. The Dispa-SET database should also comply with this convention (described in the next sections). Any non-recognized technology or fuel will be discarded in the pre-processing.

Dispa-SET configuration
-----------------------

The Dispa-SET model uses YAML configuration files to define simulation parameters. You can create and edit these configuration files using the built-in configuration editor, which provides a user-friendly interface for managing all simulation parameters.

To launch the configuration editor, use the following command in your terminal:

.. code-block:: bash

   python -m dispaset.preprocessing.config_editor

This will start a local web server and open the configuration editor in your default browser. The editor organizes parameters into logical sections for easier navigation and editing.

The table below provides a comprehensive list of all configuration parameters available in Dispa-SET, organized by their type, section in the configuration editor, whether they have default values, and if they're required.

.. list-table:: Dispa-SET Configuration Parameters
   :header-rows: 1
   :widths: 20 10 30 15 10 10

   * - Parameter
     - Type
     - Description
     - Section
     - Default value?
     - Required?
   * - AllowCurtailment
     - Boolean
     - Authorize renewable technologies to operate below their Availability Factor (AF).
     - General
     - No
     - Yes
   * - BoundarySectorData
     - File path
     - Table with the list of boundary sectors and their data.
     - Sector Coupling
     - No
     - No
   * - BoundarySectorInterconnections
     - File path
     - Historical flows between the boundary sectors.
     - Sector Coupling
     - No
     - No
   * - BoundarySectorMaxSpillage
     - File path
     - Maximum spillage in a boundary sector (in the units of the boundary sector).
     - Sector Coupling
     - No
     - No
   * - BoundarySectorNTC
     - File path
     - Table with the interconnection capacities between boundary sectors.
     - Sector Coupling
     - No
     - No
   * - CostCurtailment
     - File path
     - Cost of curtailment in EUR/MWh
     - Cost and Fuel Price
     - Yes
     - No
   * - CostLoadShedding
     - File path
     - Cost of Load shedding in EUR/MWh
     - Cost and Fuel Price
     - Yes
     - No
   * - CostOfSpillage
     - File path
     - Cost of spillage in Eur/MWh
     - Cost and Fuel Price
     - Yes
     - No
   * - CostXNotServed
     - File path
     - Cost of Userved Energy in each boundary Sector
     - Cost and Fuel Price
     - No
     - No
   * - OptimalityGap
     - Numerical
     - Relative optimality gap for the MILP solver (e.g., 0.0006 for 0.06%).
     - General
     - No
     - Yes
   * - DataTimeStep
     - Numerical
     - Time step resolution of the input time series data, in hours (e.g., 1.0, 0.5).
     - General
     - No
     - Yes
   * - Demand
     - File path
     - Demand data file with each zone as column headers. Path relative to Database directory or absolute. Use '##' for zone placeholder.
     - Time Series Data
     - No
     - Yes
   * - DemandFlexibility
     - Numerical
     - Number of hours of equivalent demand shifting storage capacity (used with ShareOfFlexibleDemand).
     - Time Series Data
     - Yes
     - No
   * - Description
     - String
     - A brief description of the simulation case.
     - General
     - No
     - No
   * - FFRGainLimit
     - File path
     - Max gain for aFFR
     - Reserve Parameters
     - Yes
     - No
   * - FFRLimit
     - File path
     - ????
     - Reserve Parameters
     - Yes
     - No
   * - FrequencyStability
     - String
     - Method or parameters for frequency stability constraints (currently placeholder).
     - Reserve Parameters
     - No
     - No
   * - GAMS_folder
     - File path
     - Path to the GAMS installation folder, if not set in the system's GAMSPATH environment variable.
     - General
     - No
     - No
   * - GeoData
     - File path
     - Geographical coordinates or data for zones/units (e.g., for plotting). Path relative to Database directory or absolute.
     - Spatial Data
     - No
     - No
   * - GridData
     - File path
     - Path the csv table containing the Power Transfer Distribution Factor (PTDF) matrix for DC-Power Flow simulations.
     - Spatial Data
     - No
     - No
   * - HorizonLength
     - Numerical (integer)
     - Length of the optimization horizon for each rolling window, in days.
     - General
     - No
     - Yes
   * - HydroScheduling
     - Off/Zonal/Regional
     - Type of mid-term hydro scheduling optimization (Off=none, Zonal=per zone, Regional=per defined region).
     - Hydro Parameters
     - No
     - Yes
   * - HydroSchedulingHorizon
     - Annual/Stop-date driven
     - Optimization horizon for the mid-term hydro scheduling (determines target reservoir levels).
     - Hydro Parameters
     - No
     - Yes
   * - InertiaLimit
     - Numerical
     - ???
     - Reserve Parameters
     - No
     - No
   * - InitialFinalReservoirLevel
     - Boolean
     - If 1 (True), enforces initial and final reservoir levels based on ReservoirLevelInitial/Final defaults or historical ReservoirLevels file. If 0 (False), these constraints are relaxed.
     - Hydro Parameters
     - No
     - Yes
   * - Interconnections
     - File path
     - Table with the historical physical flows between zones (positive flow from first zone to second). Used for flows to/from non-simulated zones (Rest of World). Path relative to Database directory or absolute.
     - Spatial Data
     - No
     - Yes
   * - LoadShedding
     - File path
     - Time series with the fraction of allowed load shedding for each zone
     - Time Series Data
     - Yes
     - No
   * - LookAhead
     - Numerical (integer)
     - Length of the perfect foresight period within the rolling horizon, in days (must be <= HorizonLength).
     - General
     - No
     - Yes
   * - NTC
     - File path
     - Table with the Net Transfer Capacities (NTCs) between simulated zones (e.g., 'Z1 -> Z2'). Path relative to Database directory or absolute.
     - Spatial Data
     - No
     - Yes
   * - Outages
     - File path
     - Table with the outage time series for each unit (0=available, 1=full outage, 0.5=half outage). Path relative to Database directory or absolute. Use '##' for unit name placeholder.
     - Time Series Data
     - No
     - Yes
   * - PowerPlantData
     - File path
     - Main table defining power plant units and their technical/economic parameters. Path relative to Database directory or absolute.
     - Unit Data
     - No
     - Yes
   * - PriceOfAmmonia
     - File path
     - Price of Ammonia fuel in EUR/MWh. Can be a file path (time series per zone) or a single value if default is used.
     - Cost and Fuel Price
     - Yes
     - No
   * - PriceOfBiomass
     - File path
     - Table with the cost time series
     - Cost and Fuel Price
     - Yes
     - No
   * - PriceOfBlackCoal
     - File path
     - Table with the cost time series
     - Cost and Fuel Price
     - Yes
     - No
   * - PriceOfCO2
     - File path
     - Table with the cost time series
     - Cost and Fuel Price
     - Yes
     - No
   * - PriceOfFuelOil
     - File path
     - Table with the cost time series
     - Cost and Fuel Price
     - Yes
     - No
   * - PriceOfGas
     - File path
     - Table with the cost time series
     - Cost and Fuel Price
     - Yes
     - No
   * - PriceOfLignite
     - File path
     - Table with the cost time series
     - Cost and Fuel Price
     - Yes
     - No
   * - PriceOfNuclear
     - File path
     - Table with the cost time series
     - Cost and Fuel Price
     - Yes
     - No
   * - PriceOfPeat
     - File path
     - Price of Peat fuel in EUR/MWh. Can be a file path (time series per zone) or a single value if default is used.
     - Cost and Fuel Price
     - Yes
     - No
   * - PriceTransmission
     - File path
     - Cost associated with power transmission (e.g., losses, fees). Can be a file path or a single value if default is used.
     - Cost and Fuel Price
     - Yes
     - No
   * - PrimaryReserveLimit
     - File path
     - Primary reserve requirement (e.g., FCR). Can be a file path (time series per zone) or a single value if default is used.
     - Reserve Parameters
     - No
     - No
   * - RenewablesAF
     - File path
     - Table with the availability factors (AF) of all renewable units (value between 0 and 1). Path relative to Database directory or absolute. Use '##' for zone/unit placeholder.
     - Time Series Data
     - No
     - Yes
   * - Reserve2D
     - File path
     - Table with the downwards secondary reserve requirements (e.g., aFRR down) for each zone. Used if ReserveCalculation=Exogenous. Path relative to Database directory or absolute.
     - Reserve Parameters
     - No
     - No
   * - Reserve2U
     - File path
     - Table with the upwards secondary reserve requirements (e.g., aFRR up) for each zone. Used if ReserveCalculation=Exogenous. Path relative to Database directory or absolute.
     - Reserve Parameters
     - No
     - No
   * - ReserveCalculation
     - Generic/Percentage/Probabilistic/Exogenous
     - Method for the reserve calculation (exogenous requires Reserve2D and Reserve2U to be provided)
     - Reserve Parameters
     - No
     - Yes
   * - ReserveParticipation
     - String List
     - List of conventional/storage unit technologies participating in secondary reserves.
     - Reserve Parameters
     - No
     - Yes
   * - ReserveParticipation_CHP
     - String List
     - List of CHP unit technologies participating in secondary reserves.
     - Reserve Parameters
     - No
     - No
   * - ReservoirLevelFinal
     - Numerical
     - Target final reservoir level as a fraction of capacity (0 to 1). Used if InitialFinalReservoirLevel=1 and ReservoirLevels file is not provided.
     - Hydro Parameters
     - Yes
     - No
   * - ReservoirLevelInitial
     - Numerical
     - Target initial reservoir level as a fraction of capacity (0 to 1). Used if InitialFinalReservoirLevel=1 and ReservoirLevels file is not provided.
     - Hydro Parameters
     - Yes
     - No
   * - ReservoirLevels
     - File path
     - Time series of target reservoir levels (fraction of capacity). Used as minimum target at the end of each rolling horizon optimization if InitialFinalReservoirLevel=1. Path relative to Database directory or absolute.
     - Hydro Parameters
     - No
     - No
   * - ReservoirScaledInflows
     - File path
     - Natural inflows to hydro reservoirs, scaled by turbine capacity (MWh of inflow per hour / MW of turbine capacity). Path relative to Database directory or absolute. Use '##' for unit placeholder.
     - Hydro Parameters
     - No
     - No
   * - SectorXDemand
     - File path
     - Demand time series for each boundary sector. Path relative to Database directory or absolute.
     - Sector Coupling
     - No
     - No
   * - SectorXFlexibleDemand
     - File path
     - Flexible demand time series for each boundary sector. Path relative to Database directory or absolute.
     - Sector Coupling
     - No
     - No
   * - SectorXFlexibleSupply
     - File path
     - Flexible (or volume-based) supply to a boundary sector. Path relative to Database directory or absolute.
     - Sector Coupling
     - No
     - No
   * - SectorXFloodControl
     - File path
     - Level of the storage in the boundary sector before flood control is activated. Path relative to Database directory or absolute.
     - Sector Coupling
     - No
     - No
   * - SectorXReservoirLevels
     - File path
     - Time series with the state of charge of the storage for each boundary sector. Path relative to Database directory or absolute.
     - Sector Coupling
     - No
     - No
   * - ShareOfFlexibleDemand
     - File path
     - Fraction of the demand that can be time-shifted (used with DemandFlexibility). Can be a file path (time series per zone) or a single value if default is used.
     - Time Series Data
     - Yes
     - No
   * - ShareOfQuickStartUnits
     - Numerical
     - Share of units considered "quick start" for reserve calculations or other rules.
     - Reserve Parameters
     - Yes
     - No
   * - SimulationDirectory
     - File path
     - Folder where the temporary simulation data and model files are copied before solving. Path relative to project root or absolute.
     - General
     - No
     - Yes
   * - SimulationTimeStep
     - Numerical
     - Time step resolution of the optimization model, in hours (e.g., 1.0). Should be >= DataTimeStep.
     - General
     - No
     - Yes
   * - SimulationType
     - LP/LP clustered/Integer clustering/Standard/No clustering
     - Formulation of the optimization problem (Mixed Integer Linear Program, Linear Program, variations with clustering).
     - General
     - No
     - Yes
   * - StartDate
     - Date tuple (YYYY, M, D, H, MIN, S)
     - Starting date and time of the simulation period. Must match the start of the input time series data.
     - General
     - No
     - Yes
   * - StopDate
     - Date tuple (YYYY, M, D, H, MIN, S)
     - Stopping date and time of the simulation period. Must match the end of the input time series data.
     - General
     - No
     - Yes
   * - StorageAlertLevels
     - File path
     - Time series with the alert levels (e.g., minimum operational level as fraction of capacity) for each storage unit. Path relative to Database directory or absolute.
     - Hydro Parameters
     - No
     - No
   * - StorageFloodControl
     - File path
     - Time series with the flood control levels (e.g., maximum operational level as fraction of capacity) for each storage unit. Path relative to Database directory or absolute.
     - Hydro Parameters
     - No
     - No
   * - SystemGainLimit
     - File path
     - Limit on system gain (related to frequency control/reserves). Can be a file path or a single value if default is used.
     - Reserve Parameters
     - Yes
     - No
   * - TransmissionGridType
     - NTC/DC-Power Flow
     - Method for modeling the transmission grid constraints (Net Transfer Capacity or DC Power Flow approximation using PTDF matrix).
     - Spatial Data
     - No
     - No
   * - ValueOfLostLoad
     - Numerical
     - Value of Lost Load (VOLL) in EUR/MWh, used for economic calculations or potentially as an alternative to CostLoadShedding.
     - Cost and Fuel Price
     - Yes
     - No
   * - WaterValue
     - Numerical
     - Shadow price or value of water stored in reservoirs (EUR/MWh), used in hydro scheduling or as an opportunity cost.
     - Hydro Parameters
     - Yes
     - No
   * - modifiers.Demand
     - Numerical
     - Multiplier applied to the demand time series (e.g., 1.1 for +10%).
     - Data Adjustment
     - No
     - No
   * - modifiers.Solar
     - Numerical
     - Multiplier applied to the solar availability factors.
     - Data Adjustment
     - No
     - No
   * - modifiers.Storage
     - Numerical
     - Multiplier applied to storage capacities.
     - Data Adjustment
     - No
     - No
   * - modifiers.Wind
     - Numerical
     - Multiplier applied to the wind availability factors.
     - Data Adjustment
     - No
     - No
   * - mts_zones
     - String List
     - List of zones included in the mid-term hydro scheduling optimization. Must be a subset of 'zones'.
     - Zones
     - No
     - Yes
   * - zones
     - String List
     - List of all power zones included in the simulation.
     - Zones
     - No
     - Yes

Understanding the Configuration Editor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The configuration editor is organized into several tabs or sections to help you manage different aspects of your simulation:

1.  **General**: Basic simulation parameters including time horizon, description, solver settings, and simulation type.
2.  **Zones**: Definition and configuration of geographical zones included in the simulation and mid-term scheduling.
3.  **Data Paths**: Paths to input files for time series data (demand, renewables, outages), unit data, spatial data (NTC, grid), hydro data, and sector coupling data.
4.  **Parameters**: Economic parameters like costs (curtailment, load shedding) and fuel prices. Includes default value settings.
5.  **Hydro**: Specific parameters for hydro power modeling, including scheduling methods, reservoir level targets, and water values.
6.  **Reserves**: Settings for system reserves (primary, secondary), calculation methods, and unit participation rules.
7.  **Sector Coupling**: Configuration for interactions between the power sector and other energy sectors (paths to boundary sector data).
8.  **Advanced**: Additional settings like data adjustment multipliers and potentially less commonly used parameters.

Parameters marked as "Required" must have a value for the simulation to run properly. Parameters with default values will use the specified default (shown in parentheses in the table) if no explicit value or file path is provided in the main configuration section. File paths should typically be relative to the main Dispa-SET Database directory (using '##' as placeholders where applicable) or absolute paths. Date tuples should be provided as (Year, Month, Day, Hour, Minute, Second).


Technologies
------------

The Dispa-SET input distinguishes between the technologies defined in the table below. The VRES column indicates the variable renewable technologies (set "tr" in the optimisation) and the Storage column indicates the technologies which can accumulate energy. 

.. table:: Dispa-SET technologies

   +---------------+-------------------------------------------+-------+--------+
   |Technology     |Description                                |VRES   |Storage |
   +===============+===========================================+=======+========+
   |**Power only**                                                              |
   +---------------+-------------------------------------------+-------+--------+
   |HDAM           |Conventional hydro dam                     |N      |Y       |
   +---------------+-------------------------------------------+-------+--------+
   |HROR           |Hydro run-of-river                         |Y      |N       |
   +---------------+-------------------------------------------+-------+--------+
   |HPHS           |Pumped hydro storage                       |N      |Y       |
   +---------------+-------------------------------------------+-------+--------+
   |PHOT           |Solar photovoltaic                         |Y      |N       |
   +---------------+-------------------------------------------+-------+--------+
   |WAVE           |Wave energy                                |Y      |N       |
   +---------------+-------------------------------------------+-------+--------+
   |WHEN           |Waste heat engine                          |N      |N       |
   +---------------+-------------------------------------------+-------+--------+
   |WTOF           |Offshore wind turbine                      |Y      |N       |
   +---------------+-------------------------------------------+-------+--------+
   |WTON           |Onshore wind turbine                       |Y      |N       |
   +---------------+-------------------------------------------+-------+--------+
   |**Combined heat and power**                                                 |
   +---------------+-------------------------------------------+-------+--------+
   |COMC           |Combined cycle                             |N      |N       |
   +---------------+-------------------------------------------+-------+--------+
   |GTUR           |Gas turbine                                |N      |N       |
   +---------------+-------------------------------------------+-------+--------+
   |ICEN           |Internal combustion engine                 |N      |N       |
   +---------------+-------------------------------------------+-------+--------+
   |SCSP           |Concentrated Solar Power                   |Y      |Y       |
   +---------------+-------------------------------------------+-------+--------+
   |STUR           |Steam turbine                              |N      |N       |
   +---------------+-------------------------------------------+-------+--------+
   |**Storage**                                                                 |
   +---------------+-------------------------------------------+-------+--------+
   |BATS           |Stationary batteries                       |N      |Y       |
   +---------------+-------------------------------------------+-------+--------+
   |BEVS           |Battery-powered electric vehicles          |N      |Y       |
   +---------------+-------------------------------------------+-------+--------+
   |CAES           |Compressed air energy storage              |N      |Y       |
   +---------------+-------------------------------------------+-------+--------+
   |P2GS           |Power-to-gas storage                       |N      |Y       |
   +---------------+-------------------------------------------+-------+--------+
   |THMS           |Thermal storage                            |N      |Y       |
   +---------------+-------------------------------------------+-------+--------+
   |**Heat only**                                                               |
   +---------------+-------------------------------------------+-------+--------+
   |GETH           |Geothermal district heating                |Y      |N       |
   +---------------+-------------------------------------------+-------+--------+
   |HOBO           |Heat only boiler                           |N      |N       |
   +---------------+-------------------------------------------+-------+--------+
   |SOTH           |Solar thermal district heating             |Y      |N       |
   +---------------+-------------------------------------------+-------+--------+
   |**Power to heat**                                                           |
   +---------------+-------------------------------------------+-------+--------+
   |ABHP           |Absorption heat pump (solar/geothermal/gas)|Y/N    |N       |
   +---------------+-------------------------------------------+-------+--------+
   |ASHP           |Air source heat pump                       |Y/N    |N       |
   +---------------+-------------------------------------------+-------+--------+
   |GSHP           |Ground source heat pump                    |Y/N    |N       |
   +---------------+-------------------------------------------+-------+--------+
   |HYHP           |Hybrid heat pump (Ground/air & HP/GAS-OIL  |Y/N    |N       |
   +---------------+-------------------------------------------+-------+--------+
   |P2HT           |Power-to-heat                              |Y/N    |N       |
   +---------------+-------------------------------------------+-------+--------+
   |REHE           |Resistive heater                           |Y/N    |N       |
   +---------------+-------------------------------------------+-------+--------+
   |WSHP           |Water source heat pump                     |Y/N    |N       |
   +---------------+-------------------------------------------+-------+--------+

Fuels
-----

Dispa-SET only considers a limited number of different fuel types. They are summarised in the following table, together with some examples.

.. table:: Dispa-SET fuels

	======= =============
	Fuel	Examples
	======= =============
	AIR     Air energy from the surrounding environment (used by heat pumps and other heat generation technologies)
	AMO	    Ammonia
	BIO	    Bagasse, Biodiesel, Gas From Biomass, Gasification, Biomass, Briquettes, Cattle Residues, Rice Hulls Or Padi Husk, Straw, Wood Gas (From Wood Gasification), Wood Waste Liquids Excl Blk Liq (Incl Red Liquor, Sludge, Wood,Spent Sulfite Liquor And Oth Liquids, Wood And Wood Waste
	GAS	    Blast Furnace Gas, Boiler Natural Gas, Butane, Coal Bed Methane, Coke Oven Gas, Flare Gas, Gas (Generic), Methane, Mine Gas, Natural Gas, Propane, Refinery Gas, Sour Gas, Synthetic Natural Gas, Top Gas, Voc Gas & Vapor, Waste Gas, Wellhead Gas
	GEO	    Geothermal steam
	HRD	    Anthracite, Other Anthracite, Bituminous Coal, Coker By-Product, Coal Gas (From Coal Gasification), Coke, Coal (Generic), Coal-Oil Mixture, Other Coal, Coal And Pet Coke Mi, Coal Tar Oil, Anthracite Coal Waste, Coal-Water Mixture, Gob, Hard Coal / Anthracite, Imported Coal, Other Solids, Soft Coal, Anthracite Silt, Steam Coal, Subbituminous, Pelletized Synthetic Fuel From Coal, Bituminous Coal Waste)
	HYD	    Hydrogen
	LIG	    Lignite black, Lignite brown, lignite
	NUC	    U (Uranium), Pu (Plutonium)
	OIL	    Crude Oil, Distillate Oil, Diesel Fuel, No. 1 Fuel Oil, No. 2 Fuel Oil, No. 3 Fuel Oil, No. 4 Fuel Oil, No. 5 Fuel Oil, No. 6 Fuel Oil, Furnace Fuel, Gas Oil, Gasoline, Heavy Oil Mixture, Jet Fuel, Kerosene, Light Fuel Oil, Liquefied Propane Gas, Methanol, Naphtha, ,Gas From Fuel Oil Gasification, Fuel Oil, Other Liquid, Orimulsion, Petroleum Coke, Petroleum Coke Synthetic Gas, Black Liquor, Residual Oils, Re-Refined Motor Oil, Oil Shale, Tar, Topped Crude Oil, Waste Oil
	OTH     All other energy carriers
	PEA	    Peat Moss
	SUN	    Solar energy
	WAT	    Hydro energy
	WIN	    Wind energy
	WST	    Digester Gas (Sewage Sludge Gas), Gas From Refuse Gasification, Hazardous Waste, Industrial Waste, Landfill Gas, Poultry Litter, Manure, Medical Waste, Refused Derived Fuel, Refuse, Waste Paper And Waste Plastic, Refinery Waste, Tires, Agricultural Waste, Waste Coal, Waste Water Sludge, Waste
	WHT     Waste heat, Excess heat
	======= =============

Different fuels may be used to power a given technology, e.g. steam turbines may be fired with almost any fuel type. In Dispa-SET, each unit must be defined with the pair of values (technology,fuel). The next tables is derived from a commercial power plant database and indicates the number of occurences of each combination. It appears clearly that, even through some combinations are irrelevant, both characteristics are needed to define a power plant type.

======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ==========
f/t	COMC	GTUR	HDAM	HPHS	HROR	ICEN	PHOT	STUR	WTOF	WTON	Total
======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ==========
AMO	1	1
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

	================================================ =============== ===========
	Description                                      Field name      Units
	================================================ =============== ===========
	Unit name                                        Unit            n.a.
	Installed Power or Heat Capacity (for one unit)  PowerCapacity   MW		
	Number of thermal blocks belonging to one unit   Nunits          n.a.
	Technology                                       Technology      n.a.	
	Primary fuel                                     Fuel            n.a.		
	Zone (Power)                                     Zone            n.a.		
	Zone (Heat)                                      Zone_th         n.a.
	Efficiency                                       Efficiency      %
	Efficiency at minimum load                       MinEfficiency   %
	CO2 intensity                                    CO2Intensity    TCO2/MWh
	Minimum load                                     PartLoadMin     %
	Ramp up rate                                     RampUpRate      %/min
	Ramp down rate                                   RampDownRate    %/min)
	Start-up time                                    StartUPTime     h
	Minimum up time                                  MinUpTime       h
	Minimum down time                                MinDownTime     h
	No load cost                                     NoLoadCost      EUR/h
	Start-up cost                                    StartUpCost     EUR
	Ramping cost                                     RampingCost     EUR/MW
	================================================ =============== ===========


NB: the fields indicated with % as unit must be entered in a non-dimensional way (i.e. 90% should be written 0.9).

Storage units
^^^^^^^^^^^^^

Some parameters must only be defined for the units equipped with storage. They can be left blank for all other units.

.. table:: Specific fields for storage units

	=============================== =======================	===========
	Description			Field name		Units
	=============================== =======================	===========
	Storage capacity 		STOCapacity		MWh
	Self-discharge rate		STOSelfDischarge	%/d
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
    % of storage heat losses per day          STOSelfDischarge   %/d
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
    % of storage heat losses per day	      STOSelfDischarge   %/d
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

Storage units are an extension of the regular units, including additional constraints and parameters. In the power plant table, four additional parameters are required: storage capacity (in MWh), self-discharge (in %/d), discharge power (in MW) and discharge efficiency (in %). 

Some other parameters must be introduced in the form of time series in the "HydroData" section of the Dispa-SET database. There are described hereunder.

It should be noted that the nomenclature adopted for the modeling of storage units refers to the characteristics of hydro units with water reservoirs. However, these parameters (e.g. inflows, level) can easily be transposed to the case of alternative storage units such as batteries or CAES.

Inflows
^^^^^^^
The Inflows are defined as the contribution of exogenous sources to the level (or state of charge) or the reservoir. They are expressed in MWh of potential energy. If the inflows are provided as m³/h, they must be converted.

The input to dispaset is defined as "StorageInflows". It is the normalized values of the inflow with respect to the nominal power of the storage unit (in discharge mode). As an example, if the inflow value at a certain time is 100MWh/h and if the turbining capacity of the hydro plant is 200 MW, the scaled inflow value must be defined as 0.5.

Scaled inflows should be provided in the form of time series with the unit name or the technology as columns header.


Storage level
^^^^^^^^^^^^^
Because emptying the storage has a zero marginal cost, a non-constrained optimization tends to leave the storage completely empty at the end of the optimisation horizon. For that reason, a minimum storage level is imposed at the last hour of each horizon. In Dispa-SET, a typical optimisation horizon is a few days. The model is therefore not capable of optimising the storage level e.g. for seasonal variations. The minimum storage level at the last hour is therefore an exogenous input. It can be selected from a historical level or obtained from a long-term hydro scheduling optimization.

The level input in the Dispa-SET database is normalized with respect to the storage capacity: its minimum value is zero and its maximum is one. 

Variable capacity storage
^^^^^^^^^^^^^^^^^^^^^^^^^
In special cases, it might be necessary to simulate a storage unit whose capacity varies in time. A typical example is the simulation of the storage capacity provided by electric vehicles: depending on the time of the day, the connected battery capacity varies. 

This special case can be simulated using the "AvailabilityFactor" input. In the case of a storage unit, reduces the available capacity by a factor varying from 0 to 1. 

Other storage units
-------------------
Other storage units include H2 storage, batteries (BATS) and electric vehicles (BEVS). In case of H2 storage, the parameter StorageInflow are defined null at all times whereas StorageOutflow corresponds to the hydrogen demand at each timsestep. For batteries and BEVS, both parameters are set to 0 all the time. 

Power plant outages
-------------------
In the current version, Dispa-SET does not distinguish planned outages from unplanned outages. They are characterized for each unit by the "OutageFactor" parameter. This parameter varies from 0 (no outage) to 1 (full outage). The available unit power is thus given by its nominal capacity multiplied by (1-OutageFactor). 

The outages are provided in the dedicated section of the Database for each unit. They consist of a time series with the unit name as columns header.


Interconnections
----------------

Two cases should be distinguished when considering interconnections:

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





