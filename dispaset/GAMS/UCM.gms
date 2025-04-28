$Title UCM model

$eolcom //
$offempty
Option threads=16;
Option IterLim=1000000000;
Option ResLim = 10000000000;
*Option optca=0.0;


// Reduce .lst file size

// Turn off the listing of the input file
$offlisting
$offlog

// Turn off the listing and cross-reference of the symbols used
$offsymxref
$offsymlist

option
    limrow = 0,     // equations listed per block
    limcol = 0,     // variables listed per block
    solprint = off,     // solver's solution output printed
    sysout = off;       // solver's system output printed

Option solver=gurobi;

*===============================================================================
*Definition of the dataset-related options
*===============================================================================

* Name of the input file (Ideally, stick to the default Input.gdx)
*$set InputFileName Input.gdx
$set InputFileName Inputs.gdx

* Definition of the equations that will be present in LP or MIP
* (1 for LP 0 for MIP TC)
$setglobal LPFormulation 0
$setglobal MTS 0

* Flag to retrieve status or not
* (1 to retrieve 0 to not)
$setglobal RetrieveStatus 0

* Activate the flexible demand equations
$setglobal ActivateFlexibleDemand 1

* Activate advanced reserve demand
$setglobal ActivateAdvancedReserves 0

*New
* Definition of the equations that will be present in Frequency Constrainted UC/OD
* (1 for FC-UC/OD 0 for UC/OD)
$setglobal FC 0

*===============================================================================
*Definition of   sets and parameters
*===============================================================================
SETS
*Markets
mk               Markets
*Units
au               All Units
u(au)            Generation units
chp(u)           CHP units
p2x(au)          Power to X
x2p(au)          X to Power
xu(au)           Boundary sector only units
s(au)            Storage Units (with reservoir)
th(au)           Units with thermal storage
hu(au)           Heat only units
*New
cu(au)           Conventional units only
ba(au)           Batteries only
*Fuels
f                Fuel types
*Technologies
t                Generation technologies
tr(t)            Renewable generation technologies
tc(t)
*Pollutants
p                Pollutants
*Nodes
n                Nodes
nx               Boundary sector nodes
*Lines
l                Lines
l_int(l)         Lines between internal zones
l_RoW(l)         Lines to the rest of the world
lx               Boundary sector lines
slx              Boundary sector spillage lines
*Simulations options
h                Hours
i(h)             Subset of simulated hours for one iteration
wat(au)          hydro technologies
z(h)             Subset of every simulated hour
;

Alias(mk,mkmk);
Alias(n,nn);
Alias(l,ll);
Alias(u,uu);
Alias(t,tt);
Alias(f,ff);
Alias(p,pp);
Alias(s,ss);
Alias(h,hh);
Alias(i,ii);

*Parameters as defined in the input file
* \u indicate that the value is provided for one single unit
PARAMETERS
AvailabilityFactor(au,h)                    [%]             Availability factor
CHPPowerLossFactor(u)                       [%]             Power loss when generating heat
CHPPowerToHeat(u)                           [%]             Nominal power-to-heat factor
CHPMaxHeat(chp)                             [MW\u]          Maximum heat capacity of chp plant
CHPType
CommittedInitial(au)                        [n.a.]          Initial committment status
Config
CostFixed(au)                               [EUR\h]         Fixed costs
CostRampUp(au)                              [EUR\MW]        Ramp-up costs
CostRampDown(au)                            [EUR\MW]        Ramp-down costs
CostShutDown(au)                            [EUR\u]         Shut-down costs
CostStartUp(au)                             [EUR\u]         Start-up costs
CostVariable(au,h)                          [EUR\MW]        Variable costs
CostOfSpillage(au,h)                          [EUR\MW]        Cost of spillage
CostWaterValue(au,h)                        [EUR\MW]        Cost of storage level violation for each unit
CostStorageAlert(au,h)                      [EUR\MW]        Cost of violating storage alert level
CostFloodControl(au,h)                      [EUR\MW]        Cost of violating storage flood level
CostXStorageAlert(nx,h)                     [EUR\MW]        Cost of violating storage alert level boundary sector
CostXFloodControl(nx,h)                     [EUR\MW]        Cost of violating storage flood control level boundary sector
CostXSpillage(slx,h)                        [EUR\MW]        Cost of spillage for boundary sector
CostXNotServed(nx,h)                        [EUR\MWh]       Cost of supplying energy to boundary sector via other means
CostLoadShedding(n,h)                       [EUR\MWh]       Cost of load shedding
Curtailment(n)                              [n.a]           Curtailment allowed or not {1 0} at node n
CostCurtailment(n,h)                        [EUR\MWh]       Cost of VRES curtailment
Demand(mk,n,h)                              [MW]            Demand
Efficiency(au,h)                            [%]             Efficiency
X2PowerConversionMultiplier(nx,au,h)        [%]             Discharge Efficiency boundary sector
Power2XConversionMultiplier(nx,au,h)        [%]             Charging Efficiency boundary sector
EmissionMaximum(n,p)                        [tP]            Emission limit
EmissionRate(au,p)                          [tP\MWh]        P emission rate
FlowMaximum(l,h)                            [MW]            Line limits
FlowMinimum(l,h)                            [MW]            Minimum flow
FlowXMaximum(lx,h)                          [MW]            Boundary sector line limits
FlowXMinimum(lx,h)                          [MW]            Boundary sector line minimum flow
Fuel(au,f)                                  [n.a.]          Fuel type {1 0}
SectorXDemand(nx,h)                         [MWh\nx]        Demand profile in boundary sectors
LineNode(l,n)                               [n.a.]          Incidence matrix {-1 +1}
LineXNode(lx,nx)                            [n.a.]          Incidence matrix {-1 +1}
SectorXSpillageNode(slx,nx)                 [n.a.]          Incidence matrix {-1 +1}
SectorXMaximumSpillage(slx,h)               [MW]            Maximum allowed spillage
LoadShedding(n,h)                           [MW]            Load shedding capacity
Location(au,n)                              [n.a.]          Location {1 0}
LocationX(au,nx)                            [n.a.]          Location {1 0}
OutageFactor(au,h)                          [%]             Outage Factor (100% = full outage)
PartLoadMin(au)                             [%]             Minimum part load
PowerCapacity(au)                           [MW\u]          Installed capacity
PowerInitial(au)                            [MW\u]          Power output before initial period
PowerMinStable(au)                          [MW\u]          Minimum power output
PriceTransmission(l,h)                      [EUR\MWh]       Transmission price
StorageChargingCapacity(au)                 [MW\u]          Storage capacity
StorageChargingEfficiency(au)               [%]             Charging efficiency
StorageSelfDischarge(au)                    [%\day]         Self-discharge of the storage units
RampDownMaximum(au)                         [MW\h\u]        Ramp down limit
RampShutDownMaximum(au)                     [MW\h\u]        Shut-down ramp limit
RampStartUpMaximum(au)                      [MW\h\u]        Start-up ramp limit
RampStartUpMaximumH(au,h)                   [MW\h\u]        Start-up ramp limit - Clustered formulation
RampShutDownMaximumH(au,h)                  [MW\h\u]        Shut-down ramp limit - Clustered formulation
RampUpMaximum(au)                           [MW\h\u]        Ramp up limit
Reserve(au)                                 [n.a.]          Reserve technology {1 0}
StorageCapacity(au)                         [MWh\u]         Storage capacity
StorageDischargeEfficiency(au)              [%]             Discharge efficiency
StorageOutflow(au,h)                        [MW\u]          Storage outflows
StorageInflow(au,h)                         [MW\u]          Storage inflows (potential energy)
StorageInitial(au)                          [MWh]           Storage level before initial period
StorageProfile(au,h)                        [%]             Storage level to be respected at the end of each horizon
StorageMinimum(au)                          [MWh]           Storage minimum
Technology(au,t)                            [n.a.]          Technology type {1 0}
TimeDownMinimum(au)                         [h]             Minimum down time
TimeUpMinimum(au)                           [h]             Minimum up time
SectorXFlexDemandInput(nx,h)                [MWh]           Flexible demand inside BS at each timestep (unless for MTS)
SectorXFlexDemandInputInitial(nx)           [MWh]           Cumulative flexible demand inside the loop
SectorXFlexMaxCapacity(nx)                  [MW]            Max capacity for BS Flexible demand
SectorXFlexSupplyInput(nx,h)                [MWh]           Flexible demand inside BS at each timestep (unless for MTS)
SectorXFlexSupplyInputInitial(nx)           [MWh]           Cumulative flexible demand inside the loop
SectorXFlexMaxSupply(nx)                    [MW]            Max capacity for BS Flexible demand
$If %RetrieveStatus%==1 CommittedCalc(u,z)  [n.a.]          Committment status as for the MILP
Nunits(au)                                  [n.a.]          Number of units inside the cluster (upper bound value for integer variables)
K_QuickStart(n)                             [n.a.]          Part of the reserve that can be provided by offline quickstart units
QuickStartPower(au,h)                       [MW\h\u]        Available max capacity in tertiary regulation up from fast-starting power plants - TC formulation
StorageAlertLevel(au,h)                     [MWh]           Storage alert - Will only be violated to avoid power rationing
StorageFloodControl(au,h)                   [MWh]           Storage flood control
StorageHours(au)                            [h]             Storage hours
SectorXStorageCapacity(nx)                  [MWh]           Storage capacity of the boundary sector
SectorXStorageSelfDischarge(nx)             [%]             Boundary sector storage self discharge
SectorXStorageHours(nx)                     [h]             Boundary sector storage hours
SectorXStorageMinimum(nx)                   [MWh]           Boundary sector storage minimum
SectorXStoragePowerMax(nx)                   [MW]            Maximum power capacity of the boundary sector storage
$If %MTS% == 0 SectorXStorageInitial(nx)    [MWh]           Boundary sector storage initial state of charge
SectorXStorageProfile(nx,h)                 [%]             Boundary sector storage level respected at the end of each horizon
SectorXAlertLevel(nx,h)                     [MWh]           Storage alert of the boundary sector - Will only be violated to avoid power rationing
SectorXFloodControl(nx,h)                   [MWh]           Storage flood control of the boundary sector
PTDF(l_int,n)                               [p.u.]          Power Transfer Distribution Factor Matrix

* New
$If %MTS% == 0 InertiaConstant(au)          [s]             Inertia Constant
$If %MTS% == 0 InertiaLimit(h)              [s\h]           Inertia Limit
$If %MTS% == 0 Droop(au)                    [%]             Droop
$If %MTS% == 0 SystemGainLimit(h)           [GW\Hz]         System Gain Limit
$If %MTS% == 0 FFRGainLimit(h)              [GW\Hz]         FFR Gain Limit
$If %MTS% == 0 PrimaryReserveLimit(h)       [MWh]           Primary Reserve
$If %MTS% == 0 FFRLimit(h)                  [MWh]           Fast Frequency Reserve

;



*Parameters as used within the loop
PARAMETERS
CostLoadShedding(n,h)                       [EUR\MW]        Value of lost load
LoadMaximum(au,h)                           [%]             Maximum load given AF and OF
PowerMustRun(au,h)                          [MW\u]          Minimum power output
StorageFinalMin(au)                         [MWh]           Minimum storage level at the end of the optimization horizon
$If %MTS% == 0 SectorXStorageFinalMin(nx)   [MWh]           Minimum boundary sector storage level at the end of the optimization horizon
MaxFlexDemand(n)                            [MW]            Maximum value of the flexible demand parameter
MaxOverSupply(n,h)                          [MWh]           Maximum flexible demand accumultation
AccumulatedOverSupply_inital(n)             [MWh]           Initial value of the flexible demand accumulation
;

* Scalar variables necessary to the loop:
scalar FirstHour,LastHour,LastKeptHour,day,ndays,failed,srp,nsrp;
FirstHour = 1;
scalar TimeStep;

*MADRID

scalar SystemFrequency, RoCoFMax, DeltaFrequencyMax, MaxPrimaryAllowed, ConversionFactor;
SystemFrequency = 50;
RoCoFMax = 0.5;
DeltaFrequencyMax = 0.8;
MaxPrimaryAllowed = 0.15;
ConversionFactor = 1000;

*Threshold values for p2h partecipation to reserve market as spinning/non-spinning reserves (TO BE IMPLEMENTED IN CONFIGFILE)
srp = 1;
nsrp = 3;

*===============================================================================
*Data import
*===============================================================================

$gdxin %inputfilename%

$LOAD mk
$LOAD n
$LOAD nx
$LOAD l
$LOAD l_RoW
$LOAD l_int
$LOAD lx
$LOAD slx
$LOAD au
$LOAD u
$LOAD t
$LOAD tr
$LOAD f
$LOAD p
$LOAD s
$LOAD wat
$LOAD p2x
$LOAD x2p
$LOAD chp
*New
$LOAD cu
$LOAD ba
$LOAD tc
$LOAD xu
$LOAD h
$LOAD z
$LOAD AvailabilityFactor
$LOAD CHPPowerLossFactor
$LOAD CHPPowerToHeat
$LOAD CHPMaxHeat
$LOAD CHPType
$LOAD Config
$LOAD CostFixed
$LOAD CostOfSpillage
$LOAD CostXStorageAlert
$LOAD CostXFloodControl
$LOAD CostXSpillage
$LOAD CostXNotServed
$LOAD CostLoadShedding
$LOAD CostShutDown
$LOAD CostStartUp
$LOAD CostVariable
$LOAD CostStorageAlert
$LOAD CostFloodControl
$LOAD Curtailment
$LOAD CostCurtailment
$LOAD Demand
$LOAD SectorXDemand
$LOAD StorageDischargeEfficiency
$LOAD Efficiency
$LOAD X2PowerConversionMultiplier
$LOAD Power2XConversionMultiplier
$LOAD EmissionMaximum
$LOAD EmissionRate
$LOAD FlowMaximum
$LOAD FlowMinimum
$LOAD FlowXMaximum
$LOAD FlowXMinimum
$LOAD Fuel
$LOAD LineNode
$LOAD LineXNode
$LOAD SectorXSpillageNode
$LOAD SectorXMaximumSpillage
$LOAD LoadShedding
$LOAD Location
$LOAD LocationX
$LOAD Nunits
$LOAD OutageFactor
$LOAD PowerCapacity
$LOAD PowerInitial
$LOAD PartLoadMin
$LOAD PriceTransmission
$LOAD StorageChargingCapacity
$LOAD StorageChargingEfficiency
$LOAD StorageSelfDischarge
$LOAD RampDownMaximum
$LOAD RampShutDownMaximum
$LOAD RampStartUpMaximum
$LOAD RampUpMaximum
$LOAD Reserve
$LOAD StorageCapacity
$LOAD StorageInflow
$LOAD StorageInitial
$LOAD StorageProfile
$LOAD StorageMinimum
$LOAD StorageAlertLevel
$LOAD StorageFloodControl
$LOAD StorageHours
$LOAD StorageOutflow
$LOAD Technology
$LOAD TimeDownMinimum
$LOAD TimeUpMinimum
$LOAD CostRampUp
$LOAD CostRampDown
$LOAD SectorXFlexDemandInput
$LOAD SectorXFlexDemandInputInitial
$LOAD SectorXFlexMaxCapacity
$LOAD SectorXFlexSupplyInput
$LOAD SectorXFlexSupplyInputInitial
$LOAD SectorXFlexMaxSupply
$LOAD SectorXStorageCapacity
$LOAD SectorXStoragePowerMax
$LOAD SectorXStorageSelfDischarge
$LOAD SectorXStorageHours
$LOAD SectorXStorageMinimum
$LOAD SectorXAlertLevel
$LOAD SectorXFloodControl
$If %MTS% == 0 $LOAD SectorXStorageInitial
$LOAD SectorXStorageProfile
$If %RetrieveStatus% == 1 $LOAD CommittedCalc
$LOAD PTDF

* New
$If %MTS% == 0 $LOAD InertiaConstant
$If %MTS% == 0 $LOAD InertiaLimit
$If %MTS% == 0 $LOAD Droop
$If %MTS% == 0 $LOAD SystemGainLimit
$If %MTS% == 0 $LOAD FFRGainLimit
$If %MTS% == 0 $LOAD PrimaryReserveLimit
$If %MTS% == 0 $LOAD FFRLimit
;


*===============================================================================
*Definition of variables
*===============================================================================

VARIABLES
Committed(au,h)      [n.a.]  Unit committed at hour h {1 0} or integer
StartUp(au,h)        [n.a.]  Unit start up at hour h {1 0}  or integer
ShutDown(au,h)       [n.a.]  Unit shut down at hour h {1 0} or integer
;

$If %LPFormulation% == 1 POSITIVE VARIABLES Committed (au,h) ; Committed.UP(au,h) = 1 ;
$If not %LPFormulation% == 1 INTEGER VARIABLES Committed (au,h), StartUp(au,h), ShutDown(au,h) ; Committed.UP(au,h) = Nunits(au) ; StartUp.UP(au,h) = Nunits(au) ; ShutDown.UP(au,h) = Nunits(au) ;



POSITIVE VARIABLES
AccumulatedOverSupply(n,h)              [MWh]   Accumulated oversupply due to the flexible demand
CostStartUpH(au,h)                      [EUR]   Cost of starting up
CostShutDownH(au,h)                     [EUR]   cost of shutting down
CostRampUpH(au,h)                       [EUR]   Ramping cost
CostRampDownH(au,h)                     [EUR]   Ramping cost
CurtailedPower(n,h)                     [MW]    Curtailed power at node n
CurtailmentReserve_2U(n,h)              [MW]    Curtailed power used for reserves at node n
CurtailmentReserve_3U(n,h)              [MW]    Curtailed power used for reserves at node n
FlowX(lx,h)                             [MW]    Flow through boundary sector lines
Power(au,h)                             [MW]    Power output
PowerConsumption(p2x,h)                 [MW]    Power consumption by P2X units
PowerMaximum(u,h)                       [MW]    Power output
PowerMinimum(u,h)                       [MW]    Power output
ShedLoad(n,h)                           [MW]    Shed load
StorageInput(au,h)                      [MW]   Charging input for storage units
StorageLevel(au,h)                      [MWh]   Storage level of charge
SectorXStorageLevel(nx,h)               [MWh]   Storage level of charge of the boundary sector
SectorXSpillage(slx,h)                  [MW]    Spillage from boundary sector x to boundary sector y
SectorXWaterNotWithdrawn(nx,h)          [MWh]   Available water not utilized
LL_MaxPower(n,h)                        [MW]    Deficit in terms of maximum power
LL_RampUp(au,h)                         [MW]    Deficit in terms of ramping up for each plant
LL_RampDown(au,h)                       [MW]    Deficit in terms of ramping down
LL_MinPower(n,h)                        [MW]    Power exceeding the demand
LL_2U(n,h)                              [MW]    Deficit in reserve up
LL_3U(n,h)                              [MW]    Deficit in reserve up - non spinning
LL_2D(n,h)                              [MW]    Deficit in reserve down
LL_StorageAlert(au,h)                   [MWh]   Unsatisfied water level constraint for going below alert level at each hour
LL_FloodControl(au,h)                   [MWh]   Unsatisfied water level constraint for going above flood level at each hour
SectorXStorageAlertViolation(nx,h)      [MWh]   Boundary Sector Unsatisfied water level constraint for going below alert level at each hour
SectorXFloodControlViolation(nx,h)      [MWh]   Boundary Sector Unsatisfied water level constraint for going above flood control at each hour
LL_SectorXFlexDemand(nx)                [MWh]   Deficit in flex demand
LL_SectorXFlexSupply(nx)                [MWh]   Deficit in flex supply
spillage(au,h)                          [MWh]   spillage from water reservoirs
SystemCost(h)                           [EUR]   Hourly system cost
Reserve_2U(au,h)                        [MW]    Spinning reserve up
Reserve_2D(au,h)                        [MW]    Spinning reserve down
Reserve_3U(au,h)                        [MW]    Non spinning quick start reserve up
Heat(au,h)                              [MW]    Heat output by chp plant
WaterSlack(au)                          [MWh]   Unsatisfied water level constraint at end of optimization period
StorageSlack(au,h)                      [MWh]   Unsatisfied storage level constraint at end of simulation timestep
XNotServed(nx,h)                        [MW]    Boundary sector demand satisfied by other sources
StorageLevelViolation(au)               [MWh]   Unsatisfied water level constraint at end of optimization period
StorageLevelViolation_H(au,h)           [MWh]   Unsatisfied storage level constraint at end of simulation timestep
SectorXStorageLevelViolation(nx)        [MWh]   Unsatisfied boundary sector water level constraint at end of optimization period
SectorXStorageLevelViolation_H(nx,h)    [MWh]   Unsatisfied boundary sector storage level constraint at end of simulation timestep
SectorXFlexDemand(nx,h)                 [MW]    FLexible boundary sector demand at each time step of each nx node
SectorXFlexSupply(nx,h)                 [MW]    FLexible boundary sector supply at each time step of each nx node
$If %MTS% == 1 SectorXStorageInitial(nx)
$If %MTS% == 1 SectorXStorageFinalMin(nx)

*New
$If %MTS% == 0 SysInertia(h)                       [s]         System Inertia
$If %MTS% == 0 PrimaryReserve_Available(cu,h)         [MW]        Primary Reserve available in the system
$If %MTS% == 0 SystemGain(h)                       [GW\Hz]     System Gain
$If %MTS% == 0 FFR_Available(ba,h)         [MW]        Fast Frequency Reserve available in the system
$If %MTS% == 0 FFRGain(h)                       [GW\Hz]     FFR Gain
*MADRID
$If %MTS% == 0 PowerLoss(h)                  [MW]        Primary Reserve available in the system
*$If %MTS% == 0 MaxInstantGenerator(h)                  [MW]        Primary Reserve available in the system
*New
TotalDemand_2U(h)                       [MW]    Total Spinning reserve up

;

free variable
SystemCostD                         ![EUR]  Total system cost for one optimization period
DemandModulation(n,h)               [MW]    Difference between the flexible demand and the baseline
PowerX(nx,au,h)                     [MW]    Power output in boundary sector
SectorXStorageInput(nx,h)           [MW]    Boundary Sector Storage input - output
ResidualLoad(n,h)                   [MW]    Residual Load
ObjectiveFunction(h)
OptimalityGap(h)
OptimizationError
Error
InjectedPower(h,n)                  [p.u]   Power injected to each node of the system
Flow(l,h)                           [MW]    Flow through lines

*New
*LoadRamp(n,h)                           [MW]    Load ramp
;

*Binary Variables
*LoadRampAux(n, h);

*===============================================================================
*Assignment of initial values
*===============================================================================


*Initial commitment status
CommittedInitial(au)=0;
CommittedInitial(u)$(PowerInitial(u)>0)=1;

* Definition of the minimum stable load:
PowerMinStable(au) = PartLoadMin(au)*PowerCapacity(au);

LoadMaximum(au,h)= AvailabilityFactor(au,h)*(1-OutageFactor(au,h));

* parameters for clustered formulation (quickstart is defined as the capability to go to minimum power in 15 min)
QuickStartPower(au,h) = 0;
QuickStartPower(au,h)$(RampStartUpMaximum(au)>=PowerMinStable(au)*4) = PowerCapacity(au)*LoadMaximum(au,h);
RampStartUpMaximumH(au,h) = min(PowerCapacity(au)*LoadMaximum(au,h),max(RampStartUpMaximum(au),PowerMinStable(au),QuickStartPower(au,h)));
RampShutDownMaximumH(au,h) = min(PowerCapacity(au)*LoadMaximum(au,h),max(RampShutDownMaximum(au),PowerMinStable(au)));

PowerMustRun(u,h)=PowerMinStable(u)*LoadMaximum(u,h);
PowerMustRun(u,h)$(sum(tr,Technology(u,tr))>=1 and smin(n,Location(u,n)*(1-Curtailment(n)))=1) = PowerCapacity(u)*LoadMaximum(u,h);

* Part of the reserve that can be provided by offline quickstart units:
K_QuickStart(n) = Config("QuickStartShare","val");

* Flexible Demand
MaxFlexDemand(n) = smax(h,Demand("Flex",n,h));
MaxOverSupply(n,h) = Config("DemandFlexibility","val") * Demand("Flex",n,h);
AccumulatedOverSupply_inital(n) = 0;

* Time step
TimeStep = Config("SimulationTimeStep","val");

* Display RampStartUpMaximum, RampShutDownMaximum, CommittedInitial;

$offorder

*===============================================================================
*Declaration and definition of equations
*===============================================================================
EQUATIONS
EQ_Objective_function
EQ_CHP_extraction_Pmax
EQ_CHP_extraction
EQ_CHP_backpressure
*EQ_CHP_demand_satisfaction 
EQ_BS_Demand_balance
EQ_BS_Demand_balance2
EQ_CHP_max_heat
EQ_Commitment
EQ_MinUpTime
EQ_MinDownTime
EQ_RampUp_TC
EQ_RampDown_TC
EQ_CostStartUp
EQ_CostShutDown
EQ_CostRampUp
EQ_CostRampDown
EQ_Demand_balance_DA
EQ_Demand_balance_2U
EQ_Demand_balance_3U
EQ_Demand_balance_2D
EQ_P2X_Power_Balance
EQ_X2P_Power_Balance
EQ_Max_Power_Consumption
EQ_Power_Balance_of_P2X_units
EQ_Power_Balance_of_X2P_units
EQ_Max_Power_Consumption_of_BS_units
EQ_Power_must_run
EQ_Power_available
EQ_Reserve_2U_capability
EQ_Reserve_2D_capability
EQ_Reserve_3U_capability
EQ_2U_limit_chp
EQ_2D_limit_chp
EQ_3U_limit_chp
EQ_Storage_minimum
EQ_Storage_alert
EQ_Storage_flood_control
EQ_Storage_level
EQ_Storage_input
EQ_Storage_MaxDischarge
EQ_Storage_MaxCharge
EQ_Storage_balance
EQ_Storage_Cyclic
*EQ_H2_demand
EQ_Storage_boundaries
EQ_Boundary_Sector_Storage_MaxDischarge
EQ_Boundary_Sector_Storage_MaxCharge
EQ_Boundary_Sector_Storage_PowerMax
EQ_Boundary_Sector_Storage_PowerMin
EQ_Boundary_Sector_Storage_minimum
EQ_Boundary_Sector_Storage_alert
EQ_Boundary_Sector_Flood_Control
EQ_Boundary_Sector_Storage_level
EQ_Boundary_Sector_Storage_balance
EQ_Boundary_Sector_Storage_boundaries
EQ_Boundary_Sector_Storage_Cyclic
EQ_SystemCost
EQ_Emission_limits
EQ_Flow_limits_lower
EQ_Flow_limits_upper
EQ_BS_Flow_limits_lower
EQ_BS_Flow_limits_upper
EQ_Force_Commitment
EQ_Force_DeCommitment
EQ_LoadShedding
EQ_Flexible_Demand
EQ_Flexible_Demand_Max
EQ_Flexible_Demand_Modulation_Min
EQ_Flexible_Demand_Modulation_Max
EQ_No_Flexible_Demand
EQ_Tot_DemandPtL
EQ_Max_Capacity_PtL
EQ_PtL_Demand
EQ_Tot_Flex_Supply_BS
EQ_Max_Flex_Supply_BS
EQ_Tot_Flex_Demand_BS
EQ_Max_Flex_Capacity_BS
EQ_BS_Flex_Demand
EQ_Curtailed_Power
EQ_Residual_Load
EQ_BS_Spillage_limits_upper
$If %RetrieveStatus% == 1 EQ_CommittedCalc
EQ_DC_Power_Flow
EQ_Total_Injected_Power

*New
$If %MTS% == 0 EQ_SysInertia
$If %MTS% == 0 EQ_Inertia_limit
$If %MTS% == 0 EQ_SystemGain
$If %MTS% == 0 EQ_SystemGain_limit
$If %MTS% == 0 EQ_PrimaryReserve_Available
$If %MTS% == 0 EQ_PrimaryReserve_Capability
$If %MTS% == 0 EQ_PrimaryReserve_Boundary
$If %MTS% == 0 EQ_Demand_balance_PrimaryReserve
$If %MTS% == 0 EQ_FFRGain
$If %MTS% == 0 EQ_FFRGain_limit
$If %MTS% == 0 EQ_FFR_Available
$If %MTS% == 0 EQ_FFR_Capability
$If %MTS% == 0 EQ_FFR_Boundary
$If %MTS% == 0 EQ_Demand_balance_FFR

*MADRID
$If %MTS% == 0 EQ_PowerLoss
*$If %MTS% == 0 EQ_MaxInstantGenerator

*New
EQ_Tot_Demand_2U
*EQ_Load_Ramp
*EQ_Load_Ramp_Aux
;

$If %RetrieveStatus% == 0 $goto skipequation

EQ_CommittedCalc(au,z)..
         Committed(au,z)
         =E=
         CommittedCalc(au,z)
;

$label skipequation

*Objective function
$ifthen %LPFormulation% == 1
EQ_SystemCost(i)..
         SystemCost(i)
         =E=
         sum(au,CostFixed(au)*TimeStep*Committed(au,i))
         +sum(au,CostRampUpH(au,i) + CostRampDownH(au,i))
         +sum(u,CostVariable(u,i) * Power(u,i)*TimeStep)
         +sum(l,PriceTransmission(l,i)*Flow(l,i)*TimeStep)
         +sum(n,CostLoadShedding(n,i)*ShedLoad(n,i)*TimeStep)
         +sum(nx, CostXNotServed(nx,i) * XNotServed(nx,i)*TimeStep)
         +sum(chp, CostVariable(chp,i) * CHPPowerLossFactor(chp) * Heat(chp,i)*TimeStep)
         +Config("ValueOfLostLoad","val")*(sum(n,(LL_MaxPower(n,i)+LL_MinPower(n,i))*TimeStep))
         +0.8*Config("ValueOfLostLoad","val")*(sum(n,(LL_2U(n,i)+LL_2D(n,i)+LL_3U(n,i))*TimeStep))
         +0.7*Config("ValueOfLostLoad","val")*sum(au,(LL_RampUp(au,i)+LL_RampDown(au,i))*TimeStep)
         +0.7*Config("ValueOfLostLoad","val")*(sum(nx,(SectorXWaterNotWithdrawn(nx,i))*TimeStep))                                                                             
         +sum(s,CostStorageAlert(s,i)*LL_StorageAlert(s,i)*TimeStep)
         +sum(s,CostFloodControl(s,i)*LL_FloodControl(s,i)*TimeStep)
*         +Config("CostOfSpillage","val")*(sum(au,spillage(au,i))*TimeStep
         +sum(nx,CostXStorageAlert(nx,i)*SectorXStorageAlertViolation(nx,i)*TimeStep)
         +sum(nx,CostXFloodControl(nx,i)*SectorXFloodControlViolation(nx,i)*TimeStep)
         +sum(au,CostOfSpillage(au,i)*spillage(au,i)*TimeStep)
         +sum(slx,CostXSpillage(slx,i)*SectorXSpillage(slx,i)*TimeStep)
         +sum(slx,SectorXSpillage(slx,i)*TimeStep)
         +sum(n,CurtailedPower(n,i) * CostCurtailment(n,i) * TimeStep)
;
$else

EQ_SystemCost(i)..
         SystemCost(i)
         =E=
         sum(au,CostFixed(au)*TimeStep*Committed(au,i))
         +sum(au,CostStartUpH(au,i) + CostShutDownH(au,i))
         +sum(au,CostRampUpH(au,i) + CostRampDownH(au,i))
         +sum(u,CostVariable(u,i) * Power(u,i)*TimeStep)
         +sum(l,PriceTransmission(l,i)*Flow(l,i)*TimeStep)
         +sum(n,CostLoadShedding(n,i)*ShedLoad(n,i)*TimeStep)
         +sum(nx, CostXNotServed(nx,i) * XNotServed(nx,i)*TimeStep)
         +sum(chp, CostVariable(chp,i) * CHPPowerLossFactor(chp) * Heat(chp,i)*TimeStep)
         +Config("ValueOfLostLoad","val")*(sum(n,(LL_MaxPower(n,i)+LL_MinPower(n,i))*TimeStep))
         +0.8*Config("ValueOfLostLoad","val")*(sum(n,(LL_2U(n,i)+LL_2D(n,i)+LL_3U(n,i))*TimeStep))
         +0.7*Config("ValueOfLostLoad","val")*sum(au,(LL_RampUp(au,i)+LL_RampDown(au,i))*TimeStep)
         +15*(sum(nx,(SectorXWaterNotWithdrawn(nx,i))*TimeStep))                                                                                                              
         +sum(s,CostStorageAlert(s,i)*LL_StorageAlert(s,i)*TimeStep)
         +sum(s,CostFloodControl(s,i)*LL_FloodControl(s,i)*TimeStep)
*         +Config("CostOfSpillage","val")*(sum(au,spillage(au,i))*TimeStep
         +sum(nx,CostXStorageAlert(nx,i)*SectorXStorageAlertViolation(nx,i)*TimeStep)
         +sum(nx,CostXFloodControl(nx,i)*SectorXFloodControlViolation(nx,i)*TimeStep)
         +sum(au,CostOfSpillage(au,i)*spillage(au,i)*TimeStep)
         +sum(slx,CostXSpillage(slx,i)*SectorXSpillage(slx,i)*TimeStep)
         +sum(slx,SectorXSpillage(slx,i)*TimeStep)
         +sum(n,CurtailedPower(n,i) * CostCurtailment(n,i) * TimeStep)
;

$endIf
;

EQ_Objective_function..
         SystemCostD
         =E=
         sum(i,SystemCost(i))
         +Config("WaterValue","val")*sum(au,WaterSlack(au))
         +Config("WaterValue","val")*sum(au,StorageLevelViolation(au))
         +Config("WaterValue","val")*sum(nx,SectorXStorageLevelViolation(nx))
         +Config("ValueOfLostLoad","val")*sum(nx,LL_SectorXFlexDemand(nx) + LL_SectorXFlexSupply(nx))
;

* 3 binary commitment status
EQ_Commitment(au,i)..
         Committed(au,i)-CommittedInitial(au)$(ord(i) = 1)-Committed(au,i-1)$(ord(i) > 1)
         =E=
         StartUp(au,i) - ShutDown(au,i)
;

* minimum up time
EQ_MinUpTime(au,i)$(TimeStep <= TimeUpMinimum(au))..
         sum(ii$( (ord(ii) >= ord(i) - ceil(TimeUpMinimum(au)/TimeStep)) and (ord(ii) <= ord(i)) ), StartUp(au,ii))
         + sum(h$( (ord(h) >= FirstHour + ord(i) - ceil(TimeUpMinimum(au)/TimeStep) -1) and (ord(h) < FirstHour)),StartUp.L(au,h))
         =L=
         Committed(au,i)
;

* minimum down time
EQ_MinDownTime(au,i)$(TimeStep <= TimeDownMinimum(au))..
         sum(ii$( (ord(ii) >= ord(i) - ceil(TimeDownMinimum(au)/TimeStep)) and (ord(ii) <= ord(i)) ), ShutDown(au,ii))
         + sum(h$( (ord(h) >= FirstHour + ord(i) - ceil(TimeDownMinimum(au)/TimeStep) -1) and (ord(h) < FirstHour)),ShutDown.L(au,h))
         =L=
         Nunits(au)-Committed(au,i)
;

* ramp up constraints
EQ_RampUp_TC(au,i)$(sum(tr,Technology(au,tr))=0)..
         Power(au,i) - Power(au,i-1)$(ord(i) > 1) - PowerInitial(au)$(ord(i) = 1)
         =L=
         (Committed(au,i) - StartUp(au,i)) * RampUpMaximum(au) * TimeStep
         + RampStartUpMaximumH(au,i) * TimeStep * StartUp(au,i)
         - PowerMustRun(au,i) * ShutDown(au,i)
         + LL_RampUp(au,i)
;

* ramp down constraints
EQ_RampDown_TC(au,i)$(sum(tr,Technology(au,tr))=0)..
         Power(au,i-1)$(ord(i) > 1) + PowerInitial(au)$(ord(i) = 1) - Power(au,i)
         =L=
         (Committed(au,i) - StartUp(au,i)) * RampDownMaximum(au) * TimeStep
         - PowerMustRun(au,i) * StartUp(au,i)
         + RampShutDownMaximumH(au,i) * TimeStep * ShutDown(au,i)
         + LL_RampDown(au,i)
;

* Start up cost
EQ_CostStartUp(au,i)$(sum(tr,Technology(au,tr))=0)..
         CostStartUpH(au,i)
         =E=
         CostStartUp(au)*StartUp(au,i)
;

* Shut down cost
EQ_CostShutDown(au,i)$(sum(tr,Technology(au,tr))=0)..
         CostShutDownH(au,i)
         =E=
         CostShutDown(au)*ShutDown(au,i)
;

EQ_CostRampUp(au,i)$(CostRampUp(au) <> 0)..
         CostRampUpH(au,i)
         =G=
         CostRampUp(au)*(Power(au,i)-PowerInitial(au)$(ord(i) = 1)-Power(au,i-1)$(ord(i) > 1))
;

EQ_CostRampDown(au,i)$(CostRampDown(au) <> 0)..
         CostRampDownH(au,i)
         =G=
         CostRampDown(au)*(PowerInitial(au)$(ord(i) = 1)+Power(au,i-1)$(ord(i) > 1)-Power(au,i))
;

EQ_Residual_Load(n,i)..
        ResidualLoad(n,i)
        =E=
        Demand("DA",n,i)
        + Demand("Flex",n,i)
        - DemandModulation(n,i)
        + InjectedPower(i,n)
        + sum(p2x, PowerConsumption(p2x,i) * Location(p2x,n))
        - sum(u,Power(u,i)$(sum(tr,Technology(u,tr))>=1) * Location(u,n))
        - sum(l_RoW,Flow(l_RoW,i)*LineNode(l_RoW,n))
;

*Hourly demand balance in the day-ahead market for each node

EQ_Demand_balance_DA(n,i)..
         sum(u,Power(u,i)*Location(u,n))
         +sum(x2p,Power(x2p,i)*Location(x2p,n))
         +InjectedPower(i,n)
         +sum(l_RoW,Flow(l_RoW,i)*LineNode(l_RoW,n))     
         +ShedLoad(n,i)
         +LL_MaxPower(n,i)
         =E=
         Demand("DA",n,i) + Demand("Flex",n,i)
         +DemandModulation(n,i)                                                                                                                                             
         +sum(s,StorageInput(s,i)*Location(s,n))
         +sum(p2x,PowerConsumption(p2x,i)*Location(p2x,n))
         +LL_MinPower(n,i)
;

* Energy balance at the level of the flexible demand, considered as a storage capacity
EQ_Flexible_Demand(n,i)..
         DemandModulation(n,i) * TimeStep
         =e=
         AccumulatedOverSupply(n,i)
         - AccumulatedOverSupply_inital(n)$(ord(i) = 1)
         - AccumulatedOverSupply(n,i-1)$(ord(i) > 1)
;

* The accumulated oversupply is limited by size of flexible demand (i.e. storage) capacity
EQ_Flexible_Demand_max(n,i)..
         AccumulatedOverSupply(n,i)
         =l=
         MaxOverSupply(n,i)
;

* The maximum downards demand modulation  at each time step is limited by the value of the flexible demand
EQ_Flexible_Demand_Modulation_Min(n,i)..
         DemandModulation(n,i)
         =g=
         -Demand("Flex",n,i)
;

* The upwards demand modulation is arbitrarily limited by the maximum recorded value of the flexible demand
EQ_Flexible_Demand_Modulation_max(n,i)..
         DemandModulation(n,i)
         =l=
*         MaxFlexDemand(n) - Demand("Flex",n,i)                                                                                                                             
         Demand("Flex",n,i)
;

* If flexible demand is not considered, set the accumulated over supply to zero:
EQ_No_Flexible_Demand(n,i)..
         AccumulatedOverSupply(n,i)
         =e=
         0
;

*Hourly demand balance in the upwards spinning reserve market for each node
EQ_Demand_balance_2U(n,i)..
         sum((u),Reserve_2U(u,i)*Reserve(u)*Location(u,n))
         + sum((chp),Reserve_2U(chp,i)*Reserve(chp)*Location(chp,n))
         + CurtailmentReserve_2U(n,i) + LL_2U(n,i)
         =E=
*New
         Demand("2U",n,i)
$If %ActivateAdvancedReserves% == 2 +(Demand("2U",n,i) + max(smax((u,tc),PowerCapacity(u)/Nunits(u)*Technology(u,tc)*LoadMaximum(u,i)*Location(u,n)), smax(l,FlowMaximum(l,i)*LineNode(l,n))$(card(l)>0)))*(1-K_QuickStart(n))
$If %ActivateAdvancedReserves% == 1 +(Demand("2U",n,i) + smax((u,tc),PowerCapacity(u)/Nunits(u)*Technology(u,tc)*LoadMaximum(u,i)*Location(u,n)))*(1-K_QuickStart(n))
*$If %ActivateAdvancedReserves% == 0 +(Demand("2U",n,i))*(1-K_QuickStart(n))
;

*New
*Hourly total demand in the upwards spinning reserve market for the system
EQ_Tot_Demand_2U(i)..
         sum((n),Demand("2U",n,i))
         =E=
         TotalDemand_2U(i)
;



*New NOT WORKING
* Big-M constraints to linearize the division
*EQ_Load_Ramp_Aux(n, i)$(ord(i) > 1)..
*    ResidualLoad(n, i-1) - ResidualLoad(n, i) =L= 0.0000001 * (1 - LoadRampAux(n, i));
*    ResidualLoad(n, i) - ResidualLoad(n, i-1) =L= 5 * (1 - LoadRampAux(n, i));

* Define the division using the binary variable
*EQ_Load_Ramp(n, i)..
*    LoadRamp(n, i) =E= (ResidualLoad(n, i) - ResidualLoad(n, i-1))/ResidualLoad(n, i);


*Hourly demand balance in the upwards non-spinning reserve market for each node
EQ_Demand_balance_3U(n,i)..
         sum((u),(Reserve_2U(u,i) + Reserve_3U(u,i))*Reserve(u)*Location(u,n))
         + sum((chp),(Reserve_2U(chp,i) + Reserve_3U(chp,i))*Reserve(chp)*Location(chp,n))
         + CurtailmentReserve_2U(n,i) + CurtailmentReserve_3U(n,i) + LL_3U(n,i)
         =E=
         Demand("2U",n,i)
$If %ActivateAdvancedReserves% == 2 + max(smax((u,tc),PowerCapacity(u)/Nunits(u)*Technology(u,tc)*LoadMaximum(u,i)*Location(u,n)), smax(l,FlowMaximum(l,i)*LineNode(l,n))$(card(l)>0))
$If %ActivateAdvancedReserves% == 1 + smax((u,tc),PowerCapacity(u)/Nunits(u)*Technology(u,tc)*LoadMaximum(u,i)*Location(u,n))
;

*Hourly demand balance in the downwards reserve market for each node
EQ_Demand_balance_2D(n,i)..
         sum((u),Reserve_2D(u,i)*Reserve(u)*Location(u,n))
         + sum((chp),Reserve_2D(chp,i)*Reserve(chp)*Location(chp,n))
         + LL_2D(n,i)
         =E=
         Demand("2D",n,i)
$If %ActivateAdvancedReserves% == 2 + max(smax(s,StorageChargingCapacity(s)/Nunits(s)*Location(s,n)), smax(l,-FlowMaximum(l,i)*LineNode(l,n))$(card(l)>0))
$If %ActivateAdvancedReserves% == 1 + smax(s,StorageChargingCapacity(s)/Nunits(s)*Location(s,n))
;

*Curtailed power
EQ_Curtailed_Power(n,i)..
        CurtailedPower(n,i) + CurtailmentReserve_2U(n,i) + CurtailmentReserve_3U(n,i)
        =E=
        sum(u,(Nunits(u)*PowerCapacity(u)*LoadMaximum(u,i)-Power(u,i))$(sum(tr,Technology(u,tr))>=1) * Location(u,n))
;

* Reserve capability
EQ_Reserve_2U_capability(u,i)..
         Reserve_2U(u,i)
         =L=
         PowerCapacity(u)*LoadMaximum(u,i)*Committed(u,i)
         - Power(u,i)
;

EQ_Reserve_2D_capability(u,i)..
         Reserve_2D(u,i)
         =L=
         (Power(u,i) - PowerMustRun(u,i) * Committed(u,i))
         + (StorageChargingCapacity(u)*Nunits(u)-StorageInput(u,i))$(s(u))
;

EQ_Reserve_3U_capability(u,i)$(QuickStartPower(u,i) > 0)..
         Reserve_3U(u,i)
         =L=
         (Nunits(u)-Committed(u,i))*QuickStartPower(u,i)*TimeStep
;

EQ_2U_limit_chp(chp,i)..
         Reserve_2U(chp,i)
         =l=
         (StorageCapacity(chp)*Nunits(chp)-StorageLevel(chp,i))/(0.25/CHPPowerToHeat(chp))

;

EQ_2D_limit_chp(chp,i)..
         Reserve_2D(chp,i)
         =l=
         StorageLevel(chp,i)/(0.25/CHPPowerToHeat(chp))

;

EQ_3U_limit_chp(chp,i)..
         Reserve_3U(chp,i)
         =l=
         (StorageCapacity(chp)*Nunits(chp)-StorageLevel(chp,i))*(CHPPowerToHeat(chp))

;






*MADRID

$ifthen %MTS% == 0
EQ_PowerLoss(u,i)$(ord(u))..
   PowerLoss(i) =g= Power(u, i);

*EQ_MaxInstantGenerator(u, i)$(MaxInstantPower(i) - 1 < Power(u, i))..
*   MaxInstantGenerator(i) =e= ord(u);





$endIf




*New
*MARCO PRIMARY RESERVE
$ifthen %MTS% == 0
EQ_SystemGain(i)..
         SystemGain(i)
         =l=
         0
$If %FC% == 1 +sum(cu,(PowerCapacity(cu)*Committed(cu,i))/(Droop(cu)*SystemFrequency))
;

EQ_SystemGain_limit(i)..
         0
$If %FC% == 1 +SystemGainLimit(i)
         =l=
         SystemGain(i)
;

EQ_PrimaryReserve_Available(i)..
         sum(cu,PrimaryReserve_Available(cu,i))
         =l=
         SystemGain(i)*DeltaFrequencyMax
;

*EQ_PrimaryReserve_Capability(i)..
*         PrimaryReserve_Available(i)
*         =L=
*         sum(cu,(PowerCapacity(cu)*LoadMaximum(cu,i)*Committed(cu,i)*Reserve(cu))-Power(cu,i))
*;

EQ_PrimaryReserve_Capability(cu,i)..
         PrimaryReserve_Available(cu,i)
         =L=
         (PowerCapacity(cu)*LoadMaximum(cu,i)*Committed(cu,i)-Power(cu,i))*Reserve(cu)
;

EQ_PrimaryReserve_Boundary(cu,i)..
         PrimaryReserve_Available(cu,i)
         =L=
         PowerCapacity(cu)*Committed(cu,i)*Reserve(cu)*MaxPrimaryAllowed
;

*Hourly demand balance in the Primary Reserve market
EQ_Demand_balance_PrimaryReserve(i)..
         0
$If %FC% == 1 +PrimaryReserveLimit(i)
         =l=
         sum(cu,PrimaryReserve_Available(cu,i))
;
$endIf

*New
*MARCO FFR
$ifthen %MTS% == 0
EQ_FFRGain(i)..
         FFRGain(i)
         =l=
         0
$If %FC% == 1 +sum(ba,(PowerCapacity(ba)*Committed(ba,i))/(Droop(ba)*SystemFrequency))
;

EQ_FFRGain_limit(i)..
         0
$If %FC% == 1 +FFRGainLimit(i)
         =l=
         FFRGain(i)
;

EQ_FFR_Available(i)..
         sum(ba,FFR_Available(ba,i))
         =l=
         FFRGain(i)*DeltaFrequencyMax
;

*EQ_PrimaryReserve_Capability(i)..
*         PrimaryReserve_Available(i)
*         =L=
*         sum(cu,(PowerCapacity(cu)*LoadMaximum(cu,i)*Committed(cu,i)*Reserve(cu))-Power(cu,i))
*;

EQ_FFR_Capability(ba,i)..
         FFR_Available(ba,i)
         =L=
         (PowerCapacity(ba)*LoadMaximum(ba,i)*Committed(ba,i)*Reserve(ba))-Power(ba,i)
;

EQ_FFR_Boundary(ba,i)..
         FFR_Available(ba,i)
         =L=
         PowerCapacity(ba)
;

*Hourly demand balance in the FFR market
EQ_Demand_balance_FFR(i)..
         0
$If %FC% == 1 +FFRLimit(i)
         =l=
         sum(ba,FFR_Available(ba,i))
;
$endIf


* New
*MARCO INERTIA RESERVES
$ifthen %MTS% == 0
EQ_SysInertia(i)..
         SysInertia(i)
         =E=
         sum(cu,PowerCapacity(cu)*Committed(cu,i)*InertiaConstant(cu))/ConversionFactor
;

EQ_Inertia_limit(u,i)..
         0
$If %FC% == 1 +InertiaLimit(i)
         =l=
         SysInertia(i)
;
$endIf

*Minimum power output is above the must-run output level for each unit in all periods
EQ_Power_must_run(u,i)..
         PowerMustRun(u,i) * Committed(u,i)
         - (StorageInput(u,i) * CHPPowerLossFactor(u) )$(chp(u) and (CHPType(u,'Extraction') or CHPType(u,'P2H')))
         - (Heat(u,i) * CHPPowerLossFactor(u) )$(chp(u) and (CHPType(u,'Extraction') or CHPType(u,'P2H')))                                                      
         =L=
         Power(u,i)
;

*Maximum power output is below the available capacity                                                                                                          
EQ_Power_available(au,i)..
         Power(au,i)$(u(au))
         + Power(au,i)$(x2p(au))
         =L=
         PowerCapacity(au)$(u(au))*LoadMaximum(au,i)$(u(au))*Committed(au,i)$(u(au))
         + ((PowerCapacity(au))*LoadMaximum(au,i)*Nunits(au))$(x2p(au))
         + 0
;

*Boundary Sector Storage level must be above a minimum
EQ_Boundary_Sector_Storage_minimum(nx,i)..
         SectorXStorageMinimum(nx)
         =L=
         SectorXStorageLevel(nx,i)
;
*Boundary Sector Storage level should be above alert level, going below will only be violated to avoid power rationing (110% of most expensive power plant)
EQ_Boundary_Sector_Storage_alert(nx,i)..
         SectorXStorageCapacity(nx)*SectorXAlertLevel(nx,i)
         =L=
         SectorXStorageLevel(nx,i)
         + SectorXStorageAlertViolation(nx,i)
;
*Boundary Sector Storage level should be below flood control level
EQ_Boundary_Sector_Flood_Control(nx,i)$(SectorXFloodControl(nx,i) > SectorXAlertLevel(nx,i))..
         SectorXStorageCapacity(nx)*SectorXFloodControl(nx,i)
         + SectorXFloodControlViolation(nx,i)
         =G=
         SectorXStorageLevel(nx,i)
;
*Storage level must be below storage capacity
EQ_Boundary_Sector_Storage_level(nx,i)..
         SectorXStorageLevel(nx,i)
         =L=
         SectorXStorageCapacity(nx)
;

*Discharge is limited by the storage level
EQ_Boundary_Sector_Storage_MaxDischarge(nx,i)..
         - SectorXStorageInitial(nx)$(ord(i) = 1)
         - SectorXStorageLevel(nx,i-1)$(ord(i) > 1)
         =L=
         SectorXStorageInput(nx,i)*TimeStep
;

EQ_Boundary_Sector_Storage_MaxCharge(nx,i)$(SectorXStorageCapacity(nx)>0)..
        SectorXStorageInput(nx,i)*TimeStep
        =L=
        (SectorXStorageCapacity(nx) - SectorXStorageInitial(nx))$(ord(i) = 1)
        + (SectorXStorageCapacity(nx) - SectorXStorageLevel(nx,i-1))$(ord(i) > 1)
;

* Storage input is bounded by the maximum power
EQ_Boundary_Sector_Storage_PowerMax(nx,i)..
        SectorXStorageInput(nx,i)
        =L=
        SectorXStoragePowerMax(nx)
;

* Storage output (negative input) is bounded by the maximum power
EQ_Boundary_Sector_Storage_PowerMin(nx,i)..
        SectorXStorageInput(nx,i)
        =G=
        -SectorXStoragePowerMax(nx)
;

*Storage balance
EQ_Boundary_Sector_Storage_balance(nx,i)..
         SectorXStorageInitial(nx)$(ord(i) = 1)
         + SectorXStorageLevel(nx,i-1)$(ord(i) > 1)
         + SectorXStorageInput(nx,i)*TimeStep
         =E=
         SectorXStorageLevel(nx,i)
         + SectorXStorageSelfDischarge(nx)*SectorXStorageLevel(nx,i)
;

* Minimum level at the end of the optimization horizon:
EQ_Boundary_Sector_Storage_boundaries(nx,i)$(ord(i) = card(i))..
         SectorXStorageFinalMin(nx)
         =L=
         SectorXStorageLevel(nx,i) + SectorXStorageLevelViolation(nx)
;

EQ_Boundary_Sector_Storage_Cyclic(nx)$(SectorXStorageHours(nx)>=8)..
         SectorXStorageFinalMin(nx)
         =E=
         SectorXStorageInitial(nx)
;

*Storage level must be above a minimum
EQ_Storage_minimum(au,i)..
         StorageMinimum(au)$(s(au))*Nunits(au)$(s(au))
         =L=
         StorageLevel(au,i)$(s(au))
;

*Storage level should be above alert level, going below will only be violated to avoid power rationing (110% of most expensive power plant)
EQ_Storage_alert(s,i)..
         StorageCapacity(s)*Nunits(s)*min(StorageAlertLevel(s,i),AvailabilityFactor(s,i))
         =L=
         StorageLevel(s,i)
         + LL_StorageAlert(s,i)
;
EQ_Storage_flood_control(s,i)$(StorageFloodControl(s,i) > StorageAlertLevel(s,i))..
         StorageCapacity(s)*Nunits(s)*StorageFloodControl(s,i)
         + LL_FloodControl(s,i)
         =G=
         StorageLevel(s,i)
;

*Storage level must be below storage capacity
EQ_Storage_level(au,i)..
         StorageLevel(au,i)$(s(au))
         =L=
         StorageCapacity(au)$(s(au))*AvailabilityFactor(au,i)$(s(au))*Nunits(au)$(s(au))
;

* Storage charging is bounded by the maximum capacity
EQ_Storage_input(s,i)..
         StorageInput(s,i)
         =L=
         StorageChargingCapacity(s)*(Nunits(s)-Committed(s,i))
;

*Discharge is limited by the storage level
EQ_Storage_MaxDischarge(au,i)$(StorageCapacity(au)$(s(au))>PowerCapacity(au)$(s(au))*TimeStep)..
         Power(au,i)$(s(au))*TimeStep/(max(StorageDischargeEfficiency(au)$(s(au)),0.0001))
         =L=
         StorageInitial(au)$(s(au))$(ord(i) = 1)
         + StorageLevel(au,i-1)$(s(au))$(ord(i) > 1)
         + StorageInflow(au,i)$(s(au))*Nunits(au)$(s(au))*TimeStep
;

*Charging is limited by the remaining storage capacity
EQ_Storage_MaxCharge(au,i)$(StorageCapacity(au)$(s(au))>PowerCapacity(au)$(s(au))*TimeStep)..
         StorageInput(au,i)$(s(au))*StorageChargingEfficiency(au)$(s(au))*TimeStep
         =L=
         (Nunits(au)$(s(au)) * StorageCapacity(au)$(s(au))-StorageInitial(au)$(s(au)))$(ord(i) = 1)
         + (Nunits(au)$(s(au)) * StorageCapacity(au)$(s(au))*AvailabilityFactor(au,i-1)$(s(au)) - StorageLevel(au,i-1))$(ord(i) > 1)
         + StorageOutflow(au,i)$(s(au))*Nunits(au)$(s(au))*TimeStep
;

*Storage balance
EQ_Storage_balance(au,i)..
         StorageInitial(au)$(s(au))$(ord(i) = 1)
         +StorageLevel(au,i-1)$(s(au))$(ord(i) > 1)
         +StorageInflow(au,i)$(s(au))*Nunits(au)$(s(au))*TimeStep
         +StorageInput(au,i)$(s(au))*StorageChargingEfficiency(au)$(s(au))*TimeStep
         =E=
         StorageLevel(au,i)$(s(au))
         +StorageOutflow(au,i)$(s(au))*Nunits(au)$(s(au))*TimeStep
         +spillage(au,i)$(s(au) or wat(au))
         +Power(au,i)$(s(au))*TimeStep/(max(StorageDischargeEfficiency(au)$(s(au)),0.0001))
         +StorageSelfDischarge(au)$(s(au))*StorageLevel(au,i)$(s(au))*TimeStep
;

* Minimum level at the end of the optimization horizon:
EQ_Storage_boundaries(au,i)$(ord(i) = card(i))..
         StorageFinalMin(au)$(s(au))
         =L=
         StorageLevel(au,i)$(s(au))
         + StorageLevelViolation(au)$(s(au))
         + WaterSlack(au)$(s(au))
;

EQ_Storage_Cyclic(au)$(StorageHours(au)>=8)..
         StorageFinalMin(au)
         =E=
         StorageInitial(au)
;

*Total emissions are capped
EQ_Emission_limits(i,p,n)..
         sum(u,Power(u,i)*EmissionRate(u,p)*TimeStep*Location(u,n))
         =L=
         EmissionMaximum(n,p)*TimeStep
;

*DC Power Flow
EQ_DC_Power_Flow(l_int,i)$(sum(n,abs(PTDF(l_int,n)<>0)))..
         Flow(l_int,i)
         =E=
         sum(n,-PTDF(l_int,n)*InjectedPower(i,n))
;

EQ_Total_Injected_Power(i,n)..
        InjectedPower(i,n)
         =E=
         sum(l_int,Flow(l_int,i)*LineNode(l_int,n))
;

*Flows are above minimum values
EQ_Flow_limits_lower(l,i)..
         FlowMinimum(l,i)
         =L=
         Flow(l,i)
;

*Flows are below maximum values
EQ_Flow_limits_upper(l,i)..
         Flow(l,i)
         =L=
         FlowMaximum(l,i)
;

*Boundary Sector Flows are above minimum values
EQ_BS_Flow_limits_lower(lx,i)..
         FlowXMinimum(lx,i)
         =L=
         FlowX(lx,i)
;

*Boundary Sector Flows are below maximum values
EQ_BS_Flow_limits_upper(lx,i)..
         FlowX(lx,i)
         =L=
         FlowXMaximum(lx,i)
;


*Boundary Sector Spillages are below maximum values
EQ_BS_Spillage_limits_upper(slx,i)..
         SectorXSpillage(slx,i)
         =L=
         SectorXMaximumSpillage(slx,i)
;

*Force Unit commitment/decommitment:
* E.g: renewable units with AF>0 must be committed
EQ_Force_Commitment(u,i)$((sum(tr,Technology(u,tr))>=1 and LoadMaximum(u,i)>0))..
         Committed(u,i)
         =G=
         1
;

* E.g: renewable units with AF=0 must be decommitted
EQ_Force_DeCommitment(u,i)$(LoadMaximum(u,i)=0)..
         Committed(u,i)
         =E=
         0
;

*Load shedding
EQ_LoadShedding(n,i)..
         ShedLoad(n,i)
         =L=
         LoadShedding(n,i)
;

* CHP units:
EQ_CHP_extraction(chp,i)$(CHPType(chp,'Extraction'))..
         Power(chp,i)
         =G=
         StorageInput(chp,i) * CHPPowerToHeat(chp)
         + Heat(chp,i) * CHPPowerToHeat(chp)

;

EQ_CHP_extraction_Pmax(chp,i)$(CHPType(chp,'Extraction') or CHPType(chp,'P2H'))..
         Power(chp,i)
         =L=
         PowerCapacity(chp)*Nunits(chp)
         - StorageInput(chp,i) * CHPPowerLossFactor(chp)
         - Heat(chp,i) * CHPPowerLossFactor(chp)
;

EQ_CHP_backpressure(chp,i)$(CHPType(chp,'Back-Pressure'))..
         Power(chp,i)
         =E=
         StorageInput(chp,i) * CHPPowerToHeat(chp)
         + Heat(chp,i) * CHPPowerToHeat(chp)
;

EQ_CHP_max_heat(chp,i)..
         StorageInput(chp,i)
         + Heat(chp,i)
         =L=
         CHPMaxHeat(chp)*Nunits(chp)
;

* Power to boundary sector units                                                                                                                                        
EQ_Power_Balance_of_P2X_units(nx,p2x,i)..                                                                                                                                                                                                                                                    
         PowerX(nx,p2x,i)
         =E=
         PowerConsumption(p2x,i) * Power2XConversionMultiplier(nx,p2x,i) * LocationX(p2x,nx)
;

EQ_Power_Balance_of_X2P_units(nx,x2p,i)..                                                                                                                                                                                                                                                    
         PowerX(nx,x2p,i) 
         =E=
         0 + 
         (Power(x2p,i)/X2PowerConversionMultiplier(nx,x2p,i)  * LocationX(x2p,nx))$(X2PowerConversionMultiplier(nx,x2p,i)<>0)
;

EQ_Max_Power_Consumption_of_BS_units(p2x,i)..
         PowerConsumption(p2x,i)
         =L=
         PowerCapacity(p2x) * Nunits(p2x)
;

EQ_BS_Demand_balance(nx,i)..
        sum(p2x, PowerX(nx,p2x,i))
        + sum(x2p, PowerX(nx,x2p,i))
        + sum(chp, Heat(chp,i)*LocationX(chp,nx))        
        + XNotServed(nx,i)
        + sum(lx,FlowX(lx,i)*LineXNode(lx,nx))
        + SectorXFlexSupply(nx,i)
        + sum(slx,SectorXSpillage(slx,i)*SectorXSpillageNode(slx,nx))
        =E=
        SectorXDemand(nx,i)
        + SectorXFlexDemand(nx,i)
        + SectorXStorageInput(nx,i)
        + SectorXWaterNotWithdrawn(nx,i)
;

EQ_BS_Demand_balance2(nx,i)..
        SectorXStorageInput(nx,i)
        + SectorXDemand(nx,i)$(SectorXDemand(nx,i)<0)
        =L=
        sum(p2x, PowerX(nx,p2x,i))
        + sum(x2p, PowerX(nx,x2p,i))
        + sum(lx,FlowX(lx,i)*LineXNode(lx,nx))
        + SectorXFlexSupply(nx,i)
        + sum(slx,SectorXSpillage(slx,i)*SectorXSpillageNode(slx,nx))
        + sum(chp, Heat(chp,i)*LocationX(chp,nx))
;

* Equations concerning boundary sector flexible demand
* Only valid in old formulation!!! If MTS=1: assures that total demand of boundary sector (PtL) is fulfilled
EQ_Tot_Flex_Demand_BS(nx)..
         SectorXFlexDemandInputInitial(nx)
         + sum(i,SectorXFlexDemandInput(nx,i))
         + LL_SectorXFlexDemand(nx)
         =E=
         sum(i,SectorXFlexDemand(nx,i))
;

* Capacity of boundary sector (PtL) must not be exceeded
EQ_Max_Flex_Capacity_BS(nx,i)..
         SectorXFlexDemand(nx,i)
         =L=
         SectorXFlexMaxCapacity(nx)
;

* Only valid in old formulation!!! If MTS = 0: Boundary sector (PtL) demand is not a variable anymore, but a parameter
EQ_BS_Flex_Demand(nx,i)..
         SectorXFlexDemand(nx,i)
         =E=
         SectorXFlexDemandInput(nx,i)
;


EQ_Tot_Flex_Supply_BS(nx)..
         SectorXFlexSupplyInputInitial(nx)
         + sum(i,SectorXFlexSupplyInput(nx,i))
         + LL_SectorXFlexSupply(nx)
         =G=
         sum(i,SectorXFlexSupply(nx,i))
;

* Capacity of boundary sector (PtL) must not be exceeded
EQ_Max_Flex_Supply_BS(nx,i)..
         SectorXFlexSupply(nx,i)
         =L=
         SectorXFlexMaxSupply(nx)
;

* Add missing equation definitions
EQ_P2X_Power_Balance(p2x,i)..
         Power(p2x,i)
         =E=
         PowerConsumption(p2x,i)
;

EQ_Max_Power_Consumption(p2x,i)..
         PowerConsumption(p2x,i)
         =L=
         PowerCapacity(p2x) * Nunits(p2x)
;

*===============================================================================
*Definition of models
*===============================================================================
MODEL UCM_SIMPLE /
EQ_Objective_function,
EQ_Curtailed_Power,
EQ_CHP_extraction_Pmax,
EQ_CHP_extraction,
EQ_CHP_backpressure,
*EQ_CHP_demand_satisfaction,
EQ_BS_Demand_balance,
*EQ_BS_Demand_balance2,
EQ_CHP_max_heat,
EQ_CostRampUp,
EQ_CostRampDown,
$If not %LPFormulation% == 1 EQ_CostStartUp,
$If not %LPFormulation% == 1 EQ_CostShutDown,
EQ_Commitment,
$If not %LPFormulation% == 1 EQ_MinUpTime,
$If not %LPFormulation% == 1 EQ_MinDownTime,
EQ_RampUp_TC,
EQ_RampDown_TC,
EQ_Demand_balance_DA,
EQ_Demand_balance_2U,
EQ_Demand_balance_2D,
EQ_Demand_balance_3U,
$If not %LPFormulation% == 1 EQ_Power_must_run,
EQ_P2X_Power_Balance,
EQ_Max_Power_Consumption,
EQ_Power_Balance_of_P2X_units,
EQ_Power_Balance_of_X2P_units,                                  
EQ_Max_Power_Consumption_of_BS_units,
EQ_Power_available,
EQ_Reserve_2U_capability,
EQ_Reserve_2D_capability,
EQ_Reserve_3U_capability,
EQ_2U_limit_chp,
EQ_2D_limit_chp,
EQ_3U_limit_chp,
EQ_Storage_minimum,
EQ_Storage_alert,
EQ_Storage_flood_control,
EQ_Storage_level,
EQ_Storage_input,
EQ_Storage_balance,
$If %MTS% == 1 EQ_Storage_Cyclic,
EQ_Storage_boundaries,
*EQ_H2_demand,
EQ_Boundary_Sector_Storage_MaxDischarge,
EQ_Boundary_Sector_Storage_MaxCharge,
EQ_Boundary_Sector_Storage_PowerMax,
EQ_Boundary_Sector_Storage_PowerMin,
EQ_Boundary_Sector_Storage_minimum,
EQ_Boundary_Sector_Storage_level,
EQ_Boundary_Sector_Storage_alert,
EQ_Boundary_Sector_Flood_Control,
EQ_Boundary_Sector_Storage_balance,
EQ_Boundary_Sector_Storage_boundaries,
$If %MTS% == 1 EQ_Boundary_Sector_Storage_Cyclic,
EQ_Storage_MaxCharge,
EQ_Storage_MaxDischarge ,
EQ_SystemCost,
*EQ_Emission_limits,
EQ_Flow_limits_lower,
EQ_Flow_limits_upper,
EQ_BS_Flow_limits_lower,
EQ_BS_Flow_limits_upper,
*$If not %MTS% == 1 EQ_Force_Commitment,
*$If not %MTS% == 1 EQ_Force_DeCommitment,
EQ_LoadShedding,
EQ_Residual_Load,
$If %ActivateFlexibleDemand% == 1 EQ_Flexible_Demand,
$If %ActivateFlexibleDemand% == 1 EQ_Flexible_Demand_Max,
$if not %ActivateFlexibleDemand% == 1 EQ_No_Flexible_Demand,
EQ_Flexible_Demand_Modulation_Min,
EQ_Flexible_Demand_Modulation_Max,
EQ_Tot_Flex_Demand_BS,
EQ_Max_Flex_Capacity_BS,
EQ_BS_Flex_Demand,
EQ_Tot_Flex_Supply_BS,
EQ_Max_Flex_Supply_BS,
EQ_BS_Spillage_limits_upper,
$If %RetrieveStatus% == 1 EQ_CommittedCalc,
EQ_DC_Power_Flow,
EQ_Total_Injected_Power,

*new
$If %MTS% == 0 EQ_SysInertia,
$If %MTS% == 0 EQ_Inertia_limit,
$If %MTS% == 0 EQ_SystemGain,
$If %MTS% == 0 EQ_SystemGain_limit,
$If %MTS% == 0 EQ_PrimaryReserve_Available,
$If %MTS% == 0 EQ_PrimaryReserve_Capability,
$If %MTS% == 0 EQ_PrimaryReserve_Boundary,
$If %MTS% == 0 EQ_Demand_balance_PrimaryReserve,

$If %MTS% == 0 EQ_FFRGain,
$If %MTS% == 0 EQ_FFRGain_limit,
$If %MTS% == 0 EQ_FFR_Available,
$If %MTS% == 0 EQ_FFR_Capability,
$If %MTS% == 0 EQ_FFR_Boundary,
$If %MTS% == 0 EQ_Demand_balance_FFR,
*MADRID
$If %MTS% == 0 EQ_PowerLoss,
*$If %MTS% == 0 EQ_MaxInstantGenerator,

EQ_Tot_Demand_2U,
*EQ_Load_Ramp,
/
;
UCM_SIMPLE.optcr = 0.01;
UCM_SIMPLE.optfile=1;

*===============================================================================
*Solving loop
*===============================================================================
ndays = floor(card(h)*TimeStep/24);

if (Config("RollingHorizon LookAhead","day") > ndays -1, abort "The look ahead period is longer than the simulation length";);
if (mod(Config("RollingHorizon LookAhead","day")*24,TimeStep) <> 0, abort "The look ahead period is not a multiple of TimeStep";);
if (mod(Config("RollingHorizon Length","day")*24,TimeStep) <> 0, abort "The rolling horizon length is not a multiple of TimeStep";);

* Fixing the initial guesses:
*PowerH.L(u,i)=PowerInitial(u);
*Committed.L(u,i)=CommittedInitial(u);

* Defining a parameter that records the solver status:
set  tmp   "tpm"  / "model", "solver"/  ;
parameter status(tmp,h);

scalar starttime;
set days /1,'ndays'/;
display days,ndays,TimeStep;
PARAMETER elapsed(days);

FOR(day = 1 TO ndays-Config("RollingHorizon LookAhead","day") by Config("RollingHorizon Length","day"),
         FirstHour = (day-1)*24/TimeStep+1;
         LastHour = min(card(h),FirstHour + (Config("RollingHorizon Length","day")+Config("RollingHorizon LookAhead","day")) * 24/TimeStep - 1);
         LastKeptHour = LastHour - Config("RollingHorizon LookAhead","day") * 24/TimeStep;
         i(h) = no;
         i(h)$(ord(h)>=firsthour and ord(h)<=lasthour)=yes;
         display day,FirstHour,LastHour,LastKeptHour;

*        Defining the minimum level at the end of the horizon :
         StorageFinalMin(s) =  sum(i$(ord(i)=card(i)),StorageProfile(s,i)*StorageCapacity(s)*Nunits(s)*AvailabilityFactor(s,i));
         StorageFinalMin(chp) =  sum(i$(ord(i)=card(i)),StorageProfile(chp,i)*StorageCapacity(chp)*Nunits(chp)*AvailabilityFactor(chp,i));
$If %MTS% == 0     SectorXStorageFinalMin(nx) = sum(i$(ord(i)=card(i)),SectorXStorageProfile(nx,i)*SectorXStorageCapacity(nx));
*$If %MTS% == 0     SectorXStorageInitial(nx) = sum(i$(ord(i)=1),SectorXStorageProfile(nx,i)*SectorXStorageCapacity(nx));

* Display PowerInitial,CommittedInitial,StorageFinalMin;

$If %LPFormulation% == 1          SOLVE UCM_SIMPLE USING LP MINIMIZING SystemCostD;
$If not %LPFormulation% == 1      SOLVE UCM_SIMPLE USING MIP MINIMIZING SystemCostD;

* Display EQ_Objective_function.M, EQ_CostRampUp.M, EQ_CostRampDown.M, EQ_Demand_balance_DA.M, EQ_Storage_minimum.M, EQ_Storage_level.M, EQ_Storage_input.M, EQ_Storage_balance.M, EQ_Storage_boundaries.M,  EQ_Flow_limits_lower.M ;


         status("model",i) = UCM_SIMPLE.Modelstat;
         status("solver",i) = UCM_SIMPLE.Solvestat;

         CommittedInitial(au) = sum(i$(ord(i)=LastKeptHour-FirstHour+1),Committed.L(au,i));
         PowerInitial(u) = sum(i$(ord(i)=LastKeptHour-FirstHour+1),Power.L(u,i));
         StorageInitial(s) =   sum(i$(ord(i)=LastKeptHour-FirstHour+1),StorageLevel.L(s,i));
         StorageInitial(chp) =   sum(i$(ord(i)=LastKeptHour-FirstHour+1),StorageLevel.L(chp,i));
$If %MTS% == 0   SectorXStorageInitial(nx) = sum(i$(ord(i)=LastKeptHour-FirstHour+1),SectorXStorageLevel.L(nx,i));
         SectorXFlexDemandInputInitial(nx) = sum(i$(ord(i)<LastKeptHour+1), SectorXFlexDemandInput(nx,i) - SectorXFlexDemand.L(nx,i));
         SectorXFlexSupplyInputInitial(nx) = sum(i$(ord(i)<LastKeptHour+1), SectorXFlexSupplyInput(nx,i) - SectorXFlexSupply.L(nx,i));
$If %ActivateFlexibleDemand% == 1 AccumulatedOverSupply_inital(n) = sum(i$(ord(i)=LastKeptHour-FirstHour+1),AccumulatedOverSupply.L(n,i));

* Assigning StorageLevelViolation (one value per optimization horizon) to the last element of StorageLevelViolation_H
StorageSlack.L(au,i)$(ord(i)=LastKeptHour-FirstHour+1) = Waterslack.L(au);
StorageLevelViolation_H.L(au,i)$(ord(i)=LastKeptHour-FirstHour+1) = StorageLevelViolation.L(au);
SectorXStorageLevelViolation_H.L(nx,i)$(ord(i)=LastKeptHour-FirstHour+1) = SectorXStorageLevelViolation.L(nx);
ObjectiveFunction.L(i)$(ord(i)=LastKeptHour-FirstHour+1) = SystemCostD.L;
Error.L = sum((i,n), CostLoadShedding(n,i)*ShedLoad.L(n,i)
          +Config("ValueOfLostLoad","val")*(LL_MaxPower.L(n,i)+LL_MinPower.L(n,i))
          +0.8*Config("ValueOfLostLoad","val")*(LL_2U.L(n,i)+LL_2D.L(n,i)+LL_3U.L(n,i)))
          +sum((au,i), 0.7*Config("ValueOfLostLoad","val")*(LL_RampUp.L(au,i)+LL_RampDown.L(au,i)))
;
OptimalityGap.L(i)$(ord(i)=LastKeptHour-FirstHour+1) = UCM_SIMPLE.objVal - UCM_SIMPLE.objEst;
OptimizationError.L(i)$(ord(i)=LastKeptHour-FirstHour+1) = Error.L - OptimalityGap.L(i);


*Loop variables to display after solving:
* Display LastKeptHour,PowerInitial,StorageInitial;
);

* Display PowerX.L,Flow.L,Power.L,Committed.L,ShedLoad.L,CurtailedPower.L,CurtailmentReserve_2U.L, CurtailmentReserve_3U.L,CurtailedHeat.L,StorageLevel.L,StorageInput.L,SystemCost.L,LL_MaxPower.L,LL_MinPower.L,LL_2U.L,LL_2D.L,LL_RampUp.L,LL_RampDown.L;

*===============================================================================
*Result export
*===============================================================================

PARAMETER
OutputDemand_2U(n,h)
OutputDemand_3U(n,h)
OutputDemand_2D(n,h)
OutputMaxOutageUp(n,h)
OutputMaxOutageDown(n,h)
OutputCommitted(au,h)
OutputFlow(l,h)
OutputInjectedPower(n,h)
OutputFlowX(lx,h)
OutputPower(au,h)
OutputPowerX(nx,au,h)
OutputPowerConsumption(p2x,h)
OutputResidualLoad(n,h)
OutputStorageInput(au,h)
OutputStorageLevel(au,h)
OutputStorageSlack(au,h)
OutputStorageLevelViolation_H(au,h)
OutputSectorXStorageLevel(nx,h)
OutputSectorXSelfDischarge(nx,h)
OutputSectorXStorageShadowPrice(nx,h)
OutputSectorXStorageLevelViolation_H(nx,h)
OutputSectorXStorageInput(nx,h)
$If %MTS% == 1 OutputSectorXStorageFinalMin(nx)
$If %MTS% == 1 OutputSectorXStorageInitial(nx)
OutputSectorXWaterNotWithdrawn(nx,h)
OutputSectorXSpillage(slx,h)
OutputSystemCost(h)
OutputSpillage(au,h)
OutputShedLoad(n,h)
OutputCurtailedPower(n,h)
OutputCurtailmentReserve_2U(n,h)
OutputCurtailmentReserve_3U(n,h)
OutputCurtailmentPerUnit(u,h)
$If %ActivateFlexibleDemand% == 1 OutputDemandModulation(n,h)
$If %ActivateFlexibleDemand% == 1 OutputAccumulatedOverSupply(n,h)
$If %ActivateFlexibleDemand% == 1 ShadowPriceDemandModulation(n,h)
ShadowPrice(n,h)
Demand_Balance_DA(n,h)
SectorXShadowPrice(nx,h)
LostLoad_MaxPower(n,h)
LostLoad_MinPower(n,h)
LostLoad_2D(n,h)
LostLoad_2U(n,h)
LostLoad_3U(n,h)
LostLoad_StorageAlert(au,h)
LostLoad_FloodControl(au,h)
$If not %LPFormulation% == 1 LostLoad_RampUp(n,h)
$If not %LPFormulation% == 1 LostLoad_RampDown(n,h)
$If not %LPFormulation% == 1 LostLoad_RampUp_Unit(au,z)
$If not %LPFormulation% == 1 LostLoad_RampDown_Unit(au,z)
OutputGenMargin(n,h)
OutputHeat(au,h)
LostLoad_WaterSlack(au)
OutputXNotServed(nx,h)
LostLoad_StorageLevelViolation(au)
LostLoad_SectorXStorageLevelViolation(nx)
OutputSectorXStorageAlertViolation(nx,h)
OutputSectorXFloodControlViolation(nx,h)
StorageShadowPrice(au,h)
OutputSectorXFlexDemand(nx,h)
OutputSectorXFlexSupply(nx,h)
OutputPowerMustRun(u,h)
$If not %LPFormulation% == 1 OutputCostStartUpH(au,h)
$If not %LPFormulation% == 1 OutputCostShutDownH(au,h)
$If not %LPFormulation% == 1 OutputCostRampUpH(au,h)
$If not %LPFormulation% == 1 OutputCostRampDownH(au,h)
ShadowPrice_2U(n,h)
ShadowPrice_2D(n,h)
ShadowPrice_3U(n,h)
OutputReserve_2U(au,h)
OutputReserve_2D(au,h)
OutputReserve_3U(au,h)
ShadowPrice_RampUp_TC(au,h)
ShadowPrice_RampDown_TC(au,h)
OutputRampRate(au,h)
OutputStartUp(au,h)
OutputShutDown(au,h)
OutputEmissions(n,p,z)
CapacityMargin(n,h)
OutputSystemCostD(h)
OutputOptimalityGap(h)
OutputOptimizationError(h)
OutputOptimizationCheck(h)
UnitHourlyPowerRevenue(au,h)
UnitHourly2URevenue(au,h)
UnitHourly2DRevenue(au,h)
UnitHourly3URevenue(au,h)
UnitHourlyHeatRevenue(au,h)
UnitHourlyRevenue(au,h)
UnitHourlyFixedCost(au,h)
UnitHourlyVariableCost(au,h)
UnitHourlyStartUpCost(au,h)
UnitHourlyShutDownCost(au,h)
UnitHourlyRampingCost(au,h)
UnitHourlyProductionCost(au,h)
UnitHourlyProfit(au,h)

*New
$If %MTS% == 0 OutputSysInertia(h)
$If %MTS% == 0 OutputSystemGain(h)
$If %MTS% == 0 OutputPrimaryReserve_Available(cu,h)
$If %MTS% == 0 OutputFFRGain(h)
$If %MTS% == 0 OutputFFR_Available(ba,h)
*MADRID
$If %MTS% == 0 OutputPowerLoss(h)
*$If %MTS% == 0 OutputMaxInstantGenerator(h)

OutputTotalDemand_2U(h)
*OutputLoadRamp(n,h)
;

$If %ActivateAdvancedReserves% == 2 OutputMaxOutageUp(n,z)=max(smax((au,tc),PowerCapacity(au)/Nunits(au)*Technology(au,tc)*LoadMaximum(au,z)*Location(au,n)), smax(l,FlowMaximum(l,z)*LineNode(l,n))$(card(l)>0));
$If %ActivateAdvancedReserves% == 2 OutputMaxOutageDown(n,z)=max(smax(s,StorageChargingCapacity(s)/Nunits(s)*Location(s,n)), smax(l,-FlowMaximum(l,z)*LineNode(l,n))$(card(l)>0));
$If %ActivateAdvancedReserves% == 2 OutputDemand_2U(n,z)=(Demand("2U",n,z) + OutputMaxOutageUp(n,z))*(1-K_QuickStart(n));
$If %ActivateAdvancedReserves% == 2 OutputDemand_3U(n,z)=Demand("2U",n,z) + OutputMaxOutageUp(n,z);
$If %ActivateAdvancedReserves% == 2 OutputDemand_2D(n,z)=Demand("2D",n,z) + OutputMaxOutageDown(n,z);
$If %ActivateAdvancedReserves% == 1 OutputMaxOutageUp(n,z)=smax((au,tc),PowerCapacity(au)/Nunits(au)*Technology(au,tc)*LoadMaximum(au,z)*Location(au,n));
$If %ActivateAdvancedReserves% == 1 OutputMaxOutageDown(n,z)=smax(s,StorageChargingCapacity(s)/Nunits(s)*Location(s,n));
$If %ActivateAdvancedReserves% == 1 OutputDemand_2U(n,z)=(Demand("2U",n,z) + OutputMaxOutageUp(n,z))*(1-K_QuickStart(n));
$If %ActivateAdvancedReserves% == 1 OutputDemand_3U(n,z)=Demand("2U",n,z) + OutputMaxOutageUp(n,z);
$If %ActivateAdvancedReserves% == 1 OutputDemand_2D(n,z)=Demand("2D",n,z) + OutputMaxOutageDown(n,z);

*New
*$If %ActivateAdvancedReserves% == 0 OutputDemand_2U(n,z)=(Demand("2U",n,z))*(1-K_QuickStart(n));
$If %ActivateAdvancedReserves% == 0 OutputDemand_2U(n,z)=Demand("2U",n,z);

*New
*$If %ActivateAdvancedReserves% == 0 OutputDemand_3U(n,z)=Demand("2U",n,z);
$If %ActivateAdvancedReserves% == 0 OutputDemand_3U(n,z)=Demand("2U",n,z);

$If %ActivateAdvancedReserves% == 0 OutputDemand_2D(n,z)=Demand("2D",n,z);


OutputCommitted(au,z)=Committed.L(au,z);
OutputFlow(l,z)=Flow.L(l,z);
OutputInjectedPower(n,z)=InjectedPower.L(z,n);
OutputFlowX(lx,z)=FlowX.L(lx,z);
OutputPower(au,z)=Power.L(au,z);
OutputPowerX(nx,au,z)=PowerX.L(nx,au,z);
OutputPowerConsumption(p2x,z)=PowerConsumption.L(p2x,z);
OutputResidualLoad(n,z)=ResidualLoad.L(n,z);
OutputHeat(au,z)=Heat.L(au,z);
OutputXNotServed(nx,z) = XNotServed.L(nx,z);
OutputStorageInput(s,z)=StorageInput.L(s,z);
OutputStorageLevel(s,z)=StorageLevel.L(s,z)/max(1,StorageCapacity(s)*Nunits(s)*AvailabilityFactor(s,z));
OutputStorageSlack(au,z) = StorageSlack.L(au,z);
OutputStorageLevelViolation_H(au,z) = StorageLevelViolation_H.L(au,z);
OutputSectorXStorageLevel(nx,z) = SectorXStorageLevel.L(nx,z)/max(1,SectorXStorageCapacity(nx));
*OutputSectorXStorageLevel(nx,z)$(SectorXStorageHours(nx)>=8) = SectorXStorageLevel.L(nx,z)/max(1,SectorXStorageCapacity(nx));
OutputSectorXSelfDischarge(nx,z) = SectorXStorageSelfDischarge(nx)*SectorXStorageLevel.L(nx,z)*TimeStep;
OutputSectorXStorageShadowPrice(nx,z) = EQ_Boundary_Sector_Storage_balance.m(nx,z);
OutputSectorXStorageLevelViolation_H(nx,z) = SectorXStorageLevelViolation_H.l(nx,z);
$If %MTS% == 1 OutputSectorXStorageFinalMin(nx) = SectorXStorageFinalMin.L(nx)/max(1,SectorXStorageCapacity(nx));
$If %MTS% == 1 OutputSectorXStorageInitial(nx) = SectorXStorageInitial.L(nx)/max(1,SectorXStorageCapacity(nx));
OutputSectorXStorageInput(nx,z) = SectorXStorageInput.l(nx,z);
OutputSectorXWaterNotWithdrawn(nx,z) = SectorXWaterNotWithdrawn.l(nx,z);
OutputSectorXSpillage(slx,z) = SectorXSpillage.l(slx,z);
OutputSystemCost(z)=SystemCost.L(z);
OutputSpillage(au,z)  = Spillage.L(au,z) ;
OutputShedLoad(n,z) = ShedLoad.L(n,z);
OutputCurtailedPower(n,z)=CurtailedPower.L(n,z);
OutputCurtailmentReserve_2U(n,z)=CurtailmentReserve_2U.L(n,z);
OutputCurtailmentReserve_3U(n,z)=CurtailmentReserve_3U.L(n,z);
OutputCurtailmentPerUnit(u,z)=(Nunits(u)*PowerCapacity(u)*LoadMaximum(u,z)-Power.L(u,z))$(sum(tr,Technology(u,tr))>=1);
$If %ActivateFlexibleDemand% == 1 OutputDemandModulation(n,z)=DemandModulation.L(n,z);
$If %ActivateFlexibleDemand% == 1 OutputAccumulatedOverSupply(n,z)=AccumulatedOverSupply.L(n,z);
$If %ActivateFlexibleDemand% == 1 ShadowPriceDemandModulation(n,z)=EQ_Flexible_Demand.m(n,z);
LostLoad_MaxPower(n,z)  = LL_MaxPower.L(n,z);
LostLoad_MinPower(n,z)  = LL_MinPower.L(n,z);
LostLoad_2D(n,z) = LL_2D.L(n,z);
LostLoad_2U(n,z) = LL_2U.L(n,z);
LostLoad_3U(n,z) = LL_3U.L(n,z);
LostLoad_StorageAlert(au,z) = LL_StorageAlert.L(au,z);
LostLoad_FloodControl(au,z) = LL_FloodControl.L(au,z);
OutputSectorXStorageAlertViolation(nx,z) = SectorXStorageAlertViolation.L(nx,z);
OutputSectorXFloodControlViolation(nx,z) = SectorXFloodControlViolation.L(nx,z);
$If not %LPFormulation% == 1 LostLoad_RampUp(n,z)    = sum(au,LL_RampUp.L(au,z)*Location(au,n));                                                                
$If not %LPFormulation% == 1 LostLoad_RampDown(n,z)  = sum(au,LL_RampDown.L(au,z)*Location(au,n));
$If not %LPFormulation% == 1 LostLoad_RampUp_Unit(au,z) = LL_RampUp.L(au,z);
$If not %LPFormulation% == 1 LostLoad_RampDown_Unit(au,z) = LL_RampDown.L(au,z);
ShadowPrice(n,z) = EQ_Demand_balance_DA.m(n,z);
Demand_balance_DA(n,z) = EQ_Demand_balance_DA.L(n,z);
LostLoad_WaterSlack(au) = WaterSlack.L(au);
SectorXShadowPrice(nx,z) = EQ_BS_Demand_balance.m(nx,z);
LostLoad_StorageLevelViolation(au) = StorageLevelViolation.L(au);
LostLoad_SectorXStorageLevelViolation(nx) = SectorXStorageLevelViolation.L(nx);
StorageShadowPrice(s,z) = 0 ;
OutputSectorXFlexDemand(nx,z) = SectorXFlexDemand.L(nx,z);
OutputSectorXFlexSupply(nx,z) = SectorXFlexSupply.L(nx,z);
StorageShadowPrice(s,z) = EQ_Storage_balance.m(s,z);
OutputPowerMustRun(u,z) = PowerMustRun(u,z);
$If not %LPFormulation% == 1 OutputCostStartUpH(au,z) = CostStartUpH.L(au,z);
$If not %LPFormulation% == 1 OutputCostShutDownH(au,z) = CostShutDownH.L(au,z);
$If not %LPFormulation% == 1 OutputCostRampUpH(au,z) = CostRampUpH.L(au,z);
$If not %LPFormulation% == 1 OutputCostRampDownH(au,z) = CostRampDownH.L(au,z);

ShadowPrice_2U(n,z) =  EQ_Demand_balance_2U.m(n,z);
ShadowPrice_2D(n,z) =  EQ_Demand_balance_2D.m(n,z);
ShadowPrice_3U(n,z) =  EQ_Demand_balance_3U.m(n,z);

OutputReserve_2U(au,z) = Reserve_2U.L(au,z);
OutputReserve_2D(au,z) = Reserve_2D.L(au,z);
OutputReserve_3U(au,z) = Reserve_3U.L(au,z);

ShadowPrice_RampUp_TC(u,z) = EQ_RampUp_TC.m(u,z);
ShadowPrice_RampDown_TC(u,z) = EQ_RampDown_TC.m(u,z);
OutputRampRate(u,z) = - Power.L(u,z-1)$(ord(z) > 1) - PowerInitial(u)$(ord(z) = 1) + Power.L(u,z);
OutputStartUp(au,z) = StartUp.L(au,z);
OutputShutDown(au,z) = ShutDown.L(au,z);

*New
$If %MTS%==0 OutputSysInertia(z) = SysInertia.L(z);
$If %MTS%==0 OutputSystemGain(z) = SystemGain.L(z);
$If %MTS%==0 OutputPrimaryReserve_Available(cu,z) = PrimaryReserve_Available.L(cu,z);
$If %MTS%==0 OutputFFRGain(z) = FFRGain.L(z);
$If %MTS%==0 OutputFFR_Available(ba,z) = FFR_Available.L(ba,z);
*OutputSysInertia(z) = SysInertia.L(z);
*OutputSystemGain(z) = SystemGain.L(z);
*OutputPrimaryReserve_Available(cu,z) = PrimaryReserve_Available.L(cu,z);
*MADRID
$If %MTS%==0 OutputPowerLoss(z) = PowerLoss.L(z);
*OutputPowerLoss(z) = PowerLoss.L(z);
*$If %MTS%==0 OutputMaxInstantGenerator(z) = MaxInstantGenerator.L(z);

*New
OutputTotalDemand_2U(z) = TotalDemand_2U.L(z);
*OutputLoadRamp(n,z) = LoadRamp.L(n,z);

OutputEmissions(n,p,z) = sum(u,Power.L(u,z)*EmissionRate(u,p)*Location(u,n));
                        

OutputSystemCostD(z) = ObjectiveFunction.L(z);
OutputOptimalityGap(z) = OptimalityGap.L(z);
OutputOptimizationError(z) = OptimizationError.L(z);
OutputOptimizationCheck(z) = OptimizationError.L(z) - OptimalityGap.L(z);
UnitHourlyPowerRevenue(au,z) = sum(n, EQ_Demand_balance_DA.m(n,z) * Location(au,n) * Power.L(au,z));
UnitHourly2URevenue(au,z) = sum(n, OutputReserve_2U(au,z) * ShadowPrice_2U(n,z) * Location(au,n));
UnitHourly2DRevenue(au,z) = sum(n, OutputReserve_2D(au,z) * ShadowPrice_2D(n,z) * Location(au,n));
UnitHourly3URevenue(au,z) = sum(n, OutputReserve_3U(au,z) * ShadowPrice_3U(n,z) * Location(au,n));
UnitHourlyRevenue(au,z) = UnitHourlyPowerRevenue(au,z) + UnitHourly2URevenue(au,z) + UnitHourly2DRevenue(au,z) + UnitHourly3URevenue(au,z);
UnitHourlyFixedCost(u,z) = Committed.L(u,z) * CostFixed(u);
UnitHourlyVariableCost(au,z) = Power.L(au,z) * CostVariable(au,z);
UnitHourlyStartUpCost(u,z) = StartUp.L(u,z) * CostStartUp(u);
UnitHourlyShutDownCost(u,z) = ShutDown.L(u,z) * CostShutDown(u);
UnitHourlyRampingCost(au,z) = CostRampUpH.L(au,z) + CostRampDownH.L(au,z);
UnitHourlyProductionCost(au,z) = sum(u, UnitHourlyFixedCost(u,z) + UnitHourlyStartUpCost(u,z) + UnitHourlyShutDownCost(u,z) + UnitHourlyRampingCost(u,z))
                                + UnitHourlyVariableCost(au,z);
UnitHourlyProfit(au,z) = UnitHourlyRevenue(au,z) - UnitHourlyProductionCost(au,z);


EXECUTE_UNLOAD "Results.gdx"
OutputCommitted,
OutputFlow,
OutputInjectedPower,
OutputFlowX,
OutputPower,
OutputPowerX,
OutputPowerConsumption,
OutputResidualLoad,
OutputHeat,
OutputXNotServed,
OutputStorageInput,
OutputStorageLevel,
OutputStorageSlack,
OutputStorageLevelViolation_H,
OutputSectorXStorageLevel,
OutputSectorXSelfDischarge,
OutputSectorXStorageShadowPrice,
OutputSectorXStorageLevelViolation_H,
OutputSectorXStorageInput,
OutputSectorXWaterNotWithdrawn,
OutputSectorXSpillage,
$If %MTS% == 1 OutputSectorXStorageFinalMin,
$If %MTS% == 1 OutputSectorXStorageInitial,
OutputSystemCost,
OutputSpillage,
OutputShedLoad,
OutputCurtailedPower,
OutputCurtailmentReserve_2U,
OutputCurtailmentReserve_3U,
OutputCurtailmentPerUnit,
$If %ActivateFlexibleDemand% == 1 OutputDemandModulation,
$If %ActivateFlexibleDemand% == 1 OutputAccumulatedOverSupply,
$If %ActivateFlexibleDemand% == 1 ShadowPriceDemandModulation,
OutputGenMargin,
LostLoad_MaxPower,
LostLoad_MinPower,
LostLoad_2D,
LostLoad_2U,
LostLoad_3U,
LostLoad_StorageAlert,
LostLoad_FloodControl,
OutputSectorXStorageAlertViolation,
OutputSectorXFloodControlViolation,
$If not %LPFormulation% == 1 LostLoad_RampUp,
$If not %LPFormulation% == 1 LostLoad_RampDown,
$If not %LPFormulation% == 1 LostLoad_RampUp_Unit,
$If not %LPFormulation% == 1 LostLoad_RampDown_Unit,
ShadowPrice,
Demand_balance_DA,
ShadowPrice_2U,
ShadowPrice_2D,
ShadowPrice_3U,
LostLoad_WaterSlack,
LostLoad_StorageLevelViolation,
LostLoad_SectorXStorageLevelViolation,
StorageShadowPrice,
OutputSectorXFlexDemand,
OutputSectorXFlexSupply,
SectorXShadowPrice,
OutputPowerMustRun,
$If not %LPFormulation% == 1 OutputCostStartUpH,
$If not %LPFormulation% == 1 OutputCostShutDownH,
$If not %LPFormulation% == 1 OutputCostRampUpH,
$If not %LPFormulation% == 1 OutputCostRampDownH,
ShadowPrice_2U,
ShadowPrice_2D,
ShadowPrice_3U,
OutputReserve_2U,
OutputReserve_2D,
OutputReserve_3U,
ShadowPrice_RampUp_TC,
ShadowPrice_RampDown_TC,
OutputRampRate,
OutputStartUp,
OutputShutDown,
OutputEmissions,
CapacityMargin,
OutputSystemCostD,
OutputOptimalityGap,
OutputOptimizationError,
OutputOptimizationCheck,
OutputMaxOutageUp,
OutputMaxOutageDown,
OutputDemand_2U,
OutputDemand_3U,
OutputDemand_2D,

*New
$If %MTS%==0 OutputSysInertia,
$If %MTS%==0 OutputSystemGain,
$If %MTS%==0 OutputPrimaryReserve_Available,
$If %MTS%==0 OutputFFRGain,
$If %MTS%==0 OutputFFR_Available,
*MADRID
$If %MTS%==0 OutputPowerLoss,
*$If %MTS%==0 OutputMaxInstantGenerator,

status,
UnitHourlyPowerRevenue
UnitHourly2URevenue
UnitHourly2DRevenue
UnitHourly3URevenue
UnitHourlyRevenue
UnitHourlyFixedCost
UnitHourlyVariableCost
UnitHourlyStartUpCost
UnitHourlyShutDownCost
UnitHourlyRampingCost
UnitHourlyProductionCost
UnitHourlyProfit

;

*display OutputPowerConsumption, heat.L, heatslack.L, powerconsumption.L, power.L;
*$If %MTS%==1 display OutputPowerConsumption, heat.L, powerconsumption.L, power.L, EQ_Storage_Cyclic.L, EQ_Boundary_Sector_Storage_Cyclic.L;

$onorder

$exit

