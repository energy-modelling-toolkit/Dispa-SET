$Title UCM model

$eolcom //
Option threads=12;
Option IterLim=1000000000;
Option ResLim = 10000000000;
*Option optca=0.0;


// Reduce .lst file size

// Turn off the listing of the input file
$offlisting
$offlog

// Turn off the listing and cross-reference of the symbols used
$offsymxref offsymlist

option
    limrow = 0,     // equations listed per block
    limcol = 0,     // variables listed per block
    solprint = off,     // solver's solution output printed
    sysout = off;       // solver's system output printed

$onempty

*===============================================================================
*Definition of the dataset-related options
*===============================================================================

* Print results to excel files (0 for no, 1 for yes)
$set Verbose 0

* Set debug mode. !! This breaks the loop and requires a debug.gdx file !!
* (0 for no, 1 for yes)
$set Debug 0

* Print results to excel files (0 for no, 1 for yes)
$set PrintResults 0

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

*===============================================================================
*Definition of   sets and parameters
*===============================================================================
SETS
mk               Markets
n                Nodes
n_th             Thermal nodes
n_bs             Boundary sector nodes
l                Lines
l_bs             Boundary sector lines
au               All Units
u(au)            Generation units
t                Generation technologies
tr(t)            Renewable generation technologies
tc(t)
f                Fuel types
p                Pollutants
p2bs(au)         Power to X
bsu(au)          Boundary sector only units
s(au)             Storage Units (with reservoir)
chp(u)           CHP units
p2h(au)          Power to heat units
th(au)           Units with thermal storage
hu(au)           Heat only units
thms(au)         Thermal storage units only
h                Hours
i(h)             Subset of simulated hours for one iteration
wat(au)           hydro technologies
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
*CostCurtailment(n,h)                       [EUR\MW]        Curtailment costs
CostFixed(au)                               [EUR\h]         Fixed costs
CostRampUp(u)                               [EUR\MW]        Ramp-up costs
CostRampDown(u)                             [EUR\MW]        Ramp-down costs
CostShutDown(u)                             [EUR\u]         Shut-down costs
CostStartUp(u)                              [EUR\u]         Start-up costs
CostVariable(au,h)                          [EUR\MW]        Variable costs
CostHeatSlack(n_th,h)                       [EUR\MWh]       Cost of supplying heat via other means
CostBoundarySectorSlack(n_bs,h)             [EUR\MWh]       Cost of supplying energy to boundary sector via other means
CostLoadShedding(n,h)                       [EUR\MWh]       Cost of load shedding
Curtailment(n)                              [n.a]           Curtailment allowed or not {1 0} at node n
CostCurtailment(n,h)                        [EUR\MWh]       Cost of VRES curtailment
Demand(mk,n,h)                              [MW]            Demand
Efficiency(au,h)                            [%]             Efficiency
EfficiencyBoundarySector(n_bs,au,h)         [%]             Discharge Efficiency boundary sector
ChargingEfficiencyBoundarySector(n_bs,au,h) [%]             Charging Efficiency boundary sector
EmissionMaximum(n,p)                        [tP]            Emission limit
EmissionRate(au,p)                          [tP\MWh]        P emission rate
FlowMaximum(l,h)                            [MW]            Line limits
FlowMinimum(l,h)                            [MW]            Minimum flow
FlowBSMaximum(l_bs,h)                       [MW]            Boundary sector line limits
FlowBSMinimum(l_bs,h)                       [MW]            Boundary sector line minimum flow
Fuel(u,f)                                   [n.a.]          Fuel type {1 0}
HeatDemand(n_th,h)                          [MWh\u]         Heat demand profile for chp units
BoundarySectorDemand(n_bs,h)                [MWh\n_bs]      Demand profile in boundary sectors
LineNode(l,n)                               [n.a.]          Incidence matrix {-1 +1}
BSLineNode(l_bs,n_bs)                       [n.a.]          Incidence matrix {-1 +1}
LoadShedding(n,h)                           [MW]            Load shedding capacity
Location(au,n)                              [n.a.]          Location {1 0}
Location_th(au,n_th)                        [n.a.]          Location {1 0}
Location_bs(au,n_bs)                        [n.a.]          Location {1 0}
Markup(u,h)                                 [EUR\MW]        Markup
OutageFactor(au,h)                          [%]             Outage Factor (100% = full outage)
PartLoadMin(au)                             [%]             Minimum part load
PowerCapacity(au)                           [MW\u]          Installed capacity
PowerInitial(u)                             [MW\u]          Power output before initial period
PowerMinStable(au)                          [MW\u]          Minimum power output
PriceTransmission(l,h)                      [EUR\MWh]       Transmission price
StorageChargingCapacity(au)                 [MW\u]          Storage capacity
StorageChargingEfficiency(au)               [%]             Charging efficiency
StorageSelfDischarge(au)                    [%\day]         Self-discharge of the storage units
RampDownMaximum(u)                          [MW\h\u]        Ramp down limit
RampShutDownMaximum(au)                     [MW\h\u]        Shut-down ramp limit
RampStartUpMaximum(au)                      [MW\h\u]        Start-up ramp limit
RampStartUpMaximumH(au,h)                   [MW\h\u]        Start-up ramp limit - Clustered formulation
RampShutDownMaximumH(au,h)                  [MW\h\u]        Shut-down ramp limit - Clustered formulation
RampUpMaximum(u)                            [MW\h\u]        Ramp up limit
Reserve(au)                                 [n.a.]          Reserve technology {1 0}
StorageCapacity(au)                         [MWh\u]         Storage capacity
StorageDischargeEfficiency(au)              [%]             Discharge efficiency
StorageOutflow(au,h)                        [MW\u]          Storage outflows
StorageInflow(au,h)                         [MW\u]          Storage inflows (potential energy)
StorageInitial(au)                          [MWh]           Storage level before initial period
StorageProfile(au,h)                        [%]             Storage level to be resepected at the end of each horizon
StorageMinimum(au)                          [MWh]           Storage minimum
Technology(au,t)                            [n.a.]          Technology type {1 0}
TimeDownMinimum(au)                         [h]             Minimum down time
TimeUpMinimum(au)                           [h]             Minimum up time
BSFlexDemandInput(n_bs,h)                   [MWh]           Flexible demand inside BS at each timestep (unless for MTS)
BSFlexDemandInputInitial(n_bs)              [MWh]           Cumulative flexible demand inside the loop
BSFlexMaxCapacity(n_bs)                     [MW]            Max capacity for BS Flexible demand
BSFlexSupplyInput(n_bs,h)                   [MWh]           Flexible demand inside BS at each timestep (unless for MTS)
BSFlexSupplyInputInitial(n_bs)              [MWh]           Cumulative flexible demand inside the loop
BSFlexMaxSupply(n_bs)                       [MW]            Max capacity for BS Flexible demand
$If %RetrieveStatus%==1 CommittedCalc(u,z)  [n.a.]          Committment status as for the MILP
Nunits(au)                                  [n.a.]          Number of units inside the cluster (upper bound value for integer variables)
K_QuickStart(n)                             [n.a.]          Part of the reserve that can be provided by offline quickstart units
QuickStartPower(au,h)                       [MW\h\u]        Available max capacity in tertiary regulation up from fast-starting power plants - TC formulation
BoundarySectorStorageCapacity(n_bs)         [MWh]           Storage capacity of the boundary sector
BoundarySectorStorageSelfDischarge(n_bs)    [%]             Boundary sector storage self discharge
BoundarySectorStorageMinimum(n_bs)          [MWh]           Boundary sector storage minimum
BoundarySectorStorageInitial(n_bs)          [MWh]           Boundary sector storage initial state of charge
BoundarySectorStorageProfile(n_bs,h)        [%]             Boundary sector storage level respected at the end of each horizon
;

*Parameters as used within the loop
PARAMETERS
CostLoadShedding(n,h)               [EUR\MW]        Value of lost load
LoadMaximum(au,h)                   [%]             Maximum load given AF and OF
PowerMustRun(au,h)                  [MW\u]          Minimum power output
StorageFinalMin(au)                 [MWh]           Minimum storage level at the end of the optimization horizon
BoundarySectorStorageFinalMin(n_bs) [MWh]           Minimum boundary sector storage level at the end of the optimization horizon
MaxFlexDemand(n)                    [MW]            Maximum value of the flexible demand parameter
MaxOverSupply(n,h)                  [MWh]           Maximum flexible demand accumultation
AccumulatedOverSupply_inital(n)     [MWh]           Initial value of the flexible demand accumulation
;

* Scalar variables necessary to the loop:
scalar FirstHour,LastHour,LastKeptHour,day,ndays,failed,srp,nsrp;
FirstHour = 1;
scalar TimeStep;

*Threshold values for p2h partecipation to reserve market as spinning/non-spinning reserves (TO BE IMPLEMENTED IN CONFIGFILE)
srp = 1;
nsrp = 3;

*===============================================================================
*Data import
*===============================================================================

$gdxin %inputfilename%

$LOAD mk
$LOAD n
$LOAD n_th
$LOAD n_bs
$LOAD l
$LOAD l_bs
$LOAD au
$LOAD u
$LOAD t
$LOAD tr
$LOAD f
$LOAD p
$LOAD s
$LOAD wat
$LOAD p2bs
$LOAD chp
$LOAD p2h
$LOAD th
$LOAD tc
$LOAD hu
$LOAD bsu
$LOAD thms
$LOAD h
$LOAD z
$LOAD AvailabilityFactor
$LOAD CHPPowerLossFactor
$LOAD CHPPowerToHeat
$LOAD CHPMaxHeat
$LOAD CHPType
$LOAD Config
$LOAD CostFixed
$LOAD CostHeatSlack
$LOAD CostBoundarySectorSlack
$LOAD CostLoadShedding
$LOAD CostShutDown
$LOAD CostStartUp
$LOAD CostVariable
$LOAD Curtailment
$LOAD CostCurtailment
$LOAD Demand
$LOAD BoundarySectorDemand
$LOAD StorageDischargeEfficiency
$LOAD Efficiency
$LOAD EfficiencyBoundarySector
$LOAD ChargingEfficiencyBoundarySector
$LOAD EmissionMaximum
$LOAD EmissionRate
$LOAD FlowMaximum
$LOAD FlowMinimum
$LOAD FlowBSMaximum
$LOAD FlowBSMinimum
$LOAD Fuel
$LOAD HeatDemand
$LOAD LineNode
$LOAD BSLineNode
$LOAD LoadShedding
$LOAD Location
$LOAD Location_th
$LOAD Location_bs
$LOAD Markup
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
$LOAD StorageOutflow
$LOAD Technology
$LOAD TimeDownMinimum
$LOAD TimeUpMinimum
$LOAD CostRampUp
$LOAD CostRampDown
$LOAD BSFlexDemandInput
$LOAD BSFlexDemandInputInitial
$LOAD BSFlexMaxCapacity
$LOAD BSFlexSupplyInput
$LOAD BSFlexSupplyInputInitial
$LOAD BSFlexMaxSupply
$LOAD BoundarySectorStorageCapacity
$LOAD BoundarySectorStorageSelfDischarge
$LOAD BoundarySectorStorageMinimum
$LOAD BoundarySectorStorageInitial
$LOAD BoundarySectorStorageProfile
$If %RetrieveStatus% == 1 $LOAD CommittedCalc
;

$If %Verbose% == 0 $goto skipdisplay

Display
mk,
n,
n_th,
n_bs,
l,
l_bs,
u,
t,
tr,
tc,
f,
p,
s,
p2bs,
wat,
chp,
p2h,
th,
h2,
hu,
thms,
h,
AvailabilityFactor,
CHPPowerLossFactor,
CHPPowerToHeat,
CHPMaxHeat,
CHPType,
Config,
CostFixed,
CostShutDown,
CostStartUp,
CostRampUp,
CostVariable,
Demand,
StorageDischargeEfficiency,
Efficiency,
EmissionMaximum,
EmissionRate,
FlowMaximum,
FlowMinimum,
FlowBSMaximum,
FlowBSMinimum,
Fuel,
HeatDemand,
LineNode,
BSLineNode,
Location,
Location_th,
Location_bs,
LoadShedding,
Markup,
OutageFactor,
PartLoadMin,
PowerCapacity,
PowerInitial,
PriceTransmission,
StorageChargingCapacity,
StorageChargingEfficiency,
StorageSelfDischarge,
RampDownMaximum,
RampShutDownMaximum,
RampStartUpMaximum,
RampUpMaximum,
Reserve,
StorageCapacity,
StorageInflow,
StorageInitial,
StorageProfile,
StorageMinimum,
StorageOutflow,
Technology,
TimeDownMinimum,
TimeUpMinimum,
BSFlexDemandInput,
BSFlexDemandInputInitial,
BSFlexMaxCapacity,
BSFlexSupplyInput,
BSFlexSupplyInputInitial,
BSFlexMaxSupply,
BoundarySectorStorageCapacity,
BoundarySectorStorageSelfDischarge,
BoundarySectorStorageMinimum,
BoundarySectorStorageInitial,
BoundarySectorStorageProfile,
$If %RetrieveStatus% == 1 , CommittedCalc
;

$label skipdisplay

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
AccumulatedOverSupply(n,h)          [MWh]   Accumulated oversupply due to the flexible demand
CostStartUpH(u,h)                   [EUR]   Cost of starting up
CostShutDownH(u,h)                  [EUR]   cost of shutting down
CostRampUpH(u,h)                    [EUR]   Ramping cost
CostRampDownH(u,h)                  [EUR]   Ramping cost
CurtailedPower(n,h)                 [MW]    Curtailed power at node n
CurtailmentReserve_2U(n,h) 			[MW]    Curtailed power used for reserves at node n
CurtailmentReserve_3U(n,h) 			[MW]    Curtailed power used for reserves at node n
CurtailedHeat(n_th,h)               [MW]    Curtailed heat at node n_th
Flow(l,h)                           [MW]    Flow through lines
FlowBS(l_bs,h)                      [MW]    Flow through boundary sector lines
Power(au,h)                         [MW]    Power output
PowerConsumption(au,h)              [MW]    Power consumption by P2X units
PowerMaximum(u,h)                   [MW]    Power output
PowerMinimum(u,h)                   [MW]    Power output
ShedLoad(n,h)                       [MW]    Shed load
StorageInput(au,h)                  [MWh]   Charging input for storage units
StorageLevel(au,h)                  [MWh]   Storage level of charge
BoundarySectorStorageLevel(n_bs,h)  [MWh]   Storage level of charge of the boundary sector
BoundarySectorSpillage(n_bs,h)      [MWh]   Spillage from boundary sector storage
LL_MaxPower(n,h)                    [MW]    Deficit in terms of maximum power
LL_RampUp(u,h)                      [MW]    Deficit in terms of ramping up for each plant
LL_RampDown(u,h)                    [MW]    Deficit in terms of ramping down
LL_MinPower(n,h)                    [MW]    Power exceeding the demand
LL_2U(n,h)                          [MW]    Deficit in reserve up
LL_3U(n,h)                          [MW]    Deficit in reserve up - non spinning
LL_2D(n,h)                          [MW]    Deficit in reserve down
spillage(au,h)                      [MWh]   spillage from water reservoirs
SystemCost(h)                       [EUR]   Hourly system cost
Reserve_2U(au,h)                    [MW]    Spinning reserve up
Reserve_2D(au,h)                    [MW]    Spinning reserve down
Reserve_3U(au,h)                    [MW]    Non spinning quick start reserve up
Heat(au,h)                          [MW]    Heat output by chp plant
HeatSlack(n_th,h)                   [MW]    Heat satisfied by other sources
BoundarySectorSlack(n_bs,h)         [MW]    Boundary sector demand satisfied by other sources
WaterSlack(au)                      [MWh]   Unsatisfied water level constraint at end of optimization period
StorageSlack(au,h)                  [MWh]   Unsatisfied storage level constraint at end of simulation timestep
BoundarySectorWaterSlack(n_bs)      [MWh]   Unsatisfied boundary sector water level constraint at end of optimization period
BoundarySectorStorageSlack(n_bs,h)  [MWh]   Unsatisfied boundary sector storage level constraint at end of simulation timestep
BSFlexDemand(n_bs,h)                [MW]    FLexible boundary sector demand at each time step of each n_bs node
BSFlexSupply(n_bs,h)                [MW]    FLexible boundary sector supply at each time step of each n_bs node
;

free variable
SystemCostD                         ![EUR]  Total system cost for one optimization period
DemandModulation(n,h)               [MW]    Difference between the flexible demand and the baseline
PowerBoundarySector(n_bs,au,h)      [MW]    Power output in boundary sector
BoundarySectorStorageInput(n_bs,h)  [MW]    Boundary Sector Storage input - output
ResidualLoad(n,h)                   [MW]    Residual Load
ObjectiveFunction(h)
OptimalityGap(h)
OptimizationError
Error
;

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

$If %Verbose% == 1 Display RampStartUpMaximum, RampShutDownMaximum, CommittedInitial;

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
EQ_Heat_Demand_balance
EQ_BS_Demand_balance
EQ_CHP_max_heat
EQ_Heat_Storage_balance
EQ_Heat_Storage_minimum
EQ_Heat_Storage_level
EQ_Heat_Storage_input
EQ_Heat_Storage_MaxDischarge
EQ_Heat_Storage_MaxCharge
EQ_Heat_Storage_boundaries
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
EQ_Max_Power_Consumption
EQ_Power_Balance_of_BS_units
EQ_Max_Power_Consumption_of_BS_units
EQ_Boundary_sector_only_power_available
EQ_Boundary_sector_only_power_available_min
EQ_Power_must_run
EQ_Power_available
EQ_Reserve_2U_capability
EQ_Reserve_2D_capability
EQ_Reserve_3U_capability
EQ_p2h_Reserve_2U_capability
EQ_p2h_Reserve_2D_capability
EQ_p2h_Reserve_3U_capability
EQ_2U_limit_p2h
EQ_2D_limit_p2h
EQ_3U_limit_p2h
EQ_2U_limit_chp
EQ_2D_limit_chp
EQ_3U_limit_chp
EQ_Storage_minimum
EQ_Storage_level
EQ_Storage_input
EQ_Storage_MaxDischarge
EQ_Storage_MaxCharge
EQ_Storage_balance
EQ_Storage_boundaries
EQ_Boundary_Sector_Storage_MaxDischarge
EQ_Boundary_Sector_Storage_MaxCharge
EQ_Boundary_Sector_Storage_minimum
EQ_Boundary_Sector_Storage_level
EQ_Boundary_Sector_Storage_balance
EQ_Boundary_Sector_Storage_boundaries
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
EQ_Tot_Flex_Supply_BS
EQ_Max_Flex_Supply_BS
EQ_Tot_Flex_Demand_BS
EQ_Max_Flex_Capacity_BS
EQ_BS_Flex_Demand
EQ_Curtailed_Power
EQ_Residual_Load
$If %RetrieveStatus% == 1 EQ_CommittedCalc
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
         +sum(u,CostRampUpH(u,i) + CostRampDownH(u,i))
         +sum(u,CostVariable(u,i) * Power(u,i)*TimeStep)
         +sum(hu,CostVariable(hu,i) * Heat(hu,i)*TimeStep)
         +sum(p2h,CostVariable(p2h,i) * Heat(p2h,i)*TimeStep)
         +sum((n_bs,bsu),CostVariable(bsu,i) * PowerBoundarySector(n_bs,bsu,i)*TimeStep)
         +sum(l,PriceTransmission(l,i)*Flow(l,i)*TimeStep)
         +sum(n,CostLoadShedding(n,i)*ShedLoad(n,i)*TimeStep)
         +sum(n_th, CostHeatSlack(n_th,i) * HeatSlack(n_th,i)*TimeStep)
         +sum(n_bs, CostBoundarySectorSlack(n_bs,i) * BoundarySectorSlack(n_bs,i)*TimeStep)
         +sum(chp, CostVariable(chp,i) * CHPPowerLossFactor(chp) * Heat(chp,i)*TimeStep)
         +Config("ValueOfLostLoad","val")*(sum(n,(LL_MaxPower(n,i)+LL_MinPower(n,i))*TimeStep))
         +0.8*Config("ValueOfLostLoad","val")*(sum(n,(LL_2U(n,i)+LL_2D(n,i)+LL_3U(n,i))*TimeStep))
         +0.7*Config("ValueOfLostLoad","val")*sum(u,(LL_RampUp(u,i)+LL_RampDown(u,i))*TimeStep)
         +Config("CostOfSpillage","val")*(sum(au,spillage(au,i)) + sum(n_bs,BoundarySectorSpillage(n_bs,i)))
         +sum(n,CurtailedPower(n,i) * CostCurtailment(n,i) * TimeStep)
;
$else

EQ_SystemCost(i)..
         SystemCost(i)
         =E=
         sum(au,CostFixed(au)*TimeStep*Committed(au,i))
         +sum(u,CostStartUpH(u,i) + CostShutDownH(u,i))
         +sum(u,CostRampUpH(u,i) + CostRampDownH(u,i))
         +sum(u,CostVariable(u,i) * Power(u,i)*TimeStep)
         +sum(hu,CostVariable(hu,i) * Heat(hu,i)*TimeStep)
         +sum(p2h,CostVariable(p2h,i) * Heat(p2h,i)*TimeStep)
         +sum((n_bs,bsu),CostVariable(bsu,i) * PowerBoundarySector(n_bs,bsu,i)*TimeStep)
         +sum(l,PriceTransmission(l,i)*Flow(l,i)*TimeStep)
         +sum(n,CostLoadShedding(n,i)*ShedLoad(n,i)*TimeStep)
         +sum(n_th, CostHeatSlack(n_th,i) * HeatSlack(n_th,i)*TimeStep)
         +sum(n_bs, CostBoundarySectorSlack(n_bs,i) * BoundarySectorSlack(n_bs,i)*TimeStep)
         +sum(chp, CostVariable(chp,i) * CHPPowerLossFactor(chp) * Heat(chp,i)*TimeStep)
         +Config("ValueOfLostLoad","val")*(sum(n,(LL_MaxPower(n,i)+LL_MinPower(n,i))*TimeStep))
         +0.8*Config("ValueOfLostLoad","val")*(sum(n,(LL_2U(n,i)+LL_2D(n,i)+LL_3U(n,i))*TimeStep))
         +0.7*Config("ValueOfLostLoad","val")*sum(u,(LL_RampUp(u,i)+LL_RampDown(u,i))*TimeStep)
         +Config("CostOfSpillage","val")*(sum(au,spillage(au,i)) + sum(n_bs,BoundarySectorSpillage(n_bs,i)))
         +sum(n,CurtailedPower(n,i) * CostCurtailment(n,i) * TimeStep)
;

$endIf
;

EQ_Objective_function..
         SystemCostD
         =E=
         sum(i,SystemCost(i))
         +Config("WaterValue","val")*sum(au,WaterSlack(au))
         +Config("WaterValue","val")*sum(n_bs,BoundarySectorWaterSlack(n_bs))
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
EQ_RampUp_TC(u,i)$(sum(tr,Technology(u,tr))=0)..
         Power(u,i-1)$(ord(i) > 1) + PowerInitial(u)$(ord(i) = 1) - Power(u,i)
         =L=
         (Committed(u,i) - StartUp(u,i)) * RampUpMaximum(u) * TimeStep
         + RampStartUpMaximumH(u,i) * TimeStep * StartUp(u,i)
         - PowerMustRun(u,i) * ShutDown(u,i)
         + LL_RampUp(u,i)
;

* ramp down constraints
EQ_RampDown_TC(u,i)$(sum(tr,Technology(u,tr))=0)..
         Power(u,i) - Power(u,i-1)$(ord(i) > 1) - PowerInitial(u)$(ord(i) = 1)
         =L=
         (Committed(u,i) - StartUp(u,i)) * RampDownMaximum(u) * TimeStep
         - PowerMustRun(u,i) * StartUp(u,i)
         + RampShutDownMaximumH(u,i) * TimeStep * ShutDown(u,i)
         + LL_RampDown(u,i)
;

* Start up cost
EQ_CostStartUp(u,i)$(CostStartUp(u) <> 0)..
         CostStartUpH(u,i)
         =E=
         CostStartUp(u)*StartUp(u,i)
;

* Shut down cost
EQ_CostShutDown(u,i)$(CostShutDown(u) <> 0)..
         CostShutDownH(u,i)
         =E=
         CostShutDown(u)*ShutDown(u,i)
;

EQ_CostRampUp(u,i)$(CostRampUp(u) <> 0)..
         CostRampUpH(u,i)
         =G=
         CostRampUp(u)*(Power(u,i)-PowerInitial(u)$(ord(i) = 1)-Power(u,i-1)$(ord(i) > 1))
;

EQ_CostRampDown(u,i)$(CostRampDown(u) <> 0)..
         CostRampDownH(u,i)
         =G=
         CostRampDown(u)*(PowerInitial(u)$(ord(i) = 1)+Power(u,i-1)$(ord(i) > 1)-Power(u,i))
;

EQ_Residual_Load(n,i)..
        ResidualLoad(n,i)
        =E=
        Demand("DA",n,i)
        + Demand("Flex",n,i)
        + DemandModulation(n,i)
        + sum(au, PowerConsumption(au,i) * Location(au,n))
        - sum(u,Power(u,i)$(sum(tr,Technology(u,tr))>=1) * Location(u,n))
        - sum(l,Flow(l,i)*LineNode(l,n))
;

*Hourly demand balance in the day-ahead market for each node
EQ_Demand_balance_DA(n,i)..
         sum(u,Power(u,i)*Location(u,n))
         +sum(p2bs,Power(p2bs,i)*Location(p2bs,n))
         +sum(l,Flow(l,i)*LineNode(l,n))
         +ShedLoad(n,i)
         +LL_MaxPower(n,i)
         =E=
         Demand("DA",n,i) + Demand("Flex",n,i)
         +DemandModulation(n,i)
         +sum(s,StorageInput(s,i)*Location(s,n))
         +sum(p2h,PowerConsumption(p2h,i)*Location(p2h,n))
         +sum(p2bs,PowerConsumption(p2bs,i)*Location(p2bs,n))
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
         MaxFlexDemand(n) - Demand("flex",n,i)
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
*         + sum((p2h),Reserve_2U(p2h,i)*Reserve(p2h)*Location(p2h,n))
         + sum((chp),Reserve_2U(chp,i)*Reserve(chp)*Location(chp,n))
         + CurtailmentReserve_2U(n,i) + LL_2U(n,i)
         =E=
         +(Demand("2U",n,i) + max(smax((u,tc),PowerCapacity(u)*Technology(u,tc)*LoadMaximum(u,i)*Location(u,n)), smax(l,FlowMaximum(l,i)*LineNode(l,n))$(card(l)>0)))*(1-K_QuickStart(n))
;

*Hourly demand balance in the upwards non-spinning reserve market for each node
EQ_Demand_balance_3U(n,i)..
         sum((u),(Reserve_2U(u,i) + Reserve_3U(u,i))*Reserve(u)*Location(u,n))
*         + sum((p2h),(Reserve_2U(p2h,i) + Reserve_3U(p2h,i))*Reserve(p2h)*Location(p2h,n))
         + sum((chp),(Reserve_2U(chp,i) + Reserve_3U(chp,i))*Reserve(chp)*Location(chp,n))
         + CurtailmentReserve_2U(n,i) + CurtailmentReserve_3U(n,i) + LL_3U(n,i)
         =E=
         Demand("2U",n,i)
         + max(smax((u,tc),PowerCapacity(u)*Technology(u,tc)*LoadMaximum(u,i)*Location(u,n)), smax(l,FlowMaximum(l,i)*LineNode(l,n))$(card(l)>0))
;

*Hourly demand balance in the downwards reserve market for each node
EQ_Demand_balance_2D(n,i)..
         sum((u),Reserve_2D(u,i)*Reserve(u)*Location(u,n))
*         + sum((p2h),Reserve_2D(p2h,i)*Reserve(p2h)*Location(p2h,n))
         + sum((chp),Reserve_2D(chp,i)*Reserve(chp)*Location(chp,n))
         + LL_2D(n,i)
         =E=
         Demand("2D",n,i)
         + max(smax(s,StorageChargingCapacity(s)*Location(s,n)), smax(l,-FlowMaximum(l,i)*LineNode(l,n))$(card(l)>0))
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

EQ_p2h_Reserve_2U_capability(p2h,i)$(StorageCapacity(p2h)>=srp*Efficiency(p2h,i)*PowerCapacity(p2h))..
         Reserve_2U(p2h,i)
         =l=
         PowerConsumption(p2h,i)
;

EQ_p2h_Reserve_2D_capability(p2h,i)$(StorageCapacity(p2h)>=srp*Efficiency(p2h,i)*PowerCapacity(p2h))..
         Reserve_2D(p2h,i)
         =l=
         PowerCapacity(p2h)*Nunits(p2h)*LoadMaximum(p2h,i)-PowerConsumption(p2h,i)
;

EQ_p2h_Reserve_3U_capability(p2h,i)$(QuickStartPower(p2h,i) > 0 and StorageCapacity(p2h)>=srp*Efficiency(p2h,i)*PowerCapacity(p2h))..
         Reserve_3U(p2h,i)
         =L=
         Nunits(p2h)*QuickStartPower(p2h,i)*TimeStep-PowerConsumption(p2h,i)
;

EQ_2U_limit_p2h(p2h,i)..
         Reserve_2U(p2h,i)
         =L=
         StorageLevel(p2h,i)/(0.25*Efficiency(p2h,i))

;

EQ_2D_limit_p2h(p2h,i)..
         Reserve_2D(p2h,i)
         =l=
         (StorageCapacity(p2h)*Nunits(p2h)-StorageLevel(p2h,i))/(0.25*Efficiency(p2h,i))

;

EQ_3U_limit_p2h(p2h,i)..
         Reserve_3U(p2h,i)
         =l=
         StorageLevel(p2h,i)/(Efficiency(p2h,i))

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
*Minimum power output is above the must-run output level for each unit in all periods
EQ_Power_must_run(u,i)..
         PowerMustRun(u,i) * Committed(u,i)
         - (StorageInput(u,i) * CHPPowerLossFactor(u) )$(chp(u) and (CHPType(u,'Extraction') or CHPType(u,'P2H')))
         =L=
         Power(u,i)
;

*Maximum power output is below the available capacity
EQ_Power_available(au,i)..
         Power(au,i)$(u(au))
         + Power(au,i)$(p2bs(au))
         + Heat(au,i)$(hu(au))
         + Heat(au,i)$(thms(au))
         + Power(au,i)$(p2h(au))
         =L=
         PowerCapacity(au)$(u(au))*LoadMaximum(au,i)$(u(au))*Committed(au,i)$(u(au))
         + PowerCapacity(au)$(p2bs(au))*LoadMaximum(au,i)$(p2bs(au))*Nunits(au)$(p2bs(au))
         + PowerCapacity(au)$(hu(au))*LoadMaximum(au,i)$(hu(au))*Committed(au,i)$(hu(au))
         + PowerCapacity(au)$(thms(au))*LoadMaximum(au,i)$(thms(au))*Committed(au,i)$(thms(au))
		 + 0
;

* Maximum boundary sector output is below the available capacity
EQ_Boundary_sector_only_power_available(n_bs,bsu,i)..
         PowerBoundarySector(n_bs,bsu,i)
         =L=
         PowerCapacity(bsu)*LoadMaximum(bsu,i)*Nunits(bsu)*Location_bs(bsu,n_bs)
;

EQ_Boundary_sector_only_power_available_min(n_bs,bsu,i)..
         0
         =L=
         PowerBoundarySector(n_bs,bsu,i)
;

*Boundary Sector Storage level must be above a minimum
EQ_Boundary_Sector_Storage_minimum(n_bs,i)..
         BoundarySectorStorageMinimum(n_bs)
         =L=
         BoundarySectorStorageLevel(n_bs,i)
;

*Storage level must be below storage capacity
EQ_Boundary_Sector_Storage_level(n_bs,i)..
         BoundarySectorStorageLevel(n_bs,i)
         =L=
         BoundarySectorStorageCapacity(n_bs)
;

*Discharge is limited by the storage level
EQ_Boundary_Sector_Storage_MaxDischarge(n_bs,i)..
         - BoundarySectorStorageInitial(n_bs)$(ord(i) = 1)
         - BoundarySectorStorageLevel(n_bs,i-1)$(ord(i) > 1)
         =L=
         BoundarySectorStorageInput(n_bs,i)*TimeStep
;

EQ_Boundary_Sector_Storage_MaxCharge(n_bs,i)$(BoundarySectorStorageCapacity(n_bs)>0)..
        BoundarySectorStorageInput(n_bs,i)*TimeStep
        =L=
        (BoundarySectorStorageCapacity(n_bs) - BoundarySectorStorageInitial(n_bs))$(ord(i) = 1)
        + (BoundarySectorStorageCapacity(n_bs) - BoundarySectorStorageLevel(n_bs,i-1))$(ord(i) > 1)
;


*Storage balance
EQ_Boundary_Sector_Storage_balance(n_bs,i)..
         BoundarySectorStorageInitial(n_bs)$(ord(i) = 1)
         + BoundarySectorStorageLevel(n_bs,i-1)$(ord(i) > 1)
         + BoundarySectorStorageInput(n_bs,i)*TimeStep
         =E=
         BoundarySectorStorageLevel(n_bs,i)
         + BoundarySectorStorageSelfDischarge(n_bs)*BoundarySectorStorageLevel(n_bs,i)*TimeStep
;

* Minimum level at the end of the optimization horizon:
EQ_Boundary_Sector_Storage_boundaries(n_bs,i)$(ord(i) = card(i))..
         BoundarySectorStorageFinalMin(n_bs)
         =L=
         BoundarySectorStorageLevel(n_bs,i) + BoundarySectorWaterSlack(n_bs)
;

*Storage level must be above a minimum
EQ_Storage_minimum(au,i)..
         StorageMinimum(au)$(s(au))*Nunits(au)$(s(au))
         =L=
         StorageLevel(au,i)$(s(au))
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

* The system could curtail by pumping and turbining at the same time if Nunits>1. This should be included into the curtailment equation!

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
         +spillage(au,i)$(s(au))
         +Power(au,i)$(s(au))*TimeStep/(max(StorageDischargeEfficiency(au)$(s(au)),0.0001))
         +StorageSelfDischarge(au)$(s(au))*StorageLevel(au,i)$(s(au))*TimeStep
;

* Minimum level at the end of the optimization horizon:
EQ_Storage_boundaries(au,i)$(ord(i) = card(i))..
         StorageFinalMin(au)$(s(au))
         =L=
         StorageLevel(au,i)$(s(au))
         + WaterSlack(au)$(s(au))
;

*Total emissions are capped
EQ_Emission_limits(i,p,n)..
         sum(u,Power(u,i)*EmissionRate(u,p)*TimeStep*Location(u,n))
         =L=
         EmissionMaximum(n,p)*TimeStep
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
EQ_BS_Flow_limits_lower(l_bs,i)..
         FlowBSMinimum(l_bs,i)
         =L=
         FlowBS(l_bs,i)
;

*Boundary Sector Flows are below maximum values
EQ_BS_Flow_limits_upper(l_bs,i)..
         FlowBS(l_bs,i)
         =L=
         FlowBSMaximum(l_bs,i)
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
;

EQ_CHP_extraction_Pmax(chp,i)$(CHPType(chp,'Extraction') or CHPType(chp,'P2H'))..
         Power(chp,i)
         =L=
         PowerCapacity(chp)*Nunits(chp)  - StorageInput(chp,i) * CHPPowerLossFactor(chp)
;

EQ_CHP_backpressure(chp,i)$(CHPType(chp,'Back-Pressure'))..
         Power(chp,i)
         =E=
         StorageInput(chp,i) * CHPPowerToHeat(chp)
;

EQ_CHP_max_heat(chp,i)..
         StorageInput(chp,i)
         =L=
         CHPMaxHeat(chp)*Nunits(chp)
;

* Power to X units
EQ_P2X_Power_Balance(au,i)..
         StorageInput(au,i)$(p2h(au))
         =E=
         PowerConsumption(au,i)$(p2h(au)) * Efficiency(au,i)$(p2h(au))
;

EQ_Max_Power_Consumption(au,i)..
         PowerConsumption(au,i)$(p2h(au))
         =L=
*         PowerCapacity(au)$(p2h(au)) * Committed(au,i)$(p2h(au))
         PowerCapacity(au)$(p2h(au)) * Nunits(au)$(p2h(au))
;

* Power to boundary sector units
EQ_Power_Balance_of_BS_units(n_bs,p2bs,i)..
         PowerBoundarySector(n_bs,p2bs,i)
         =E=
         (PowerConsumption(p2bs,i) * ChargingEfficiencyBoundarySector(n_bs,p2bs,i) * Location_bs(p2bs,n_bs))$(ChargingEfficiencyBoundarySector(n_bs,p2bs,i) <> 0)
         + (Power(p2bs,i) / (EfficiencyBoundarySector(n_bs,p2bs,i)+0.0001) * Location_bs(p2bs,n_bs))$(EfficiencyBoundarySector(n_bs,p2bs,i) <> 0)
;

EQ_Max_Power_Consumption_of_BS_units(p2bs,i)..
         PowerConsumption(p2bs,i)
         =L=
         StorageChargingCapacity(p2bs) * Nunits(p2bs) *TimeStep
;


EQ_Heat_Demand_balance(n_th,i)..
         sum(chp, Heat(chp,i)*Location_th(chp,n_th))
         + sum(p2h, Heat(p2h,i)*Location_th(p2h,n_th))
         + sum(hu, Heat(hu,i)*Location_th(hu,n_th))
         + sum(thms, Heat(thms,i)*Location_th(thms,n_th))
         =E=
         HeatDemand(n_th, i)
         - HeatSlack(n_th,i)
         + sum(thms, StorageInput(thms,i)*Location_th(thms,n_th))
;


EQ_BS_Demand_balance(n_bs,i)..
        sum(p2bs, PowerBoundarySector(n_bs,p2bs,i))
        + sum(bsu, PowerBoundarySector(n_bs,bsu,i))
        + BoundarySectorSlack(n_bs,i)
        + sum(l_bs,FlowBS(l_bs,i)*BSLineNode(l_bs,n_bs))
        + BSFlexSupply(n_bs,i)
        =E=
        BoundarySectorDemand(n_bs,i)
        + BSFlexDemand(n_bs,i)
        + BoundarySectorStorageInput(n_bs,i)
        + BoundarySectorSpillage(n_bs,i)
;

*Heat Storage balance
EQ_Heat_Storage_balance(th,i)..
         StorageInitial(th)$(ord(i) = 1)
         + StorageLevel(th,i-1)$(ord(i) > 1)
         + StorageInput(th,i)*TimeStep
         =E=
         StorageLevel(th,i)
         + Heat(th,i)*TimeStep
         + StorageSelfDischarge(th)*StorageLevel(th,i)*TimeStep
;
* The self-discharge proportional to the charging level is a bold hypothesis, but it avoids keeping self-discharging if the level reaches zero

* Heat Storage level must be above a minimum
EQ_Heat_Storage_minimum(th,i)..
         StorageMinimum(th)
         =L=
         StorageLevel(th,i)
;

*Storage level must be below a maximum (to make the reserves available)
EQ_Heat_Storage_level(th,i)..
         StorageLevel(th,i)
         =L=
         StorageCapacity(th)*AvailabilityFactor(th,i)*Nunits(th)
;

* Heat Storage charging is bounded by the maximum capacity
EQ_Heat_Storage_input(thms,i)..
         StorageInput(thms,i)
         =L=
         StorageChargingCapacity(thms)*(Nunits(thms)-Committed(thms,i))
;

* Heat storage discharge is limited by the remaining storage capacity
EQ_Heat_Storage_MaxDischarge(thms,i)$(StorageCapacity(thms)>PowerCapacity(thms)*TimeStep)..
         Heat(thms,i)*TimeStep/(max(StorageDischargeEfficiency(thms),0.0001))
         =L=
         StorageInitial(thms)$(ord(i) = 1)
         + StorageLevel(thms,i-1)$(ord(i) > 1)
         + StorageInflow(thms,i)*Nunits(thms)*TimeStep
;

* Heat storage charging is limited by the remaining storage capacity
EQ_Heat_Storage_MaxCharge(thms,i)$(StorageCapacity(thms)>PowerCapacity(thms)*TimeStep)..
         StorageInput(thms,i)*StorageChargingEfficiency(thms)*TimeStep
         =L=
         (Nunits(thms) * StorageCapacity(thms)-StorageInitial(thms))$(ord(i) = 1)
         + (Nunits(thms) * StorageCapacity(thms)*AvailabilityFactor(thms,i-1)
         - StorageLevel(thms,i-1))$(ord(i) > 1)
         + StorageOutflow(thms,i)*Nunits(thms)*TimeStep
;

* Minimum heat storage level at the end of the optimization horizon:
EQ_Heat_Storage_boundaries(thms,i)$(ord(i) = card(i))..
         StorageFinalMin(thms)
         =L=
         StorageLevel(thms,i)
         + WaterSlack(thms)
;

* Equations concerning boundary sector (PtL) flexible demand
* Only valid in old formulation!!! If MTS=1: assures that total demand of boundary sector (PtL) is fulfilled
EQ_Tot_Flex_Demand_BS(n_bs)..
         BSFlexDemandInputInitial(n_bs)
         + sum(i,BSFlexDemandInput(n_bs,i))
         =E=
         sum(i,BSFlexDemand(n_bs,i))
;

* Capacity of boundary sector (PtL) must not be exceeded
EQ_Max_Flex_Capacity_BS(n_bs,i)..
         BSFlexDemand(n_bs,i)
         =L=
         BSFlexMaxCapacity(n_bs)*TimeStep
;

* Only valid in old formulation!!! If MTS = 0: Boundary sector (PtL) demand is not a variable anymore, but a parameter
EQ_BS_Flex_Demand(n_bs,i)..
         BSFlexDemand(n_bs,i)
         =E=
         BSFlexDemandInput(n_bs,i)
;


EQ_Tot_Flex_Supply_BS(n_bs)..
         BSFlexSupplyInputInitial(n_bs)
         + sum(i,BSFlexSupplyInput(n_bs,i))
         =G=
         sum(i,BSFlexSupply(n_bs,i))
;

* Capacity of boundary sector (PtL) must not be exceeded
EQ_Max_Flex_Supply_BS(n_bs,i)..
         BSFlexSupply(n_bs,i)
         =L=
         BSFlexMaxSupply(n_bs)*TimeStep
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
EQ_Heat_Demand_balance,
EQ_BS_Demand_balance,
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
EQ_Power_Balance_of_BS_units,
EQ_Max_Power_Consumption_of_BS_units,
EQ_Boundary_sector_only_power_available
EQ_Boundary_sector_only_power_available_min
EQ_Power_available,
EQ_Heat_Storage_balance,
EQ_Heat_Storage_minimum,
EQ_Heat_Storage_level,
EQ_Heat_Storage_input,
EQ_Heat_Storage_MaxDischarge,
EQ_Heat_Storage_MaxCharge,
EQ_Heat_Storage_boundaries,
EQ_Reserve_2U_capability,
EQ_Reserve_2D_capability,
EQ_Reserve_3U_capability,
EQ_p2h_Reserve_2U_capability,
EQ_p2h_Reserve_2D_capability,
EQ_p2h_Reserve_3U_capability,
EQ_2U_limit_p2h,
EQ_2D_limit_p2h,
EQ_3U_limit_p2h,
EQ_2U_limit_chp,
EQ_2D_limit_chp,
EQ_3U_limit_chp,
EQ_Storage_minimum,
EQ_Storage_level,
EQ_Storage_input,
EQ_Storage_balance,
EQ_Storage_boundaries,
EQ_Boundary_Sector_Storage_MaxDischarge,
EQ_Boundary_Sector_Storage_MaxCharge,
EQ_Boundary_Sector_Storage_minimum,
EQ_Boundary_Sector_Storage_level,
EQ_Boundary_Sector_Storage_balance,
EQ_Boundary_Sector_Storage_boundaries,
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
*$If %MTS% == 1 EQ_Tot_Flex_Demand_BS,
*$If %MTS% == 1 EQ_Max_Flex_Capacity_BS,
*$If %MTS% == 0 EQ_BS_Flex_Demand,
EQ_Tot_Flex_Demand_BS,
EQ_Max_Flex_Capacity_BS,
EQ_Tot_Flex_Supply_BS,
EQ_Max_Flex_Supply_BS,
$If %RetrieveStatus% == 1 EQ_CommittedCalc
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

* Some parameters used for debugging:
failed=0;
parameter CommittedInitial_dbg(au), PowerInitial_dbg(u), StorageInitial_dbg(au), StorageFinalMin_dbg(au), AccumulatedOverSupply_inital_dbg(n);

* Fixing the initial guesses:
*PowerH.L(u,i)=PowerInitial(u);
*Committed.L(u,i)=CommittedInitial(u);

* Defining a parameter that records the solver status:
set  tmp   "tpm"  / "model", "solver"/  ;
parameter status(tmp,h);

$if %Debug% == 1 $goto DebugSection

display "OK";

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
         StorageFinalMin(thms) =  sum(i$(ord(i)=card(i)),StorageProfile(thms,i)*StorageCapacity(thms)*Nunits(thms)*AvailabilityFactor(thms,i));
		 StorageFinalMin(chp) =  sum(i$(ord(i)=card(i)),StorageProfile(chp,i)*StorageCapacity(chp)*Nunits(chp)*AvailabilityFactor(chp,i));
         BoundarySectorStorageFinalMin(n_bs) = sum(i$(ord(i)=card(i)),BoundarySectorStorageProfile(n_bs,i)*BoundarySectorStorageCapacity(n_bs));

$If %Verbose% == 1   Display PowerInitial,CommittedInitial,StorageFinalMin;
$If %Verbose% == 1   Display PowerInitial,StorageFinalMin;

$If %LPFormulation% == 1          SOLVE UCM_SIMPLE USING LP MINIMIZING SystemCostD;
$If not %LPFormulation% == 1      SOLVE UCM_SIMPLE USING MIP MINIMIZING SystemCostD;

$If %Verbose% == 0 $goto skipdisplay3
Display EQ_Objective_function.M, EQ_CostRampUp.M, EQ_CostRampDown.M, EQ_Demand_balance_DA.M, EQ_Storage_minimum.M, EQ_Storage_level.M, EQ_Storage_input.M, EQ_Storage_balance.M, EQ_Storage_boundaries.M,  EQ_Flow_limits_lower.M ;
$label skipdisplay3

         status("model",i) = UCM_SIMPLE.Modelstat;
         status("solver",i) = UCM_SIMPLE.Solvestat;

*if(UCM_SIMPLE.Modelstat <> 1 and UCM_SIMPLE.Modelstat <> 8 and not failed, CommittedInitial_dbg(au) = CommittedInitial(au); PowerInitial_dbg(u) = PowerInitial(u); StorageInitial_dbg(au) = StorageInitial(au);
*                                                                           EXECUTE_UNLOAD "debug.gdx" day, status, CommittedInitial_dbg, PowerInitial_dbg, StorageInitial_dbg;
*                                                                           failed=1;);
         CommittedInitial_dbg(au) = CommittedInitial(au);
         PowerInitial_dbg(u) = PowerInitial(u);
         StorageInitial_dbg(au) = StorageInitial(au);
         StorageFinalMin_dbg(au) = StorageFinalMin(au);
$If %ActivateFlexibleDemand% == 1 AccumulatedOverSupply_inital_dbg(n) = AccumulatedOverSupply_inital(n);
         EXECUTE_UNLOAD "debug.gdx" day, status, CommittedInitial_dbg, PowerInitial_dbg, StorageInitial_dbg, StorageFinalMin_dbg,
$If %ActivateFlexibleDemand% == 1 AccumulatedOverSupply_inital_dbg;

         CommittedInitial(au) = sum(i$(ord(i)=LastKeptHour-FirstHour+1),Committed.L(au,i));
         PowerInitial(u) = sum(i$(ord(i)=LastKeptHour-FirstHour+1),Power.L(u,i));
         StorageInitial(s) =   sum(i$(ord(i)=LastKeptHour-FirstHour+1),StorageLevel.L(s,i));
         StorageInitial(thms) =   sum(i$(ord(i)=LastKeptHour-FirstHour+1),StorageLevel.L(thms,i));
         StorageInitial(chp) =   sum(i$(ord(i)=LastKeptHour-FirstHour+1),StorageLevel.L(chp,i));
         BoundarySectorStorageInitial(n_bs) = sum(i$(ord(i)=LastKeptHour-FirstHour+1),BoundarySectorStorageLevel.L(n_bs,i));
         BSFlexDemandInputInitial(n_bs) = sum(i$(ord(i)<LastKeptHour+1), BSFlexDemandInput(n_bs,i) - BSFlexDemand.L(n_bs,i));
         BSFlexSupplyInputInitial(n_bs) = sum(i$(ord(i)<LastKeptHour+1), BSFlexSupplyInput(n_bs,i) - BSFlexSupply.L(n_bs,i));
$If %ActivateFlexibleDemand% == 1 AccumulatedOverSupply_inital(n) = sum(i$(ord(i)=LastKeptHour-FirstHour+1),AccumulatedOverSupply.L(n,i));

* Assigning waterslack (one value per optimization horizon) to the last element of storageslack
StorageSlack.L(au,i)$(ord(i)=LastKeptHour-FirstHour+1) = Waterslack.L(au);
BoundarySectorStorageSlack.L(n_bs,i)$(ord(i)=LastKeptHour-FirstHour+1) = BoundarySectorWaterslack.L(n_bs);
ObjectiveFunction.L(i)$(ord(i)=LastKeptHour-FirstHour+1) = SystemCostD.L;
Error.L = sum((i,n), CostLoadShedding(n,i)*ShedLoad.L(n,i)
                               +Config("ValueOfLostLoad","val")*(LL_MaxPower.L(n,i)+LL_MinPower.L(n,i))
                               +0.8*Config("ValueOfLostLoad","val")*(LL_2U.L(n,i)+LL_2D.L(n,i)+LL_3U.L(n,i)))
                      +sum((u,i), 0.7*Config("ValueOfLostLoad","val")*(LL_RampUp.L(u,i)+LL_RampDown.L(u,i)))
                      +sum((i,n_th), CostHeatSlack(n_th,i) * HeatSlack.L(n_th,i))
;
OptimalityGap.L(i)$(ord(i)=LastKeptHour-FirstHour+1) = UCM_SIMPLE.objVal - UCM_SIMPLE.objEst;
OptimizationError.L(i)$(ord(i)=LastKeptHour-FirstHour+1) = Error.L - OptimalityGap.L(i);


*Loop variables to display after solving:
$If %Verbose% == 1 Display LastKeptHour,PowerInitial,StorageInitial;
);

*CurtailedPower.L(n,z)=sum(u,(Nunits(u)*PowerCapacity(u)*LoadMaximum(u,z)-Power.L(u,z))$(sum(tr,Technology(u,tr))>=1) * Location(u,n)) + sum(s,spillage.L(s,z)* Location(s,n));
CurtailedHeat.L(n_th,z)=sum(hu,(Nunits(hu)*PowerCapacity(hu)*LoadMaximum(hu,z)-Heat.L(hu,z))$(sum(tr,Technology(hu,tr))>=1) * Location_th(hu,n_th));


$If %Verbose% == 1 Display Flow.L,Power.L,Committed.L,ShedLoad.L,CurtailedPower.L,CurtailmentReserve_2U.L, CurtailmentReserve_3U.L,CurtailedHeat.L,StorageLevel.L,StorageInput.L,SystemCost.L,LL_MaxPower.L,LL_MinPower.L,LL_2U.L,LL_2D.L,LL_RampUp.L,LL_RampDown.L;
$If %Verbose% == 1 Display Flow.L,Power.L,ShedLoad.L,CurtailedPower.L,CurtailmentReserve.L,CurtailmentReserve_3U.L,CurtailedHeat.L,StorageLevel.L,StorageInput.L,SystemCost.L,LL_MaxPower.L,LL_MinPower.L,LL_2U.L,LL_2D.L;

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
OutputFlowBoundarySector(l_bs,h)
OutputPower(au,h)
OutputPowerBoundarySector(n_bs,au,h)
OutputPowerConsumption(au,h)
OutputResidualLoad(n,h)
OutputStorageInput(au,h)
OutputStorageLevel(au,h)
OutputStorageSlack(au,h)
OutputBoundarySectorStorageLevel(n_bs,h)
OutputBoundarySectorStorageShadowPrice(n_bs,h)
OutputBoundarySectorStorageSlack(n_bs,h)
OutputBoundarySectorStorageInput(n_bs,h)
OutputBoundarySectorSpillage(n_bs,h)
OutputSystemCost(h)
OutputSpillage(au,h)
OutputShedLoad(n,h)
OutputCurtailedPower(n,h)
OutputCurtailmentReserve_2U(n,h)
OutputCurtailmentReserve_3U(n,h)
OutputCurtailedHeat(n_th,h)
$If %ActivateFlexibleDemand% == 1 OutputDemandModulation(n,h)
ShadowPrice(n,h)
HeatShadowPrice(n_th,h)
BoundarySectorShadowPrice(n_bs,h)
LostLoad_MaxPower(n,h)
LostLoad_MinPower(n,h)
LostLoad_2D(n,h)
LostLoad_2U(n,h)
LostLoad_3U(n,h)
$If %MTS%==0 LostLoad_RampUp(n,h)
$If %MTS%==0 LostLoad_RampDown(n,h)
OutputGenMargin(n,h)
OutputHeat(au,h)
OutputHeatSlack(n_th,h)
OutputBoundarySectorSlack(n_bs,h)
LostLoad_WaterSlack(au)
LostLoad_BoundarySectorWaterSlack(n_bs)
StorageShadowPrice(au,h)
OutputBSFlexDemand(n_bs,h)
OutputBSFlexSupply(n_bs,h)
OutputPowerMustRun(u,h)
$If %MTS%==0 OutputCostStartUpH(u,h)
$If %MTS%==0 OutputCostShutDownH(u,h)
$If %MTS%==0 OutputCostRampUpH(u,h)
$If %MTS%==0 OutputCostRampDownH(u,h)
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
;

OutputMaxOutageUp(n,z)=max(smax((au,tc),PowerCapacity(au)*Technology(au,tc)*LoadMaximum(au,z)*Location(au,n)), smax(l,FlowMaximum(l,z)*LineNode(l,n))$(card(l)>0));
OutputMaxOutageDown(n,z)=max(smax(s,StorageChargingCapacity(s)*Location(s,n)), smax(l,-FlowMaximum(l,z)*LineNode(l,n))$(card(l)>0));
OutputDemand_2U(n,z)=(Demand("2U",n,z) + OutputMaxOutageUp(n,z))*(1-K_QuickStart(n));
OutputDemand_3U(n,z)=Demand("2U",n,z) + OutputMaxOutageUp(n,z);
OutputDemand_2D(n,z)=Demand("2D",n,z) + OutputMaxOutageDown(n,z);
OutputCommitted(au,z)=Committed.L(au,z);
OutputFlow(l,z)=Flow.L(l,z);
OutputFlowBoundarySector(l_bs,z)=FlowBS.L(l_bs,z);
OutputPower(au,z)=Power.L(au,z);
OutputPowerBoundarySector(n_bs,au,z)=PowerBoundarySector.L(n_bs,au,z);
OutputPowerConsumption(au,z)=PowerConsumption.L(au,z);
OutputResidualLoad(n,z)=ResidualLoad.L(n,z);
OutputHeat(au,z)=Heat.L(au,z);
OutputHeatSlack(n_th,z)=HeatSlack.L(n_th,z);
OutputBoundarySectorSlack(n_bs,z) = BoundarySectorSlack.L(n_bs,z);
OutputStorageInput(s,z)=StorageInput.L(s,z);
OutputStorageInput(th,z)=StorageInput.L(th,z);
OutputStorageLevel(s,z)=StorageLevel.L(s,z)/max(1,StorageCapacity(s)*Nunits(s)*AvailabilityFactor(s,z));
OutputStorageLevel(th,z)=StorageLevel.L(th,z)/max(1,StorageCapacity(th)*Nunits(th));
OutputStorageSlack(au,z) = StorageSlack.L(au,z);
OutputBoundarySectorStorageLevel(n_bs,z) = BoundarySectorStorageLevel.L(n_bs,z);
OutputBoundarySectorStorageShadowPrice(n_bs,z) = EQ_Boundary_Sector_Storage_balance.m(n_bs,z);
OutputBoundarySectorStorageSlack(n_bs,z) = BoundarySectorStorageSlack.l(n_bs,z);
OutputBoundarySectorStorageInput(n_bs,z) = BoundarySectorStorageInput.l(n_bs,z);
OutputBoundarySectorSpillage(n_bs,z) = BoundarySectorSpillage.l(n_bs,z);
OutputSystemCost(z)=SystemCost.L(z);
OutputSpillage(au,z)  = Spillage.L(au,z) ;
OutputShedLoad(n,z) = ShedLoad.L(n,z);
OutputCurtailedPower(n,z)=CurtailedPower.L(n,z);
OutputCurtailmentReserve_2U(n,z)=CurtailmentReserve_2U.L(n,z);
OutputCurtailmentReserve_3U(n,z)=CurtailmentReserve_3U.L(n,z);
OutputCurtailedHeat(n_th,z)=CurtailedHeat.L(n_th,z);
$If %ActivateFlexibleDemand% == 1 OutputDemandModulation(n,z)=DemandModulation.L(n,z);
LostLoad_MaxPower(n,z)  = LL_MaxPower.L(n,z);
LostLoad_MinPower(n,z)  = LL_MinPower.L(n,z);
LostLoad_2D(n,z) = LL_2D.L(n,z);
LostLoad_2U(n,z) = LL_2U.L(n,z);
LostLoad_3U(n,z) = LL_3U.L(n,z);
$If %MTS%==0 LostLoad_RampUp(n,z)    = sum(u,LL_RampUp.L(u,z)*Location(u,n));
$If %MTS%==0 LostLoad_RampDown(n,z)  = sum(u,LL_RampDown.L(u,z)*Location(u,n));
ShadowPrice(n,z) = EQ_Demand_balance_DA.m(n,z);
HeatShadowPrice(n_th,z) = EQ_Heat_Demand_balance.m(n_th,z);
BoundarySectorShadowPrice(n_bs,z) = EQ_BS_Demand_balance.m(n_bs,z);
LostLoad_WaterSlack(au) = WaterSlack.L(au);
LostLoad_BoundarySectorWaterSlack(n_bs) = BoundarySectorWaterSlack.L(n_bs);
StorageShadowPrice(s,z) = 0 ;
OutputBSFlexDemand(n_bs,z) = BSFlexDemand.L(n_bs,z);
OutputBSFlexSupply(n_bs,z) = BSFlexSupply.L(n_bs,z);
StorageShadowPrice(s,z) = EQ_Storage_balance.m(s,z);
StorageShadowPrice(th,z) = EQ_Heat_Storage_balance.m(th,z);
OutputPowerMustRun(u,z) = PowerMustRun(u,z);
$If (%MTS%==0 or %LPFormulation% == 1) OutputCostStartUpH(u,z) = CostStartUpH.L(u,z);
$If (%MTS%==0 or %LPFormulation% == 1) OutputCostShutDownH(u,z) = CostShutDownH.L(u,z);
$If (%MTS%==0 or %LPFormulation% == 1) OutputCostRampUpH(u,z) = CostRampUpH.L(u,z);
$If (%MTS%==0 or %LPFormulation% == 1) OutputCostRampDownH(u,z) = CostRampDownH.L(u,z);

ShadowPrice_2U(n,z) =  EQ_Demand_balance_2U.m(n,z);
ShadowPrice_2D(n,z) =  EQ_Demand_balance_2D.m(n,z);
ShadowPrice_3U(n,z) =  EQ_Demand_balance_3U.m(n,z);

OutputReserve_2U(au,z) = Reserve_2U.L(au,z);
OutputReserve_2D(au,z) = Reserve_2D.L(au,z);
OutputReserve_3U(au,z) = Reserve_3U.L(au,z);

ShadowPrice_RampUp_TC(u,z) = EQ_RampUp_TC.m(u,z);
ShadowPrice_RampDown_TC(u,z) = EQ_RampDown_TC.m(u,z);
OutputRampRate(u,z) = - Power.L(u,z-1)$(ord(z) > 1) - PowerInitial(u)$(ord(z) = 1) + Power.L(u,z);
OutputRampRate(hu,z) = - Heat.L(hu,z-1)$(ord(z) > 1) + Heat.L(hu,z);
OutputStartUp(au,z) = StartUp.L(au,z);
OutputShutDown(au,z) = ShutDown.L(au,z);

OutputEmissions(n,p,z) = (sum(u,Power.L(u,z)*EmissionRate(u,p)*Location(u,n))
                        + sum(hu,Heat.L(hu,z)*EmissionRate(hu,p)*Location(hu,n)))
                        / (sum(u,Power.L(u,z)*Location(u,n)) + sum(hu,Heat.L(hu,z)*Location(hu,n)));

CapacityMargin(n,z) = (sum(u, Nunits(u)*PowerCapacity(u)$(not s(u))*LoadMaximum(u,z)*Location(u,n))
                      + min(sum(s, Nunits(s)*PowerCapacity(s)*LoadMaximum(s,z)*Location(s,n)), sum(s, StorageLevel.L(s,z)*StorageCapacity(s)))
                      + sum(l, Flow.L(l,z)*LineNode(l,n))
                      + CurtailedPower.L(n,z)
*                      + sum(l,(FlowMaximum(l,z)-Flow.L(l,z))*LineNode(l,n))$(card(l)>0)
                      - Demand("DA",n,z)
                      - DemandModulation.L(n,z)
                      - sum(p2h,PowerConsumption.L(p2h,z)*Location(p2h,n))
                      - sum(s, StorageInput.L(s,z)*Location(s,n))
                      - sum(au, (Reserve_2U.L(au,z) + Reserve_3U.L(au,z))*Location(au,n))
                      - CurtailmentReserve_2U.L(n,z)
                      - CurtailmentReserve_3U.L(n,z)
);
OutputSystemCostD(z) = ObjectiveFunction.L(z);
OutputOptimalityGap(z) = OptimalityGap.L(z);
OutputOptimizationError(z) = OptimizationError.L(z);
OutputOptimizationCheck(z) = OptimizationError.L(z) - OptimalityGap.L(z);

EXECUTE_UNLOAD "Results.gdx"
OutputCommitted,
OutputFlow,
OutputFlowBoundarySector,
OutputPower,
OutputPowerBoundarySector,
OutputPowerConsumption,
OutputResidualLoad,
OutputHeat,
OutputHeatSlack,
OutputBoundarySectorSlack,
OutputStorageInput,
OutputStorageLevel,
OutputStorageSlack,
OutputBoundarySectorStorageLevel,
OutputBoundarySectorStorageShadowPrice,
OutputBoundarySectorStorageSlack,
OutputBoundarySectorStorageInput,
OutputBoundarySectorSpillage,
OutputSystemCost,
OutputSpillage,
OutputShedLoad,
OutputCurtailedPower,
OutputCurtailmentReserve_2U,
OutputCurtailmentReserve_3U,
OutputCurtailedHeat,
$If %ActivateFlexibleDemand% == 1 OutputDemandModulation,
OutputGenMargin,
LostLoad_MaxPower,
LostLoad_MinPower,
LostLoad_2D,
LostLoad_2U,
LostLoad_3U,
$If %MTS%==0 LostLoad_RampUp,
$If %MTS%==0 LostLoad_RampDown,
ShadowPrice,
ShadowPrice_2U,
ShadowPrice_2D,
ShadowPrice_3U,
HeatShadowPrice,
LostLoad_WaterSlack,
LostLoad_BoundarySectorWaterSlack,
StorageShadowPrice,
OutputBSFlexDemand,
OutputBSFlexSupply,
BoundarySectorShadowPrice,
OutputPowerMustRun,
$If %MTS%==0 OutputCostStartUpH,
$If %MTS%==0 OutputCostShutDownH,
$If %MTS%==0 OutputCostRampUpH,
$If %MTS%==0 OutputCostRampDownH,
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
status
;

display OutputPowerConsumption, heat.L, heatslack.L, powerconsumption.L, power.L;

$onorder
* Exit here if the PrintResult option is set to 0:
$if not %PrintResults%==1 $exit

EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=Technology rng=Technology!A1 rdim=2 cdim=0'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=PowerCapacity rng=PowerCapacity!A1 rdim=1 cdim=0'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=PowerInitial rng=PowerInitialA1 rdim=1 cdim=0'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=RampDownMaximum rng=RampDownMaximum!A1 rdim=1 cdim=0'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=RampShutDownMaximum rng=RampShutDownMaximum!A1 rdim=1 cdim=0'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=RampStartUpMaximum rng=RampStartUpMaximum!A1 rdim=1 cdim=0'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=RampUpMaximum rng=RampUpMaximum!A1 rdim=1 cdim=0'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=TimeUpMinimum rng=TimeUpMinimum!A1 rdim=1 cdim=0'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=TimeDownMinimum rng=TimeDownMinimum!A1 rdim=1 cdim=0'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=Reserve rng=Reserve!A1 rdim=1 cdim=0'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=StorageCapacity rng=StorageCapacity!A1 rdim=1 cdim=0'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=Y par=StorageInflow rng=StorageInflow!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=StorageCapacity rng=StorageCapacity!A1 rdim=1 cdim=0'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=LoadShedding rng=LoadShedding!A1 rdim=1 cdim=0'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=FlowMaximum rng=FlowMaximum!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=AvailabilityFactor rng=AvailabilityFactor!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=Y par=OutageFactor rng=OutageFactor!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=Demand rng=Demand!A1 rdim=2 cdim=1'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=PartLoadMin rng=PartLoadMin!A1 rdim=1 cdim=0'

EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N var=CurtailedPower rng=CurtailedPower!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N var=ShedLoad rng=ShedLoad!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=OutputCommitted rng=Committed!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=OutputFlow rng=Flow!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=OutputPower rng=Power!A5 epsout=0 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=OutputPowerConsumption rng=Power!A5 epsout=0 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=OutputStorageInput rng=StorageInput!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=OutputStorageLevel rng=StorageLevel!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=OutputSystemCost rng=SystemCost!A1 rdim=1 cdim=0'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=Y var=LostLoad_MaxPower rng=LostLoad_MaxPower!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=Y var=LostLoad_MinPower rng=LostLoad_MinPower!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=Y var=LostLoad_2D rng=LostLoad_2D!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=Y var=LostLoad_2U rng=LostLoad_2U!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=Y var=OutputStorageSlack rng=OutputStorageSlack!A1 rdim=1 cdim=1'
$If %MTS%==0 EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=Y var=LostLoad_RampUp rng=LostLoad_RampUp!A1 rdim=1 cdim=1'
$If %MTS%==0 EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=Y var=LostLoad_RampDown rng=LostLoad_RampDown!A1 rdim=1 cdim=1'


$exit

$Label DebugSection

$gdxin debug.gdx
$LOAD day
$LOAD PowerInitial_dbg
$If %MTS% == 0 $LOAD CommittedInitial_dbg
$LOAD StorageInitial_dbg
$LOAD StorageFinalMin_dbg
$If %ActivateFlexibleDemand% == 1 $LOAD AccumulatedOverSupply_inital_dbg
;
PowerInitial(u) = PowerInitial_dbg(u);
$If %MTS%==0 CommittedInitial(au) = CommittedInitial_dbg(au);
StorageInitial(au) = StorageInitial_dbg(au);
$If %ActivateFlexibleDemand% == 1 AccumulatedOverSupply_inital(n) = AccumulatedOverSupply_inital_dbg(n);
FirstHour = (day-1)*24/TimeStep+1;
LastHour = min(card(h),FirstHour + (Config("RollingHorizon Length","day")+Config("RollingHorizon LookAhead","day")) * 24/TimeStep - 1);
LastKeptHour = LastHour - Config("RollingHorizon LookAhead","day") * 24/TimeStep;
i(h) = no;
i(h)$(ord(h)>=firsthour and ord(h)<=lasthour)=yes;
*StorageFinalMin(s) =  sum(i$(ord(i)=card(i)),StorageProfile(s,i)*StorageCapacity(s)*Nunits(s)*AvailabilityFactor(s,i));
StorageFinalMin(au) =  sum(i$(ord(i)=card(i)),StorageProfile(au,i)*StorageCapacity(au)*Nunits(au)*AvailabilityFactor(au,i));
*StorageFinalMin(au) = StorageFinalMin_dbg(au);

display day,FirstHour,LastHour,LastKeptHour;
Display PowerInitial,CommittedInitial,StorageInitial,StorageFinalMin;
$If %LPFormulation% == 1          SOLVE UCM_SIMPLE USING LP MINIMIZING SystemCostD;
$If not %LPFormulation% == 1      SOLVE UCM_SIMPLE USING MIP MINIMIZING SystemCostD;
$If %LPFormulation% == 1          Display EQ_Objective_function.M, EQ_CostRampUp.M, EQ_CostRampDown.M, EQ_Demand_balance_DA.M, EQ_Storage_minimum.M, EQ_Storage_level.M, EQ_Storage_input.M, EQ_Storage_balance.M, EQ_Storage_boundaries.M, EQ_Storage_MaxCharge.M, EQ_Storage_MaxDischarge.M, EQ_Flow_limits_lower.M ;
$If not %LPFormulation% == 1      Display EQ_Objective_function.M, EQ_CostStartUp.M, EQ_CostShutDown.M, EQ_Demand_balance_DA.M, EQ_Storage_minimum.M, EQ_Storage_level.M, EQ_Storage_input.M, EQ_Storage_balance.M, EQ_Storage_boundaries.M, EQ_Storage_MaxCharge.M, EQ_Storage_MaxDischarge.M, EQ_Flow_limits_lower.M ;

Display Flow.L,Power.L,Committed.L,ShedLoad.L,StorageLevel.L,StorageInput.L,SystemCost.L,Spillage.L,StorageLevel.L,StorageInput.L,LL_MaxPower.L,LL_MinPower.L,LL_2U.L,LL_2D.L,LL_RampUp.L,LL_RampDown.L,WaterSlack.L;
