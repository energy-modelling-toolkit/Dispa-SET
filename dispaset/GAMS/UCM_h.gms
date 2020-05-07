$Title UCM model

$eolcom //
Option threads=8;
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
l                Lines
au               All Units
u(au)           Generation units
t                Generation technologies
tr(t)            Renewable generation technologies
f                Fuel types
p                Pollutants
s(u)             Storage Units (with reservoir)
chp(u)           CHP units
p2h(au)          Power to heat units
th(au)           Units with thermal storage
h                Hours
i(h)             Subset of simulated hours for one iteration
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
AvailabilityFactor(au,h)          [%]      Availability factor
CHPPowerLossFactor(u)            [%]      Power loss when generating heat
CHPPowerToHeat(u)                [%]      Nominal power-to-heat factor
CHPMaxHeat(chp)                  [MW\u]     Maximum heat capacity of chp plant
CHPType
CommittedInitial(u)              [n.a.]   Initial committment status
Config
*CostCurtailment(n,h)             [EUR\MW]  Curtailment costs
CostFixed(u)                     [EUR\h]    Fixed costs
CostRampUp(u)                    [EUR\MW] Ramp-up costs
CostRampDown(u)                  [EUR\MW] Ramp-down costs
CostShutDown(u)                  [EUR\u]    Shut-down costs
CostStartUp(u)                   [EUR\u]    Start-up costs
CostVariable(u,h)                [EUR\MW]   Variable costs
CostHeatSlack(th,h)             [EUR\MWh]  Cost of supplying heat via other means
CostLoadShedding(n,h)            [EUR\MWh] Cost of load shedding
Curtailment(n)                   [n.a]    Curtailment allowed or not {1 0} at node n
Demand(mk,n,h)                   [MW]     Demand
Efficiency(p2h,h)                  [%]      Efficiency
EmissionMaximum(n,p)             [tP]     Emission limit
EmissionRate(u,p)                [tP\MWh] P emission rate
FlowMaximum(l,h)                 [MW]     Line limits
FlowMinimum(l,h)                 [MW]     Minimum flow
Fuel(u,f)                        [n.a.]   Fuel type {1 0}
HeatDemand(au,h)                [MWh\u]  Heat demand profile for chp units
LineNode(l,n)                    [n.a.]   Incidence matrix {-1 +1}
LoadShedding(n,h)                [MW]   Load shedding capacity
Location(au,n)                    [n.a.]   Location {1 0}
Markup(u,h)                      [EUR\MW]   Markup
OutageFactor(au,h)                [%]      Outage Factor (100% = full outage)
PartLoadMin(au)                 [%]      Minimum part load
PowerCapacity(au)               [MW\u]     Installed capacity
PowerInitial(u)                  [MW\u]     Power output before initial period
PowerMinStable(au)              [MW\u]     Minimum power output
PriceTransmission(l,h)           [EUR\MWh]  Transmission price
StorageChargingCapacity(au)       [MW\u]     Storage capacity
StorageChargingEfficiency(au)   [%]      Charging efficiency
StorageSelfDischarge(au)        [%\day]  Self-discharge of the storage units
RampDownMaximum(u)               [MW\h\u]   Ramp down limit
RampShutDownMaximum(au)           [MW\h\u]   Shut-down ramp limit
RampStartUpMaximum(au)            [MW\h\u]   Start-up ramp limit
RampStartUpMaximumH(u,h)         [MW\h\u]   Start-up ramp limit - Clustered formulation
RampShutDownMaximumH(u,h)        [MW\h\u]   Shut-down ramp limit - Clustered formulation
RampUpMaximum(u)                 [MW\h\u]   Ramp up limit
Reserve(au)                       [n.a.]   Reserve technology {1 0}
StorageCapacity(au)             [MWh\u]    Storage capacity
StorageDischargeEfficiency(au)  [%]      Discharge efficiency
StorageOutflow(u,h)              [MWh\u]    Storage outflows
StorageInflow(u,h)               [MWh\u]    Storage inflows (potential energy)
StorageInitial(au)                [MWh]    Storage level before initial period
StorageProfile(u,h)              [MWh]    Storage level to be resepected at the end of each horizon
StorageMinimum(au)                [MWh\u]    Storage minimum
Technology(au,t)                  [n.a.]   Technology type {1 0}
TimeDownMinimum(u)               [h]      Minimum down time
TimeUpMinimum(u)                 [h]      Minimum up time
$If %RetrieveStatus% == 1 CommittedCalc(u,z)               [n.a.]   Committment status as for the MILP
Nunits(au)                      [n.a.]   Number of units inside the cluster (upper bound value for integer variables)
K_QuickStart(n)                  [n.a.]   Part of the reserve that can be provided by offline quickstart units
QuickStartPower(au,h)            [MW\h\u]   Available max capacity in tertiary regulation up from fast-starting power plants - TC formulation
;

*Parameters as used within the loop
PARAMETERS
CostLoadShedding(n,h)            [EUR\MW]  Value of lost load
LoadMaximum(au,h)                 [%]     Maximum load given AF and OF
PowerMustRun(u,h)                [MW\u]    Minimum power output
StorageFinalMin(s)               [MWh]   Minimum storage level at the end of the optimization horizon
MaxFlexDemand(n)                 [MW]    Maximum value of the flexible demand parameter
MaxOverSupply(n,h)               [MWh]   Maximum flexible demand accumultation
AccumulatedOverSupply_inital(n)  [MWh]   Initial value of the flexible demand accumulation
;

* Scalar variables necessary to the loop:
scalar FirstHour,LastHour,LastKeptHour,day,ndays,failed;
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
$LOAD l
$LOAD au
$LOAD u
$LOAD t
$LOAD tr
$LOAD f
$LOAD p
$LOAD s
$LOAD chp
$LOAD p2h
$LOAD th
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
$LOAD CostLoadShedding
$LOAD CostShutDown
$LOAD CostStartUp
$LOAD CostVariable
$LOAD Curtailment
$LOAD Demand
$LOAD StorageDischargeEfficiency
$LOAD Efficiency
$LOAD EmissionMaximum
$LOAD EmissionRate
$LOAD FlowMaximum
$LOAD FlowMinimum
$LOAD Fuel
$LOAD HeatDemand
$LOAD LineNode
$LOAD LoadShedding
$LOAD Location
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
$If %RetrieveStatus% == 1 $LOAD CommittedCalc
;


$If %Verbose% == 0 $goto skipdisplay

Display
mk,
n,
l,
u,
t,
tr,
f,
p,
s,
chp,
p2h,
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
Fuel,
HeatDemand,
LineNode,
Location,
LoadShedding
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
TimeUpMinimum
$If %RetrieveStatus% == 1 , CommittedCalc
;

$label skipdisplay

*===============================================================================
*Definition of variables
*===============================================================================

VARIABLES
Committed(u,h)      [n.a.]  Unit committed at hour h {1 0} or integer
StartUp(u,h)        [n.a.]  Unit start up at hour h {1 0}  or integer
ShutDown(u,h)       [n.a.]  Unit shut down at hour h {1 0} or integer
;

$If %LPFormulation% == 1 POSITIVE VARIABLES Committed (u,h) ; Committed.UP(u,h) = 1 ;
$If not %LPFormulation% == 1 INTEGER VARIABLES Committed (u,h), StartUp(u,h), ShutDown(u,h) ; Committed.UP(u,h) = Nunits(u) ; StartUp.UP(u,h) = Nunits(u) ; ShutDown.UP(u,h) = Nunits(u) ;



POSITIVE VARIABLES
AccumulatedOverSupply(n,h) [MWh]   Accumulated oversupply due to the flexible demand
CostStartUpH(u,h)          [EUR]   Cost of starting up
CostShutDownH(u,h)         [EUR]   cost of shutting down
CostRampUpH(u,h)           [EUR]   Ramping cost
CostRampDownH(u,h)         [EUR]   Ramping cost
CurtailedPower(n,h)        [MW]    Curtailed power at node n
Flow(l,h)                  [MW]    Flow through lines
Power(u,h)                 [MW]    Power output
PowerConsumption(p2h,h)    [MW]    Power consumption by P2H units
PowerMaximum(u,h)          [MW]    Power output
PowerMinimum(u,h)          [MW]    Power output
ShedLoad(n,h)              [MW]    Shed load
StorageInput(au,h)        [MWh]   Charging input for storage units
StorageLevel(au,h)        [MWh]   Storage level of charge
LL_MaxPower(n,h)           [MW]    Deficit in terms of maximum power
LL_RampUp(u,h)             [MW]    Deficit in terms of ramping up for each plant
LL_RampDown(u,h)           [MW]    Deficit in terms of ramping down
LL_MinPower(n,h)           [MW]    Power exceeding the demand
LL_2U(n,h)                 [MW]    Deficit in reserve up
LL_3U(n,h)                 [MW]    Deficit in reserve up - non spinning
LL_2D(n,h)                 [MW]    Deficit in reserve down
spillage(s,h)              [MWh]   spillage from water reservoirs
SystemCost(h)              [EUR]   Hourly system cost
Reserve_2U(au,h)            [MW]    Spinning reserve up
Reserve_2D(au,h)            [MW]    Spinning reserve down
Reserve_3U(au,h)            [MW]    Non spinning quick start reserve up
Heat(au,h)                [MW]    Heat output by chp plant
HeatSlack(au,h)           [MW]    Heat satisfied by other sources
WaterSlack(s)              [MWh]   Unsatisfied water level constraint
;

free variable
SystemCostD                ![EUR]   Total system cost for one optimization period
DemandModulation(n,h)      [MW] Difference between the flexible demand and the baseline
;

*===============================================================================
*Assignment of initial values
*===============================================================================


*Initial commitment status
CommittedInitial(u)=0;
CommittedInitial(u)$(PowerInitial(u)>0)=1;

* Definition of the minimum stable load:
PowerMinStable(au) = PartLoadMin(au)*PowerCapacity(au);

LoadMaximum(au,h)= AvailabilityFactor(au,h)*(1-OutageFactor(au,h));

* parameters for clustered formulation (quickstart is defined as the capability to go to minimum power in 15 min)
QuickStartPower(au,h) = 0;
QuickStartPower(au,h)$(RampStartUpMaximum(au)>=PowerMinStable(au)*4) = PowerCapacity(au)*LoadMaximum(au,h);

RampStartUpMaximumH(u,h) = min(PowerCapacity(u)*LoadMaximum(u,h),max(RampStartUpMaximum(u),PowerMinStable(u),QuickStartPower(u,h)));
RampShutDownMaximumH(u,h) = min(PowerCapacity(u)*LoadMaximum(u,h),max(RampShutDownMaximum(u),PowerMinStable(u)));

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
EQ_CHP_demand_satisfaction
EQ_CHP_max_heat
EQ_Heat_Storage_balance
EQ_Heat_Storage_minimum
EQ_Heat_Storage_level
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
EQ_P2H
EQ_Max_P2H
EQ_Power_must_run
EQ_Power_available
EQ_Reserve_2U_capability
EQ_Reserve_2D_capability
EQ_Reserve_3U_capability
EQ_p2h_Reserve_2U_capability
EQ_p2h_Reserve_2D_capability
EQ_p2h_Reserve_3U_capability
EQ_Storage_minimum
EQ_Storage_level
EQ_Storage_input
EQ_Storage_MaxDischarge
EQ_Storage_MaxCharge
EQ_Storage_balance
EQ_Storage_boundaries
EQ_SystemCost
EQ_Emission_limits
EQ_Flow_limits_lower
EQ_Flow_limits_upper
EQ_Force_Commitment
EQ_Force_DeCommitment
EQ_LoadShedding
EQ_Flexible_Demand
EQ_Flexible_Demand_Max
EQ_Flexible_Demand_Modulation_Min
EQ_Flexible_Demand_Modulation_Max
EQ_No_Flexible_Demand
$If %RetrieveStatus% == 1 EQ_CommittedCalc
;

$If %RetrieveStatus% == 0 $goto skipequation

EQ_CommittedCalc(u,z)..
         Committed(u,z)
         =E=
         CommittedCalc(u,z)
;

$label skipequation

*Objective function
$IfThen %LPFormulation% == 1

EQ_SystemCost(i)..
         SystemCost(i)
         =E=
         sum(u,CostFixed(u)*TimeStep*Committed(u,i))
         +sum(u,CostRampUpH(u,i) + CostRampDownH(u,i))
         +sum(u,CostVariable(u,i) * Power(u,i)*TimeStep)
         +sum(l,PriceTransmission(l,i)*Flow(l,i)*TimeStep)
         +sum(n,CostLoadShedding(n,i)*ShedLoad(n,i)*TimeStep)
         +sum(th, CostHeatSlack(th,i) * HeatSlack(th,i)*TimeStep)
         +sum(chp, CostVariable(chp,i) * CHPPowerLossFactor(chp) * Heat(chp,i)*TimeStep)
         +Config("ValueOfLostLoad","val")*(sum(n,(LL_MaxPower(n,i)+LL_MinPower(n,i))*TimeStep))
         +0.8*Config("ValueOfLostLoad","val")*(sum(n,(LL_2U(n,i)+LL_2D(n,i)+LL_3U(n,i))*TimeStep))
         +0.7*Config("ValueOfLostLoad","val")*sum(u,(LL_RampUp(u,i)+LL_RampDown(u,i))*TimeStep)
         +Config("CostOfSpillage","val")*sum(s,spillage(s,i));
$else

EQ_SystemCost(i)..
         SystemCost(i)
         =E=
         sum(u,CostFixed(u)*TimeStep*Committed(u,i))
         +sum(u,CostStartUpH(u,i) + CostShutDownH(u,i))
         +sum(u,CostRampUpH(u,i) + CostRampDownH(u,i))
         +sum(u,CostVariable(u,i) * Power(u,i)*TimeStep)
         +sum(l,PriceTransmission(l,i)*Flow(l,i)*TimeStep)
         +sum(n,CostLoadShedding(n,i)*ShedLoad(n,i)*TimeStep)
         +sum(th, CostHeatSlack(th,i) * HeatSlack(th,i)*TimeStep)
         +sum(chp, CostVariable(chp,i) * CHPPowerLossFactor(chp) * Heat(chp,i)*TimeStep)
         +Config("ValueOfLostLoad","val")*(sum(n,(LL_MaxPower(n,i)+LL_MinPower(n,i))*TimeStep))
         +0.8*Config("ValueOfLostLoad","val")*(sum(n,(LL_2U(n,i)+LL_2D(n,i)+LL_3U(n,i))*TimeStep))
         +0.7*Config("ValueOfLostLoad","val")*sum(u,(LL_RampUp(u,i)+LL_RampDown(u,i))*TimeStep)
         +Config("CostOfSpillage","val")*sum(s,spillage(s,i));

$endIf
;



EQ_Objective_function..
         SystemCostD
         =E=
         sum(i,SystemCost(i))
         +Config("WaterValue","val")*sum(s,WaterSlack(s))
;

* 3 binary commitment status
EQ_Commitment(u,i)..
         Committed(u,i)-CommittedInitial(u)$(ord(i) = 1)-Committed(u,i-1)$(ord(i) > 1)
         =E=
         StartUp(u,i) - ShutDown(u,i)
;

* minimum up time
EQ_MinUpTime(u,i)$(TimeStep <= TimeUpMinimum(u))..
         sum(ii$( (ord(ii) >= ord(i) - ceil(TimeUpMinimum(u)/TimeStep)) and (ord(ii) <= ord(i)) ), StartUp(u,ii))
         + sum(h$( (ord(h) >= FirstHour + ord(i) - ceil(TimeUpMinimum(u)/TimeStep) -1) and (ord(h) < FirstHour)),StartUp.L(u,h))
         =L=
         Committed(u,i)
;

* minimum down time
EQ_MinDownTime(u,i)$(TimeStep <= TimeDownMinimum(u))..
         sum(ii$( (ord(ii) >= ord(i) - ceil(TimeDownMinimum(u)/TimeStep)) and (ord(ii) <= ord(i)) ), ShutDown(u,ii))
         + sum(h$( (ord(h) >= FirstHour + ord(i) - ceil(TimeDownMinimum(u)/TimeStep) -1) and (ord(h) < FirstHour)),ShutDown.L(u,h))
         =L=
         Nunits(u)-Committed(u,i)
;

* ramp up constraints
EQ_RampUp_TC(u,i)$(sum(tr,Technology(u,tr))=0)..
         - Power(u,i-1)$(ord(i) > 1) - PowerInitial(u)$(ord(i) = 1) + Power(u,i)
         =L=
         (Committed(u,i) - StartUp(u,i)) * RampUpMaximum(u) * TimeStep + RampStartUpMaximumH(u,i) * TimeStep *  StartUp(u,i) - PowerMustRun(u,i) * ShutDown(u,i) + LL_RampUp(u,i)
;

* ramp down constraints
EQ_RampDown_TC(u,i)$(sum(tr,Technology(u,tr))=0)..
         Power(u,i-1)$(ord(i) > 1) + PowerInitial(u)$(ord(i) = 1) - Power(u,i)
         =L=
         (Committed(u,i) - StartUp(u,i)) * RampDownMaximum(u) * TimeStep + RampShutDownMaximumH(u,i) * TimeStep * ShutDown(u,i) - PowerMustRun(u,i) * StartUp(u,i) + LL_RampDown(u,i)
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

*Hourly demand balance in the day-ahead market for each node
EQ_Demand_balance_DA(n,i)..
         sum(u,Power(u,i)*Location(u,n))
          +sum(l,Flow(l,i)*LineNode(l,n))
         =E=
         Demand("DA",n,i) + Demand("Flex",n,i)
         +DemandModulation(n,i)
         +sum(s,StorageInput(s,i)*Location(s,n))
         -ShedLoad(n,i)
         +sum(p2h,PowerConsumption(p2h,i)*Location(p2h,n))
         -LL_MaxPower(n,i)
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
         sum((au),Reserve_2U(au,i)*Reserve(au)*Location(au,n))
         =G=
         +Demand("2U",n,i)*(1-K_QuickStart(n))
         -LL_2U(n,i)
;

*Hourly demand balance in the upwards non-spinning reserve market for each node
EQ_Demand_balance_3U(n,i)..
         sum((au),(Reserve_2U(au,i) + Reserve_3U(au,i))*Reserve(au)*Location(au,n))
         =G=
         +Demand("2U",n,i)
         -LL_3U(n,i)
;

*Hourly demand balance in the downwards reserve market for each node
EQ_Demand_balance_2D(n,i)..
         sum((au),Reserve_2D(au,i)*Reserve(au)*Location(au,n))
         =G=
         Demand("2D",n,i)
         -LL_2D(n,i)
;

EQ_Reserve_2U_capability(u,i)..
         Reserve_2U(u,i)
         =L=
         PowerCapacity(u)*LoadMaximum(u,i)*Committed(u,i) - Power(u,i)
;

EQ_Reserve_2D_capability(u,i)..
         Reserve_2D(u,i)
         =L=
         (Power(u,i) - PowerMustRun(u,i) * Committed(u,i)) + (StorageChargingCapacity(u)*Nunits(u)-StorageInput(u,i))$(s(u))
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

EQ_p2h_Reserve_3U_capability(p2h,i)$(QuickStartPower(p2h,i) > 0 and StorageCapacity(p2h)>=trp*Efficiency(p2h,i)*PowerCapacity(p2h))..
         Reserve_3U(p2h,i)
         =L=
         Nunits(p2h)*QuickStartPower(p2h,i)*TimeStep-PowerConsumption(p2h,i)
;

*Minimum power output is above the must-run output level for each unit in all periods
EQ_Power_must_run(u,i)..
         PowerMustRun(u,i) * Committed(u,i) - (StorageInput(u,i) * CHPPowerLossFactor(u) )$(chp(u) and (CHPType(u,'Extraction') or CHPType(u,'P2H')))
         =L=
         Power(u,i)
;

*Maximum power output is below the available capacity
EQ_Power_available(u,i)..
         Power(u,i)
         =L=
         PowerCapacity(u)
                 *LoadMaximum(u,i)
                        *Committed(u,i)
;

*Storage level must be above a minimum
EQ_Storage_minimum(s,i)..
         StorageMinimum(s)* Nunits(s)
         =L=
         StorageLevel(s,i)
;

*Storage level must be below storage capacity
EQ_Storage_level(s,i)..
         StorageLevel(s,i)
         =L=
         StorageCapacity(s)*AvailabilityFactor(s,i)*Nunits(s)
;

* Storage charging is bounded by the maximum capacity
EQ_Storage_input(s,i)..
         StorageInput(s,i)
         =L=
         StorageChargingCapacity(s)*(Nunits(s)-Committed(s,i))
;

* The system could curtail by pumping and turbining at the same time if Nunits>1. This should be included into the curtailment equation!

*Discharge is limited by the storage level
EQ_Storage_MaxDischarge(s,i)..
         Power(s,i)*TimeStep/(max(StorageDischargeEfficiency(s),0.0001))
         +StorageOutflow(s,i)*Nunits(s)*TimeStep +Spillage(s,i) - StorageInflow(s,i)*Nunits(s)*TimeStep
         =L=
         StorageLevel(s,i)
;

*Charging is limited by the remaining storage capacity
EQ_Storage_MaxCharge(s,i)..
         StorageInput(s,i)*StorageChargingEfficiency(s)*TimeStep
         -StorageOutflow(s,i)*Nunits(s)*TimeStep -spillage(s,i) + StorageInflow(s,i)*Nunits(s)*TimeStep
         =L=
         StorageCapacity(s)*AvailabilityFactor(s,i)*Nunits(s) - StorageLevel(s,i)
;

*Storage balance
EQ_Storage_balance(s,i)..
         StorageInitial(s)$(ord(i) = 1)
         +StorageLevel(s,i-1)$(ord(i) > 1)
*         +StorageLevelH(h--1,s)
         +StorageInflow(s,i)*Nunits(s)*TimeStep
         +StorageInput(s,i)*StorageChargingEfficiency(s)*TimeStep
         =E=
         StorageLevel(s,i)
         +StorageOutflow(s,i)*Nunits(s)*TimeStep
         +spillage(s,i)
         +Power(s,i)*TimeStep/(max(StorageDischargeEfficiency(s),0.0001))
;

* Minimum level at the end of the optimization horizon:
EQ_Storage_boundaries(s,i)$(ord(i) = card(i))..
         StorageFinalMin(s)
         =L=
         StorageLevel(s,i) + WaterSlack(s)
;

*Total emissions are capped
EQ_Emission_limits(n,i,p)..
         sum(u,Power(u,i)*EmissionRate(u,p)*TimeStep*Location(u,n))
         =L=
         EmissionMaximum(n,p)
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

* Power to heat units
EQ_P2H(p2h,i)..
         StorageInput(p2h,i)
         =E=
         PowerConsumption(p2h,i) * Efficiency(p2h,i)
;


EQ_Max_P2H(p2h,i)..
         PowerConsumption(p2h,i)
         =L=
         PowerCapacity(p2h) * Nunits(p2h)
;

EQ_CHP_demand_satisfaction(th,i)..
         Heat(th,i) + HeatSlack(th,i)
         =E=
         HeatDemand(th,i)
;

*Heat Storage balance
EQ_Heat_Storage_balance(th,i)..
          StorageInitial(th)$(ord(i) = 1)
         +StorageLevel(th,i-1)$(ord(i) > 1)
         +StorageInput(th,i)*TimeStep
         =E=
         StorageLevel(th,i)
         +Heat(th,i)*TimeStep + StorageSelfDischarge(th) * StorageLevel(th,i) * TimeStep/24
;
* The self-discharge proportional to the charging level is a bold hypothesis, but it avoids keeping self-discharging if the level reaches zero

*Storage level must be above a minimum
EQ_Heat_Storage_minimum(th,i)..
         StorageMinimum(th)*Nunits(th)
         =L=
         StorageLevel(th,i)
;

*Storage level must be below storage capacity
EQ_Heat_Storage_level(th,i)..
         StorageLevel(th,i)
         =L=
         StorageCapacity(th)*Nunits(th)
;

* Minimum level at the end of the optimization horizon:
*EQ_Heat_Storage_boundaries(chp,i)$(ord(i) = card(i))..
*         StorageFinalMin(chp)
*         =L=
*         StorageLevel(chp,i)
*;
*===============================================================================
*Definition of models
*===============================================================================
MODEL UCM_SIMPLE /
EQ_Objective_function,
EQ_CHP_extraction_Pmax,
EQ_CHP_extraction,
EQ_CHP_backpressure,
EQ_CHP_demand_satisfaction,
EQ_CHP_max_heat,
EQ_CostRampUp,
EQ_CostRampDown,
$If not %LPFormulation% == 1 EQ_CostStartUp,
$If not %LPFormulation% == 1 EQ_CostShutDown,
$If not %LPFormulation% == 1 EQ_Commitment,
$If not %LPFormulation% == 1 EQ_MinUpTime,
$If not %LPFormulation% == 1 EQ_MinDownTime,
EQ_RampUp_TC,
EQ_RampDown_TC,
EQ_Demand_balance_DA,
EQ_Demand_balance_2U,
EQ_Demand_balance_2D,
EQ_Demand_balance_3U,
$If not %LPFormulation% == 1 EQ_Power_must_run,
EQ_P2H,
EQ_Max_P2H,
EQ_Power_available,
EQ_Heat_Storage_balance,
EQ_Heat_Storage_minimum,
EQ_Heat_Storage_level,
EQ_Reserve_2U_capability,
EQ_Reserve_2D_capability,
EQ_Reserve_3U_capability,
EQ_p2h_Reserve_2U_capability,
EQ_p2h_Reserve_2D_capability,
EQ_p2h_Reserve_3U_capability,
EQ_Storage_minimum,
EQ_Storage_level,
EQ_Storage_input,
EQ_Storage_balance,
EQ_Storage_boundaries,
EQ_Storage_MaxCharge
EQ_Storage_MaxDischarge
EQ_SystemCost
*EQ_Emission_limits,
EQ_Flow_limits_lower,
EQ_Flow_limits_upper,
$If not %LPFormulation% == 1 EQ_Force_Commitment,
$If not %LPFormulation% == 1 EQ_Force_DeCommitment,
EQ_LoadShedding,
$If %ActivateFlexibleDemand% == 1 EQ_Flexible_Demand,
$If %ActivateFlexibleDemand% == 1 EQ_Flexible_Demand_Max,
$if not %ActivateFlexibleDemand% == 1 EQ_No_Flexible_Demand,
EQ_Flexible_Demand_Modulation_Min,
EQ_Flexible_Demand_Modulation_Max,
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
parameter CommittedInitial_dbg(u), PowerInitial_dbg(u), StorageInitial_dbg(s);

* Fixing the initial guesses:
*PowerH.L(u,i)=PowerInitial(u);
*Committed.L(u,i)=CommittedInitial(u);

* Defining a parameter that records the solver status:
set  tmp   "tpm"  / "model", "solver" /  ;
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

*        Defining the minimum level at the end of the horizon, ensuring that it is feasible with the provided inflows:
         StorageFinalMin(s) =  min(StorageInitial(s) + (sum(i,StorageInflow(s,i)*TimeStep) - sum(i,StorageOutflow(s,i)*TimeStep))*Nunits(s), sum(i$(ord(i)=card(i)),StorageProfile(s,i)*Nunits(s)*StorageCapacity(s)*AvailabilityFactor(s,i)));
*        Correcting the minimum level to avoid the infeasibility in case it is too close to the StorageCapacity:
         StorageFinalMin(s) = min(StorageFinalMin(s),Nunits(s)*StorageCapacity(s) - Nunits(s)*smax(i,StorageInflow(s,i)*TimeStep));

$If %Verbose% == 1   Display PowerInitial,CommittedInitial,StorageFinalMin;
$If %Verbose% == 1   Display PowerInitial,StorageFinalMin;

$If %LPFormulation% == 1          SOLVE UCM_SIMPLE USING LP MINIMIZING SystemCostD;
$If not %LPFormulation% == 1      SOLVE UCM_SIMPLE USING MIP MINIMIZING SystemCostD;

$If %Verbose% == 0 $goto skipdisplay3
$Ifthen %LPFormulation% == 1             Display EQ_Objective_function.M, EQ_CostRampUp.M, EQ_CostRampDown.M, EQ_Demand_balance_DA.M, EQ_Storage_minimum.M, EQ_Storage_level.M, EQ_Storage_input.M, EQ_Storage_balance.M, EQ_Storage_boundaries.M, EQ_Storage_MaxCharge.M, EQ_Storage_MaxDischarge.M, EQ_Flow_limits_lower.M ;
$else not %LPFormulation% == 1           Display EQ_Objective_function.M, EQ_CostStartUp.M, EQ_CostShutDown.M, EQ_Commitment.M, EQ_MinUpTime.M, EQ_MinDownTime.M, EQ_RampUp_TC.M, EQ_RampDown_TC.M, EQ_Demand_balance_DA.M, EQ_Demand_balance_2U.M, EQ_Demand_balance_2D.M, EQ_Demand_balance_3U.M, EQ_Reserve_2U_capability.M, EQ_Reserve_2D_capability.M, EQ_Reserve_3U_capability.M, EQ_Power_must_run.M, EQ_Power_available.M, EQ_Storage_minimum.M, EQ_Storage_level.M, EQ_Storage_input.M, EQ_Storage_balance.M, EQ_Storage_boundaries.M, EQ_Storage_MaxCharge.M, EQ_Storage_MaxDischarge.M, EQ_SystemCost.M, EQ_Flow_limits_lower.M, EQ_Flow_limits_upper.M, EQ_Force_Commitment.M, EQ_Force_DeCommitment.M, EQ_LoadShedding.M ;
$endIf;
$label skipdisplay3

         status("model",i) = UCM_SIMPLE.Modelstat;
         status("solver",i) = UCM_SIMPLE.Solvestat;

if(UCM_SIMPLE.Modelstat <> 1 and UCM_SIMPLE.Modelstat <> 8 and not failed, CommittedInitial_dbg(u) = CommittedInitial(u); PowerInitial_dbg(u) = PowerInitial(u); StorageInitial_dbg(s) = StorageInitial(s);
                                                                           EXECUTE_UNLOAD "debug.gdx" day, status, CommittedInitial_dbg, PowerInitial_dbg, StorageInitial_dbg;
                                                                           failed=1;);

         CommittedInitial(u)=sum(i$(ord(i)=LastKeptHour-FirstHour+1),Committed.L(u,i));
         PowerInitial(u) = sum(i$(ord(i)=LastKeptHour-FirstHour+1),Power.L(u,i));

         StorageInitial(s) =   sum(i$(ord(i)=LastKeptHour-FirstHour+1),StorageLevel.L(s,i));
         StorageInitial(chp) =   sum(i$(ord(i)=LastKeptHour-FirstHour+1),StorageLevel.L(chp,i));

$If %ActivateFlexibleDemand% == 1 AccumulatedOverSupply_inital(n) = sum(i$(ord(i)=LastKeptHour-FirstHour+1),AccumulatedOverSupply.L(n,i));


*Loop variables to display after solving:
$If %Verbose% == 1 Display LastKeptHour,PowerInitial,StorageInitial;
);

CurtailedPower.L(n,z)=sum(u,(Nunits(u)*PowerCapacity(u)*LoadMaximum(u,z)-Power.L(u,z))$(sum(tr,Technology(u,tr))>=1) * Location(u,n));

$If %Verbose% == 1 Display Flow.L,Power.L,Committed.L,ShedLoad.L,CurtailedPower.L,StorageLevel.L,StorageInput.L,SystemCost.L,LL_MaxPower.L,LL_MinPower.L,LL_2U.L,LL_2D.L,LL_RampUp.L,LL_RampDown.L;

*===============================================================================
*Result export
*===============================================================================

PARAMETER
OutputCommitted(u,h)
OutputFlow(l,h)
OutputPower(u,h)
OutputPowerConsumption(p2h,h)
OutputStorageInput(au,h)
OutputStorageLevel(au,h)
OutputSystemCost(h)
OutputSpillage(s,h)
OutputShedLoad(n,h)
OutputCurtailedPower(n,h)
$If %ActivateFlexibleDemand% == 1 OutputDemandModulation(n,h)
ShadowPrice(n,h)
ReserveUpShadowPrice(n,h)
ReserveDownShadowPrice(n,h)
HeatShadowPrice(au,h)
LostLoad_MaxPower(n,h)
LostLoad_MinPower(n,h)
LostLoad_2D(n,h)
LostLoad_2U(n,h)
LostLoad_3U(n,h)
LostLoad_RampUp(n,h)
LostLoad_RampDown(n,h)
OutputGenMargin(n,h)
OutputHeat(au,h)
OutputHeatSlack(au,h)
LostLoad_WaterSlack(s)
StorageShadowPrice(au,h)
OutputReserve_2U(au,h)
OutputReserve_2D(au,h)
OutputReserve_3U(au,h)

;

OutputCommitted(u,z)=Committed.L(u,z);
OutputFlow(l,z)=Flow.L(l,z);
OutputPower(u,z)=Power.L(u,z);
OutputPowerConsumption(p2h,z)=PowerConsumption.L(p2h,z);
OutputHeat(au,z)=Heat.L(au,z);
OutputHeatSlack(au,z)=HeatSlack.L(au,z);
OutputStorageInput(s,z)=StorageInput.L(s,z);
OutputStorageInput(th,z)=StorageInput.L(th,z);
OutputStorageLevel(s,z)=StorageLevel.L(s,z)/StorageCapacity(s);
OutputStorageLevel(th,z)=StorageLevel.L(th,z)/StorageCapacity(th);
OutputSystemCost(z)=SystemCost.L(z);
OutputSpillage(s,z)  = Spillage.L(s,z) ;
OutputShedLoad(n,z) = ShedLoad.L(n,z);
OutputCurtailedPower(n,z)=CurtailedPower.L(n,z);
$If %ActivateFlexibleDemand% == 1 OutputDemandModulation(n,z)=DemandModulation.L(n,z);
LostLoad_MaxPower(n,z)  = LL_MaxPower.L(n,z);
LostLoad_MinPower(n,z)  = LL_MinPower.L(n,z);
LostLoad_2D(n,z) = LL_2D.L(n,z);
LostLoad_2U(n,z) = LL_2U.L(n,z);
LostLoad_3U(n,z) = LL_3U.L(n,z);
LostLoad_RampUp(n,z)    = sum(u,LL_RampUp.L(u,z)*Location(u,n));
LostLoad_RampDown(n,z)  = sum(u,LL_RampDown.L(u,z)*Location(u,n));
ShadowPrice(n,z) = EQ_Demand_balance_DA.m(n,z);
ReserveUpShadowPrice(n,z) = EQ_Demand_balance_2U.m(n,z);
ReserveDownShadowPrice(n,z) = EQ_Demand_balance_2D.m(n,z);
HeatShadowPrice(au,z) = EQ_CHP_demand_satisfaction.m(au,z);
LostLoad_WaterSlack(s) = WaterSlack.L(s);
StorageShadowPrice(s,z) = EQ_Storage_balance.m(s,z);
StorageShadowPrice(th,z) = EQ_Heat_Storage_balance.m(th,z);
OutputReserve_2U(au,z)=Reserve_2U.L(au,z);
OutputReserve_2D(au,z)=Reserve_2D.L(au,z);
OutputReserve_3U(au,z)=Reserve_3U.L(au,z);

EXECUTE_UNLOAD "Results.gdx"
$If %MTS%==0 OutputCommitted,
OutputFlow,
OutputPower,
OutputPowerConsumption,
OutputHeat,
OutputHeatSlack,
OutputStorageInput,
OutputStorageLevel,
OutputSystemCost,
OutputSpillage,
OutputShedLoad,
OutputCurtailedPower,
$If %ActivateFlexibleDemand% == 1 OutputDemandModulation,
OutputGenMargin,
LostLoad_MaxPower,
LostLoad_MinPower,
LostLoad_2D,
LostLoad_2U,
LostLoad_3U,
LostLoad_RampUp,
LostLoad_RampDown,
ShadowPrice,
ReserveUpShadowPrice,
ReserveDownShadowPrice,
HeatShadowPrice,
LostLoad_WaterSlack,
StorageShadowPrice,
OutputReserve_2U,
OutputReserve_2D,
OutputReserve_3U,
status
;

display OutputPowerConsumption, heat.L, heatslack.L, powerconsumption.L;

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
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=Y var=LostLoad_RampUp rng=LostLoad_RampUp!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=Y var=LostLoad_RampDown rng=LostLoad_RampDown!A1 rdim=1 cdim=1'


$exit

$Label DebugSection

$gdxin debug.gdx
$LOAD day
$LOAD PowerInitial_dbg
$LOAD CommittedInitial_dbg
$LOAD StorageInitial_dbg
;
PowerInitial(u) = PowerInitial_dbg(u);
CommittedInitial(u) = CommittedInitial_dbg(u);
StorageInitial(s) = StorageInitial_dbg(s);
FirstHour = (day-1)*24/TimeStep+1;
LastHour = min(card(h),FirstHour + (Config("RollingHorizon Length","day")+Config("RollingHorizon LookAhead","day")) * 24/TimeStep - 1);
LastKeptHour = LastHour - Config("RollingHorizon LookAhead","day") * 24/TimeStep;
i(h) = no;
i(h)$(ord(h)>=firsthour and ord(h)<=lasthour)=yes;
StorageFinalMin(s) =  min(StorageInitial(s) + sum(i,StorageInflow(s,i)*TimeStep) - sum(i,StorageOutflow(s,i)*TimeStep) , sum(i$(ord(i)=card(i)),StorageProfile(s,i)*StorageCapacity(s)*AvailabilityFactor(s,i)));
StorageFinalMin(s) = min(StorageFinalMin(s),StorageCapacity(s) - smax(i,StorageInflow(s,i)*TimeStep));

$If %Verbose% == 1   Display TimeUpLeft_initial,TimeUpLeft_JustStarted,PowerInitial,CommittedInitial,StorageFinalMin;
$If %LPFormulation% == 1          SOLVE UCM_SIMPLE USING LP MINIMIZING SystemCostD;
$If not %LPFormulation% == 1      SOLVE UCM_SIMPLE USING MIP MINIMIZING SystemCostD;
$If %LPFormulation% == 1          Display EQ_Objective_function.M, EQ_CostRampUp.M, EQ_CostRampDown.M, EQ_Demand_balance_DA.M, EQ_Demand_balance_2U.M, EQ_Demand_balance_2D.M, EQ_Demand_balance_3U.M, EQ_Power_available.M, EQ_Ramp_up.M, EQ_Ramp_down.M, EQ_Storage_minimum.M, EQ_Storage_level.M, EQ_Storage_input.M, EQ_Storage_balance.M, EQ_Storage_boundaries.M, EQ_Storage_MaxCharge.M, EQ_Storage_MaxDischarge.M, EQ_Flow_limits_lower.M ;
$If not %LPFormulation% == 1      Display EQ_Objective_function.M, EQ_CostStartUp.M, EQ_CostShutDown.M, EQ_Demand_balance_DA.M, EQ_Demand_balance_2U.M, EQ_Demand_balance_2D.M, EQ_Demand_balance_3U.M, EQ_Power_must_run.M, EQ_Power_available.M, EQ_Ramp_up.M, EQ_Ramp_down.M, EQ_MaxShutDowns.M, EQ_MaxShutDowns_JustStarted.M, EQ_MaxStartUps.M, EQ_MaxStartUps_JustStopped.M, EQ_Storage_minimum.M, EQ_Storage_level.M, EQ_Storage_input.M, EQ_Storage_balance.M, EQ_Storage_boundaries.M, EQ_Storage_MaxCharge.M, EQ_Storage_MaxDischarge.M, EQ_Flow_limits_lower.M ;

display day,FirstHour,LastHour,LastKeptHour;
Display StorageFinalMin,PowerInitial,CommittedInitial,StorageFinalMin;
Display Flow.L,Power.L,Committed.L,ShedLoad.L,StorageLevel.L,StorageInput.L,SystemCost.L,Spillage.L,StorageLevel.L,StorageInput.L,LL_MaxPower.L,LL_MinPower.L,LL_2U.L,LL_2D.L,LL_RampUp.L,LL_RampDown.L;
