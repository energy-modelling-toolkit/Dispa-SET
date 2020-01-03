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
$set InputFileName Inputs.gdx

* Definition of the equations in LP
$setglobal LPFormulation 1
* Flag to retrieve status or not
* (1 to retrieve 0 to not)
$setglobal RetrieveStatus 0

*===============================================================================
*Definition of   sets and parameters
*===============================================================================
SETS
mk               Markets
n                Nodes
l                Lines
u                Units
t                Generation technologies
tr(t)            Renewable generation technologies
f                Fuel types
p                Pollutants
s(u)             Storage Units (with reservoir)
h                Hours
i(h)             Subset of simulated hours for one iteration
z(h)             Subset of all simulated hours
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
AvailabilityFactor(u,h)          [%]      Availability factor
Config
*CostCurtailment(n,h)             [EUR\MW]  Curtailment costs
CostFixed(u)                     [EUR\h]    Fixed costs
CostRampUp(u)                    [EUR\MW\h] Ramp-up costs
CostRampDown(u)                  [EUR\MW\h] Ramp-down costs
CostVariable(u,h)                [EUR\MW]   Variable costs
CostLoadShedding(n,h)            [EUR\MWh] Cost of load shedding
Curtailment(n)                   [n.a]    Curtailment allowed or not {1 0} at node n
Demand(mk,n,h)                   [MW]     Demand
*Efficiency(u)                    [%]      Efficiency
EmissionMaximum(n,p)             [tP]     Emission limit
EmissionRate(u,p)                [tP\MWh] P emission rate
FlowMaximum(l,h)                 [MW]     Line limits
FlowMinimum(l,h)                 [MW]     Minimum flow
Fuel(u,f)                        [n.a.]   Fuel type {1 0}
LineNode(l,n)                    [n.a.]   Incidence matrix {-1 +1}
LoadShedding(n,h)                [MW]   Load shedding capacity
Location(u,n)                    [n.a.]   Location {1 0}
Markup(u,h)                      [EUR\MW]   Markup
OutageFactor(u,h)                [%]      Outage Factor (100% = full outage)
PartLoadMin(u)                   [%]      Minimum part load
PowerCapacity(u)                 [MW\u]     Installed capacity
PowerInitial(u)                  [MW\u]     Power output before initial period
PowerMinStable(u)                [MW\u]     Minimum power output
PriceTransmission(l,h)           [EUR\MWh]  Transmission price
StorageChargingCapacity(u)       [MW\u]     Storage capacity
StorageChargingEfficiency(u)     [%]      Charging efficiency
StorageSelfDischarge(u)          [%\day]  Self-discharge of the storage units
Reserve(t)                       [n.a.]   Reserve technology {1 0}
StorageCapacity(u)               [MWh\u]    Storage capacity
StorageDischargeEfficiency(u)    [%]      Discharge efficiency
StorageOutflow(u,h)              [MWh\u]    Storage outflows
StorageInflow(u,h)               [MWh\u]    Storage inflows (potential energy)
StorageInitial(u)                [MWh]    Storage level before initial period
StorageProfile(u,h)              [MWh]    Storage level to be resepected at the end of each horizon
StorageMinimum(u)                [MWh\u]    Storage minimum
Technology(u,t)                  [n.a.]   Technology type {1 0}
Nunits(u)                        [n.a.]   Number of units inside the cluster (upper bound value for integer variables)
K_QuickStart(n)                      [n.a.]   Part of the reserve that can be provided by offline quickstart units
QuickStartPower(u,h)            [MW\h\u]   Available max capacity in tertiary regulation up from fast-starting power plants - TC formulation
;

*Parameters as used within the loop
PARAMETERS
CostLoadShedding(n,h)            [EUR\MW]  Value of lost load
LoadMaximum(u,h)                 [%]     Maximum load given AF and OF
StorageFinalMin(s)               [MWh]   Minimum storage level at the end of the optimization horizon
;

* Scalar variables necessary to the loop:
scalar FirstHour,LastHour,LastKeptHour,day,ndays,failed;
FirstHour = 1;

*===============================================================================
*Data import
*===============================================================================

$gdxin %inputfilename%

$LOAD mk
$LOAD n
$LOAD l
$LOAD u
$LOAD t
$LOAD tr
$LOAD f
$LOAD p
$LOAD s
$LOAD h
$LOAD z
$LOAD AvailabilityFactor
$LOAD Config
$LOAD CostFixed
$LOAD CostLoadShedding
$LOAD CostVariable
$LOAD Curtailment
$LOAD Demand
$LOAD StorageDischargeEfficiency
*$LOAD Efficiency
$LOAD EmissionMaximum
$LOAD EmissionRate
$LOAD FlowMaximum
$LOAD FlowMinimum
$LOAD Fuel
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
$LOAD Reserve
$LOAD StorageCapacity
$LOAD StorageInflow
$LOAD StorageInitial
$LOAD StorageProfile
$LOAD StorageMinimum
$LOAD StorageOutflow
$LOAD Technology
$LOAD CostRampUp
$LOAD CostRampDown
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
h,
AvailabilityFactor,
Config,
CostFixed,
CostRampUp,
CostVariable,
Demand,
StorageDischargeEfficiency,
*Efficiency,
EmissionMaximum,
EmissionRate,
FlowMaximum,
FlowMinimum,
Fuel,
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
Reserve,
StorageCapacity,
StorageInflow,
StorageInitial,
StorageProfile,
StorageMinimum,
StorageOutflow,
Technology,
;

$label skipdisplay

*===============================================================================
*Definition of variables
*===============================================================================
POSITIVE VARIABLES
CostRampUpH(u,h)           [EUR]   Ramping cost
CostRampDownH(u,h)         [EUR]   Ramping cost
CurtailedPower(n,h)        [MW]    Curtailed power at node n
Flow(l,h)                  [MW]    Flow through lines
Power(u,h)                 [MW]    Power output
PowerMaximum(u,h)          [MW]    Power output
PowerMinimum(u,h)          [MW]    Power output
ShedLoad(n,h)              [MW]    Shed load
StorageInput(u,h)          [MWh]   Charging input for storage units
StorageLevel(u,h)          [MWh]   Storage level of charge
LL_MaxPower(n,h)           [MW]    Deficit in terms of maximum power
LL_MinPower(n,h)           [MW]    Power exceeding the demand
LL_2U(n,h)                 [MW]    Deficit in reserve up
LL_3U(n,h)                 [MW]    Deficit in reserve up - non spinning
LL_2D(n,h)                 [MW]    Deficit in reserve down
spillage(s,h)              [MWh]   spillage from water reservoirs
SystemCost(h)              [EUR]   Hourly system cost
Reserve_2U(u,h)            [MW]    Spinning reserve up
Reserve_2D(u,h)            [MW]    Spinning reserve down
Reserve_3U(u,h)            [MW]    Non spinning quick start reserve up
WaterSlack(s)              [MWh]   Unsatisfied water level constraint
;

free variable
SystemCostD                ![EUR]   Total system cost for one optimization period
;

*===============================================================================
*Assignment of initial values
*===============================================================================
* Definition of the minimum stable load:
PowerMinStable(u) = PartLoadMin(u)*PowerCapacity(u);

LoadMaximum(u,h)= AvailabilityFactor(u,h)*(1-OutageFactor(u,h))*Nunits(u);

* parameters for clustered formulation (quickstart is defined as the capability to go to minimum power in 15 min)
QuickStartPower(u,h) = 0;
QuickStartPower(u,h) = PowerCapacity(u)*LoadMaximum(u,h);

* Part of the reserve that can be provided by offline quickstart units:
K_QuickStart(n) = Config("QuickStartShare","val");

$offorder

*===============================================================================
*Declaration and definition of equations
*===============================================================================
EQUATIONS
EQ_Objective_function
EQ_CostRampUp
EQ_CostRampDown
EQ_Demand_balance_DA
EQ_Demand_balance_2U
EQ_Demand_balance_3U
EQ_Demand_balance_2D
EQ_Power_available
EQ_Reserve_2U_capability
EQ_Reserve_2D_capability
EQ_Reserve_3U_capability
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
EQ_LoadShedding
;

$If %RetrieveStatus% == 0 $goto skipequation

$label skipequation

*Objective function

EQ_SystemCost(i)..
         SystemCost(i)
         =E=
         sum(u,CostFixed(u))
         +sum(u,CostRampUpH(u,i) + CostRampDownH(u,i))
         +sum(u,CostVariable(u,i) * Power(u,i))
         +sum(l,PriceTransmission(l,i)*Flow(l,i))
         +sum(n,CostLoadShedding(n,i)*ShedLoad(n,i))
         +Config("ValueOfLostLoad","val")*(sum(n,LL_MaxPower(n,i)+LL_MinPower(n,i)))
         +0.8*Config("ValueOfLostLoad","val")*(sum(n,LL_2U(n,i)+LL_2D(n,i)+LL_3U(n,i)))
         +Config("CostOfSpillage","val")*sum(s,spillage(s,i));

EQ_Objective_function..
         SystemCostD
         =E=
         sum(i,SystemCost(i))
         +Config("WaterValue","val")*sum(s,WaterSlack(s))
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
         Demand("DA",n,i)
         +sum(s,StorageInput(s,i)*Location(s,n))
         -ShedLoad(n,i)
         -LL_MaxPower(n,i)
         +LL_MinPower(n,i)
;

*Hourly demand balance in the upwards spinning reserve market for each node
EQ_Demand_balance_2U(n,i)..
         sum((u,t),Reserve_2U(u,i)*Technology(u,t)*Reserve(t)*Location(u,n))
         =G=
         +Demand("2U",n,i)*(1-K_QuickStart(n))
         -LL_2U(n,i)
;

*Hourly demand balance in the upwards non-spinning reserve market for each node
EQ_Demand_balance_3U(n,i)..
         sum((u,t),(Reserve_2U(u,i) + Reserve_3U(u,i))*Technology(u,t)*Reserve(t)*Location(u,n))
         =G=
         +Demand("2U",n,i)
         -LL_3U(n,i)
;

*Hourly demand balance in the downwards reserve market for each node
EQ_Demand_balance_2D(n,i)..
         sum((u,t),Reserve_2D(u,i)*Technology(u,t)*Reserve(t)*Location(u,n))
         =G=
         Demand("2D",n,i)
         -LL_2D(n,i)
;

EQ_Reserve_2U_capability(u,i)..
         Reserve_2U(u,i)
         =L=
         PowerCapacity(u)*LoadMaximum(u,i) - Power(u,i)
;

* New equation
EQ_Reserve_2D_capability(u,i)..
         Reserve_2D(u,i)
         =L=
         Power(u,i) + (StorageChargingCapacity(u)*Nunits(u)-StorageInput(u,i))$(s(u))
;

EQ_Reserve_3U_capability(u,i)$(QuickStartPower(u,i) > 0)..
         Reserve_3U(u,i)
         =L=
         Nunits(u)*QuickStartPower(u,i)
;

*Maximum power output is below the available capacity
EQ_Power_available(u,i)..
         Power(u,i)
         =L=
         PowerCapacity(u)
                 *LoadMaximum(u,i)
;

*Storage level must be above a minimum
EQ_Storage_minimum(s,i)..
         StorageMinimum(s)* Nunits(s)
         =L=
         StorageLevel(s,i)
;

*Storage level must below storage capacity
EQ_Storage_level(s,i)..
         StorageLevel(s,i)
         =L=
         StorageCapacity(s)*AvailabilityFactor(s,i)*Nunits(s)
;

* Storage charging is bounded by the maximum capacity
EQ_Storage_input(s,i)..
         StorageInput(s,i)
         =L=
         StorageChargingCapacity(s)*Nunits(s)
;
* The system could curtail by pumping and turbining at the same time if Nunits>1. This should be included into the curtailment equation!

*Discharge is limited by the storage level
EQ_Storage_MaxDischarge(s,i)..
         Power(s,i)/(max(StorageDischargeEfficiency(s),0.0001))
         +StorageOutflow(s,i)*Nunits(s) +Spillage(s,i) - StorageInflow(s,i)*Nunits(s)
         =L=
         StorageLevel(s,i)
;

*Charging is limited by the remaining storage capacity
EQ_Storage_MaxCharge(s,i)..
         StorageInput(s,i) * StorageChargingEfficiency(s)
         -StorageOutflow(s,i)*Nunits(s) -spillage(s,i) + StorageInflow(s,i)*Nunits(s)
         =L=
         StorageCapacity(s)*AvailabilityFactor(s,i)*Nunits(s) - StorageLevel(s,i)
;

*Storage balance
EQ_Storage_balance(s,i)..
         StorageInitial(s)$(ord(i) = 1)
         +StorageLevel(s,i-1)$(ord(i) > 1)
*         +StorageLevelH(h--1,s)
         +StorageInflow(s,i)*Nunits(s)
         +StorageInput(s,i)*StorageChargingEfficiency(s)
         =E=
         StorageLevel(s,i)
         +StorageOutflow(s,i)*Nunits(s)
         +spillage(s,i)
         +Power(s,i)/(max(StorageDischargeEfficiency(s),0.0001))
;

* Minimum level at the end of the optimization horizon:
EQ_Storage_boundaries(s,i)$(ord(i) = card(i))..
         StorageFinalMin(s)
         =L=
         StorageLevel(s,i) + WaterSlack(s)
;

*Total emissions are capped
EQ_Emission_limits(n,i,p)..
         sum(u,Power(u,i)*EmissionRate(u,p)*Location(u,n))
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

*Load shedding
EQ_LoadShedding(n,i)..
         ShedLoad(n,i)
         =L=
         LoadShedding(n,i)
;

*===============================================================================
*Definition of models
*===============================================================================
MODEL UCM_SIMPLE /
EQ_Objective_function,
EQ_CostRampUp,
EQ_CostRampDown,
EQ_Demand_balance_DA,
EQ_Demand_balance_2U,
EQ_Demand_balance_2D,
EQ_Demand_balance_3U,
EQ_Power_available,
EQ_Reserve_2U_capability,
EQ_Reserve_2D_capability,
EQ_Reserve_3U_capability,
EQ_Storage_minimum,
EQ_Storage_level,
EQ_Storage_input,
EQ_Storage_balance,
EQ_Storage_boundaries,
EQ_Storage_MaxCharge
EQ_Storage_MaxDischarge
EQ_SystemCost
EQ_Emission_limits,
EQ_Flow_limits_lower,
EQ_Flow_limits_upper,
EQ_LoadShedding,
/
;
UCM_SIMPLE.optcr = 0.01;
UCM_SIMPLE.optfile=1;

*===============================================================================
*Solving loop
*===============================================================================

ndays = floor(card(h)/24);
Config("RollingHorizon LookAhead","day")=0;
Config("RollingHorizon Length","day")=ndays;

if (Config("RollingHorizon LookAhead","day") > ndays, abort "The look ahead period is longer than the simulation length";);

* Some parameters used for debugging:
failed=0;
parameter PowerInitial_dbg(u), StorageInitial_dbg(s);

* Defining a parameter that records the solver status:
set  tmp   "tpm"  / "model", "solver" /  ;
parameter status(tmp,h);

$if %Debug% == 1 $goto DebugSection

display "OK";

scalar starttime;
set days /1,'ndays'/;
display days;
PARAMETER elapsed(days);

*FOR(day = 1 TO ndays-Config("RollingHorizon LookAhead","day") by Config("RollingHorizon Length","day"),
FOR(day = 1 TO ndays by ndays,
         FirstHour = (day-1)*24+1;
*         LastHour = min(card(h),FirstHour + (Config("RollingHorizon Length","day")+Config("RollingHorizon LookAhead","day")) * 24 - 1);
         LastHour = min(card(h),FirstHour + ndays * 24 - 1);
*         LastKeptHour = LastHour - Config("RollingHorizon LookAhead","day") * 24;
         LastKeptHour = LastHour;
         i(h) = no;
         i(h)$(ord(h)>=firsthour and ord(h)<=lasthour)=yes;
         display day,FirstHour,LastHour,LastKeptHour;

*        Defining the minimum level at the end of the horizon, ensuring that it is feasible with the provided inflows:
         StorageFinalMin(s) =  min(StorageInitial(s) + (sum(i,StorageInflow(s,i)) - sum(i,StorageOutflow(s,i)))*Nunits(s), sum(i$(ord(i)=card(i)),StorageProfile(s,i)*Nunits(s)*StorageCapacity(s)*AvailabilityFactor(s,i)));
*        Correcting the minimum level to avoid the infeasibility in case it is too close to the StorageCapacity:
         StorageFinalMin(s) = min(StorageFinalMin(s),Nunits(s)*StorageCapacity(s) - Nunits(s)*smax(i,StorageInflow(s,i)));

$If %Verbose% == 1   Display PowerInitial,StorageFinalMin;

$If %LPFormulation% == 1          SOLVE UCM_SIMPLE USING LP MINIMIZING SystemCostD;

$If %Verbose% == 0 $goto skipdisplay2
$If %LPFormulation% == 1          Display EQ_Objective_function.M, EQ_CostRampUp.M, EQ_CostRampDown.M, EQ_Demand_balance_DA.M, EQ_Storage_minimum.M, EQ_Storage_level.M, EQ_Storage_input.M, EQ_Storage_balance.M, EQ_Storage_boundaries.M, EQ_Storage_MaxCharge.M, EQ_Storage_MaxDischarge.M, EQ_Flow_limits_lower.M ;
$label skipdisplay2

         status("model",i) = UCM_SIMPLE.Modelstat;
         status("solver",i) = UCM_SIMPLE.Solvestat;

if(UCM_SIMPLE.Modelstat <> 1 and UCM_SIMPLE.Modelstat <> 8 and not failed, PowerInitial_dbg(u) = PowerInitial(u); StorageInitial_dbg(s) = StorageInitial(s);
                                                                           EXECUTE_UNLOAD "debug.gdx" day, status, PowerInitial_dbg, StorageInitial_dbg;
                                                                           failed=1;);
* Initial power output
         PowerInitial(u) = sum(i$(ord(i)=LastKeptHour-FirstHour+1),Power.L(u,i));
* Initial storage level
         StorageInitial(s) =   sum(i$(ord(i)=LastKeptHour-FirstHour+1),StorageLevel.L(s,i));

*Loop variables to display after solving:
$If %Verbose% == 1 Display LastKeptHour,PowerInitial;

);

CurtailedPower.L(n,z)=sum(u,(Nunits(u)*PowerCapacity(u)*LoadMaximum(u,z)-Power.L(u,z))$(sum(tr,Technology(u,tr))>=1) * Location(u,n));

$If %Verbose% == 1 Display Flow.L,Power.L,ShedLoad.L,CurtailedPower.L,StorageLevel.L,StorageInput.L,SystemCost.L,LL_MaxPower.L,LL_MinPower.L,LL_2U.L,LL_2D.L;


*===============================================================================
*Result export
*===============================================================================

PARAMETER
OutputFlow(l,h)
OutputPower(u,h)
OutputStorageInput(u,h)
OutputStorageLevel(u,h)
OutputSystemCost(h)
OutputSpillage(s,h)
OutputShedLoad(n,h)
OutputCurtailedPower(n,h)
ShadowPrice(n,h)
LostLoad_MaxPower(n,h)
LostLoad_MinPower(n,h)
LostLoad_2D(n,h)
LostLoad_2U(n,h)
LostLoad_3U(n,h)
OutputGenMargin(n,h)
LostLoad_WaterSlack(s)
StorageShadowPrice(s,h)
;

OutputFlow(l,z) = Flow.L(l,z);
OutputPower(u,z) = Power.L(u,z);
OutputStorageInput(s,z) = StorageInput.L(s,z);
OutputStorageLevel(s,z) = StorageLevel.L(s,z)/StorageCapacity(s);
OutputSystemCost(z) = SystemCost.L(z);
OutputSpillage(s,z)  = Spillage.L(s,z) ;
OutputShedLoad(n,z) = ShedLoad.L(n,z);
OutputCurtailedPower(n,z) = CurtailedPower.L(n,z);
LostLoad_MaxPower(n,z)  = LL_MaxPower.L(n,z);
LostLoad_MinPower(n,z)  = LL_MinPower.L(n,z);
LostLoad_2D(n,z) = LL_2D.L(n,z);
LostLoad_2U(n,z) = LL_2U.L(n,z);
LostLoad_3U(n,z) = LL_3U.L(n,z);
ShadowPrice(n,z) = EQ_Demand_balance_DA.m(n,z);
LostLoad_WaterSlack(s) = WaterSlack.L(s);
StorageShadowPrice(s,z) = EQ_Storage_balance.m(s,z);


EXECUTE_UNLOAD "Results_simple.gdx"
OutputFlow,
OutputPower,
OutputStorageInput,
OutputStorageLevel,
OutputSystemCost,
OutputSpillage,
OutputShedLoad,
OutputCurtailedPower,
OutputGenMargin,
LostLoad_MaxPower,
LostLoad_MinPower,
LostLoad_2D,
LostLoad_2U,
LostLoad_3U,
ShadowPrice,
LostLoad_WaterSlack,
StorageShadowPrice,
status
;

$onorder
* Exit here if the PrintResult option is set to 0:
$if not %PrintResults%==1 $exit

EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=Technology rng=Technology!A1 rdim=2 cdim=0'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=PowerCapacity rng=PowerCapacity!A1 rdim=1 cdim=0'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=PowerInitial rng=PowerInitialA1 rdim=1 cdim=0'
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
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=OutputFlow rng=Flow!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=OutputPower rng=Power!A5 epsout=0 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=OutputStorageInput rng=StorageInput!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=OutputStorageLevel rng=StorageLevel!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=OutputSystemCost rng=SystemCost!A1 rdim=1 cdim=0'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=Y var=LostLoad_MaxPower rng=LostLoad_MaxPower!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=Y var=LostLoad_MinPower rng=LostLoad_MinPower!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=Y var=LostLoad_2D rng=LostLoad_2D!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=Y var=LostLoad_2U rng=LostLoad_2U!A1 rdim=1 cdim=1'

$exit

$Label DebugSection

$gdxin debug.gdx
$LOAD day
$LOAD PowerInitial_dbg
$LOAD StorageInitial_dbg
;
PowerInitial(u) = PowerInitial_dbg(u); StorageInitial(s) = StorageInitial_dbg(s);
FirstHour = (day-1)*24+1;
LastHour = min(card(h),FirstHour + (Config("RollingHorizon Length","day")+Config("RollingHorizon LookAhead","day")) * 24 - 1);
LastKeptHour = LastHour - Config("RollingHorizon LookAhead","day") * 24;
i(h) = no;
i(h)$(ord(h)>=firsthour and ord(h)<=lasthour)=yes;
StorageFinalMin(s) =  min(StorageInitial(s) + sum(i,StorageInflow(s,i)) - sum(i,StorageOutflow(s,i)) , sum(i$(ord(i)=card(i)),StorageProfile(s,i)*StorageCapacity(s)*AvailabilityFactor(s,i)));
StorageFinalMin(s) = min(StorageFinalMin(s),StorageCapacity(s) - smax(i,StorageInflow(s,i)));
$If %Verbose% == 1   Display PowerInitial,StorageFinalMin;
$If %LPFormulation% == 1          SOLVE UCM_SIMPLE USING LP MINIMIZING SystemCostD;
$If %LPFormulation% == 1          Display EQ_Objective_function.M, EQ_CostRampUp.M, EQ_CostRampDown.M, EQ_Demand_balance_DA.M, EQ_Power_available.M, EQ_Storage_minimum.M, EQ_Storage_level.M, EQ_Storage_input.M, EQ_Storage_balance.M, EQ_Storage_boundaries.M, EQ_Storage_MaxCharge.M, EQ_Storage_MaxDischarge.M, EQ_Flow_limits_lower.M ;

display day,FirstHour,LastHour,LastKeptHour;
Display StorageFinalMin,PowerInitial,StorageFinalMin;
Display Flow.L,Power.L,ShedLoad.L,StorageLevel.L,StorageInput.L,SystemCost.L,Spillage.L,StorageLevel.L,StorageInput.L,LL_MaxPower.L,LL_MinPower.L,LL_2U.L,LL_2D.L;
