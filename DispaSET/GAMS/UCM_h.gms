$Title UCM model

$eolcom //
Option threads=2;
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
* (1 for LP 0 for MIP)
$setglobal LPFormulation 0
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
s(u)             Hydro Storage Units (with reservoir)
chp(u)           CHP units
h                Hours
i(h)             Subset of simulated hours for one iteration
z(h)             Subset of all simulated hours
;

*$if %LookAhead%==1 SET i   TimeStep    /1*48/ ;
*$if %LookAhead%==0 SET i   TimeStep    /1*24/ ;


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
PARAMETERS
AvailabilityFactor(u,h)          [%]      Availability factor
CHPPowerLossFactor(u)          [%]      Power loss when generating heat
CHPPowerToHeat(u)              [%]      Nominal power-to-heat factor
CHPMaxHeat(chp)                [MW]     Maximum heat capacity of chp plant
CHPType
CommittedInitial(u)              [n.a.]   Initial committment status
Config
*CostCurtailment(n,h)             [EUR\MW]  Curtailment costs
CostFixed(u)                     [EUR\h]    Fixed costs
CostRampUp(u)                    [EUR\MW\h] Ramp-up costs
CostRampDown(u)                  [EUR\MW\h] Ramp-down costs
CostShutDown(u)                  [EUR]      Shut-down costs
CostStartUp(u)                   [EUR]      Start-up costs
CostVariable(u,h)                [EUR\MW]   Variable costs
CostHeatSlack(chp,h)             [EUR\MWh]  Cost of supplying heat via other means
CostLoadShedding(n,h)            [EUR\MWh] Cost of load shedding
Curtailment(n)                   [n.a]    Curtailment allowed or not {1 0} at node n
Demand(mk,n,h)                   [MW]     Demand
Efficiency(u)                    [%]      Efficiency
EmissionMaximum(n,p)             [tP]     Emission limit
EmissionRate(u,p)                [tP\MWh] P emission rate
FlowMaximum(l,h)                 [MW]     Line limits
FlowMinimum(l,h)                 [MW]     Minimum flow
FuelPrice(n,f,h)                 [EUR\F]    Fuel price
Fuel(u,f)                        [n.a.]   Fuel type {1 0}
HeatDemand(chp,h)                [MWh]    Heat demand profile for chp units
LineNode(l,n)                    [n.a.]   Incidence matrix {-1 +1}
LoadShedding(n,h)                [MW]   Load shedding capacity
Location(u,n)                    [n.a.]   Location {1 0}
Markup(u,h)                      [EUR\MW]   Markup
OutageFactor(u,h)                [%]      Outage Factor (100% = full outage)
PartLoadMin(u)                   [%]      Minimum part load
PowerCapacity(u)                 [MW]     Installed capacity
PowerInitial(u)                  [MW]     Power output before initial period
PowerMinStable(u)                [MW]     Minimum power output
PriceTransmission(l,h)           [EUR\MWh]  Transmission price
StorageChargingCapacity(u)       [MW]     Storage capacity
StorageChargingEfficiency(u)     [%]      Charging efficiency
StorageSelfDischarge(u)          [%\day]  Self-discharge of the storage units
RampDownMaximum(u)               [MW\h]   Ramp down limit
RampShutDownMaximum(u)           [MW\h]   Shut-down ramp limit
RampStartUpMaximum(u)            [MW\h]   Start-up ramp limit
RampUpMaximum(u)                 [MW\h]   Ramp up limit
Reserve(t)                       [n.a.]   Reserve technology {1 0}
StorageCapacity(u)               [MWh]    Storage capacity
StorageDischargeEfficiency(u)    [%]      Discharge efficiency
StorageOutflow(u,h)              [MWh]    Storage outflows
StorageInflow(u,h)               [MWh]    Storage inflows (potential energy)
StorageInitial(u)                [MWh]    Storage level before initial period
StorageProfile(u,h)              [MWh]    Storage level to be resepected at the end of each horizon
StorageMinimum(u)                [MWh]    Storage minimum
Technology(u,t)                  [n.a.]   Technology type {1 0}
TimeDown(u,h)                    [h]      Hours down
TimeDownLeft_initial(u)          [h]      Required time down left at the beginning of the simulated time period
TimeDownLeft_JustStopped(u,h)    [h]      Required time down left at hour h if the Unit has just been stopped
TimeDownInitial(u)               [h]      Hours down before initial period
TimeDownMinimum(u)               [h]      Minimum down time
TimeUpLeft_initial(u)            [h]      Required time up left at the beginning of the simulated time period
TimeUpInitial(u)                 [h]      Hours on before initial period
TimeUpMinimum(u)                 [h]      Minimum up time
FlexibilityUp(u)                 [MW\h]   Flexibility (up) of fast-starting power plants
FlexibilityDown(u)               [MW\h]   Flexibility (down) of a committed power plant
$If %RetrieveStatus% == 1 CommittedCalc(u,z)               [n.a.]   Committment status as for the MILP
;

*Parameters as used within the loop
PARAMETERS
TimeUpLeft_JustStarted(u,h)      [h]     Required time up left at hour h if the Unit has just been started
CostLoadShedding(n,h)            [EUR\MW]  Value of lost load
TimeUp(u,h)                      [h]     Hours up
LoadMaximum(u,h)                 [%]     Maximum load given AF and OF
PowerMustRun(u,h)                [MW]    Minimum power output
StorageFinalMin(s)             [MWh]   Minimum storage level at the end of the optimization horizon
;

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
$LOAD chp
$LOAD h
$LOAD z
*$LOAD d
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
$LOAD FuelPrice
$LOAD Fuel
$LOAD HeatDemand
$LOAD LineNode
$LOAD LoadShedding
$LOAD Location
$LOAD Markup
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
$LOAD TimeDownInitial
$LOAD TimeDownMinimum
$LOAD TimeUpInitial
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
FuelPrice,
Fuel,
HeatDemand,
HeatSlack,
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
TimeDownInitial,
TimeDownMinimum,
TimeUpInitial,
TimeUpMinimum
$If %RetrieveStatus% == 1 , CommittedCalc
;

$label skipdisplay

*===============================================================================
*Definition of variables
*===============================================================================
VARIABLES
Committed(u,h)             [n.a.]  Unit committed at hour h {1 0}
;

$If %LPFormulation% == 1 POSITIVE VARIABLES Committed (u,h) ; Committed.UP(u,h) = 1 ;
$If not %LPFormulation% == 1 BINARY VARIABLES Committed (u,h) ;

POSITIVE VARIABLES
CostStartUpH(u,h)          [EUR]   Cost of starting up
CostShutDownH(u,h)         [EUR]   cost of shutting down
CostRampUpH(u,h)           [EUR]   Ramping cost
CostRampDownH(u,h)         [EUR]   Ramping cost
CurtailedPower(n,h)        [MW]    Curtailed power at node n
Flow(l,h)                  [MW]    Flow through lines
MaxRamp2U(u,h)             [MW\h]  Maximum 15-min Ramp-up capbility
MaxRamp2D(u,h)             [MW\h]  Maximum 15-min Ramp-down capbility
Power(u,h)                 [MW]    Power output
PowerMaximum(u,h)          [MW]    Power output
PowerMinimum(u,h)          [MW]    Power output
ShedLoad(n,h)              [MW]    Shed load
StorageInput(u,h)          [MWh]   Charging input for storage units
StorageLevel(u,h)          [MWh]   Storage level of charge
LostLoad_MaxPower(n,h)     [MW]    Deficit in terms of maximum power
LostLoad_RampUp(u,h)       [MW]    Deficit in terms of ramping up for each plant
LostLoad_RampDown(u,h)     [MW]    Deficit in terms of ramping down
LostLoad_MinPower(n,h)     [MW]    Power exceeding the demand
LostLoad_Reserve2U(n,h)    [MW]    Deficit in reserve up
LostLoad_Reserve2D(n,h)    [MW]    Deficit in reserve down
spillage(s,h)              [MWh]   spillage from water reservoirs
SystemCost(h)              [EUR]   Hourly system cost
Heat(chp,h)                [MW]    Heat output by chp plant
HeatSlack(chp,h)           [MW]    Heat satisfied by other sources
;

free variable
SystemCostD               ![EUR]   Total system cost for one optimization period
;

*===============================================================================
*Assignment of initial values
*===============================================================================

*Forecasted upwards reserve margin (UCTE). Only if not provided in the parameters
Demand("2U",n,h)$(Demand("2U",n,h)=0)=sqrt(10*smax(hh,Demand("DA",n,hh))+150**2)-150;
*Forecasted downwards reserve margin (UCTE)
Demand("2D",n,h)$(Demand("2D",n,h)=0)=0.5*Demand("2U",n,h);

*Initial commitment status
CommittedInitial(u)=0;
CommittedInitial(u)$(PowerInitial(u)>0)=1;

* Definition of the minimum stable load:
PowerMinStable(u) = PartLoadMin(u)*PowerCapacity(u);

* Start-up and Shutdown ramping constraints. This remains to be solved
RampStartUpMaximum(u) = max(RampStartUpMaximum(u),PowerMinStable(u));
RampShutDownMaximum(u) = max(RampShutDownMaximum(u),PowerMinStable(u));

* If the plant is stopped, its 15-min ramp-up capability is RampStartUpMaximum if it can start in this timeframe:
FlexibilityUp(u) = RampStartUpMaximum(u)$(RampStartUpMaximum(u)>=PowerMinStable(u)*4);

* If the plant is started, its 15-min ramp-down capability is either RampShutDownMaximum if it is fast enough, or RampDownMaximum otherwise
*  RampDownMaximum(u)$(RampShutDownMaximum(u)<PowerMinStable(u)*4)
FlexibilityDown(u) = RampShutDownMaximum(u)$(RampShutDownMaximum(u)>=PowerMinStable(u)*4);

LoadMaximum(u,h)= AvailabilityFactor(u,h)*(1-OutageFactor(u,h));

PowerMustRun(u,h)=PowerMinStable(u)*LoadMaximum(u,h);
PowerMustRun(u,h)$(sum(tr,Technology(u,tr))>=1 and smin(n,Location(u,n)*(1-Curtailment(n)))=1) = PowerCapacity(u)*LoadMaximum(u,h);

$If %Verbose% == 1 Display RampStartUpMaximum, RampShutDownMaximum, CommittedInitial, FlexibilityUp, FlexibilityDown;

$offorder

*===============================================================================
*Declaration and definition of equations
*===============================================================================
EQUATIONS
EQ_Objective_function
EQ_CHP_extraction_Pmax
EQ_CHP_extraction
EQ_CHP_backpressure
EQ_CHP_P2H
EQ_CHP_demand_satisfaction
EQ_CHP_max_heat
EQ_Heat_Storage_balance
EQ_Heat_Storage_minimum
EQ_Heat_Storage_level
EQ_CostStartUp
EQ_CostShutDown
EQ_CostRampUp
EQ_CostRampDown
*EQ_CostRamping2
EQ_Demand_balance_DA
EQ_Demand_balance_2U
EQ_Demand_balance_2D
EQ_Power_must_run
EQ_Power_bound_lower
EQ_Power_bound_upper
EQ_Power_available
EQ_Ramp_up
EQ_Ramp_down
EQ_Max_RampUp1
EQ_Max_RampUp2
EQ_Max_RampDown1
EQ_Max_RampDown2
EQ_Minimum_time_up_A
EQ_Minimum_time_up_B
EQ_Minimum_time_up_C
EQ_Minimum_time_up_JustStarted
EQ_Minimum_time_down_A
EQ_Minimum_time_down_B
EQ_Minimum_time_down_C
EQ_Minimum_time_down_JustStopped
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
*EQ_Curtailment
EQ_LoadShedding
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
EQ_SystemCost(i)..
         SystemCost(i)
         =E=
         sum(u,CostFixed(u)*Committed(u,i))
         +sum(u,CostStartUpH(u,i) + CostShutDownH(u,i))
         +sum(u,CostRampUpH(u,i) + CostRampDownH(u,i))
         +sum(u,CostVariable(u,i) * Power(u,i))
         +sum(l,PriceTransmission(l,i)*Flow(l,i))
         +sum(n,CostLoadShedding(n,i)*ShedLoad(n,i))
         +sum(chp, CostHeatSlack(chp,i) * HeatSlack(chp,i))
         +sum(chp, CostVariable(chp,i) * CHPPowerLossFactor(chp) * Heat(chp,i))
         +100E3*(sum(n,LostLoad_MaxPower(n,i)+LostLoad_MinPower(n,i)))
         +80E3*(sum(n,LostLoad_Reserve2U(n,i)+LostLoad_Reserve2D(n,i)))
         +70E3*sum(u,LostLoad_RampUp(u,i)+LostLoad_RampDown(u,i))
         +1*sum(s,spillage(s,i))
;

EQ_Objective_function..
         SystemCostD
         =E=
         sum(i,SystemCost(i))
;

EQ_CostStartUp(u,i)$(CostStartUp(u) <> 0)..
         CostStartUpH(u,i)
         =G=
         CostStartUp(u)*(Committed(u,i)-CommittedInitial(u)$(ord(i) = 1)-Committed(u,i-1)$(ord(i) > 1))
;

EQ_CostShutDown(u,i)$(CostShutDown(u) <> 0)..
         CostShutDownH(u,i)
         =G=
         CostShutDown(u)*(CommittedInitial(u)$(ord(i) = 1)+Committed(u,i-1)$(ord(i) > 1)-Committed(u,i))
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

*EQ_CostRamping2(u,i)$(CostRampUp(u)=0 and CostRampDown(u)=0)..
*         CostRamping(u,i)
*         =e=
*         0
*;

*Hourly demand balance in the day-ahead market for each node
EQ_Demand_balance_DA(n,i)..
         sum(u,Power(u,i)*Location(u,n))
*         +sum(s,StorageOutputH(s,i)*Location(s,n))
          +sum(l,Flow(l,i)*LineNode(l,n))
         =E=
         Demand("DA",n,i)
         +sum(s,StorageInput(s,i)*Location(s,n))
         -ShedLoad(n,i)
         -LostLoad_MaxPower(n,i)
         +LostLoad_MinPower(n,i)
;

* Maximum 15-min ramping up, in MW/h:
Eq_Max_RampUp1(u,i)$(sum(tr,Technology(u,tr))=0)..
         MaxRamp2U(u,i)
         =L=
         RampUpMaximum(u)*Committed(u,i)
         + FlexibilityUp(u)*(1-Committed(u,i))
*         +RampStartUpMaximum(u)$(RampStartUpMaximum(u)>=PowerMinStable(u)*4)*(1-Committed(u,i))
;

* Maximum 15-min ramping up, in MW/h:
Eq_Max_RampUp2(u,i)$(sum(tr,Technology(u,tr))=0)..
         MaxRamp2U(u,i)
         =L=
         (PowerCapacity(u)*LoadMaximum(u,i) - Power(u,i))*4
;

* Maximum 15-min shutting down, in MW/h:
Eq_Max_RampDown1(u,i)$(sum(tr,Technology(u,tr))=0)..
         MaxRamp2D(u,i)
         =L=
         max(RampDownMaximum(u),FlexibilityDown(u))*Committed(u,i)
;

* Maximum 15-min ramping down, in MW/h:
Eq_Max_RampDown2(u,i)$(sum(tr,Technology(u,tr))=0)..
         MaxRamp2D(u,i)
         =L=
         (Power(u,i) - PowerMinStable(u)$(RampShutDownMaximum(u)<PowerMinStable(u)*4)*Committed(u,i))*4
;

EQ_Demand_balance_2U(n,i)..
         sum((u,t),MaxRamp2U(u,i)*Technology(u,t)*Reserve(t)*Location(u,n))
*         +CurtailedPowerH(n,i)*Curtailment(n,i)
         =G=
         +Demand("2U",n,i)
         -LostLoad_reserve2U(n,i)
;

*Hourly demand balance in the downwards reserve market for each node
EQ_Demand_balance_2D(n,i)..
         sum((u,t),MaxRamp2D(u,i)*Technology(u,t)*Reserve(t)*Location(u,n))
         =G=
         Demand("2D",n,i)
         -sum(s,(StorageChargingCapacity(s)-StorageInput(s,i)) )*4
         -LostLoad_reserve2D(n,i)
;


*Minimum power output is above the must-run output level for each unit in all periods
EQ_Power_must_run(u,i)..
         PowerMustRun(u,i) * Committed(u,i) - (StorageInput(u,i) * CHPPowerLossFactor(u) )$(chp(u) and CHPType(u,'Extraction'))
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

*Maximum power output with respect to power output in the previous period (ramping up constraint).
EQ_Ramp_up(u,i)$(sum(tr,Technology(u,tr))=0)..
         Power(u,i)
         =L=
         (PowerInitial(u)
         +RampUpMaximum(u)
                 *CommittedInitial(u)
         +RampStartUpMaximum(u)
                 *(1-CommittedInitial(u)))$(ord(i) = 1)
         +(Power(u,i-1)
         +RampUpMaximum(u)
                 *Committed(u,i-1)
         +RampStartUpMaximum(u)
                 *(1-Committed(u,i-1)))$(ord(i) > 1)
         +LostLoad_RampUp(u,i)
;

*If the unit keeps committed the reduction in power output is lower than the
*ramp-down limit. If the unit is de-committed the reduction is lower than the
*shut-down ramp limit
EQ_Ramp_down(u,i)$(sum(tr,Technology(u,tr))=0)..
         (PowerInitial(u)-Power(u,i))$(ord(i) = 1)
         +(Power(u,i-1)-Power(u,i))$(ord(i) > 1)
         =L=
         RampDownMaximum(u) * Committed(u,i)
         + RampShutDownMaximum(u) * (1-Committed(u,i))
         + LostLoad_RampDown(u,i)
;

*Minimum time up constraints
*EQ_Minimum_time_up(u)..
*         Committed(u,i)
*         =G=
*         sum(ii$((ord(ii) >= ord(i)+1-TimeUpMinimum(u)) and (ord(ii))<=ord(i))),Committed(u,i)-CommittedH(h-1,u))
*;

EQ_Minimum_time_up_A(u)..
         sum(i$(ord(i) <= TimeUpLeft_initial(u)),1-Committed(u,i))
         =E=
         0
;

EQ_Minimum_time_up_B(u,i)$((TimeUpLeft_initial(u)+1 <= ord(i)) and (ord(i) <= card(i)-TimeUpMinimum(u)+1))..
         sum(ii$((ord(i) <= ord(ii)) and (ord(ii) <= (ord(i)+TimeUpMinimum(u)-1))),Committed(u,ii))
         =G=
         TimeUpMinimum(u)
                 *(Committed(u,i)-CommittedInitial(u)$(ord(i) = 1)-Committed(u,i-1)$(ord(i) > 1))
;

EQ_Minimum_time_up_C(u,i)$((card(i)-TimeUpMinimum(u)+2 <= ord(i)) and (ord(i)<=card(i)))..
         sum(ii$((ord(i) <= ord(ii)) and (ord(ii) <= card(i))),Committed(u,ii)-(Committed(u,i)-CommittedInitial(u)$(ord(i) = 1)-Committed(u,i-1)$(ord(i) > 1)))
         =G=
         0
;

EQ_Minimum_time_up_JustStarted(u,i)$(ord(i) > 1)..
         sum(ii$((ord(i) <= ord(ii)) and (ord(ii) <= (ord(i)+TimeUpLeft_JustStarted(u,i)-1))),Committed(u,ii))
         =G=
         TimeUpLeft_JustStarted(u,i)
                 *(Committed(u,i)-CommittedInitial(u)$(ord(i) = 1)-Committed(u,i-1)$(ord(i) > 1))
;

*Minimum time down constraints
*EQ_Minimum_time_down(u)..
*         1-Committed(u,i)
*         =G=
*         sum(ii$((ord(ii) >= ord(h)+1-TimeDownMinimum(u)) and (ord(ii))<=ord(h))),CommittedH(h-1,u)-Committed(u,i))
*;

EQ_Minimum_time_down_A(u)..
         sum(i$(ord(i) <= TimeDownLeft_initial(u)),Committed(u,i))
         =E=
         0
;

EQ_Minimum_time_down_B(u,i)$((TimeDownLeft_initial(u)+1 <= ord(i)) and (ord(i) <= card(i)-TimeDownMinimum(u)+1))..
         sum(ii$((ord(i) <= ord(ii)) and (ord(ii) <= (ord(i)+TimeDownMinimum(u)-1))),1-Committed(u,ii))
         =G=
         TimeDownMinimum(u)
                 *(CommittedInitial(u)$(ord(i) = 1)+Committed(u,i-1)$(ord(i) > 1)-Committed(u,i))
;

EQ_Minimum_time_down_C(u,i)$((card(i)-TimeDownMinimum(u)+2 <= ord(i)) and (ord(i)<=card(i)))..
         sum(ii$((ord(i) <= ord(ii)) and (ord(ii) <= card(i))),1-Committed(u,ii)-(CommittedInitial(u)$(ord(i) = 1)+Committed(u,i-1)$(ord(i) > 1)-Committed(u,i)))
         =G=
         0
;

*IH: why do we need this equation?, to replace the two previous?
EQ_Minimum_time_down_JustStopped(u,i)$(TimeDownLeft_initial(u)+1 <= ord(i))..
         sum(ii$((ord(i) <= ord(ii)) and (ord(ii) <= (ord(i)+TimeDownLeft_JustStopped(u,i)-1))),1-Committed(u,ii))
         =G=
         TimeDownLeft_JustStopped(u,i)
                 *(CommittedInitial(u)$(ord(i) = 1)+Committed(u,i-1)$(ord(i) > 1)-Committed(u,i))
;

*Storage level must be above a minimum
EQ_Storage_minimum(s,i)..
         StorageMinimum(s)
         =L=
         StorageLevel(s,i)
;

*Storage level must below storage capacity
EQ_Storage_level(s,i)..
         StorageLevel(s,i)
         =L=
         StorageCapacity(s)*AvailabilityFactor(s,i)
;

* Storage charging is bounded by the maximum capacity
EQ_Storage_input(s,i)..
         StorageInput(s,i)
         =L=
         StorageChargingCapacity(s)*(1-Committed(s,i))
;

*Discharge is limited by the storage level
EQ_Storage_MaxDischarge(s,i)..
         Power(s,i)/(max(StorageDischargeEfficiency(s),0.0001))
         +StorageOutflow(s,i) +Spillage(s,i) - StorageInflow(s,i)
         =L=
         StorageLevel(s,i)
;

*Charging is limited by the remaining storage capacity
EQ_Storage_MaxCharge(s,i)..
         StorageInput(s,i) * StorageChargingEfficiency(s)
         -StorageOutflow(s,i) -spillage(s,i) + StorageInflow(s,i)
         =L=
         StorageCapacity(s)*AvailabilityFactor(s,i) - StorageLevel(s,i)
;

*Storage balance
EQ_Storage_balance(s,i)..
         StorageInitial(s)$(ord(i) = 1)
         +StorageLevel(s,i-1)$(ord(i) > 1)
*         +StorageLevelH(h--1,s)
         +StorageInflow(s,i)
         +StorageInput(s,i)*StorageChargingEfficiency(s)
         =E=
         StorageLevel(s,i)
         +StorageOutflow(s,i)
         +spillage(s,i)
         +Power(s,i)/(max(StorageDischargeEfficiency(s),0.0001))
;

* Minimum level at the end of the optimization horizon:
EQ_Storage_boundaries(s,i)$(ord(i) = card(i))..
         StorageFinalMin(s)
         =L=
         StorageLevel(s,i)
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

*Force Unit commitment/decommitment:
* E.g: renewable units with AF>0 must be committed
EQ_Force_Commitment(u,i)$((sum(tr,Technology(u,tr))>=1 and LoadMaximum(u,i)>0) or (ord(i)=4 and ord(u)=129))..
         Committed(u,i)
         =E=
         1;

* E.g: renewable units with AF=0 must be decommitted
EQ_Force_DeCommitment(u,i)$(LoadMaximum(u,i)=0 or ord(u)=200)..
         Committed(u,i)
         =E=
         0;

*Curtailment allowed or not
*EQ_Curtailment(n,i)..
*         CurtailedPowerH(n,i)
*         =L=
*         sum(u,Power(u,i)*sum(tr,Technology(u,tr))*Location(u,n))
*            *Curtailment(n)
*;

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

EQ_CHP_extraction_Pmax(chp,i)$(CHPType(chp,'Extraction'))..
         Power(chp,i)
         =L=
         PowerCapacity(chp)  - StorageInput(chp,i) * CHPPowerLossFactor(chp)
;

EQ_CHP_backpressure(chp,i)$(CHPType(chp,'Back-Pressure'))..
         Power(chp,i)
         =E=
         StorageInput(chp,i) * CHPPowerToHeat(chp)
;

EQ_CHP_P2H(chp,i)$(CHPType(chp,'P2H'))..
         Power(chp,i)
         =E=
         PowerCapacity(chp) - StorageInput(chp,i) * CHPPowerLossFactor(chp)
;

EQ_CHP_max_heat(chp,i)..
         StorageInput(chp,i)
         =L=
         CHPMaxHeat(chp)
;

EQ_CHP_demand_satisfaction(chp,i)..
         Heat(chp,i) + HeatSlack(chp,i)
         =E=
         HeatDemand(chp,i)
;

*Heat Storage balance
EQ_Heat_Storage_balance(chp,i)..
          StorageInitial(chp)$(ord(i) = 1)
         +StorageLevel(chp,i-1)$(ord(i) > 1)
         +StorageInput(chp,i)
         =E=
         StorageLevel(chp,i)
         +Heat(chp,i) + StorageSelfDischarge(chp) * StorageLevel(chp,i)/24
;
* The self-discharge proportional to the charging level is a bold hypothesis, but it avoids keeping self-discharging if the level reaches zero

*Storage level must be above a minimum
EQ_Heat_Storage_minimum(chp,i)..
         StorageMinimum(chp)
         =L=
         StorageLevel(chp,i)
;

*Storage level must below storage capacity
EQ_Heat_Storage_level(chp,i)..
         StorageLevel(chp,i)
         =L=
         StorageCapacity(chp)
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
EQ_CHP_P2H,
EQ_CHP_demand_satisfaction,
EQ_CHP_max_heat,
$If not %LPFormulation% == 1 EQ_CostStartUp,
$If not %LPFormulation% == 1 EQ_CostShutDown,
$If %LPFormulation% == 1 EQ_CostRampUp,
$If %LPFormulation% == 1 EQ_CostRampDown,
EQ_Demand_balance_DA,
EQ_Demand_balance_2U,
EQ_Demand_balance_2D,
$If not %LPFormulation% == 1 EQ_Power_must_run,
EQ_Power_available,
EQ_Heat_Storage_balance,
EQ_Heat_Storage_minimum,
EQ_Heat_Storage_level,
EQ_Ramp_up,
EQ_Ramp_down,
*EQ_Minimum_time_up,
$If not %LPFormulation% == 1 EQ_Minimum_time_up_A,
*EQ_Minimum_time_up_B,
*EQ_Minimum_time_up_C,
$If not %LPFormulation% == 1 EQ_Minimum_time_up_JustStarted,
*EQ_Minimum_time_down,
$If not %LPFormulation% == 1 EQ_Minimum_time_down_A,
*EQ_Minimum_time_down_B,
*EQ_Minimum_time_down_C,
$If not %LPFormulation% == 1 EQ_Minimum_time_down_JustStopped,
EQ_Max_RampUp1,
EQ_Max_RampUp2,
EQ_Max_RampDown1,
EQ_Max_RampDown2,
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
EQ_Force_Commitment,
EQ_Force_DeCommitment,
*EQ_Curtailment,
EQ_LoadShedding,
$If %RetrieveStatus% == 1 EQ_CommittedCalc
/
;
UCM_SIMPLE.optcr = 0.01;
*UCM_SIMPLE.epgap = 0.005
*UCM_SIMPLE.probe = 3
*UCM_SIMPLE.optfile=1;

*===============================================================================
*Solving loop
*===============================================================================

* Scalar variables necessary to the loop:
scalar FirstHour,LastHour,LastKeptHour,day,ndays,failed;
ndays = floor(card(h)/24);
if (Config("RollingHorizon LookAhead","day") > ndays -1, abort "The look ahead period is longer than the simulation length";);

* Some parameters used for debugging:
failed=0;
parameter TimeUpInitial_dbg(u), TimeDownInitial_dbg(u), CommittedInitial_dbg(u), PowerInitial_dbg(u), StorageInitial_dbg(s);

* Fixing the initial guesses:
*PowerH.L(u,i)=PowerInitial(u);
*Committed.L(u,i)=CommittedInitial(u);

* Defining a parameter that records the solver status:
set  tmp   "tpm"  / "model", "solver" /  ;
parameter status(tmp,h);

$if %Debug% == 1 $goto DebugSection

FOR(day = 1 TO ndays-Config("RollingHorizon LookAhead","day") by Config("RollingHorizon Length","day"),
         FirstHour = (day-1)*24+1;
         LastHour = min(card(h),FirstHour + (Config("RollingHorizon Length","day")+Config("RollingHorizon LookAhead","day")) * 24 - 1);
         LastKeptHour = LastHour - Config("RollingHorizon LookAhead","day") * 24;
         i(h) = no;
         i(h)$(ord(h)>=firsthour and ord(h)<=lasthour)=yes;
         display day,FirstHour,LastHour,LastKeptHour;

         TimeUpLeft_initial(u)=min(card(i),(TimeUpMinimum(u)-TimeUpInitial(u))*CommittedInitial(u));
         TimeUpLeft_JustStarted(u,i) = min(card(i)-ord(i)+1,TimeUpMinimum(u));
         TimeDownLeft_initial(u)=min(card(i),(TimeDownMinimum(u)-TimeDownInitial(u))*(1-CommittedInitial(u)));
         TimeDownLeft_JustStopped(u,i) = min(card(i)-ord(i)+1,TimeDownMinimum(u));

*        Defining the minimum level at the end of the horizon, ensuring that it is feasible with the provided inflows:
         StorageFinalMin(s) =  min(StorageInitial(s) + sum(i,StorageInflow(s,i)) - sum(i,StorageOutflow(s,i)) , sum(i$(ord(i)=card(i)),StorageProfile(s,i)*StorageCapacity(s)*AvailabilityFactor(s,i)));
*        Correcting the minimum level to avoid the infeasibility in case it is too close to the StorageCapacity:
         StorageFinalMin(s) = min(StorageFinalMin(s),StorageCapacity(s) - smax(i,StorageInflow(s,i)));

$If %Verbose% == 1   Display TimeUpLeft_initial,TimeUpLeft_JustStarted,TimeDownLeft_initial,TimeDownLeft_JustStopped,TimeUpInitial,TimeDownInitial,PowerInitial,CommittedInitial,StorageFinalMin;

$If %LPFormulation% == 1          SOLVE UCM_SIMPLE USING LP MINIMIZING SystemCostD;
$If not %LPFormulation% == 1      SOLVE UCM_SIMPLE USING MIP MINIMIZING SystemCostD;

$If %LPFormulation% == 1          Display EQ_Objective_function.M, EQ_CostRampUp.M, EQ_CostRampDown.M, EQ_Demand_balance_DA.M, EQ_Power_available.M, EQ_Ramp_up.M, EQ_Ramp_down.M, EQ_Max_RampUp1.M, EQ_Max_RampUp2.M,EQ_Max_RampDown1.M, EQ_Max_RampDown2.M, EQ_Storage_minimum.M, EQ_Storage_level.M, EQ_Storage_input.M, EQ_Storage_balance.M, EQ_Storage_boundaries.M, EQ_Storage_MaxCharge.M, EQ_Storage_MaxDischarge.M, EQ_Flow_limits_lower.M ;
$If not %LPFormulation% == 1      Display EQ_Objective_function.M, EQ_CostStartUp.M, EQ_CostShutDown.M, EQ_Demand_balance_DA.M, EQ_Power_must_run.M, EQ_Power_available.M, EQ_Ramp_up.M, EQ_Ramp_down.M, EQ_Minimum_time_up_A.M, EQ_Minimum_time_up_JustStarted.M, EQ_Minimum_time_down_A.M, EQ_Minimum_time_down_JustStopped.M, EQ_Max_RampUp1.M, EQ_Max_RampUp2.M, EQ_Max_RampDown1.M, EQ_Max_RampDown2.M, EQ_Storage_minimum.M, EQ_Storage_level.M, EQ_Storage_input.M, EQ_Storage_balance.M, EQ_Storage_boundaries.M, EQ_Storage_MaxCharge.M, EQ_Storage_MaxDischarge.M, EQ_Flow_limits_lower.M ;

         status("model",i) = UCM_SIMPLE.Modelstat;
         status("solver",i) = UCM_SIMPLE.Solvestat;

if(UCM_SIMPLE.Modelstat <> 1 and UCM_SIMPLE.Modelstat <> 8 and not failed, TimeUpInitial_dbg(u) = TimeUpInitial(u); TimeDownInitial_dbg(u) = TimeDownInitial(u); CommittedInitial_dbg(u) = CommittedInitial(u); PowerInitial_dbg(u) = PowerInitial(u); StorageInitial_dbg(s) = StorageInitial(s);
                                                                           EXECUTE_UNLOAD "debug.gdx" day, status, TimeUpInitial_dbg, TimeDownInitial_dbg, CommittedInitial_dbg, PowerInitial_dbg, StorageInitial_dbg;
                                                                           failed=1;);

*Time counters
         Loop(i,
              TimeUp(u,i)$(ord(i) = 1 and Committed.L(u,i) = 1)=TimeUpInitial(u)+1;
              TimeUp(u,i)$(ord(i) = 1 and Committed.L(u,i) = 0)=0;
              TimeUp(u,i)$(ord(i) > 1 and Committed.L(u,i) = 1) = TimeUp(u,i-1)+1;
              TimeUp(u,i)$(ord(i) > 1 and Committed.L(u,i) = 0) = 0;

              TimeDown(u,i)$(ord(i) = 1 and Committed.L(u,i) = 0) = TimeDownInitial(u)+1;
              TimeDown(u,i)$(ord(i) = 1 and Committed.L(u,i) = 1) = 0;
              TimeDown(u,i)$(ord(i) > 1 and Committed.L(u,i) = 0) = TimeDown(u,i-1)+1;
              TimeDown(u,i)$(ord(i) > 1 and Committed.L(u,i) = 1) = 0;
              );

         TimeUpInitial(u)=sum(i$(ord(i)=LastKeptHour-FirstHour+1),TimeUp(u,i));
         TimeDownInitial(u)=sum(i$(ord(i)=LastKeptHour-FirstHour+1),TimeDown(u,i));
         CommittedInitial(u)=sum(i$(ord(i)=LastKeptHour-FirstHour+1),Committed.L(u,i));
         PowerInitial(u) = sum(i$(ord(i)=LastKeptHour-FirstHour+1),Power.L(u,i));

         StorageInitial(s) =   sum(i$(ord(i)=LastKeptHour-FirstHour+1),StorageLevel.L(s,i));
         StorageInitial(chp) =   sum(i$(ord(i)=LastKeptHour-FirstHour+1),StorageLevel.L(chp,i));


*Loop variables to display after solving:
$If %Verbose% == 1 Display LastKeptHour,PowerInitial,TimeUp,TimeDown,MaxRamp2D.L,MaxRamp2U.L,CostStartUpH.L,CostShutDownH.L,CostRampUpH.L;

);

CurtailedPower.L(n,z)=sum(u,(PowerCapacity(u)*LoadMaximum(u,z)-Power.L(u,z))$(sum(tr,Technology(u,tr))>=1) * Location(u,n));

$If %Verbose% == 1 Display Flow.L,Power.L,Committed.L,ShedLoad.L,CurtailedPower.L,StorageLevel.L,StorageInput.L,SystemCost.L,MaxRamp2U.L,MaxRamp2D.L,LostLoad_MaxPower.L,LostLoad_MinPower.L,LostLoad_reserve2U.L,LostLoad_reserve2D.L,LostLoad_RampUP.L,LostLoad_RampDown.L;


*===============================================================================
*Result export
*===============================================================================

PARAMETER
OutputCommitted(u,h)
OutputHeat(chp,h)
OutputFlow(l,h)
OutputPower(u,h)
OutputStorageInput(u,h)
OutputStorageLevel(u,h)
*OutputTimeDown(u,i)
*OutputTimeUp(u,i)
OutputSystemCost(h)
OutputSpillage(s,h)
OutputShedLoad(n,h)
OutputCurtailedPower(n,h)
OutputHeat(chp,h)
OutputHeatSlack(chp,h)
ShadowPrice(n,h)
;

OutputCommitted(u,z)=Committed.L(u,z);
OutputFlow(l,z)=Flow.L(l,z);
OutputPower(u,z)=Power.L(u,z);
OutputHeat(chp,z)=Heat.L(chp,z);
OutputHeatSlack(chp,z)=HeatSlack.L(chp,z);
OutputStorageInput(s,z)=StorageInput.L(s,z);
OutputStorageInput(chp,z)=StorageInput.L(chp,z);
OutputStorageLevel(s,z)=StorageLevel.L(s,z);
OutputStorageLevel(chp,z)=StorageLevel.L(chp,z);
*OutputTimeDown(u,i)=TimeDown(u,i);
*OutputTimeUp(u,h)=TimeUp(u,i);
OutputSystemCost(z)=SystemCost.L(z);
OutputSpillage(s,z)  = Spillage.L(s,z) ;
OutputShedLoad(n,z) = ShedLoad.L(n,z);
OutputCurtailedPower(n,z)=CurtailedPower.L(n,z);
ShadowPrice(n,z) = EQ_Demand_balance_DA.m(n,z);

EXECUTE_UNLOAD "Results.gdx"
OutputCommitted,
OutputFlow,
OutputPower,
OutputHeat,
OutputHeatSlack,
OutputStorageInput,
OutputStorageLevel,
*OutputTimeDown,
*OutputTimeUp,
OutputSystemCost,
OutputSpillage,
OutputShedLoad,
OutputCurtailedPower,
LostLoad_MaxPower,
LostLoad_MinPower,
LostLoad_Reserve2D,
LostLoad_Reserve2U,
LostLoad_RampUp,
LostLoad_RampDown,
ShadowPrice,
status
;

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
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=TimeUpInitial rng=TimeUpInitial!A1 rdim=1 cdim=0'
EXECUTE 'GDXXRW.EXE "%inputfilename%" O="Results.xlsx" Squeeze=N par=TimeDownInitial rng=TimeDownInitial!A1 rdim=1 cdim=0'
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

*EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=DayOutput rng=Power!B3 cdim=1'
*EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=MonthOutput rng=Power!B2 cdim=1'
*EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=YearOutput rng=Power!B1 cdim=1'

EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N var=CurtailedPower rng=CurtailedPower!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N var=ShedLoad rng=ShedLoad!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=OutputCommitted rng=Committed!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=OutputFlow rng=Flow!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=OutputPower rng=Power!A5 epsout=0 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=OutputStorageInput rng=StorageInput!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=OutputStorageLevel rng=StorageLevel!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=OutputSystemCost rng=SystemCost!A1 rdim=1 cdim=0'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=Y var=LostLoad_MaxPower rng=LostLoad_MaxPower!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=Y var=LostLoad_MinPower rng=LostLoad_MinPower!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=Y var=LostLoad_Reserve2D rng=LostLoad_Reserve2D!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=Y var=LostLoad_Reserve2U rng=LostLoad_Reserve2U!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=Y var=LostLoad_RampUp rng=LostLoad_RampUp!A1 rdim=1 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=Y var=LostLoad_RampDown rng=LostLoad_RampDown!A1 rdim=1 cdim=1'



$exit

$Label DebugSection

$gdxin debug.gdx
$LOAD day
$LOAD PowerInitial_dbg
$LOAD CommittedInitial_dbg
$LOAD StorageInitial_dbg
$LOAD TimeDownInitial_dbg
$LOAD TimeUpInitial_dbg
;
PowerInitial(u) = PowerInitial_dbg(u); CommittedInitial(u) = CommittedInitial_dbg(u); StorageInitial(s) = StorageInitial_dbg(s); TimeDownInitial(u) = TimeDownInitial_dbg(u); TimeUpInitial(u) = TimeUpInitial_dbg(u);
FirstHour = (day-1)*24+1;
LastHour = min(card(h),FirstHour + (Config("RollingHorizon Length","day")+Config("RollingHorizon LookAhead","day")) * 24 - 1);
LastKeptHour = LastHour - Config("RollingHorizon LookAhead","day") * 24;
i(h) = no;
i(h)$(ord(h)>=firsthour and ord(h)<=lasthour)=yes;
TimeUpLeft_initial(u)=min(card(i),(TimeUpMinimum(u)-TimeUpInitial(u))*CommittedInitial(u));
TimeUpLeft_JustStarted(u,i) = min(card(i)-ord(i)+1,TimeUpMinimum(u));
TimeDownLeft_initial(u)=min(card(i),(TimeDownMinimum(u)-TimeDownInitial(u))*(1-CommittedInitial(u)));
TimeDownLeft_JustStopped(u,i) = min(card(i)-ord(i)+1,TimeDownMinimum(u));
StorageFinalMin(s) =  min(StorageInitial(s) + sum(i,StorageInflow(s,i)) - sum(i,StorageOutflow(s,i)) , sum(i$(ord(i)=card(i)),StorageProfile(s,i)*StorageCapacity(s)*AvailabilityFactor(s,i)));
StorageFinalMin(s) = min(StorageFinalMin(s),StorageCapacity(s) - smax(i,StorageInflow(s,i)));
$If %Verbose% == 1   Display TimeUpLeft_initial,TimeUpLeft_JustStarted,TimeDownLeft_initial,TimeDownLeft_JustStopped,TimeUpInitial,TimeDownInitial,PowerInitial,CommittedInitial,StorageFinalMin;
$If %LPFormulation% == 1          SOLVE UCM_SIMPLE USING LP MINIMIZING SystemCostD;
$If not %LPFormulation% == 1      SOLVE UCM_SIMPLE USING MIP MINIMIZING SystemCostD;
$If %LPFormulation% == 1          Display EQ_Objective_function.M, EQ_CostRampUp.M, EQ_CostRampDown.M, EQ_Demand_balance_DA.M, EQ_Power_available.M, EQ_Ramp_up.M, EQ_Ramp_down.M, EQ_Max_RampUp1.M, EQ_Max_RampUp2.M,EQ_Max_RampDown1.M, EQ_Max_RampDown2.M, EQ_Storage_minimum.M, EQ_Storage_level.M, EQ_Storage_input.M, EQ_Storage_balance.M, EQ_Storage_boundaries.M, EQ_Storage_MaxCharge.M, EQ_Storage_MaxDischarge.M, EQ_Flow_limits_lower.M ;
$If not %LPFormulation% == 1      Display EQ_Objective_function.M, EQ_CostStartUp.M, EQ_CostShutDown.M, EQ_Demand_balance_DA.M, EQ_Power_must_run.M, EQ_Power_available.M, EQ_Ramp_up.M, EQ_Ramp_down.M, EQ_Minimum_time_up_A.M, EQ_Minimum_time_up_JustStarted.M, EQ_Minimum_time_down_A.M, EQ_Minimum_time_down_JustStopped.M, EQ_Max_RampUp1.M, EQ_Max_RampUp2.M, EQ_Max_RampDown1.M, EQ_Max_RampDown2.M, EQ_Storage_minimum.M, EQ_Storage_level.M, EQ_Storage_input.M, EQ_Storage_balance.M, EQ_Storage_boundaries.M, EQ_Storage_MaxCharge.M, EQ_Storage_MaxDischarge.M, EQ_Flow_limits_lower.M ;

display day,FirstHour,LastHour,LastKeptHour;
Display StorageFinalMin,TimeUpLeft_initial,TimeUpLeft_JustStarted,TimeDownLeft_initial,TimeDownLeft_JustStopped,TimeUpInitial,TimeDownInitial,PowerInitial,CommittedInitial,StorageFinalMin;
Display Flow.L,Power.L,Committed.L,ShedLoad.L,StorageLevel.L,StorageInput.L,SystemCost.L,MaxRamp2U.L,MaxRamp2D.L,Spillage.L,StorageLevel.L,StorageInput.L,LostLoad_MaxPower.L,LostLoad_MinPower.L,LostLoad_reserve2U.L,LostLoad_reserve2D.L,LostLoad_RampUP.L,LostLoad_RampDown.L;
