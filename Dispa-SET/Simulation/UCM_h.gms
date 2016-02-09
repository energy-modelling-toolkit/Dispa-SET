$Title UCM model
*Option LimRow=100000000;
Option IterLim=100000000;
Option ResLim = 10000000000;
*Option SolPrint=On;
Option SolPrint=Silent;



*===============================================================================
*Definition of the dataset-related options
*===============================================================================

* Print results to excel files (0 for no, 1 for yes)
$set PrintResults 1

* Name of the input file (Ideally, stick to the default Input.gdx)
*$set InputFileName Input.gdx
$set InputFileName Inputs.gdx

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
AvailabilityFactor(u,h)          [%]     Availability factor
CommittedInitial(u)              [n.a.]  Initial committment status
Config
*CostCurtailment(n,h)             [€\MW]  Curtailment costs
CostFixed(u)                     [€\h]   Fixed costs
CostRampUp(u)                    [€\MW\h]Ramp-up costs
CostRampDown(u)                  [€\MW\h]Ramp-down costs
CostShutDown(u)                  [€]     Shut-down costs
CostStartUp(u)                   [€]     Start-up costs
CostVariable(u,h)              [€\MW]  Variable costs
Curtailment(n)                   [n.a]   Curtailment allowed or not {1 0} at node n
Demand(mk,n,h)                   [MW]    Demand
Efficiency(u)                    [%]     Efficiency
EmissionMaximum(n,p)             [tP]    Emission limit
EmissionRate(u,p)                [tP\MW] P emission rate
FlowMaximum(l,h)                 [MW]    Line limits
FlowMinimum(l,h)                 [MW]    Minimum flow
FuelPrice(n,f,h)                 [€\F]   Fuel price
Fuel(u,f)                        [n.a.]  Fuel type {1 0}
LineNode(l,n)                    [n.a.]  Incidence matrix {-1 +1}
LoadShedding(n)                  [n.a.]  Load shedding capacity
Location(u,n)                    [n.a.]  Location {1 0}
Markup(u,h)                      [€\MW]  Markup
OutageFactor(u,h)                [%]     Outage Factor (100% = full outage)
PartLoadMin(u)                   [%]     Minimum part load
PermitPrice(p)                   [€\tP]  Permit price
PowerCapacity(u)                 [MW]    Installed capacity
PowerInitial(u)                  [MW]    Power output before initial period
PowerMinStable(u)                [MW]    Minimum power output
PriceTransmission(l,h)           [€\MWh] Transmission price
StorageChargingCapacity(u)        [MW]   Storage capacity
StorageChargingEfficiency(u)      [%]    Charging efficiency
RampDownMaximum(u)               [MW\h]  Ramp down limit
RampShutDownMaximum(u)           [MW\h]  Shut-down ramp limit
RampStartUpMaximum(u)            [MW\h]  Start-up ramp limit
RampUpMaximum(u)                 [MW\h]  Ramp up limit
Reserve(t)                       [n.a.]  Reserve technology {1 0}
StorageCapacity(u)               [MWh] Storage capacity
StorageDischargeEfficiency(u)    [%]     Discharge efficiency
StorageOutflow(u,h)              [MWh]  Storage outflows
StorageInflow(u,h)               [MWh]  Storage inflows (potential energy)
StorageInitial(u)                [MWh] Storage level before initial period
StorageMinimum(u)                [MWh] Storage minimum
Technology(u,t)                  [n.a.]  Technology type {1 0}
TimeDown(u,h)                    [h]     Hours down
TimeDownLeft_initial(u)          [h]     Required time down left at the beginning of the simulated time period
TimeDownLeft_JustStopped(u,h)    [h]     Required time down left at hour h if the Unit has just been stopped
TimeDownInitial(u)               [h]     Hours down before initial period
TimeDownMinimum(u)               [h]     Minimum down time
TimeUpLeft_initial(u)            [h]     Required time up left at the beginning of the simulated time period
TimeUpInitial(u)                 [h]     Hours on before initial period
TimeUpMinimum(u)                 [h]     Minimum up time
FlexibilityUp(u)                [MW\h]  Flexibility (up) of fast-starting power plants
FlexibilityDown(u)              [MW\h]  Flexibility (down) of a committed power plant
;

*Parameters as used within the loop
PARAMETERS
TimeUpLeft_JustStarted(u,h)      [h]     Required time up left at hour h if the Unit has just been started
CostLoadShedding(n,h)            [€\MW]  Value of lost load
TimeUp(u,h)                      [h]     Hours up
LoadMaximum(u,h)                 [%]     Maximum load given AF and OF
PowerMustRun(u,h)                [MW]    Minimum power output
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
$LOAD h
*$LOAD d
$LOAD AvailabilityFactor
$LOAD Config
$LOAD CostFixed
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
$LOAD LineNode
$LOAD LoadShedding
$LOAD Location
$LOAD Markup
$LOAD OutageFactor
$LOAD PermitPrice
$LOAD PowerCapacity
$LOAD PowerInitial
$LOAD PartLoadMin
$LOAD PriceTransmission
$LOAD StorageChargingCapacity
$LOAD StorageChargingEfficiency
$LOAD RampDownMaximum
$LOAD RampShutDownMaximum
$LOAD RampStartUpMaximum
$LOAD RampUpMaximum
$LOAD Reserve
$LOAD StorageCapacity
$LOAD StorageInflow
$LOAD StorageInitial
$LOAD StorageMinimum
$LOAD StorageOutflow
$LOAD Technology
$LOAD TimeDownInitial
$LOAD TimeDownMinimum
$LOAD TimeUpInitial
$LOAD TimeUpMinimum
;

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
CostShutDown,
CostStartUp,
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
LineNode,
Location,
LoadShedding
Markup,
OutageFactor,
PartLoadMin,
PermitPrice,
PowerCapacity,
PowerInitial,
PriceTransmission,
StorageChargingCapacity,
StorageChargingEfficiency,
RampDownMaximum,
RampShutDownMaximum,
RampStartUpMaximum,
RampUpMaximum,
Reserve,
StorageCapacity,
StorageInflow,
StorageInitial,
StorageMinimum,
StorageOutflow,
Technology,
TimeDownInitial,
TimeDownMinimum,
TimeUpInitial
TimeUpMinimum
;

*===============================================================================
*Definition of variables
*===============================================================================
BINARY VARIABLES
Committed(u,h)             [n.a.]  Unit committed at hour h {1 0}
;

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
StorageInput(s,h)          [MWh]   Charging input for storage units
StorageLevel(s,h)          [MWh]   Storage level of charge
LostLoad_MaxPower(n,h)     [MW]    Deficit in terms of maximum power
LostLoad_RampUp(u,h)       [MW]    Deficit in terms of ramping up for each plant
LostLoad_RampDown(u,h)     [MW]    Deficit in terms of ramping down
LostLoad_MinPower(n,h)     [MW]    Power exceeding the demand
LostLoad_Reserve2U(n,h)    [MW]    Deficit in reserve up
LostLoad_Reserve2D(n,h)    [MW]    Deficit in reserve down
SystemCost(h)              [EUR]   Hourly system cost
;

free variable
SystemCostD               ![EUR]   Total system cost  for one optimization period
;

*===============================================================================
*Assignment of initial values
*===============================================================================
*CostVariable(u,h,b)=Markup(u,h,b)+sum((n,f),(Fuel(u,f)*FuelPrice(n,h,f)*Location(u,n))/Efficiency(u))+sum(p,EmissionRate(u,p)*PermitPrice(p));

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

CostLoadShedding(n,h)=1000;

Display RampStartUpMaximum, RampShutDownMaximum, CommittedInitial, FlexibilityUp, FlexibilityDown;

CostRampup(u)=0;
CostRampDown(u)=0;

$offorder

*===============================================================================
*Declaration and definition of equations
*===============================================================================
EQUATIONS
EQ_Objective_function
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
EQ_PowerMaximum_previous
EQ_PowerMaximum_following
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
;

*Objective function

EQ_SystemCost(i)..
         SystemCost(i)
         =E=
         sum((u),CostFixed(u)*Committed(u,i))
         +sum((u),CostStartUpH(u,i) + CostShutDownH(u,i))
         +sum((u),CostRampUpH(u,i) + CostRampDownH(u,i))
         +sum((u),CostVariable(u,i)*Power(u,i))
         +sum((l),PriceTransmission(l,i)*Flow(l,i))
         +sum((n),CostLoadShedding(n,i)*ShedLoad(n,i))
         +30E3*(sum((n),LostLoad_MaxPower(n,i)+LostLoad_MinPower(n,i)))
         +20E3*(sum((n),LostLoad_Reserve2U(n,i)+LostLoad_Reserve2D(n,i)))
         +10E3*sum((u),LostLoad_RampUp(u,i)+LostLoad_RampDown(u,i))
;

EQ_Objective_function..
         SystemCostD
         =E=
         sum(i,SystemCost(i))
;

EQ_CostStartUp(u,i)$(CostStartUp(u) <> 0)..
         CostStartUpH(u,i)
         =g=
         CostStartUp(u)*(Committed(u,i)-CommittedInitial(u)$(ord(i) = 1)-Committed(u,i-1)$(ord(i) > 1))
;

EQ_CostShutDown(u,i)$(CostShutDown(u) <> 0)..
         CostShutDownH(u,i)
         =g=
         CostShutDown(u)*(CommittedInitial(u)$(ord(i) = 1)+Committed(u,i-1)$(ord(i) > 1)-Committed(u,i))
;

EQ_CostRampUp(u,i)$(CostRampUp(u) <> 0)..
         CostRampUpH(u,i)
         =g=
         CostRampUp(u)*(Power(u,i)-PowerInitial(u)$(ord(i) = 1)-Power(u,i-1)$(ord(i) > 1))
;

EQ_CostRampDown(u,i)$(CostRampDown(u) <> 0)..
         CostRampDownH(u,i)
         =g=
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


*Hourly demand balance in the upwards reserve market for each node
*EQ_Demand_balance_2U(n,i)..
*         sum((u,t),PowerMaximumH(u,i)*Technology(u,t)*Reserve(t)*Location(u,n))
*         +CurtailedPowerH(n,i)*Curtailment(n,i)
*         =G=
*         DemandH(n,i,"DA")
*         +DemandH(n,i,"2U")
*         -LostLoad_reserve2U(n,i)
*;

*Hourly demand balance in the downwards reserve market for each node
*EQ_Demand_balance_2D(n,i)..
*         sum((u,t),PowerMinimumH(u,i)*Technology(u,t)*Reserve(t)*Location(u,n))
*         =L=
*         DemandH(n,i,"DA")
*         -DemandH(n,i,"2D")
*         +LostLoad_reserve2DH(n,i)
*;

*Minimum power output is above the must-run output level for each unit in all periods
EQ_Power_must_run(u,i)..
         PowerMustRun(u,i)
         *Committed(u,i)
*         +(PowerCapacity(u)
*                 *LoadMaximum(u,i)
*                 -PowerMustRunH(u,i))
*                         *Committed(u,i)
*                         *(1-Curtailment(u,i))
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

*Maximum power output with respect to power output in the previous period.
EQ_PowerMaximum_previous(u,i)$(sum(tr,Technology(u,tr))=0)..
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

*Maximum power output with respect to power output in the following period
EQ_PowerMaximum_following(u,i)$(sum(tr,Technology(u,tr))=0)..
         Power(u,i)$(ord(i) < card(i))
         =L=
         (PowerCapacity(u)
                 *LoadMaximum(u,i+1)
                         *Committed(u,i+1)
         +RampShutDownMaximum(u)
                 *(1-Committed(u,i+1)))$(ord(i) < card(i))
         +LostLoad_RampDown(u,i)
;

*If the unit keeps committed the reduction in power output is lower than the
*ramp-down limit. If the unit is de-committed the reduction is lower than the
*shut-down ramp limit
EQ_Ramp_down(u,i)$(sum(tr,Technology(u,tr))=0)..
         (PowerInitial(u)-Power(u,i))$(ord(i) = 1)
         +(Power(u,i-1)-Power(u,i))$(ord(i) > 1)
         =L=
         PowerCapacity(u)
                         *(1-Committed(u,i))
         +(RampDownMaximum(u)
                 *CommittedInitial(u))$(ord(i) = 1)
         +(RampDownMaximum(u)
                 *Committed(u,i-1))$(ord(i) > 1)
         +LostLoad_RampDown(u,i)
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

*IH: why do we need this equation?, to replace the two previous?
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
         StorageCapacity(s)
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
         +StorageOutflow(s,i) - StorageInflow(s,i)
         =L=
         StorageLevel(s,i)
;

*Charging is limited by the remaining storage capacity
EQ_Storage_MaxCharge(s,i)..
         StorageInput(s,i) * StorageChargingEfficiency(s)
         -StorageOutflow(s,i) + StorageInflow(s,i)
         =L=
         StorageCapacity(s) - StorageLevel(s,i)
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
         +Power(s,i)/(max(StorageDischargeEfficiency(s),0.0001))
;

* Assuming cyclic boundary counditions:
EQ_Storage_boundaries(s,i)$(ord(i) = card(i))..
         StorageInitial(s)
         =E=
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
         LoadShedding(n)
;


*===============================================================================
*Definition of models
*===============================================================================
MODEL UCM_SIMPLE /
EQ_Objective_function,
EQ_CostStartUp,
EQ_CostShutDown,
*EQ_CostRampUp,
*EQ_CostRampDown,
EQ_Demand_balance_DA,
EQ_Demand_balance_2U,
EQ_Demand_balance_2D,
EQ_Power_must_run,
EQ_Power_available,
EQ_PowerMaximum_previous,
EQ_PowerMaximum_following,
EQ_Ramp_down,
*EQ_Minimum_time_up,
EQ_Minimum_time_up_A,
*EQ_Minimum_time_up_B,
*EQ_Minimum_time_up_C,
EQ_Minimum_time_up_JustStarted,
*EQ_Minimum_time_down,
EQ_Minimum_time_down_A,
*EQ_Minimum_time_down_B,
*EQ_Minimum_time_down_C,
EQ_Minimum_time_down_JustStopped,
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
EQ_LoadShedding
/
;
UCM_SIMPLE.optcr = 0.04;
*UCM_SIMPLE.epgap = 0.005
*UCM_SIMPLE.probe = 3
*UCM_SIMPLE.optfile=1;

*===============================================================================
*Solving loop
*===============================================================================

* Define the general index of the first and last days and of the start and stop days:
scalar index_first,index_last,index_start,index_stop;
index_first=jdate(Config("FirstDay","year"),Config("FirstDay","month"),Config("FirstDay","day")) ;
index_last=jdate(Config("LastDay","year"),Config("LastDay","month"),Config("LastDay","day")) ;
index_start=jdate(Config("DayStart","year"),Config("DayStart","month"),Config("DayStart","day"))  ;
index_stop=jdate(Config("DayStop","year"),Config("DayStop","month"),Config("DayStop","day"));

* Check that the length of index h corresponds to the provided first and last days:
if (card(h) <> (index_last-index_first+1)*24, abort "The number of time indexes does not correspond to the provided first and last days";);
* Check that the dates make sense
if (index_start < index_first or index_start > index_last or index_stop < index_start or index_stop > index_last, abort "The start and stop dates must be comprised between the first and last dates";);
* Check the rolling horizon length
if (Config("RollingHorizon Length","day") + Config("RollingHorizon LookAhead","day") > index_stop - index_start + 1, abort "The rolling horizon is longer than the simulation length";);




display index_first,index_last,index_start,index_stop;

* Scalar variables necessary to the loop:
scalar FirstHour,LastHour,LastKeptHour,day;

* Fixing the initial guesses:
*PowerH.L(u,i)=PowerInitial(u);
*Committed.L(u,i)=CommittedInitial(u);

FOR(day=index_start-index_first+1 TO index_stop-index_first-(Config("RollingHorizon Length","day")+Config("RollingHorizon LookAhead","day"))+2 by Config("RollingHorizon Length","day"),
         FirstHour = (day-1)*24+1;
         LastHour = FirstHour + (Config("RollingHorizon Length","day")+Config("RollingHorizon LookAhead","day")) * 24 - 1;
         LastKeptHour = FirstHour + Config("RollingHorizon Length","day") * 24 - 1;
         i(h) = no;
         i(h)$(ord(h)>=firsthour and ord(h)<=lasthour)=yes;
         display FirstHour,LastHour,LastKeptHour;

*        Update the subset of all simulated hours:
         z(h)$(ord(h)>=firsthour and ord(h)<=lastkepthour)=yes;

         TimeUpLeft_initial(u)=min(card(i),(TimeUpMinimum(u)-TimeUpInitial(u))*CommittedInitial(u));
         TimeUpLeft_JustStarted(u,i) = min(card(i)-ord(i)+1,TimeUpMinimum(u));
         TimeDownLeft_initial(u)=min(card(i),(TimeDownMinimum(u)-TimeDownInitial(u))*(1-CommittedInitial(u)));
         TimeDownLeft_JustStopped(u,i) = min(card(i)-ord(i)+1,TimeDownMinimum(u));

         Display TimeUpLeft_initial,TimeUpLeft_JustStarted,TimeDownLeft_initial,TimeDownLeft_JustStopped,TimeUpInitial,TimeDownInitial,PowerInitial,CommittedInitial;

         SOLVE UCM_SIMPLE USING MIP MINIMIZING SystemCostD;

         Display EQ_Objective_function.M,EQ_CostStartUp.M,EQ_CostShutDown.M,EQ_Demand_balance_DA.M,EQ_Power_must_run.M,EQ_Power_available.M,EQ_PowerMaximum_previous.M,EQ_PowerMaximum_following.M,EQ_Ramp_down.M,EQ_Minimum_time_up_A.M,EQ_Minimum_time_up_JustStarted.M,EQ_Minimum_time_down_A.M,EQ_Minimum_time_down_JustStopped.M,   EQ_Max_RampUp1.M,    EQ_Max_RampUp2.M,EQ_Max_RampDown1.M,  EQ_Max_RampDown2.M, EQ_Storage_minimum.M, EQ_Storage_level.M, EQ_Storage_input.M, EQ_Storage_balance.M,EQ_Storage_boundaries.M,EQ_Storage_MaxCharge.m,EQ_Storage_MaxDischarge.m,EQ_Flow_limits_lower.M;

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

*Loop variables to display after solving:
         Display LastKeptHour,PowerInitial,TimeUp,TimeDown,MaxRamp2D.L,MaxRamp2U.L,CostStartUpH.L,CostShutDownH.L;

);

CurtailedPower.L(n,z)=sum(u,(PowerCapacity(u)*LoadMaximum(u,z)-Power.L(u,z))$(sum(tr,Technology(u,tr))>=1) * Location(u,n));

Display Flow.L,Power.L,Committed.L,ShedLoad.L,CurtailedPower.L,StorageLevel.L,StorageInput.L,SystemCost.L,MaxRamp2U.L,MaxRamp2D.L,LostLoad_MaxPower.L,LostLoad_MinPower.L,LostLoad_reserve2U.L,LostLoad_reserve2D.L,LostLoad_RampUP.L,LostLoad_RampDown.L;


*===============================================================================
*Result export
*===============================================================================

* Parameters that contain the date information for the output:
parameter YearOutput(h),MonthOutput(h),DayOutput(h);
YearOutput(z) = gyear(index_start+floor((ord(z)-1)/24));
MonthOutput(z) = gmonth(index_start+floor((ord(z)-1)/24));
DayOutput(z) = gday(index_start+floor((ord(z)-1)/24));

PARAMETER
OutputCommitted(u,h)
*OutputDown(u,h)
*OutputUp(u,h)
OutputFlow(l,h)
OutputPower(u,h)
OutputStorageInput(s,h)
OutputStorageLevel(s,h)
*OutputTimeDown(u,i)
*OutputTimeUp(u,i)
OutputSystemCost(h)
;

OutputCommitted(u,z)=Committed.L(u,z);
OutputFlow(l,z)=Flow.L(l,z);
OutputPower(u,z)=Power.L(u,z);
OutputStorageInput(s,z)=StorageInput.L(s,z);
OutputStorageLevel(s,z)=StorageLevel.L(s,z);
*OutputTimeDown(u,i)=TimeDown(u,i);
*OutputTimeUp(u,h)=TimeUp(u,i);
OutputSystemCost(z)=SystemCost.L(z);

*OutputFlow(l,h)$(OutputFlow(l,h) = 0)=eps;
*OutputPower(u,z)$(OutputPower(u,z) = 0)=eps;
*OutputStorageInput(s,h)$(OutputStorageInput(s,h) = 0)=eps;
*OutputStorageLevel(s,h)$(OutputStorageLevel(s,h) = 0)=eps;

EXECUTE_UNLOAD "Results.gdx"
YearOutput ,
MonthOutput,
DayOutput,
OutputCommitted,
OutputFlow,
OutputPower,
OutputStorageInput,
OutputStorageLevel,
*OutputTimeDown,
*OutputTimeUp,
OutputSystemCost,
LostLoad_MaxPower,
LostLoad_MinPower,
LostLoad_Reserve2D,
LostLoad_Reserve2U,
LostLoad_RampUp,
LostLoad_RampDown,
ShedLoad,
CurtailedPower
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

EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=DayOutput rng=Power!B3 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=MonthOutput rng=Power!B2 cdim=1'
EXECUTE 'GDXXRW.EXE "Results.gdx" O="Results.xlsx" Squeeze=N par=YearOutput rng=Power!B1 cdim=1'

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


