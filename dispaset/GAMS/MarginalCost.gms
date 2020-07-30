If(%RunHeat%,
Demand("DA", n, h) = Demand("DA",n,h)+ 1E-6 * P_diff_b.l(h) + 1E-6 * P_base_add(h);
P_diff_b.fx(h) = 0  ;
T_relax_tot_b.fx(h)=0  ;
*Committed.fx(u,z) = OutputCommitted(u,z)  ;


$gdxin Inputs.gdx

$LOADR PowerInitial
$LOADR StorageInitial
$LOADR TimeDownInitial
$LOADR TimeUpInitial
;

CommittedInitial(u)=0;
CommittedInitial(u)$(PowerInitial(u)>0)=1;

MODEL UCM_COST /
EQ_Objective_function,
EQ_CHP_extraction_Pmax,
EQ_CHP_extraction,
EQ_CHP_backpressure,
EQ_CHP_demand_satisfaction,
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
UCM_COST.optcr = 0.01;

$offorder
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

         status("model",i) = UCM.Modelstat;
         status("solver",i) = UCM.Solvestat;

if(UCM.Modelstat <> 1 and UCM.Modelstat <> 8 and not failed, TimeUpInitial_dbg(u) = TimeUpInitial(u); TimeDownInitial_dbg(u) = TimeDownInitial(u); CommittedInitial_dbg(u) = CommittedInitial(u); PowerInitial_dbg(u) = PowerInitial(u); StorageInitial_dbg(s) = StorageInitial(s);
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

ShadowPrice(n,z) = EQ_Demand_balance_DA.m(n,z);

$onorder

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
);

