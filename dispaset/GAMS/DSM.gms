$Title DSM model

$ontext
DSM heat model to be used alone or to be coupled with the UCM.gms file
$offtext


$setglobal RunHeat 0

*===============================================================================
* DEFINITION OF SETS
*===============================================================================

SETS
ind              Houses
h
xx               Air Wall Floor
yy
i(h)
k
*type /HP, WH/
;

alias(xx,xxx);

*===============================================================================
* DEFINITION OF PARAMETERS
*===============================================================================

PARAMETERS
P_base(h)

cst      n_cas_agg(house cases) K1-K2(HP parameters) T_low-high_sp_b(limits on DHW temperature) T_DHW_0(initial DHW temperature) Q_dot_WH_n (maximum WH power)
a_SH(xx,xx,ind)          []      SH state space parameters
b_SH(xx,yy,ind)          []      SH state space parameters
ratio(ind)               []
occ_type(ind)            []      Part of houses typologies
Q_dot_WH_n(ind)          [W]     Nominal (maximum) water heater power
P_heater_n(ind)          [W]     Nominal (maximal) additional heater power (HP)
T_SH_0(xx,ind)           [K]     Initial house temperatures
T_DHW_0(ind)             [K]     Initial domestic hot water temperatures

* Full year parameters
a_DHW_f(*,ind)           []      DHW state space parameters
b_DHW_f(xx,*,ind)        []      DHW state space parameters
COP_fl_DHW_f(*,ind)      []      Heat pump COP for DHW at full load
I_glob_h_f(*)            [W\m2]  Global horizontal irradiation
P_base_HP_mean_f(*)      [W]     Base electric consumption of heat pumps
P_base_WH_mean_f(*)      [W]     Base electric consumption of water heaters
Q_dot_fl_DHW_f(*,ind)    [W]     Heat generation of HP at full load for DHW
Q_dot_fl_SH_f(*,ind)     [W]     Heat generation of HP at full load for space heating
Q_dot_int_f(xx,*,ind)    [W]     Internal gains
T_amb_f(*)               [K]     Ambiant temperature
T_low_sp_b_f(ind,*)      [K]     Minimum space heating temperature
T_high_sp_b_f(ind,*)     [K]     Maximum space heating temperature
T_mains_c_f(*)           [K]     City water temperature
W_dot_fl_SH_f(*,ind)     [W]     Heat pump consumption at full load for DHW

* Configuration data
Config                           Give the start and end date
NumberHouses                     Give the number of houses ("SH" and "WH")        /HP.Base  0, HP.Simulated 0, WH.Base 0, WH.Simulated 0/

* Initial temperatures for iterations
T_SH_initial(xx,ind)     [K]     House temperatures
T_DHW_HP_initial(ind)    [K]     DHW tank temperature (HP)
T_DHW_WH_initial(ind)    [K]     DHW tank temperature (WH)

* Sliced parameters
a_DHW(h,ind)
b_DHW(xx,h,ind)
COP_fl_DHW(h,ind)
I_glob_h(h)
P_base_HP_mean(h)
P_base_WH_mean(h)
Q_dot_fl_DHW(h,ind)
Q_dot_fl_SH(h,ind)
Q_dot_int(xx,h,ind)
T_amb(h)
T_low_sp_b(ind,h)
T_high_sp_b(ind,h)
T_mains_c(h)
W_dot_fl_SH(h,ind)
;

*===============================================================================
* DATA IMPORT AND SLICING
*===============================================================================

$gdxin Inputs.gdx
$LOAD h
$LOAD Config
$LOADR NumberHouses
;

$gdxin InputsHeat.gdx
$LOAD ind
$LOAD xx
$LOAD yy
$LOAD cst
$LOAD k
$LOAD a_DHW_f
$LOAD a_SH
$LOAD b_DHW_f
$LOAD b_SH
$LOAD COP_fl_DHW_f
$LOAD I_glob_h_f
$LOAD occ_type
$LOAD P_base_HP_mean_f
$LOAD P_base_WH_mean_f
$LOAD P_heater_n
$LOAD Q_dot_fl_DHW_f
$LOAD Q_dot_fl_SH_f
$LOAD Q_dot_int_f
$LOAD Q_dot_WH_n
$LOAD ratio
$LOAD T_amb_f
$LOAD T_low_sp_b_f
$LOAD T_high_sp_b_f
$LOAD T_mains_c_f
$LOAD W_dot_fl_SH_f
$LOAD T_SH_0
$LOAD T_DHW_0
;

*Slicing parameters
Parameter
HourFirst        First hour of the simulation
DayFirst         First day of the simulation
DayOne           First day of the annual data
;

DayOne = jdate(Config("FirstDay","year"),1,1);
DayFirst=jdate(Config("FirstDay","year"),Config("FirstDay","month"),Config("FirstDay","day"));
HourFirst = 1 + (DayFirst - DayOne) * 24;

display DayOne, DayFirst, HourFirst;

If(HourFirst = 1,
         a_DHW(h,ind) = a_DHW_f(h,ind);
         b_DHW(xx,h,ind) = b_DHW_f(xx,h,ind);
         COP_fl_DHW(h,ind) =  COP_fl_DHW_f(h,ind);
         I_glob_h(h) = I_glob_h_f(h);
         P_base_HP_mean(h)= P_base_HP_mean_f(h);
         P_base_WH_mean(h) = P_base_WH_mean_f(h);
         Q_dot_fl_DHW(h,ind) = Q_dot_fl_DHW_f(h,ind);
         Q_dot_fl_SH(h,ind)  = Q_dot_fl_SH_f(h,ind);
         Q_dot_int(xx,h,ind) = Q_dot_int_f(xx,h,ind);
         T_amb(h) = T_amb_f(h);
         T_low_sp_b(ind,h) = T_low_sp_b_f(ind,h);
         T_high_sp_b(ind,h) = T_high_sp_b_f(ind,h);
         T_mains_c(h) = T_mains_c_f(h);
         W_dot_fl_SH(h,ind)= W_dot_fl_SH_f(h,ind);
else
         Loop((h,k)$(ord(k)=ord(h)+(HourFirst-1)),
                 a_DHW(h,ind) = a_DHW_f(k,ind);
                 b_DHW(xx,h,ind) = b_DHW_f(xx,k,ind);
                 COP_fl_DHW(h,ind) =  COP_fl_DHW_f(k,ind);
                 I_glob_h(h) = I_glob_h_f(k);
                 P_base_HP_mean(h)= P_base_HP_mean_f(k);
                 P_base_WH_mean(h) = P_base_WH_mean_f(k);
                 Q_dot_fl_DHW(h,ind) = Q_dot_fl_DHW_f(k,ind);
                 Q_dot_fl_SH(h,ind)  = Q_dot_fl_SH_f(k,ind);
                 Q_dot_int(xx,h,ind) = Q_dot_int_f(xx,k,ind);
                 T_amb(h) = T_amb_f(k);
                 T_low_sp_b(ind,h) = T_low_sp_b_f(ind,k);
                 T_high_sp_b(ind,h) = T_high_sp_b_f(ind,k);
                 T_mains_c(h) = T_mains_c_f(k);
                 W_dot_fl_SH(h,ind)= W_dot_fl_SH_f(k, ind);
         ;)
);

Parameter N_HP_add, N_WH_add, P_base_add(h), ModelHeat;


N_HP_add = max(0, NumberHouses("HP","Simulated") - NumberHouses("HP","Base"));
N_WH_add = max(0, NumberHouses("WH","Simulated") - NumberHouses("WH","Base"));
P_base_add(h) = N_HP_add * P_base_HP_mean(h) + N_WH_add * P_base_WH_mean(h);

$If %RunHeat% == 1 ModelHeat = 1;
$If not %RunHeat% == 1 ModelHeat = 0;

P_base(h) = NumberHouses("HP","Simulated") * P_base_HP_mean(h) + NumberHouses("WH","Simulated") * P_base_WH_mean(h);


*===============================================================================
* DEFINITION OF VARIABLES
*===============================================================================

POSITIVE VARIABLES
* HEAT PUMP
* Space heating
T_SH_b(xx,ind,h)         [K]     House temperatures
x_SH_b(ind,h)            [K]     Air temperature within limits
T_SH_cooler_b(ind,h)     [K]     Relaxation variable - cooler than minimum air temperature
T_SH_heater_b(ind,h)     [K]     Relaxation variable - hotter than maximum air temperature
MV_SH_b(ind,h)           [W]     Total heat power released for space heating
P_heater_b(ind,h)        [W]     Power of the additional heater
Q_dot_hp_SH_b(ind,h)     [W]     Heat power released by the heat pump for space heating
W_dot_hp_SH_b(ind,h)     [W]     Electric power used by the heat pump for space heating
* Domestic hot water
T_DHW_HP_b(ind,h)        [K]     DHW temperature (HP)
T_DHW_cooler_HP_b(ind,h) [K]     Relaxation variable - cooler than minimum DHW temperature (HP)
x_DHW_HP_b(ind,h)        [K]     DHW tank temperature within limits (HP)
MV_DHW_b(ind,h)          [W]     Heat power released for DHW (HP)
W_dot_hp_DHW_b(ind,h)    [W]     Electric power used for DHW (HP)
* Other
P_HP_b(ind,h)            [W]     Total electric power used for heating with HP
P_HP_mean_b(h)           [W]     Mean electric power used for heating with HP
y(ind,h)                 []      Fraction of the heat pump functionning in space heating

* WATER HEATER (only DHW)
T_DHW_WH_b(ind,h)        [K]     DHW temperature (WH)
T_DHW_cooler_WH_b(ind,h) [K]     Relaxation variable - cooler than minimum DHW temperature (WH)
x_DHW_WH_b(ind,h)        [K]     DHW tank temperature within limits (WH)
P_WH_b(ind,h)            [W]     Total electric and heat power used for heating with WH
P_WH_mean_b(h)           [W]     Mean electric power used for heating with WH

* TOTAL
P_tot_b(h)               [W]     Total electric consumption at step h
T_relax_tot_b(h)         [K]     Total deviation of the temperature limits (sum of relaxation variables)
;

FREE VARIABLE
P_diff_b(h)              [W]     Difference in electric consumption with the base consumption at step h
;


* Variable boundaries:
y.up(ind,h) = 1;
x_DHW_HP_b.lo(ind,h) = cst("T_tank_sp_low");
x_DHW_HP_b.up(ind,h) = cst("T_tank_sp_high");
x_DHW_WH_b.lo(ind,h) = cst("T_tank_sp_low");
x_DHW_WH_b.up(ind,h) = cst("T_tank_sp_high");

* Initial temperatures
T_SH_initial(xx,ind) = T_SH_0(xx,ind);
T_DHW_HP_initial(ind) = T_DHW_0(ind);
T_DHW_WH_initial(ind) = T_DHW_0(ind);

$$offorder

*===============================================================================
* DEFINITION OF EQUATIONS
*===============================================================================

EQUATIONS
* HEAT PUMP
*Space heating
EQ_T_SH
EQ_T_SH_initial
EQ_x_SH
EQ_x_SH_min
EQ_x_SH_max
EQ_HP_DHW
EQ_HP_SH1
EQ_HP_SH2
EQ_HP_SH3
EQ_heater
EQ_Qdot_HP_SH_max
EQ_heater_max
*Domestic hot water
EQ_T_DHW_HP
EQ_T_DHW_HP_initial
EQ_x_DHW_HP
EQ_Qdot_HP_DHW_max
*Total
EQ_P_HP
EQ_P_HP_mean

*WATER HEATER
EQ_T_DHW_WH
EQ_T_DHW_WH_initial
EQ_x_DHW_WH
EQ_P_WH_max
EQ_P_WH_mean

*TOTAL
EQ_P_tot
EQ_P_diff
EQ_T_relax_tot
*EQ_Obj
;

* Initial (guess) values:
x_SH_b.l(ind,h) = 293;
T_SH_b.l(xx,ind,h) = 293;
x_DHW_HP_b.l(ind,h) = 330;
x_DHW_WH_b.l(ind,h) = 330;
T_DHW_HP_b.l(ind,h) = 330;
T_DHW_WH_b.l(ind,h) = 330;

*--------HEAT PUMP-----------
* Space heating


EQ_T_SH(xx,ind,i)$(ord(i)>1)..
         T_SH_b(xx,ind,i)
         =e=
         sum(xxx,a_SH(xx,xxx,ind)*T_SH_b(xxx,ind,i-1))
         + b_SH(xx,'1',ind) * T_amb(i)
         + b_SH(xx,'2',ind) * I_glob_h(i)
         + b_SH(xx,'3',ind) * MV_SH_b(ind,i)
         + b_SH(xx,'4',ind) * Q_dot_int(xx,i,ind);

EQ_T_SH_initial(xx,ind,i)$(ord(i)=1)..
         T_SH_b(xx,ind,i)
         =e=
         sum(xxx,a_SH(xx,xxx,ind)*T_SH_initial(xxx,ind))
         + b_SH(xx,'1',ind) * T_amb(i)
         + b_SH(xx,'2',ind) * I_glob_h(i)
         + b_SH(xx,'3',ind) * MV_SH_b(ind,i)
         + b_SH(xx,'4',ind) * Q_dot_int(xx,i,ind);

EQ_x_SH(ind,i)..
         T_SH_b('1',ind,i)
         =e=
         x_SH_b(ind,i)
         - T_SH_cooler_b(ind,i)
         + T_SH_heater_b(ind,i);

EQ_x_SH_min(ind,i)..
         x_SH_b(ind,i) =g= T_low_sp_b(ind,i);

EQ_x_SH_max(ind,i)..
         T_high_sp_b(ind,i) =g= x_SH_b(ind,i);

EQ_HP_SH1(ind,i)..
         W_dot_hp_SH_b(ind,i) =g= 0.77 * Q_dot_hp_SH_b(ind,i)/Q_dot_fl_SH(i,ind)*W_dot_fl_SH(i,ind);

EQ_HP_SH2(ind,i)..
         W_dot_hp_SH_b(ind,i) =g= (0.6881+((cst("K_2")-cst("K_1"))+2*(1-cst("K_2"))*0.75)*(Q_dot_hp_SH_b(ind,i)/Q_dot_fl_SH(i,ind)-0.75))*W_dot_fl_SH(i,ind);

EQ_HP_SH3(ind,i)..
         W_dot_hp_SH_b(ind,i) =l= Q_dot_hp_SH_b(ind,i) * W_dot_fl_SH(i,ind) /Q_dot_fl_SH(i,ind);

EQ_heater(ind,i)..
        P_heater_b(ind,i) =e= MV_SH_b(ind,i)- Q_dot_hp_SH_b(ind,i);

EQ_Qdot_HP_SH_max(ind,i)..
         Q_dot_hp_SH_b(ind,i) =l= Q_dot_fl_SH(i,ind)*y(ind,i);

EQ_heater_max(ind,i)..
        P_heater_b(ind,i) =l= P_heater_n(ind)*y(ind,i);

* Domestic hot water

EQ_T_DHW_HP(ind,i)$(ord(i)>1)..
         T_DHW_HP_b(ind,i)
         =e=
         a_DHW(i,ind)*T_DHW_HP_b(ind,i-1)
         + b_DHW('1',i,ind)*T_amb(i)
         + b_DHW('2',i,ind)*T_mains_c(i)
         + b_DHW('3',i,ind)*MV_DHW_b(ind,i);

EQ_T_DHW_HP_initial(ind,i)$(ord(i)=1)..
         T_DHW_HP_b(ind,i)
         =e=
         a_DHW(i,ind)*T_DHW_HP_initial(ind)
         + b_DHW('1',i,ind)*T_amb(i)
         + b_DHW('2',i,ind)*T_mains_c(i)
         + b_DHW('3',i,ind)*MV_DHW_b(ind,i);

EQ_x_DHW_HP(ind,i)..
         T_DHW_HP_b(ind,i)
         =e=
         x_DHW_HP_b(ind,i)
        - T_DHW_cooler_HP_b(ind,i);

EQ_HP_DHW(ind,i)..
        W_dot_hp_DHW_b(ind,i) =e= MV_DHW_b(ind,i)/COP_fl_DHW(i,ind);

EQ_Qdot_HP_DHW_max(ind,i)..
         MV_DHW_b(ind,i) =l= Q_dot_fl_DHW(i,ind)*(1-y(ind,i));

* Total

EQ_P_HP(ind,i)..
        P_HP_b(ind,i) =e= W_dot_hp_SH_b(ind,i)+W_dot_hp_DHW_b(ind,i)+P_heater_b(ind,i);

EQ_P_HP_mean(i)..
         P_HP_mean_b(i) =e= sum(ind,occ_type(ind) * P_HP_b(ind,i))
;


*-----------WATER HEATER-------------

EQ_T_DHW_WH(ind,i)$(ord(i)>1)..
         T_DHW_WH_b(ind,i)
         =e=
         a_DHW(i,ind)*(T_DHW_WH_b(ind,i-1))
         + b_DHW('1',i,ind)*T_amb(i)
         + b_DHW('2',i,ind)*T_mains_c(i)
         + b_DHW('3',i,ind)*P_WH_b(ind,i);

EQ_T_DHW_WH_initial(ind,i)$(ord(i)=1)..
         T_DHW_WH_b(ind,i)
         =e=
         a_DHW(i,ind)*T_DHW_WH_initial(ind)
         + b_DHW('1',i,ind)*T_amb(i)
         + b_DHW('2',i,ind)*T_mains_c(i)
         + b_DHW('3',i,ind)*P_WH_b(ind,i);

EQ_x_DHW_WH(ind,i)..
         T_DHW_WH_b(ind,i)
         =e=
         x_DHW_WH_b(ind,i)
        - T_DHW_cooler_WH_b(ind,i);


EQ_P_WH_max(ind,i)..
         P_WH_b(ind,i) =l= Q_dot_WH_n(ind);

EQ_P_WH_mean(i)..
         P_WH_mean_b(i) =e= sum(ind,occ_type(ind) * P_WH_b(ind,i))
;


*------------TOTAL-------------

EQ_P_tot(i)..
         P_tot_b(i) =e= NumberHouses('HP','Simulated') * P_HP_mean_b(i)
                        + NumberHouses('WH','Simulated') * P_WH_mean_b(i);

EQ_P_diff(i)..
         P_diff_b(i) =e= P_tot_b(i)
                         - P_base(i);

EQ_T_relax_tot(i)..
         T_relax_tot_b(i) =e= sum((ind), occ_type(ind) * (T_SH_cooler_b(ind,i) + T_SH_heater_b(ind,i)+ T_DHW_cooler_HP_b(ind,i) + T_DHW_cooler_WH_b(ind,i)));

*EQ_obj..
*         Obj =e= sum(i,P_tot_b(i)+100E9*T_relax_tot_b(i));



*===============================================================================
* MODEL
*===============================================================================

MODEL DSM /
EQ_T_DHW_HP
EQ_T_DHW_HP_initial
EQ_T_SH
EQ_T_SH_initial
EQ_x_DHW_HP
EQ_x_SH
EQ_x_SH_min
EQ_x_SH_max
EQ_HP_DHW
EQ_HP_SH1
EQ_HP_SH2
EQ_HP_SH3
EQ_heater
EQ_Qdot_HP_DHW_max
EQ_Qdot_HP_SH_max
EQ_heater_max
EQ_P_HP
EQ_P_HP_mean
EQ_T_DHW_WH
EQ_T_DHW_WH_initial
EQ_x_DHW_WH
EQ_P_WH_max
EQ_P_WH_mean
EQ_P_tot
EQ_P_diff
EQ_T_relax_tot
*EQ_obj
/
;

* -- Test
*T_SH_b.fx(xx,ind,'1')=T_SH_initial(xx,ind);
*T_DHW_HP_b.fx(ind,'1')=T_DHW_HP_initial(ind) ;
*T_DHW_WH_b.fx(ind,'1')=T_DHW_HP_initial(ind) ;
*i(h)=yes$(ord(h)<10);
*SOLVE DSM USING LP MINIMIZING Obj;



