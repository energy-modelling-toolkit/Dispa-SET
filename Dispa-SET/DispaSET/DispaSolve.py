# -*- coding: utf-8 -*-
"""
This worksheet contains the two main functions to solve the Dispa-SET optimization problem using PYOMO.

Functions: 
- DispaSolve: Implements the rolling horizon
- DispOptim: Performs the optimization for each horizon

@author: 'Sylvain Quoilin'
"""

#######################################################################################################################
############################################ Dispa-SET: main model ####################################################
#######################################################################################################################


import sys
from pyomo.environ import *
from pyomo.opt import TerminationCondition

import numpy as np
import pandas as pd
import logging

from DispaTools import pyomo_format, pyomo_to_pandas


def  DispOptim(sets, parameters, Mixed_Integer_LP=False):
    '''
    This is the main optimization function of Dispaset. 
    Two operation are performed:
    1. Translation of the dispaset data format into the pyomo format
    2. Definition of the Pyomo optimization model as a ConcreteModel
    
    :param sets: Dictionary containing the sets (defined as a list of strings or integegers)
    :param parameters: Dictionary containing the parameters of the optimization problem (in the DispaSET 2.1.1 format)
    
    :return: The Pyomo optimization model instance
    '''
    
    # Definition of model:
    model = ConcreteModel('DispaSET')

#######################################################################################################################
############################################ Definition of the sets ###################################################
#######################################################################################################################

    # Assign all sets to the pyomo model (to be changed using real Pyomo sets!!)

    model.h = Set(initialize=sets['h'])                    # Hours
    model.d = Set(initialize=sets['d'])             # Day
    model.mk = Set(initialize=sets['mk'])   # Market
    model.n = Set(initialize=sets['n'])    # Nodes
    model.p = Set(initialize=sets['p'])    # Pollutant
    model.l = Set(initialize=sets['l'])    # Lines
    model.f = Set(initialize=sets['f'])    # Fuel types
    model.u = Set(initialize=sets['u'])    #units
    model.t = Set(initialize=sets['t'])     # Generation technologies
    model.s = Set(initialize=sets['s'])                       # Storage Unit (with reservoir)
    model.tr = Set(initialize=sets['tr'])                      # Renewable generation technologies

        
    # Transform the parameters into the pyomo format:
    params = {}
    for key in parameters.keys():
        params[key] = pyomo_format(sets, parameters[key])
    #params = pyomo_format(sets,parameters)


########################################################################################################################
######################################### Definition of the Parameters #################################################
########################################################################################################################

    model.LoadMaximum = Param(sets['u'],sets['h'],initialize=params['LoadMaximum'])          # [%] Availability factor
    model.CostFixed = Param(sets['u'],initialize=params['CostFixed'])                                   # [€/h] Fixed costs
    model.CommittedInitial = Param(sets['u'],initialize=params['CommittedInitial'])                                   # [€/h] Fixed costs
    model.CostLoadShedding = Param(model.n,model.h,initialize=params['CostLoadShedding'])                                         # [€/MW] Value of lost load (Initialize after)
    model.CostRampDown = Param(sets['u'],initialize=params['CostRampDown'])                         # [€/MW/h] Ramp-up costs 
    model.CostRampUp = Param(sets['u'],initialize=params['CostRampUp'])                                                    # [€/MW/h] Ramp-down costs 
    model.CostShutDown = Param(sets['u'],initialize=params['CostShutDown'])                              # [€] Shut-down costs
    model.CostStartUp = Param(sets['u'],initialize=params['CostStartUp'])                                # [€] Start-up costs
    model.CostVariable = Param(sets['u'],sets['h'],initialize=params['CostVariable'])             # [€/MW] Variable costs
    model.Demand = Param(sets['mk'],sets['n'],sets['h'],initialize=params['Demand'],mutable=True)                        # [MW] Demand
    model.Efficiency = Param(sets['u'],initialize=params['Efficiency'])                                  # [%] Efficiency
    model.EmissionMaximum = Param(sets['n'],sets['p'],initialize=params['EmissionMaximum'])                # [tP] Emission limit
    model.EmissionRate = Param(sets['u'],sets['p'],initialize=params['EmissionRate'])                      # [tP/MW] P emission rate
    model.FlowMaximum = Param(sets['l'],sets['h'],initialize=params['FlowMaximum'])                        # [MW] Line limits
    model.FlowMinimum = Param(sets['l'],sets['h'],initialize=params['FlowMinimum'])                        # [MW] Minimum flow
    model.FuelPrice = Param(sets['n'],sets['f'],sets['h'],initialize=params['FuelPrice'])                    # [€/F] Fuel price
    model.Fuel = Param(sets['u'],sets['f'],initialize=params['Fuel'])                                      # [n.a.] Fuel type {1 0}
    model.LineNode = Param(sets['l'],sets['n'],initialize=params['LineNode'])                              # [n.a.] Incidence matrix {-1 +1}
    model.LoadShedding = Param(sets['n'],initialize=params['LoadShedding'])                              # [n.a.] Load shedding capacity
    model.Location = Param(sets['u'],sets['n'],initialize=params['Location'])                              # [n.a.] Location {1 0}
    model.OutageFactor = Param(sets['u'],sets['h'],initialize=params['OutageFactor'])                      # [%] Outage factor (100% = full outage)
    model.PartLoadMin = Param(sets['u'],initialize=params['PartLoadMin'])                                # [%] Minimum part load
    model.PermitPrice = Param(sets['p'],initialize=params['PermitPrice'])                                # €/tP] Permit price
    model.PowerCapacity = Param(sets['u'],initialize=params['PowerCapacity'])                            # [MW] Installed capacity
    model.PowerInitial = Param(sets['u'],initialize=params['PowerInitial'])                              # [MW] Power output before initial period
    model.PowerMustRun = Param(sets['u'],sets['h'],initialize=params['PowerMustRun'])          # [MW] must run power
    model.PowerMinStable = Param(sets['u'],initialize=params['PowerMinStable'])                                                   # [MW] Minimum power output (Initialize after)
    model.PriceTransmission = Param(sets['l'],sets['h'],initialize=params['PriceTransmission'])           # [€/Mwh] Transmission price
    model.StorageChargingCapacity = Param(sets['u'],initialize=params['StorageChargingCapacity'])        # [MW] Storage capacity
    model.StorageChargingEfficiency = Param(sets['u'],initialize=params['StorageChargingEfficiency'])    # [%] Charging efficiency
    model.RampDownMaximum = Param(sets['u'],initialize=params['RampDownMaximum'])                       # [MW/h] Ramp down limit
    model.RampShutDownMaximum = Param(sets['u'],initialize=params['RampShutDownMaximum'])   # [MW/h] Shut-down ramp limit
    model.RampStartUpMaximum = Param(sets['u'],initialize=params['RampStartUpMaximum'])     # [MW/h] Start-up ramp limit
    model.RampUpMaximum = Param(sets['u'],initialize=params['RampUpMaximum'])                            # [MW/h] Ramp up limit
    model.Reserve = Param(sets['t'],initialize=params['Reserve'])                                        # [n.a.] Reserve technology {1 0}
    model.StorageCapacity = Param(sets['u'],initialize=params['StorageCapacity'])                        # [MWh] Storage capacity
    model.StorageDischargeEfficiency = Param(sets['u'],initialize=params['StorageDischargeEfficiency'])  # [%] Discharge efficiency
    model.StorageInflow = Param(sets['u'],sets['h'],initialize=params['StorageInflow'])                    # [MWh] Storage inflows (potential energy)
    model.StorageInitial = Param(sets['u'],initialize=params['StorageInitial'])                          # [MWh] Storage level before initial period
    model.StorageMinimum = Param(sets['u'],initialize=params['StorageMinimum'])                          # [MWh] Storage minimum
    model.StorageOutflow = Param(sets['u'],sets['h'],initialize=params['StorageOutflow'])                  # [MWh] Storage outflows
    model.Technology = Param(sets['u'],sets['t'],initialize=params['Technology'])                         # [n.a] Technology type {1 0}
    model.TimeDown = Param(sets['u'],sets['h'])                                                  # [h] Hours down
    model.TimeDownLeft_initial = Param(sets['u'],initialize=params['TimeDownLeft_initial'])                                             # [h] Required time down left at the beginning of the simulated time period (Initialize after)
    model.TimeDownLeft_JustStopped = Param(sets['u'],sets['h'],initialize=params['TimeDownLeft_JustStopped'])                               # [h] Required time down left at hour h if the unit has just been stopped (Initialize after)
    model.TimeUpLeft_initial = Param(sets['u'],initialize=params['TimeUpLeft_initial'])                                               # [h] Required time up left at the beginning of the simulated period (Initialize after)
    model.TimeUpLeft_JustStarted = Param(sets['u'],sets['h'],initialize=params['TimeUpLeft_JustStarted'])                                # [h] Required time up left at hour h if the unit has just been started (Initialize after)
    model.FlexibilityUp = Param(sets['u'],initialize=params['FlexibilityUp'])                                                    # [MW/h] Flexibility (up) of fast-starting power plants (Initialize after)
    model.FlexibilityDown = Param(sets['u'],initialize=params['FlexibilityDown'])                                                  # [MW/h] Flexibility (down) of a committed power plant (Initialize after)

########################################################################################################################
##########################################  Definition of variables ####################################################
########################################################################################################################

    # a) Binary Variables:
    #model.Committed = Var(sets['h'],sets['u'],within = PositiveReals,bounds=(0,1))                   # [n.a.] Unit committed at hour h {1 0}
    if Mixed_Integer_LP:
        model.Committed = Var(model.u,model.h, within = Binary)
    else:
        model.Committed = Var(model.u,model.h, within = PositiveReals, bounds=(0,1))
        
    # b) Positive Variables:
    model.CostStartUpH=Var(model.u,model.h,within = PositiveReals)             #  [€]Cost of start-up
    model.CostShutDownH=Var(model.u,model.h,within = PositiveReals)            #[ €]Cost ofshut-down
    model.CostRampUpH=Var(model.u,model.h,within = PositiveReals)              # [€]Cost of Ramping up
    model.CostRampDownH=Var(model.u,model.h,within = PositiveReals)            # [€] Cost of Ramping down
    model.Flow = Var(model.l,model.h,within = PositiveReals)                   # [MW] Flow through lines
    model.Power = Var(model.u,model.h,within = PositiveReals)                  # [MW] Power output
    model.PowerMaximum = Var(model.u,model.h,within = PositiveReals)           # [MW] Power output
    model.PowerMinimum = Var(model.u,model.h,within = PositiveReals)           # [MW] Power output
    model.ShedLoad = Var(model.n,model.h,within = PositiveReals)               # [MW] Shed load
    model.StorageInput = Var(model.s,model.h,within = PositiveReals)           # [MWh] Charging input for storage units
    model.StorageLevel = Var(model.s,model.h,within = PositiveReals)           # [MWh] Storage level of charge
    model.SystemCost = Var(model.h,within = PositiveReals) 
    model.LostLoad_MaxPower = Var(model.n,model.h,within = PositiveReals)      # [MW] Deficit in terms of maximum power
    model.LostLoad_RampUp = Var(model.u,model.h,within = PositiveReals)        # [MW] Deficit in terms of ramping up for each plant
    model.LostLoad_RampDown = Var(model.u,model.h,within = PositiveReals)      # [MW] Deficit in terms of ramping down
    model.LostLoad_MinPower = Var(model.n,model.h,within = PositiveReals)      # [MW] Power exceeding the demand
    model.LostLoad_Reserve2U = Var(model.n,model.h,within = PositiveReals)     # Deficit in reserve up
    model.LostLoad_Reserve2D = Var(model.n,model.h, within = PositiveReals)    # Deficit in reserve Down

    model.MaxRamp2U = Var(model.u,model.h,within = PositiveReals)
    model.MaxRamp2D = Var(model.u,model.h,within = PositiveReals)



########################################################################################################################
############################################### EQUATIONS ##############################################################
########################################################################################################################

# 0) Objective function:
###################################################################################################################################################################################################################################################################

    def EQ_Objective_function(model):
                return sum(model.SystemCost[h] for h in model.h)

###################################################################################################################################################################################################################################################################

# 1) Systeam cost at each hour:
###################################################################################################################################################################################################################################################################

    def EQ_SystemCost(model,h):
                return model.SystemCost[h] == (sum(model.CostFixed[u] * model.Committed[u,h] for u in model.u)
                    + sum((model.CostStartUpH[u,h] + model.CostShutDownH[u,h]) for u in model.u)
                    + sum((model.CostRampUpH[u,h] + model.CostRampDownH[u,h]) for u in model.u)
                    + sum(model.CostVariable[u,h] * model.Power[u,h] for u in model.u)
                    + sum(model.PriceTransmission[l,h] * model.Flow[l,h] for l in model.l)
                    + sum(model.CostLoadShedding[n,h] * model.ShedLoad[n,h] for n in model.n)
                    + 30000 * sum(model.LostLoad_MaxPower[n,h] + model.LostLoad_MinPower[n,h] for n in model.n)
                    + 20000 * sum(model.LostLoad_Reserve2U[n,h] + model.LostLoad_Reserve2D[n,h] for n in model.n)
                    + 10000 * sum(model.LostLoad_RampUp[u,h] + model.LostLoad_RampDown[u,h] for u in model.u))



###################################################################################################################################################################################################################################################################

# 2) CostStartUp:
###################################################################################################################################################################################################################################################################

    def EQ_CostStartup(model,u,h):
            if (model.CostStartUp[u]!=0):
                    if (h==1):
                        return model.CostStartUpH[u,h]>= model.CostStartUp[u]*(model.Committed[u,h]-model.CommittedInitial[u])
                    elif(h>1):
                        return model.CostStartUpH[u,h]>= model.CostStartUp[u]*(model.Committed[u,h]-model.Committed[u,h-1])
            else:
                return Constraint.Skip

###################################################################################################################################################################################################################################################################

# 3) CostShutDown:
###################################################################################################################################################################################################################################################################

    def EQ_CostShutDown(model,u,h):
            if (model.CostShutDown[u]!=0):
                    if (int(h)==1):
                        return model.CostShutDownH[u,h]>=model.CostShutDown[u]*(model.CommittedInitial[u]-model.Committed[u,h])
                    elif(int(h)>1):
                        return model.CostShutDownH[u,h]>=model.CostShutDown[u]*(model.Committed[u,h-1]-model.Committed[u,h])
            else:
                return Constraint.Skip

###################################################################################################################################################################################################################################################################

# 2) Ramp up cost:
###################################################################################################################################################################################################################################################################

    def EQ_CostRampUp(model,u,h):
            if (model.CostStartUp[u]!=0):
                    if (h==1):
                        return model.CostRampUpH[u,h] >= model.CostRampUp[u]*(model.Power[u,h]-model.PowerInitial[u])
                    elif(h>1):
                        return model.CostRampUpH[u,h] >= model.CostRampUp[u]*(model.Power[u,h]-model.Power[u,h-1])
            else:
                return Constraint.Skip

###################################################################################################################################################################################################################################################################

# 3) Ramp Down cost:
###################################################################################################################################################################################################################################################################

    def EQ_CostRampDown(model,u,h):
            if (model.CostRampDown[u]!=0):
                    if (int(h)==1):
                        return model.CostRampDownH[u,h]>=model.CostRampDown[u]*(model.PowerInitial[u]-model.Power[u,h])
                    elif(int(h)>1):
                        return model.CostRampDownH[u,h]>=model.CostRampDown[u]*(model.Power[u,h-1]-model.Power[u,h])
            else:
                return Constraint.Skip


###################################################################################################################################################################################################################################################################

# 4) Hourly demand balance in the day-ahead market for each node:
###################################################################################################################################################################################################################################################################

    def EQ_Demand_balance_DA(model,n,h):
        return (sum(model.Power[u,h]*model.Location[u,n] for u in model.u)
        +sum(model.Flow[l,h]*model.LineNode[l,n] for l in model.l)
        == model.Demand['DA',n,h]
        +sum(model.StorageInput[s,h]*model.Location[s,n] for s in model.s)
        -model.ShedLoad[n,h]
        -model.LostLoad_MaxPower[n,h]
        +model.LostLoad_MinPower[n,h])

###################################################################################################################################################################################################################################################################

# 5) Demand balance 2U
###################################################################################################################################################################################################################################################################

    def EQ_Demand_balance_2U(model,n,h):
       return sum(model.MaxRamp2U[u,h]*model.Technology[u,t]*model.Reserve[t]*model.Location[u,n] for u in model.u for t in model.t) >= model.Demand['2U',n,h]-model.LostLoad_Reserve2U[n,h]

###################################################################################################################################################################################################################################################################

# 6) Hourly demand balance in the downwards reserve market for each node
###################################################################################################################################################################################################################################################################

    def EQ_Demand_balance_2D(model,n,h):
        return sum(model.MaxRamp2D[u,h]*model.Technology [u,t]*model.Reserve[t]*model.Location[u,n] for u in model.u for t in model.t) >= model.Demand['2D',n,h]-sum(model.StorageChargingCapacity[s]-model.StorageInput[s,h] for s in model.s)*4 - model.LostLoad_Reserve2D[n,h]

###################################################################################################################################################################################################################################################################

# 7) Minimum power output is above the must-run output level for each unit in all periods:
###################################################################################################################################################################################################################################################################

    def EQ_Power_must_run(model,u,h):
        return model.PowerMustRun[u,h]*model.Committed[u,h]<= model.Power[u,h]

###################################################################################################################################################################################################################################################################

# 8) Maximum power output is below the available capacity:
###################################################################################################################################################################################################################################################################

    def EQ_Power_available(model,u,h):
        return model.Power[u,h] <= model.PowerCapacity[u]*model.LoadMaximum[u,h]*model.Committed[u,h]

###################################################################################################################################################################################################################################################################

# 9) Maximum power output with respect to power output in the previous period:
###################################################################################################################################################################################################################################################################

    def EQ_PowerMaximum_previous(model,u,h):
        if all(model.Technology[u,tr] ==0 for tr in model.tr):
            if (h==0):
                return model.Power[u,h] <= model.PowerInitial[u]+ model.RampUpMaximum[u]*model.CommittedInitial[u]+model.RampStartUpMaximum[u]*(1-model.CommittedInitial[u])+model.LostLoad_RampUp[u,h]

            elif(h>0):
                return model.Power[u,h] <= model.Power[u,h-1]+model.RampUpMaximum[u]*model.Committed[u,h-1]+model.RampStartUpMaximum[u]*(1-model.Committed[u,h-1])+model.LostLoad_RampUp[u,h]
        else:
            return Constraint.Skip

###################################################################################################################################################################################################################################################################

# 10) Maximum power output with respect to power output in the following period:
###################################################################################################################################################################################################################################################################

    def EQ_PowerMaximum_following(model,u,h):
          if all(model.Technology[u,tr] ==0 for tr in model.tr):
            if (h<model.h.__len__()-1):
                return model.Power[u,h] <= model.PowerCapacity[u]*model.LoadMaximum[u,h+1]*model.Committed[u,h+1]+ model.RampShutDownMaximum[u]*(1-model.Committed[u,h+1])+model.LostLoad_RampDown[u,h]
            else:
                return Constraint.Skip
          else:
            return Constraint.Skip

###################################################################################################################################################################################################################################################################

# 11) If the unit keeps committed the reduction in power output is lower than the ramp-down limit. If the unit is de-committed the reduction is lower than the shut-down ramp limit:
###################################################################################################################################################################################################################################################################

    def EQ_Ramp_Down(model,u,h):
        if all(model.Technology[u,tr]==0 for tr in model.tr):
            if (h==0):
               return  model.PowerInitial[u]-model.Power[u,h]<= model.PowerCapacity[u]*(1-model.Committed[u,h])+model.RampDownMaximum[u]*model.CommittedInitial[u]+model.LostLoad_RampDown[u,h]

            elif(h>0):
                return model.Power[u,h-1]-model.Power[u,h] <= model.PowerCapacity[u]*(1-model.Committed[u,h])+model.RampDownMaximum[u]*model.Committed[u,h-1]+model.LostLoad_RampDown[u,h]
        else:
            return Constraint.Skip

###################################################################################################################################################################################################################################################################

# 12) EQ_Minimum_time_up_A: maybe error with $
###################################################################################################################################################################################################################################################################

    def EQ_Minimum_time_up_A(model,u):
        for h in model.h:
            if (h<=model.TimeUpLeft_initial[u]):
                return sum(1-model.Committed[u,h] for h in model.h) >= 0
            else:
                return Constraint.Skip


###################################################################################################################################################################################################################################################################

# 13) EQ_Minimum_Time_up_JustStarted
###################################################################################################################################################################################################################################################################

    def EQ_Minimum_Time_up_JustStarted(model,u,h):
        idx_h = range(h,int(model.TimeUpLeft_JustStarted[u,h])+h,1)
        if(model.TimeUpLeft_JustStarted[u,h] > 0 and h>0):
            return sum(model.Committed[u,h] for h in idx_h) >= model.TimeUpLeft_JustStarted[u,h]*model.Committed[u,h] - model.TimeUpLeft_JustStarted[u,h]*model.Committed[u,h-1]
        elif(model.TimeUpLeft_JustStarted[u,h] > 0 and h==0):
            #return sum(model.Committed[u,h] for h in idx_h) >= model.TimeUpLeft_JustStarted[u,h]*(model.Committed[u,h]-model.CommittedInitial[u])
            return sum(model.Committed[u,h] for h in idx_h) >= model.TimeUpLeft_JustStarted[u,h]*model.Committed[u,h]  - model.TimeUpLeft_JustStarted[u,h]*model.CommittedInitial[u]
        else:
            return Constraint.Skip

###################################################################################################################################################################################################################################################################

# 14) EQ_Minimum_time_down_A:
###################################################################################################################################################################################################################################################################

    def EQ_Minimum_time_down_A(model,u):
        for h in model.h:
            if (h<=model.TimeDownLeft_initial[u]):
                return sum(model.Committed[u,h] for h in model.h) >= 0
            else:
                return Constraint.Skip

###################################################################################################################################################################################################################################################################

# 15) EQ_Minimum_Time_down_JustStopped
###################################################################################################################################################################################################################################################################

    def EQ_Minimum_Time_down_JustStopped(model,u,h):
        
        if(model.TimeDownLeft_initial[u]<=h):
          if (model.TimeDownLeft_JustStopped[u,h] > 0 and h==0):
              return sum(1-model.Committed[u,h] for h in range(h,model.TimeDownLeft_JustStopped[u,h]+h,1)) >= model.TimeDownLeft_JustStopped[u,h]*model.CommittedInitial[u]-model.TimeDownLeft_JustStopped[u,h]*model.Committed[u,h]

          elif(model.TimeDownLeft_JustStopped[u,h] > 0 and h>0):
              return sum(1-model.Committed[u,h] for h in range(h,model.TimeDownLeft_JustStopped[u,h]+h,1) ) >= model.TimeDownLeft_JustStopped[u,h]*model.Committed[u,h-1]-model.TimeDownLeft_JustStopped[u,h]*model.Committed[u,h]

          else: return Constraint.Skip
        else:
            return Constraint.Skip

###################################################################################################################################################################################################################################################################

# 16) Maximum 15-min ramping up, in MW/h
###################################################################################################################################################################################################################################################################

    def EQ_Max_RampUp1(model,u,h):
        if all(model.Technology[u,tr]==0 for tr in model.tr):
           return model.MaxRamp2U[u,h]<= model.RampUpMaximum[u]*model.Committed[u,h]+model.FlexibilityUp[u]*(1-model.Committed[u,h])
        else:
            return Constraint.Skip

###################################################################################################################################################################################################################################################################

# 17) Maximum 15-min ramping up, in MW/h:
###################################################################################################################################################################################################################################################################

    def EQ_Max_RampUp2(model,u,h):
        if all(model.Technology[u,tr]==0 for tr in model.tr):
            return model.MaxRamp2U[u,h]<= (model.PowerCapacity[u]*model.LoadMaximum[u,h]-model.Power[u,h])*4
        else:
            return Constraint.Skip

###################################################################################################################################################################################################################################################################

# 18) Maximum 15-min shutting down, in MW/h:
###################################################################################################################################################################################################################################################################

    def EQ_Max_RampDown1(model,u,h):
        if all(model.Technology[u,tr]==0 for tr in model.tr):
            return model.MaxRamp2D[u,h]<= max(float(model.RampDownMaximum[u]),float(model.FlexibilityDown[u]))*model.Committed[u,h]
        else:
            return Constraint.Skip

###################################################################################################################################################################################################################################################################

# 19) Maximum 15-min ramping down, in MW/h:
###################################################################################################################################################################################################################################################################

    def EQ_Max_RampDown2(model,u,h):
        if all(model.Technology[u,tr]==0 for tr in model.tr):
            if (float(model.RampShutDownMaximum[u])< float(model.PowerMinStable[u])*4):
                return model.MaxRamp2D[u,h]<= (model.Power[u,h]- model.PowerMinStable[u]*model.Committed[u,h])*4
            else:
                return model.MaxRamp2D[u,h]<= model.Power[u,h]*4
        else:
            return Constraint.Skip

###################################################################################################################################################################################################################################################################

# 20) Storage level must be above a minimum:
###################################################################################################################################################################################################################################################################

    def EQ_Storage_minimum(model,s,h):
        return model.StorageLevel[s,h] >= model.StorageMinimum[s]

###################################################################################################################################################################################################################################################################

# 21) Storage level must below storage capacity:
###################################################################################################################################################################################################################################################################

    def EQ_Storage_level(model,s,h):
        return model.StorageLevel[s,h] <= model.StorageCapacity[s]

###################################################################################################################################################################################################################################################################

# 22) Storage charging is bounded by the maximum capacity:
###################################################################################################################################################################################################################################################################

    def EQ_Storage_input(model,s,h):
        return model.StorageInput[s,h] <= model.StorageChargingCapacity[s]*(1-model.Committed[s,h])

###################################################################################################################################################################################################################################################################

# 23) Storage balance:
###################################################################################################################################################################################################################################################################

    def EQ_Storage_balance(model,s,h):
        if (h==0):
            return model.StorageInitial[s]+model.StorageInflow[s,h]+model.StorageInput[s,h]*model.StorageChargingEfficiency[s] == model.StorageLevel[s,h]+model.StorageOutflow[s,h]+model.Power[s,h]/max(float(model.StorageDischargeEfficiency[s]),0.0001)

        elif(h>0):
            return model.StorageLevel[s,h-1]+model.StorageInflow[s,h]+model.StorageInput[s,h]*model.StorageChargingEfficiency[s] == model.StorageLevel[s,h]+model.StorageOutflow[s,h]+model.Power[s,h]/max(float(model.StorageDischargeEfficiency[s]),0.0001)

###################################################################################################################################################################################################################################################################

# 24) Assuming cyclic boundary conditions:
###################################################################################################################################################################################################################################################################

    def EQ_Storage_boundaries(model,s,h):
        if (h == len(model.h)-1):
            return model.StorageLevel[s,h]== model.StorageInitial[s]
        else:
            return Constraint.Skip

###################################################################################################################################################################################################################################################################

# 25) Charging is limited by the remaining storage capacity:
###################################################################################################################################################################################################################################################################

    def EQ_Storage_MaxCharge(model,s,h):
        return model.StorageInput[s,h] * model.StorageChargingEfficiency[s]-model.StorageOutflow[s,h]+model.StorageInflow[s,h] <= model.StorageCapacity[s]-model.StorageLevel[s,h]

###################################################################################################################################################################################################################################################################

# 26) Discharge is limited by the storage level:
###################################################################################################################################################################################################################################################################

    def EQ_Storage_MaxDischarge(model,s,h):
        return model.Power[s,h]/max(float(model.StorageDischargeEfficiency[s]),0.0001)+ model.StorageOutflow[s,h]-model.StorageInflow[s,h] <= model.StorageLevel[s,h]

###################################################################################################################################################################################################################################################################

# 27) Flows are above minimum values:
###################################################################################################################################################################################################################################################################

    def EQ_Flow_limits_lower(model,l,h):
        return model.Flow[l,h] >= model.FlowMinimum[l,h]

###################################################################################################################################################################################################################################################################

# 28) Flows are below maximum values:
###################################################################################################################################################################################################################################################################

    def EQ_Flow_limits_upper(model,l,h):
        return model.Flow[l,h] <= model.FlowMaximum[l,h]

###################################################################################################################################################################################################################################################################

# 29) Force Unit Commitment/decommitment:
    # E.g: renewable unit with AF>0 must be committed:
###################################################################################################################################################################################################################################################################

    def EQ_Force_Commitment(model,u,h):
        if any(model.Technology[u,tr]>=1 for tr in model.tr) and model.LoadMaximum[u,h]>0:
            return model.Committed[u,h] == 1
        else:
            return Constraint.Skip

###################################################################################################################################################################################################################################################################

# 30) E.g: renewable unit with AF = 0 must be decommitted:
###################################################################################################################################################################################################################################################################

    def EQ_Force_DeCommitment(model,u,h):
        if(model.LoadMaximum[u,h]==0):
            return model.Committed[u,h] == 0
        else:
            return Constraint.Skip

###################################################################################################################################################################################################################################################################

# 31) Load shedding:
###################################################################################################################################################################################################################################################################
    def EQ_LoadShedding(model,n,h):
            return model.ShedLoad[n,h] <= model.LoadShedding[n]

###################################################################################################################################################################################################################################################################

#######################################################################################################################
######################################## Definition of model ##########################################################
#######################################################################################################################


    model.EQ_Objective_Function = Objective(rule=EQ_Objective_function, sense = minimize)
    
    model.EQ_SystemCost = Constraint(model.h,rule=EQ_SystemCost)
    
    if not Mixed_Integer_LP:   
        model.EQ_CostStartup = Constraint(model.u,model.h,rule=EQ_CostStartup)
        model.EQ_CostShutDown = Constraint(model.u,model.h,rule=EQ_CostShutDown)
        model.EQ_CostRampUp = Constraint(model.u,model.h,rule=EQ_CostRampUp)
        model.EQ_CostRampDown = Constraint(model.u,model.h,rule=EQ_CostRampDown)


    model.EQ_Demand_balance_DA = Constraint(model.n,model.h,rule = EQ_Demand_balance_DA)
    model.EQ_Demand_balance_2U = Constraint(model.n,model.h,rule = EQ_Demand_balance_2U)
    model.EQ_Demand_balance_2D = Constraint(model.n,model.h,rule = EQ_Demand_balance_2D)



    if not Mixed_Integer_LP:
        model.EQ_Power_must_run = Constraint(model.u,model.h,rule=EQ_Power_must_run)

    model.EQ_Power_available=Constraint(model.u,model.h,rule=EQ_Power_available)

    model.EQ_PowerMaximum_previous=Constraint(model.u,model.h,rule=EQ_PowerMaximum_previous)
    model.EQ_PowerMaximum_following=Constraint(model.u,model.h,rule=EQ_PowerMaximum_following)

    model.EQ_Ramp_Down = Constraint(model.u,model.h,rule=EQ_Ramp_Down)

    if not Mixed_Integer_LP:
        model.EQ_Minimum_time_up_A = Constraint(model.u,rule=EQ_Minimum_time_up_A)
        model.EQ_Minimum_Time_up_JustStarted=Constraint(model.u,model.h,rule=EQ_Minimum_Time_up_JustStarted)

        model.EQ_Minimum_time_down_A = Constraint(model.u,rule=EQ_Minimum_time_down_A)
        model.EQ_Minimum_Time_down_JustStopped=Constraint(model.u,model.h,rule=EQ_Minimum_Time_down_JustStopped)

    model.EQ_Max_RampUp1 = Constraint(model.u,model.h,rule=EQ_Max_RampUp1)
    model.EQ_Max_RampUp2 = Constraint(model.u,model.h,rule=EQ_Max_RampUp2)

    model.EQ_Max_RampDown1 = Constraint(model.u,model.h,rule=EQ_Max_RampDown1)
    model.EQ_Max_RampDown2 = Constraint(model.u,model.h,rule=EQ_Max_RampDown2)

    model.EQ_Storage_minimum = Constraint(model.s,model.h,rule=EQ_Storage_minimum)
    model.EQ_Storage_level= Constraint(model.s,model.h,rule=EQ_Storage_level)
    model.EQ_Storage_input=Constraint(model.s,model.h,rule=EQ_Storage_input)
    model.EQ_Storage_balance=Constraint(model.s,model.h,rule=EQ_Storage_balance)
    model.EQ_Storage_boundaries=Constraint(model.s,model.h,rule=EQ_Storage_boundaries)
    model.EQ_Storage_MaxCharge = Constraint(model.s,model.h,rule=EQ_Storage_MaxCharge)
    model.EQ_Storage_MaxDischarge = Constraint(model.s,model.h,rule=EQ_Storage_MaxDischarge)
    model.EQ_Flow_limits_upper = Constraint(model.l,model.h,rule=EQ_Flow_limits_upper)
    model.EQ_Flow_limits_lower = Constraint(model.l,model.h,rule=EQ_Flow_limits_lower)
    model.EQ_Force_Commitment = Constraint(model.u,model.h,rule=EQ_Force_Commitment)
    model.EQ_Force_DeCommitment = Constraint(model.u,model.h,rule=EQ_Force_DeCommitment)
    model.EQ_LoadShedding = Constraint(model.n,model.h,rule =EQ_LoadShedding)

    return model




def run_solver(instance, solver="cplex", solver_manager="serial", options_string=""):
    # initialize the solver / solver manager.
    solver = SolverFactory(solver)
    if solver is None:
        logging.critical( "Solver %s is not available on this machine." % solver)
        sys.exit(1)
    solver_manager = SolverManagerFactory(solver_manager) #serial or pyro

    results = solver.solve(instance, options_string=options_string, tee=True)
    if results.solver.termination_condition != TerminationCondition.optimal:
    # something went wrong
        logging.warn("Solver: %s" % results.solver.termination_condition)
        logging.debug(results.solver)                                                                   
    else:
        logging.info("Solver: %s" % results.solver.termination_condition)
        instance.solutions.load_from(results)                                                                      
        return instance
   
def after_solver(results):
    pass


#######################################################################################################################
############################### Main wrapper (rolling horizon optimization) ###########################################
#######################################################################################################################

def DispaSolve(sets, parameters, Mixed_Integer_LP=False):
    '''
    The DispaSolve function defines the rolling horizon optimization and saves each result variable in a pandas dataframe
    The definition of the rolling horizon must be included into the DispaSET Config parameter'
    
    :param sets: Dictionary containing the sets (defined as a list of strings or integers)
    :param parameters: Dictionary containing the parameters of the optimization problem (in the DispaSET 2.1.1 format)
    
    :return: Dictionary of pandas dataframes with the optimization variables
    '''    

    # Initialize the results dictionnary:

    # Load the config parameter in the pyomo format (easier to read):
    config = pyomo_format(sets,parameters['Config'])
    
    # Time parameters:
    Nhours = len(sets['h'])
    
    #Build pandas indexes based on the config variables:
    first = pd.datetime(config['FirstDay','year'],config['FirstDay','month'],config['FirstDay','day'],0,0,0)
    last = pd.datetime(config['LastDay','year'],config['LastDay','month'],config['LastDay','day'],23,59,59)
    start = pd.datetime(config['DayStart','year'],config['DayStart','month'],config['DayStart','day'],0,0,0)
    stop = pd.datetime(config['DayStop','year'],config['DayStop','month'],config['DayStop','day'],23,59,59)
    
    # Index corresponding to the data:
    index_all = pd.DatetimeIndex(start=first,end=last,freq='h')
    if len(index_all) != Nhours:
        sys.exit('The interval length in the Config parameter (' + str(len(index_all) + ' does not correspond to the number of hours in the data (' + str(Nhours)))
    
    # Index of the selected slice for simulation:
    index_sim = pd.DatetimeIndex(start=start,end=stop,freq='h')
    Nhours_sim = len(index_sim)
    Ndays = Nhours_sim/24
    
    Nunits = len(sets['u'])
    parameters['Demand']['val'][1:2,:,:] = np.zeros(parameters['Demand']['val'][1:2,:,:].shape)
    
    # Pre-processing of model parameters:
        # Forecasted upwards reserve margin (UCTE). Only if not provided in the parameters
        #Forecasted downwards reserve margin (UCTE)
    if (parameters['Demand']['val'][1:2,:,:] == 0).all():   # only applys if all upward/downard reserve requirement are zero!!!
        shape = parameters['Demand']['val'][1,:,:].shape
        maximum = np.max(parameters['Demand']['val'])
        parameters['Demand']['val'][1,:,:] = sqrt(10*maximum+ pow(150,2))-150 * np.ones(shape)     # reserve up
        parameters['Demand']['val'][2,:,:] = 0.5 * (sqrt(10*maximum+ pow(150,2))-150 * np.ones(shape))   # reserve down
    
    # Definition of the minimum stable load
    parameters['PowerMinStable'] = {'sets':['u'],'val':parameters['PowerCapacity']['val'] * parameters['PartLoadMin']['val']}
    
    # Update of the Availability factor, taking into account the outage factor:
    parameters['LoadMaximum'] = {'sets':['u','h'],'val':parameters['AvailabilityFactor']['val'] * (1-parameters['OutageFactor']['val'])}
    
    # Start-up and Shutdown ramping constraints (minimum value is the PowerMinStable)
    parameters['RampShutDownMaximum']['val'] = np.maximum(parameters['RampShutDownMaximum']['val'], parameters['PowerMinStable']['val'])
    parameters['RampStartUpMaximum']['val'] = np.maximum(parameters['RampStartUpMaximum']['val'], parameters['PowerMinStable']['val'])
    
    # If the plant is stopped, its 15-min ramp-up capability is RampStartUpMaximum if it can start in this timeframe
    parameters['FlexibilityUp'] = {'sets':parameters['RampShutDownMaximum']['sets'],'val':np.zeros(parameters['RampStartUpMaximum']['val'].shape)}
    flexible_plants_up = parameters['RampStartUpMaximum']['val'] >= 4 * parameters['PowerMinStable']['val']
    parameters['FlexibilityUp']['val'][flexible_plants_up] = parameters['RampStartUpMaximum']['val'][flexible_plants_up]
    
    # If the plant is started, its 15-min ramp-down capability is either RampShutDownMaximum if it is fast enough, or RampDownMaximum otherwise
    parameters['FlexibilityDown'] = {'sets':parameters['RampShutDownMaximum']['sets'],'val':np.zeros(parameters['RampShutDownMaximum']['val'].shape)}
    flexible_plants_down = parameters['RampShutDownMaximum']['val'] >= 4 * parameters['PowerMinStable']['val']
    parameters['FlexibilityDown']['val'][flexible_plants_down] = parameters['RampShutDownMaximum']['val'][flexible_plants_up]
    
    # Set the ramping costs to zero if not defined:
    if not 'CostRampUp' in parameters:
       parameters['CostRampUp'] = {'sets':parameters['PowerCapacity']['sets'],'val':np.zeros(parameters['PowerCapacity']['val'].shape)}
    if not 'CostRampDown' in parameters:
       parameters['CostRampDown'] = {'sets':parameters['PowerCapacity']['sets'],'val':np.zeros(parameters['PowerCapacity']['val'].shape)}
     
    # Initialize CostLoadShedding:
    parameters['CostLoadShedding'] = {'sets':['n','h'],'val':1000 * np.ones([len(sets['n']),Nhours])}
    
    #calculate the minimum must run power for renewable technologies
    parameters['PowerMustRun'] = {'sets':['u','h'],'val':np.zeros([Nunits,Nhours])}
    for u in range(Nunits):
        # find technology:
        tech = sets['t'][np.where(parameters['Technology']['val'][u,:] == 1)[0][0]]
        loc = np.where(parameters['Location']['val'][u,:] == 1)[0][0]
        curt = parameters['Curtailment']['val'][loc]
        if tech in sets['tr'] and curt != 1:
            parameters['PowerMustRun']['val'][u,:] = parameters['PowerCapacity']['val'][u] * parameters['AvailabilityFactor']['val'][u,:]
    
    # Converting boolean array to integers (pyomo does not like booleans)
    for p in parameters:
        if 'val' in parameters[p]:
            if parameters[p]['val'].dtype == np.dtype('bool'):
                parameters[p]['val'] = np.array(parameters[p]['val'],dtype='int')
    
    range_start = index_all.get_loc(start)
    days = range(range_start/24,Nhours_sim/24-1 - config['RollingHorizon LookAhead','day'],config['RollingHorizon Length','day'])
    
    results = {}

    for d in days:
        range_start = d*24
        range_stop = np.minimum(range_start + config['RollingHorizon Length','day']*24 + config['RollingHorizon LookAhead','day']*24,Nhours)
        h_range = range(range_start,range_stop,1)       # Time indexes of for the current horizon
        h_kept = range(range_start,range_stop - config['RollingHorizon LookAhead','day']*24)        # Time indexes of the values that will be kept and stored
        index_range = index_all[h_range]
        index_kept = index_all[h_kept]
        logging.info('Optimizing time interval ' + str(index_range[0]) + ' to ' + str(index_range[-1]))
        logging.info('Conserving only the interval ' + str(index_kept[0]) + ' to ' + str(index_kept[-1]))
        # Slice the time-dependent parameters to the right horizon: 
        parameters_sliced = {}
        for var in parameters:
            if parameters[var]['sets'][-1] == 'h':      # if the last set of the parameter is time
                dim = len(parameters[var]['sets'])         
                var_sliced = {}
                var_sliced['sets'] = parameters[var]['sets']
                if dim == 1:
                    var_sliced['val'] = parameters[var]['val'][h_range]
                elif dim == 2:
                    var_sliced['val'] = parameters[var]['val'][:,h_range]
                elif dim == 3:
                    var_sliced['val'] = parameters[var]['val'][:,:,h_range]
                else:
                    sys.exit('Variables with more than 3 dimensions are not allowed. ' + parameters[var]['name'] + ' has ' + str(dim) + ' Dimensions')
                parameters_sliced[var] = var_sliced
            else:
                parameters_sliced[var] = parameters[var]
            
        # Copy the sets dictionnary and slice the 'h' set:
        sets_sliced = sets.copy()
        sets_sliced['h'] = sets_sliced['h'] =  range(len(h_range))
            
        # Add the parameters specific to the optimization of a single time horizon:
        #parameters_sliced['PowerInitial'] = 
        parameters_sliced['CommittedInitial'] = {'sets':['u'],'val': parameters_sliced['PowerInitial']['val'] > 0}
            
        parameters_sliced['TimeUpLeft_initial'] = {'sets':['u'],'val':parameters_sliced['TimeUpMinimum']['val'] - parameters_sliced['TimeUpInitial']['val']*parameters_sliced['CommittedInitial']['val']}
        parameters_sliced['TimeDownLeft_initial'] = {'sets':['u'],'val':parameters_sliced['TimeDownMinimum']['val'] - parameters_sliced['TimeDownInitial']['val']*(1 - parameters_sliced['CommittedInitial']['val'])}
        
        parameters_sliced['TimeUpLeft_JustStarted'] = {'sets':['u','h'],'val':np.zeros([Nunits,len(h_range)],dtype='int')}
        parameters_sliced['TimeDownLeft_JustStopped'] = {'sets':['u','h'],'val':np.zeros([Nunits,len(h_range)],dtype='int')}
        for u in range(Nunits):
            parameters_sliced['TimeUpLeft_JustStarted']['val'][u,:] = np.minimum(len(h_range) - 1 - np.arange(len(h_range)),parameters_sliced['TimeUpMinimum']['val'][u]*np.ones(len(h_range))).astype('int')
            parameters_sliced['TimeDownLeft_JustStopped']['val'][u,:] = np.minimum(len(h_range) - 1 - np.arange(len(h_range)),parameters_sliced['TimeDownMinimum']['val'][u]*np.ones(len(h_range))).astype('int')
    
        # Optimize: 
            
        instance = DispOptim(sets_sliced, parameters_sliced, Mixed_Integer_LP)
        if Mixed_Integer_LP:
            opt = run_solver(instance, options_string="mipgap=0.0")
        else:
            opt = run_solver(instance)
           
        
        results_sliced = {}
        # TDO Iterate all VARs instead of listing everything. Can we?
        #for v in model.component_objects(Var):
        #    results_sliced[v] = pyomo_to_pandas(opt, v)

        results_sliced['Committed'] = pyomo_to_pandas(opt,'Committed')
        results_sliced['CostStartUpH'] = pyomo_to_pandas(opt,'CostStartUpH')
        results_sliced['CostShutDownH'] = pyomo_to_pandas(opt,'CostShutDownH')
        results_sliced['CostRampUpH'] = pyomo_to_pandas(opt,'CostRampUpH')
        results_sliced['CostRampDownH'] = pyomo_to_pandas(opt,'CostRampDownH')
        results_sliced['Flow'] = pyomo_to_pandas(opt,'Flow')
        results_sliced['Power'] = pyomo_to_pandas(opt,'Power')
        results_sliced['ShedLoad'] = pyomo_to_pandas(opt,'ShedLoad')
        results_sliced['StorageInput'] = pyomo_to_pandas(opt,'StorageInput')
        results_sliced['StorageLevel'] = pyomo_to_pandas(opt,'StorageLevel')
        results_sliced['SystemCost'] = pyomo_to_pandas(opt,'SystemCost')
        results_sliced['LostLoad_MaxPower'] = pyomo_to_pandas(opt,'LostLoad_MaxPower')
        results_sliced['LostLoad_RampUp'] = pyomo_to_pandas(opt,'LostLoad_RampUp')
        results_sliced['LostLoad_RampDown'] = pyomo_to_pandas(opt,'LostLoad_RampDown')
        results_sliced['LostLoad_MinPower'] = pyomo_to_pandas(opt,'LostLoad_MinPower')
        results_sliced['LostLoad_Reserve2U'] = pyomo_to_pandas(opt,'LostLoad_Reserve2U')
        results_sliced['LostLoad_Reserve2D'] = pyomo_to_pandas(opt,'LostLoad_Reserve2D')
        # Defining the main results dictionnary:
        if len(results) == 0:
            for r in results_sliced:
                results[r] = pd.DataFrame(index=index_sim, columns=results_sliced[r].columns)

        # Adding the sliced results to the main results dictionnary:
        for r in results_sliced:
            results_sliced[r].index=index_range
            results[r].loc[index_kept, :] = results_sliced[r].loc[index_kept, :] 
            
        # Calculating the times up and down of each power plant
        TimeUp = np.zeros(Nunits)
        TimeDown = np.zeros(Nunits)
        for u in range(Nunits):
            for i in index_kept.order(ascending=False):   # take the inverse #FIXME new pandas : index.kept.sort_values() 
                if results['Power'].loc[i, sets['u'][u]] > 0:
                    TimeUp[u] += 1
                else:
                    break
            for i in index_kept.order(ascending=False):   # take the inverse #FIXME new pandas : index.kept.sort_values() 
                if results['Committed'].loc[i, sets['u'][u]] == 0:
                    TimeDown[u] += 1
                else:
                    break    
            if TimeUp[u] == len(h_kept):
                TimeUp[u] += parameters['TimeUpInitial']['val'][u]
            if TimeDown[u] == len(h_kept):
                TimeDown[u] += parameters['TimeDownInitial']['val'][u]        
            
        # Updating the initial values for the next optimization:
        parameters['StorageInitial']['val'] = results['StorageLevel'].loc[index_kept[-1], :].values.astype('float')
        parameters['PowerInitial']['val'] = results['Power'].loc[index_kept[-1], :].values.astype('float')
        parameters['TimeUpInitial']['val'] = TimeUp
        parameters['TimeDownInitial']['val'] = TimeDown
    
    #    for v in opt.component_objects(Var, active=True):
    #        if v.name not in results:
    #            results[v.name] = pd.DataFrame(index=index_all)
    #        
    
        # Clearing the optimization variables for this time horizon
        del opt,parameters_sliced,sets_sliced
    
    return results
