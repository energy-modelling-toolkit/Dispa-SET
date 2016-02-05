########################################################################################################################
############################################ Pre-processing ############################################################
########################################################################################################################
__author__ = 'Sylvain Quoilin & Valentin Kaulman'


# This script is the Python pre-processing tool for the JRC Dispa-SET 2.0 model
# The output is a gdx file readable in GAMS
# It requires the wgdx.py utility, developed by S. Quoilin and the gdxx.py utility, provided in GAMS
# To install the gdxx library (for Python 2.7):
# C:\GAMS\win64\24.0\apifiles\Python\api>python gdxsetup.py install

import time
tic = time.clock()

import sys

sys.path.append('../python-files')

import numpy as np
import os
import shutil
import matplotlib.pyplot as plt
from DispaSET_io_data import *
from DispaTools import *

###################################################################################################################
#####################################   Main Inputs    ############################################################
###################################################################################################################

# Paths:
gams_dir = '/home/sylvain/progs/GAMS/matlab'
#gams_dir = 'C:\\GAMS\\win32\\24.0'

path_data_csv = 'BE-2014/csv'
path_data_pandas = 'BE-2014/pandas'

#gdx_out = "Inputs.gdx"

# Create simulation test case environment:
# Provide the path of the folder to be created. Leave empty if no simulation environment should be created.
description = 'Scenario-8-2014'

# Residual binary variable:
# if True, no power is assigned to VRE, VRE production is deducted from the load
# if False, the total load is provided to Dispaset, VRE are considered as must-run power plants (but with some curtailment allowed)
use_residual_load =False

# additional VRE, storage and demand (set Pcap to -1 to use historical values):
Pcap_sol = 5000                 # power of PV
Pcap_wind = 5000              # power of wind
Pcap_sto = -1                      # Power of pumped hydro storage
hours_sto =5                   # Number of hours of storage at full capacity
coeff_charge = 1                        # Load scaling factor

# Day/hour corresponding to the first line of the data:
y_start = 2014
m_start = 1
d_start = 1
h_start = 0    # Hours are defined from 0 to 23

# Number of hours to be simulated:
Ndays =10
Nhours = 24*Ndays

# boolean option to merge plants or not
merge_plants = True

# Allow Curtailment:
curtailment=True


###################################################################################################################
#####################################   Data Loading    ###########################################################
###################################################################################################################

# Load plant data:
plants = load_csv_to_pd(path_data_csv,'production_park_data.csv',path_data_pandas,'production_park_data')

# Load demand data:
load_data_pd = load_csv_to_pd(path_data_csv,'load_data.csv',path_data_pandas,'load_data')
load_data= np.array(load_data_pd['Vertical load [kW]'])

# Outage data:
outages_pd = load_csv_to_pd(path_data_csv,'outages.csv',path_data_pandas,'outages')
outages= np.array(outages_pd)[:,1:]

# Load the timestamps in string format:
dates_pd = load_csv_to_pd(path_data_csv,'daysofyear.csv',path_data_pandas,'daysofyear')
dates= dates_pd.Date.tolist()

# Load the historical generation data:
generation_data_pd = load_csv_to_pd(path_data_csv,'generation_data.csv',path_data_pandas,'generation_data')
generation_data= np.array(generation_data_pd)

# Load load_forecast
load_forecast_pd = load_csv_to_pd(path_data_csv,'load_forecast.csv',path_data_pandas,'load_forecast')
load_forecast= np.array(load_forecast_pd)

# Load interconnections
inter_fr_pd = load_csv_to_pd(path_data_csv,'interconnections_fr.csv',path_data_pandas,'interconnections_fr')
inter_nl_pd = load_csv_to_pd(path_data_csv,'interconnections_nl.csv',path_data_pandas,'interconnections_nl')
interconnections = {'fr':np.array(inter_fr_pd['Power exchange ']),'nl':np.array(inter_nl_pd['Power exchange'])}

# Load solar data:
solar_data_pd = load_csv_to_pd(path_data_csv,'solar_data.csv',path_data_pandas,'solar_data')
solar_data= np.array(solar_data_pd)

# Load wind data:
wind_data_pd = load_csv_to_pd(path_data_csv,'WindData.csv',path_data_pandas,'WindData')
wind_data= np.array(wind_data_pd)

# Fuel price data from IEA, monthly, in eur/MWhth
# From october 2012 to september 2013
price = {'coal_month':np.array([12.2430352051,12.2430352051,12.2430352051,10.2058479078,10.2058479078,10.2058479078,10.5107415473,10.5107415473,10.5107415473,9.8076194807,9.8076194807,9.8076194807])}
price['oil_month'] = np.array([73.7754656768,73.7754656768,73.7754656768,69.6732200066,69.6732200066,69.6732200066,67.2718323463,67.2718323463,67.2718323463,79.5365975377,79.5365975377,79.5365975377])
price['gas_month'] = np.array([28.85,28.85,28.85,30.05,30.05,30.05,29.9,29.9,29.9,30.2,30.2,30.2])



###################################################################################################################
#####################################   Data Formatting   #########################################################
###################################################################################################################


Nunits = len(plants['unit'])

# Calculating the starting row and hour:
# Find the starting row:
str_query = '(year == ' + str(y_start) + ') & (month == ' + str(m_start) + ') & (day == ' + str(d_start) + ') & (hour == ' + str(h_start) + ')'
find_rows = load_data_pd.query(str_query)
if len(find_rows)==0:
    print 'The specified starting date could not be found in the data'
else:
    row_start = find_rows.index.values[0]

hour_start = (row_start)/4 + 1
day_start = (hour_start - 1)/24 + 1
row_end = row_start + Nhours * 4

# Shrink outages to the considered days and turn into a list:
outages = outages[:,day_start-1:day_start + Ndays-1].tolist()

# Shrink the 15-min data to the considered period, convert kW to MW when necessary:
dates = dates[row_start:row_end]
generation_data = generation_data[row_start:row_end,:]/1000
load_data = load_data[row_start:row_end]/1000
load_forecast = load_forecast[row_start:row_end,:]/1000
interconnections['fr'] = interconnections['fr'][row_start:row_end]/1000
interconnections['nl'] = interconnections['nl'][row_start:row_end]/1000
solar_data = solar_data[row_start:row_end,:]
wind_data = wind_data[row_start:row_end,:]/1000

# Transforming the fuel price (monthly data) into hourly data:
x = np.array(range(1,13))
xvals = np.linspace(1,12,8760)
price['oil_hour'] = np.interp(xvals,x,price['oil_month'])
price['gas_hour'] = np.interp(xvals,x,price['gas_month'])
price['coal_hour'] = np.interp(xvals,x,price['coal_month'])
# shrink to the considered hours:
price['oil'] = price['oil_hour'][hour_start-1:hour_start+Nhours-1]
price['gas'] = price['gas_hour'][hour_start-1:hour_start+Nhours-1]
price['coal'] = price['coal_hour'][hour_start-1:hour_start+Nhours-1]

###################################################################################################################
#####################################   Merging Power plants    ###################################################
###################################################################################################################
# Merge excessively disaggregated power units.

# Definition of the simulation timestep (in hours). All effects (e.g. ramping, min
# times, etc) faster than the time step can be neglected.
timestep = 1

# Defining slices in the power plant data to categorize them:
N=20        # number of slices to consider (a lower value means that the plants will be aggregated more easily)
bounds={'part-load':np.linspace(0,1,N),'rampup':np.linspace(0,1,N),'rampdown':np.linspace(0,1,N),'startup':mylogspace(0,36,N),'minuptime':mylogspace(0,168,N),'mindowntime':mylogspace(0,168,N),'noloadcost':np.linspace(0,50,N),'startupcost':np.linspace(0,500,N),'efficiency':np.linspace(0,1,N)}

# Definition of the mapping variable, from the old power plant list the new (merged) one:
map_old_new = np.zeros(Nunits)
map_plant_orig = []

# Definition of the fingerprint value of each power plant, which is the pattern of the slices number in which each of
# its characteristics falls:
fingerprints = []
fingerprints_merged = []
for i in range(Nunits):
    fingerprints.append([find_nearest(bounds['part-load'],plants['part-load'][i]),find_nearest(bounds['rampup'],plants['rampup'][i]),find_nearest(bounds['rampdown'],plants['rampdown'][i]),find_nearest(bounds['startup'],plants['startup'][i]),find_nearest(bounds['minuptime'],plants['minuptime'][i]),find_nearest(bounds['mindowntime'],plants['mindowntime'][i]),find_nearest(bounds['noloadcost'],plants['noloadcost'][i]),find_nearest(bounds['startupcost'],plants['startupcost'][i]),find_nearest(bounds['efficiency'],plants['efficiency'][i])])

# Definition of the merged power plants dicitonnary:
plants_merged = {'unit':[],'type':[],'fuel':[],'pmax':[],'part-load':[],'rampup':[],'rampdown':[],'startup':[],'minuptime':[],'mindowntime':[],'heatrate':[],'noloadcost':[],'startupcost':[],'variablecost':[],'efficiency':[]}
outages_merged = []

for i in range(Nunits):
    merged=False
    for j in range(len(plants_merged['unit'])):
        same_type = ((plants['type'][i]==plants_merged['type'][j]) and (plants['fuel'][i]==plants_merged['fuel'][j]))
        same_fingerprint = (fingerprints[i]==fingerprints_merged[j])
        low_pmin = (plants['part-load'][i] <= 0.01 )
        highly_flexible = (plants['rampup'][i] > 1/60 and (plants['rampdown'][i] > 1/60) and (plants['startup'][i]<1))
        if same_type and same_fingerprint and low_pmin and merge_plants or same_type and highly_flexible and merge_plants:     # merge the two plants in plants_merged:
            P_old = plants_merged['pmax'][j]                    # Old power in plants_merged
            P_add = plants['pmax'][i]                           # Additional power to be added
            for key in plants_merged:
                if key in ['rampup','rampdown','minuptime','mindowntime','heatrate','noloadcost','startupcost','variablecost','efficiency']:
                    # Do a weighted average:
                    plants_merged[key][j] = (plants_merged[key][j] * P_old + plants[key][i] * P_add)/(P_add + P_old)
                elif key in ['pmax']:
                    # Do a sum:
                    plants_merged[key][j] = plants_merged[key][j] + plants[key][i]
                elif key in ['part_load']:
                    #The minimum load is the minimum load of the smallest plant
                    plants_merged[key][j] = 0
            map_old_new[i]=j
            map_plant_orig[j].append(i)
            outages_merged[j]=[outages[i][k]*P_add/(P_old+P_add) + outages_merged[j][k]*P_old/(P_old+P_add) for k in range(len(outages[i]))]
            merged = True
            break
    if not merged:       # Add a new plant in plants_merged:
        for key in plants_merged:
            plants_merged[key].append(plants[key][i])
        map_plant_orig.append([i])
        map_old_new[i]=len(map_plant_orig)-1
        outages_merged.append(outages[i])
        fingerprints_merged.append(fingerprints[i])

# Updating Nunits to the new value:
Nunits = len(plants_merged['unit'])

# Modify the unit names with the original index number. In case of merged plants, indicate all indexes + the plant type and fuel
for j in range(Nunits):
    if len(map_plant_orig[j])==1:                                     # The plant has not been merged
        plants_merged['unit'][j] = str(map_plant_orig[j]) + ' - ' + plants_merged['unit'][j]
    else:
        plants_merged['unit'][j] = str(map_plant_orig[j]) + ' - ' + plants_merged['type'][j] + ' - ' + plants_merged['fuel'][j]



###################################################################################################################
#####################################   Residual load     #########################################################
###################################################################################################################
# This section computes the Day-ahead net load and the actual net load
# The amount of VRE present in the system can be increased using by defining an additional installed power.

# Load
# The load can be defined in different ways, depending on the inclusion or not of VRE, interconnections, losses, etc.
# In this program, the following definitions is applied:
# load_tot: total Vertical load, i.e. the total consumption seen on the transmission network
# load_int: vertical load minus interconnections
# load_int_sol: vertical load minus interconnections and PV
# load_int_sol_wind: vertical load minus interc., PV and wind
# In the Elia load data, the non-monitored VRE is already substracted. It must therefore be added back.
# in the following, residual load is used for the load with additional VRE substracted.

# Historical installed capacity of VRE:
Pmax_sol = solar_data[:,4]
Pmax_wind = wind_data[:,2]

# PV generation:
solar_data[0:4451,0] = solar_data[0:4451,1]             # Day ahead values are not available for the first month, we take the morning predictions
AF_sol_15min = solar_data[:,2]/solar_data[:,4]          # Availability factor is obtained by dividing by the historical installed capacity of VRE:
AF_sol_15min_forecast = solar_data[:,0]/solar_data[:,4]

if Pcap_sol < 0:                                        # If requested PV capacity is negative, use historical values
    Psol_forecast = solar_data[:,0]
    Psol = solar_data[:,2]
else:                                                   # Else scale the Availability factor
    Psol_forecast = Pcap_sol*AF_sol_15min_forecast
    Psol = Pcap_sol*AF_sol_15min

error_sol = Psol - Psol_forecast

# Wind generation
AF_wind_15min = wind_data[:,1]/wind_data[:,2]
AF_wind_15min_forecast = wind_data[:,0]/wind_data[:,2]


if Pcap_wind < 0:
    Pwind_forecast = wind_data[:,0]                     # Forecast
    Pwind = wind_data[:,1]                              # Historical value
else:                                                   # if Pcap is defined, scale AF up
    Pwind_forecast = Pcap_wind*AF_wind_15min_forecast
    Pwind = Pcap_wind*AF_wind_15min_forecast

error_wind = Pwind - Pwind_forecast


# Water generation
# The water data include the hydro storage in turbine mode.
# To find the run-of-water production, the minimum production is taken for each day.
Pwater = generation_data[:,5]

for i in np.r_[0:len(Pwater):96]:
    P = np.sort(Pwater[i:i+95])
    P = np.minimum(P,80)            # P cannot be higher than the installed power
    P_mean = np.mean(P[1:4*4])    # taking the average of the 12 hours with less water generation
    Pwater[i:i+95]=P_mean
AF_water_temp = Pwater.astype(float)/80

# Non dispatchable plants:
P_nondispatch = 0
for i in range(len(plants_merged['type'])):
    if plants_merged['type'][i] in ['IS']:
        P_nondispatch = P_nondispatch + plants_merged['pmax'][i]
AF_nondispatch = 0.5 * np.ones(Nhours*4)                            # Hypothesis, in the absence of information...

# Interconnections:
export_net = interconnections['fr'][:] + interconnections['nl'][:]

load_tot = load_data[:]-generation_data[:,6] + solar_data[:,2] + wind_data[:,1]       # DSO-connected VRE is added to the load
load_int = load_tot + export_net
load_int_VRE = load_int - Psol - Pwind - Pwater - P_nondispatch*AF_nondispatch

net_load_forecast = (load_forecast[:,5]-generation_data[:,6])*coeff_charge
net_load = (load_data[:]-generation_data[:,6])*coeff_charge

error_load = (net_load - net_load_forecast)
error_rel_load = error_load/load_forecast[:,5]

generation = generation_data[:,0] - generation_data[:,6]    # wind is removed from the total generation


## Residual load
# Computation of a new net load, as if there was no TSO or DSO-connected VRE
# Avoiding negative loads:
residual_load = np.maximum(0,load_int_VRE)

# Total production by VRE
VRE_production = np.sum(Pwind + Psol)/4

# Transform 15-min data into hourly data:
demand_res = np.zeros([Nhours,1])
demand_tot = np.zeros([Nhours,1])
demand_wo_int = np.zeros([Nhours,1])
AF_water = np.zeros(Nhours)
AF_wind = np.zeros(Nhours)
AF_sol = np.zeros(Nhours)
dates_h = []
dates_h_ymd = []
for i in range(Nhours):
    demand_res[i] = np.mean(residual_load[1 + i*4 : (i+1)*4])
    demand_tot[i] = np.mean(load_tot[1 + i*4 : (i+1)*4])
    demand_wo_int[i] = np.mean(load_int[1 + i*4 : (i+1)*4])
    AF_water[i] = np.mean(AF_water_temp[1 + i*4 : (i+1)*4])
    AF_wind[i] = np.mean(AF_wind_15min[1 + i*4 : (i+1)*4])
    AF_sol[i] = np.mean(AF_sol_15min[1 + i*4 : (i+1)*4])
    dates_h.append(dates[i*4].encode('ascii','ignore'))
    dates_h_ymd.append(dates_h[i][:11])


###################################################################################################################
#####################################   Reserve needs     #########################################################
###################################################################################################################

# Since the time step is one hour, we assume that the required 15-min
# ramping needs are integrated into the reserves.

rampup15 = np.zeros(Nhours)
rampdown15 = np.zeros(Nhours)

rampup15[0] = max(residual_load[1:3] - residual_load[0:2])*4
rampdown15[0] = max(-residual_load[1:3] + residual_load[0:2])*4

for i in range(Nhours-1):
    rampup15[i+1] = max(residual_load[i*4+1:i*4+4] - residual_load[i*4:i*4+3])*4
    rampdown15[i+1] = max(-residual_load[i*4+1:i*4+4] + residual_load[i*4:i*4+3])*4

rampup15 = np.maximum(0,rampup15)
rampdown15 = np.maximum(0,rampdown15)

# Reserve needs according the "standard" formula:
reserve_2U = np.sqrt(10*np.max(demand_tot)+150**2)-150
reserve_2D = reserve_2U/2

# Reserve needs according to Elia:
# Evaluation from their 2018 reserve study and the VRE capaciy:
# Taking the medium scenario (p.10):
# In 2013, considering 2.2 GW PV and 1 GW wind

# Tertiary reserves:
# [year   VRE_capacity    FRRm_down    FRRm_up]
# FRRm_Elia = [ [2013,   3.2,        1120,     695 ],
#              [2018,   3.2+4.8,    1138,     1238] ]
Pcap_sol2 = Pcap_sol
Pcap_wind2 = Pcap_wind
if Pcap_sol < 0:
    Pcap_sol2 = np.mean(Pmax_sol)
if Pcap_wind < 0:
    Pcap_wind2 = np.mean(Pmax_wind)
P_VRE_GW = (Pcap_sol2+Pcap_wind2)/1000
FRRm_down = 695 + (1238 - 695)/(8 - 3.2)*(P_VRE_GW - 3.2)
FRRm_up = 1120 + (1138 - 1120)/(8 - 3.2)*(P_VRE_GW - 3.2)

FRRm_down = FRRm_down/4     # going from MW to MW/h (15 time frame)
FRRm_up = FRRm_up/4     # going from MW to MW/h (15 time frame)

#  Secondary reserves:
#  [year   VRE_capacity    FRRa]
# FRRa_Elia = [ 2013   3.2        140
#              2018   3.2+4.8    172]
FRRa = 140 + (172 - 140)/(8 - 3.2)*(P_VRE_GW - 3.2)

reserve_2U_tot = FRRa + np.maximum(rampup15,FRRm_up)
reserve_2D_tot = FRRa + np.maximum(rampdown15,FRRm_down)


###################################################################################################################
############################################   Sets    ############################################################
###################################################################################################################


# In the case of a set, the uels variable must be a list containing the different set values
# In the case of a parameter, the 'sets' variable must be a list containing the different set dictionnaries

# Simulation hours:
set_h = {'name':'h','type':'set','uels':[str(x+1) for x in range(Nhours)]}
# Simulation days:
set_d = {'name':'d','type':'set','uels':[str(x/24 + 1) for x in range(Nhours)]}
# Markets:
set_mk = {'name':'mk','type':'set','uels':['DA','2U','2D']}
# Nodes:
set_n = {'name':'n','type':'set','uels':['BE']}
# Units:
set_u = {'name':'u','type':'set','uels':plants_merged['unit']}
# Bidding blocks:
set_b = {'name':'b','type':'set','uels':['B01']}
# Countries:
set_c = {'name':'c','type':'set','uels':['BE']}
# Generation companies:
set_g = {'name':'g','type':'set','uels':['dummy']}
# Pollutants:
set_p = {'name':'p','type':'set','uels':['CO2']}
# Lines:
set_l = {'name':'l','type':'set','uels':['Dummy']}
# Fuels:
set_f = {'name':'f','type':'set','uels':['Coal','NG','Fuel','Nuclear','Bio','VRE','Elec','Other']}
# Storage Units (subset of u):
sto=[]
for j in range(Nunits):
    if plants_merged['type'][j]=='HUSTO':
        sto.append(plants_merged['unit'][j])
#set_s = {'name':'s','type':'set','uels':sto}
# Technologies:
set_t = {'name':'t','type':'set','uels':['CCGT','CL','D','GT','HU','HUSTO','IS','NU','TJ','WKK','WT','PV']}
# Renewable Technologies (subset of t):
set_tr = {'name':'tr','type':'set','uels':['HU','IS','WT','PV'],'index':[4,6,10,11]}



###################################################################################################################
############################################   Parameters    ######################################################
###################################################################################################################

list_par = []   # List of the parameters to be written to the gdx
list_par2 = [] # List of the parameters (no function of the time)


# Availability Factor:
# The AvailabilityFactor parameter is used for the Variable Renewable generation (it is the fraction of the nominal
# power available at each hour). It is only relevant if data is available for the generation of individual VRE
# or CHP power plants. In this analysis, a value of 1 is imposed to all power plants.
# The data is provided on a daily basis, transforming into hours and building the gdx structure:

values = np.zeros([Nunits*Nhours,3])
values[:,0] = np.repeat(np.arange(0,Nunits,1),Nhours)   # Repeat Nunits vector Nhours times
values[:,1] = np.tile(np.arange(0,Nhours,1),Nunits)     # superpose Nhours vector Nunits times
for i in range(Nunits):
        if plants_merged['type'][i] == 'WT':
            values[i*Nhours:(i+1)*Nhours,2]= AF_wind
        elif plants_merged['type'][i] == 'PV':
            values[i*Nhours:(i+1)*Nhours,2]= AF_sol
        elif plants_merged['type'][i] == 'HU':
            values[i*Nhours:(i+1)*Nhours,2] =  AF_water
        elif plants_merged['type'][i] == 'IS':
            values[i*Nhours:(i+1)*Nhours,2] =  AF_nondispatch[::4]
        else:
            values[i*Nhours:(i+1)*Nhours,2]= np.ones(Nhours)

sets = [set_u, set_h]
list_par.append({'name':'AvailabilityFactor','type':'parameter','sets':sets,'val':values})

AvailabilityFactor_array = np.zeros((set_u['uels'].__len__(),set_h['uels'].__len__()))

j=0
k=0
for i in range(0,values.shape[0],1):
    AvailabilityFactor_array[j][k]=values[i][2]
    if (values[i][1]==Nhours-1):
        j=j+1
        k=0
    if(values[i][1]!=Nhours-1):
        k=k+1


CostFixed_array_pyomo={}
CostShutDown_array_pyomo = {}
CostStartUp_array_pyomo= {}
Efficiency_array_pyomo={}
PartLoadMin_array_pyomo={}
value_0= np.array(plants_merged['efficiency'])
value_1 = np.array(plants_merged['part-load'])

for i in range(0,Nunits,1):
    CostFixed_array_pyomo[i]=0
    CostShutDown_array_pyomo[i]=0
    CostStartUp_array_pyomo[i]= 0
    Efficiency_array_pyomo[i]= value_0[i]
    PartLoadMin_array_pyomo[i]= value_1[i]

# Fixed cost parameter for each unit:
sets = [set_u]
list_par.append({'name':'CostFixed','type':'parameter','sets':sets,'val':np.array([[ ]])})
list_par2.append({'name':'CostFixed','type':'parameter','sets':sets,'val':CostFixed_array_pyomo})


# Shutdown cost parameter for each unit:
sets = [set_u]
list_par.append({'name':'CostShutDown','type':'parameter','sets':sets,'val':np.array([[ ]])})
list_par2.append({'name':'CostShutDown','type':'parameter','sets':sets,'val':CostShutDown_array_pyomo})

# Curtailment:
# For each node
# Set to 0 to impeed curtailment
# Set to 1 to allow. Max curtailement is the sum of renewable production at each time step
# if set higher than 1, multiplies renewable production
sets = [set_n]
list_par.append({'name':'Curtailment','type':'parameter','sets':sets,'val':np.array([[0,curtailment]])})
Curtailment_array_pyomo = {}
Curtailment_array_pyomo[0]=curtailment
list_par2.append({'name':'Curtailment','type':'parameter','sets':sets,'val':Curtailment_array_pyomo})
# Demand:
sets = [set_h, set_n, set_mk]
if use_residual_load:
    value = np.array([np.r_[0:Nhours], np.zeros(Nhours), np.zeros(Nhours), np.squeeze(demand_res)])
else:
    value = np.array([np.r_[0:Nhours], np.zeros(Nhours), np.zeros(Nhours), np.squeeze(demand_wo_int)])
# Reducing the down reserve requirement at low residual load:
reserve_2D_tot_trunc = np.minimum(np.maximum(0,demand_res.transpose()),reserve_2D_tot)
value = np.concatenate((value,np.array([np.r_[0:Nhours], np.zeros(Nhours), np.ones(Nhours), np.squeeze(reserve_2U_tot)])),axis=1)
value = np.concatenate((value,np.array([np.r_[0:Nhours], np.zeros(Nhours), np.ones(Nhours)*2, np.squeeze(reserve_2D_tot_trunc)])),axis=1)
list_par.append({'name':'Demand','type':'parameter','sets':sets,'val':value.T})
Demand_array = np.zeros((set_n['uels'].__len__(),set_mk['uels'].__len__(),set_h['uels'].__len__()))
l=0
for i in range(0,set_n['uels'].__len__(),1):
    for j in range(0,set_mk['uels'].__len__(),1):
        for k in range(0,set_h['uels'].__len__(),1):
            Demand_array [i][j][k]=value[3][l]
            l=l+1
# Plant Efficiency:
sets = [set_u]
value = np.array([range(0,Nunits),plants_merged['efficiency']])
list_par.append({'name':'Efficiency','type':'parameter','sets':sets,'val':value.T})
list_par2.append({'name':'Efficiency','type':'parameter','sets':sets,'val':Efficiency_array_pyomo})
# Fuels:
sets = [set_u, set_f]
value_f = np.zeros(Nunits)
for i in range(Nunits):
    if (plants_merged['fuel'][i] == 'NG') or (plants_merged['fuel'][i] == 'NG/BF'):
        value_f[i] = 1
    else:
        value_f[i] = 0
value = np.array([range(Nunits),value_f, np.ones(Nunits)])
list_par.append({'name':'Fuel','type':'parameter','sets':sets,'val':value.T})
Fuel_array_pyomo={}
Fuel_array = np.zeros((set_u['uels'].__len__(),2))
for i in range(0,value[0].size,1):
    Fuel_array[i][0]= value[1][i]
    Fuel_array[i][1]= value[2][i]
    for j in range(0,2,1):
        Fuel_array_pyomo[i,j]= 1
        if (value[1][i]==1):
            del Fuel_array_pyomo[i,0]
            Fuel_array_pyomo[i,1]= 1
            break
        if (value[1][i]==0):
            break
list_par2.append({'name':'Fuel','type':'parameter','sets':sets,'val':Fuel_array_pyomo})
# Variable cost:
# Note that Dispaset allows defining the price in bidding blocks. This is not used here => only the first block is considered

# Biomass price
price_bio = 250                            # in EUR/T, source: http://www.michamps4b.be/evolution-prix-du-pellet-en-belgique.php
LHV = 5.28                                 # in MWh_th/T   (19 Gj/T)
price_bio = price_bio/LHV                  # in EUR/MWh_th
price_bio = max(0,price_bio - 60/0.4)      # Including green certificates into the fuel price (eta=0.4)

values = []
for i in range(Nunits):

    if plants_merged['fuel'][i] in ['NG','NG/BF','CG']:
        value = price['gas']/plants_merged['efficiency'][i]
    elif plants_merged['fuel'][i] in ['CL','CP']:
        value = price['coal']/plants_merged['efficiency'][i]
    elif plants_merged['fuel'][i] in ['LF','LV','GO']:
        value = price['oil']/plants_merged['efficiency'][i]
    elif plants_merged['fuel'][i] in ['NU']:
        value = 6 * np.ones(Nhours)
    elif plants_merged['fuel'][i] in ['WI','WA','WR','SOL']:
        value = np.zeros(Nhours)
    elif plants_merged['fuel'][i] in ['WP','BF']:
        value = price_bio * np.ones(Nhours)/plants_merged['efficiency'][i]
    else:
        print 'No fuel price value has been found for fuel ' + plants_merged['fuel'][i] + ' in unit ' +  plants_merged['fuel'][i] + '. A null variable cost has been assigned'
        value = np.zeros(Nhours)
    if i == 0:
        values= [range(Nhours)]
        values.append([i for x in range(Nhours)])
        values.append([0]*Nhours)
        values.append(value.tolist())
    else:
        values[0] = values[0] + range(Nhours)
        values[1] = values[1] + [i for x in range(Nhours)]
        values[2] = values[2] + [0]*Nhours
        values[3] = values[3] + value.tolist()
sets = [set_h, set_u, set_b]
list_par.append({'name':'CostVariable','type':'parameter','sets':sets,'val':np.array(values).T})
CostVariable_array = np.zeros((set_b['uels'].__len__(),set_u['uels'].__len__(),set_h['uels'].__len__()))
l=0
for i in range(0,set_b['uels'].__len__(),1):
    for j in range(0,set_u['uels'].__len__(),1):
        for k in range(0,set_h['uels'].__len__(),1):
            CostVariable_array [i][j][k]=np.array(values).T[l][3]
            l=l+1
# Fuel prices (not used in this case):
sets = [set_h,set_n,set_f]
list_par.append({'name':'FuelPrice','type':'parameter','sets':sets,'val':np.array([[ ]])})
Fuel_price_array = np.zeros((set_n['uels'].__len__(),set_f['uels'].__len__(),set_h['uels'].__len__()))
# Markup prices (not used in this case):
sets = [set_h,set_n,set_f]
list_par.append({'name':'Markup','type':'parameter','sets':sets,'val':np.array([[ ]])})
Markup_array = np.zeros((set_n['uels'].__len__(),set_f['uels'].__len__(),set_h['uels'].__len__()))
# Start-up cost parameter:
sets = [set_u]
list_par.append({'name':'CostStartUp','type':'parameter','sets':sets,'val':np.array([[ ]])})
list_par2.append({'name':'CostStartUp','type':'parameter','sets':sets,'val':CostStartUp_array_pyomo})
# Maximum pollutant emissions:
# For each country and each emission type:
# In this case, we consider no limit => very high number
sets = [set_c,set_p]
list_par.append({'name':'EmissionMaximum','type':'parameter','sets':sets,'val':np.array([[0, 0, 1E15]])})
e = {}
e[0,0]=1E15
list_par2.append({'name':'EmissionMaximum','type':'parameter','sets':sets,'val':e})
# Emission rate parameter (t/MWh):
# for each unit and each emission type:
# No emission considered for the time being
sets = [set_u,set_p]
list_par.append({'name':'EmissionRate','type':'parameter','sets':sets,'val':np.array([[ ]])})
EmissionRate_array_pyomo={}
for i in range(0,set_u['uels'].__len__(),):
    for j in range(0,set_p['uels'].__len__(),1):
        EmissionRate_array_pyomo[i,j]=0
list_par2.append({'name':'EmissionRate','type':'parameter','sets':sets,'val':EmissionRate_array_pyomo})
## Maximum line capacity:
# for each hour and each line:
Nlines = 1
value = []
for i in range(Nlines):
    if i == 0:
        values= [range(Nhours)]
        values.append([i for x in range(Nhours)])
        values.append([1E15]*Nhours)
    else:
        values[0] = values[0] + range(Nhours)
        values[1] = values[1] + [i for x in range(Nhours)]
        values[2] = values[2] + [1E15]*Nhours
sets = [set_h,set_l]
list_par.append({'name':'FlowMaximum','type':'parameter','sets':sets,'val':np.array(values).T})
FlowMaximum_array = np.zeros((set_l['uels'].__len__(),set_h['uels'].__len__()))
l=0
for i in range(0,set_l['uels'].__len__(),1):
    for j in range(0,set_h['uels'].__len__(),1):
            FlowMaximum_array [i][j]=np.array(values).T[l][2]
            l=l+1
# Minimum line capacity (flows are always positive in Dispaset): 0
sets = [set_h,set_l]
list_par.append({'name':'FlowMinimum','type':'parameter','sets':sets,'val':np.array([[ ]])})
FlowMinimum_array = np.zeros((set_l['uels'].__len__(),set_h['uels'].__len__()))
# Nodes in lines - Incidence matrix:
# for each line, define the two boundary nodes:
# Assign -1 if 1 MW in the line removes power from node n
# Assign +1 if 1 MW in the line adds power to node n
sets = [set_l,set_n]
list_par.append({'name':'LineNode','type':'parameter','sets':sets,'val':np.array([[ ]])})
LineNode_array_pyomo={}
for i in range(0,set_l['uels'].__len__(),1):
    for j in range(0,set_n['uels'].__len__(),1):
        LineNode_array_pyomo[i,j]=0

list_par2.append({'name':'LineNode','type':'parameter','sets':sets,'val':LineNode_array_pyomo})
# Location of each unit:
# For each unit, assign a 1 for the pair unit/node that corresponds
values=[]
values.append(range(Nunits))
values.append([0]*Nunits)
values.append([1]*Nunits)
sets = [set_u,set_n]
list_par.append({'name':'Location','type':'parameter','sets':sets,'val':np.array(values).T})
Location_array_pyomo ={}
l=0
for i in range(0,set_u['uels'].__len__(),1):
    for j in range(0,set_n['uels'].__len__(),1):
        Location_array_pyomo[i,j]=np.array(values).T[l][2]
        l=l+1
list_par2.append({'name':'Location','type':'parameter','sets':sets,'val':Location_array_pyomo})
# Load Shedding capacity:
# For each node
# 331 MW available in the winter 2012-2013, in the form of interruptible contracts
sets = [set_n]
list_par.append({'name':'LoadShedding','type':'parameter','sets':sets,'val':np.array([[0, 331]])})
LoadShedding_array_pyomo={}
LoadShedding_array_pyomo[0]=331
list_par2.append({'name':'LoadShedding','type':'parameter','sets':sets,'val':LoadShedding_array_pyomo})
# OutageFactor parameter:
# In this analysis, we take as availability factor the outage data provided by Elia:
# The data is provided on a daily basis, transforming into hours and

values = np.zeros([Nunits*Nhours,3])
temp = 1-np.array(outages_merged)
values[:,0] = np.repeat(np.arange(0,Nunits,1),Nhours)   # Repeat Nunits vector Nhours times
values[:,1] = np.tile(np.arange(0,Nhours,1),Nunits)     # superpose Nhours vector Nunits times
for i in range(Nunits):
    values[i*Nhours:(i+1)*Nhours,2]= np.repeat(temp[i,:],24)

sets = [set_u,set_h]
list_par.append({'name':'OutageFactor','type':'parameter','sets':sets,'val':values})
OutageFactor_array = np.zeros((set_u['uels'].__len__(),set_h['uels'].__len__()))
j=0
k=0
for i in range(0,values.shape[0],1):
    if (values [i][2]<0.1):
        values[i][2]=0
    OutageFactor_array[j][k]=values[i][2]
    if (values[i][1]==Nhours-1):
        j=j+1
        k=0
    if(values[i][1]!=Nhours-1):
        k=k+1

# Ownership
# No ownership-specific data is considered in this work, adding a dummy owner for each plant:
sets = [set_u,set_g]
list_par.append({'name':'Ownership','type':'parameter','sets':sets,'val':np.array([[ ]])})
Ownership_array_pyomo ={}
for i in range(0,set_u['uels'].__len__(),1):
    for j in range(0,set_g['uels'].__len__(),1):
        Ownership_array_pyomo[i,j]=0
list_par2.append({'name':'Ownership','type':'parameter','sets':sets,'val':Ownership_array_pyomo})
# Price of pollutants (e.g. ETS system)
# Defining a cost of 5 EUR/Tco2
sets = [set_p]
list_par.append({'name':'PermitPrice','type':'parameter','sets':sets,'val':np.array([[0, 5]])})
PermitPrice_array_pyomo={}
PermitPrice_array_pyomo[0]=5
list_par2.append({'name':'PermitPrice','type':'parameter','sets':sets,'val':PermitPrice_array_pyomo})
# Transmission price
# for each line and each hour.
# Assuming 0 => empty array
sets = [set_h,set_l]
list_par.append({'name':'PriceTransmission','type':'parameter','sets':sets,'val':np.array([[ ]])})
PriceTransmission_array = np.zeros((set_l['uels'].__len__(),set_h['uels'].__len__()))
for i in range(Nunits):
    if plants_merged['type'][i] in ['HUSTO']:
        if Pcap_sto <0:
            plants_merged['pmax'][i] = plants_merged['pmax'][i]
        else:
            plants_merged['pmax'][i] = Pcap_sto
            Pcap_sto = 0
    if plants_merged['fuel'][i] in ['SOL','WI']:
        plants_merged['pmax'][i] = 0                # Sets all renewable capacities to 0
if not use_residual_load:                           # if we dont use the residual load approach, assign a power to the VRE plants
    for i in range(Nunits):                             # Only valid if there is only one of each in the unit list
        if plants_merged['fuel'][i] == 'SOL':
            plants_merged['pmax'][i] = Pcap_sol2
        if plants_merged['fuel'][i] == 'WI':
            plants_merged['pmax'][i] = Pcap_wind2
else:
    for i in range(Nunits):
        if plants_merged['type'][i] in ['IS','HU']:
            plants_merged['pmax'][i] = 0;
sets = [set_u]
value = np.array([range(0,Nunits),plants_merged['pmax']])
list_par.append({'name':'PowerCapacity','type':'parameter','sets':sets,'val':value.T})
PowerCapacity_array_pyomo ={}
for i in range(0,Nunits,1):
    PowerCapacity_array_pyomo[i]= plants_merged['pmax'][i]

list_par2.append({'name':'PowerCapacity','type':'parameter','sets':sets,'val':PowerCapacity_array_pyomo})
# Block capacity:
sets = [set_u,set_b]
value = np.array([range(0,Nunits),[0]*Nunits,plants_merged['pmax']])
list_par.append({'name':'BlockCapacity','type':'parameter','sets':sets,'val':value.T})
BlockCapacity_array_pyomo = {}
l=0
for i in range(0,set_u['uels'].__len__(),1):
    for j in range(0,set_b['uels'].__len__(),1):
        BlockCapacity_array_pyomo[i,j]=value.T[l][2]
        l=l+1

list_par2.append({'name':'BlockCapacity','type':'parameter','sets':sets,'val':BlockCapacity_array_pyomo})
# Minimum / must run capacity
sets = [set_u]
value = np.array([range(0,Nunits),plants_merged['part-load']])
list_par.append({'name':'PartLoadMin','type':'parameter','sets':sets,'val':value.T})
list_par2.append({'name':'PartLoadMin','type':'parameter','sets':sets,'val':PartLoadMin_array_pyomo})
# Ramping up/down
# for each unit, in MW/h
# PowerCapacity:
sets = [set_u]
value = [plants_merged['pmax'][x]*plants_merged['rampup'][x]*60 for x in range(Nunits)]
values = np.array([range(0,Nunits),value])
list_par.append({'name':'RampUpMaximum','type':'parameter','sets':sets,'val':values.T})
RampUpMaximum_array_pyomo ={}
for i in range(0,Nunits,1):
    RampUpMaximum_array_pyomo[i]= value[i]

list_par2.append({'name':'RampUpMaximum','type':'parameter','sets':sets,'val':RampUpMaximum_array_pyomo})
value = [plants_merged['pmax'][x]*plants_merged['rampdown'][x]*60 for x in range(Nunits)]
values = np.array([range(0,Nunits),value])
list_par.append({'name':'RampDownMaximum','type':'parameter','sets':sets,'val':values.T})
RampDownMaximum_array_pyomo ={}
for i in range(0,Nunits,1):
    RampDownMaximum_array_pyomo[i]= value[i]
list_par2.append({'name':'RampDownMaximum','type':'parameter','sets':sets,'val':RampDownMaximum_array_pyomo})
# Start up / Shut down
# for each unit, in MW/h
# In the data, the start_up is provided in hours to full load.
# start up ramping rate is given by the capacity divided by start-up time
# if start-up time is zero, it is assumed that full load can be reached in 15 minutes
sets = [set_u]
value = [plants_merged['pmax'][x]/max(plants_merged['startup'][x],0.25) for x in range(Nunits)]
values = np.array([range(0,Nunits),value])
list_par.append({'name':'RampStartUpMaximum','type':'parameter','sets':sets,'val':values.T})

# No data for shutdown times, taking the same time as start-up.
list_par.append({'name':'RampShutDownMaximum','type':'parameter','sets':sets,'val':values.T})
RampStartUpMaximum_array_pyomo ={}
for i in range(0,Nunits,1):
    RampStartUpMaximum_array_pyomo[i]= value[i]
list_par2.append({'name':'RampStartUpMaximum','type':'parameter','sets':sets,'val':RampStartUpMaximum_array_pyomo})
list_par2.append({'name':'RampShutDownMaximum','type':'parameter','sets':sets,'val':RampStartUpMaximum_array_pyomo})
# Participation to the reserve market
# Binary variable depending on the technology
# {'CCGT','CL','D','GT','HU','HUSTO','IS','NU','TJ','WKK','WT','PV'}
sets = [set_t]
Ntech = len(set_t['uels'])
value = np.array([range(Ntech),[ 1,1,1,1,0,1,0,0,1,0,0,0 ]])
list_par.append({'name':'Reserve','type':'parameter','sets':sets,'val':value.T})
value = np.array([ 1,1,1,1,0,1,0,0,1,0,0,0 ])
Reserve_array_pyomo ={}
for i in range(0,Ntech,1):
    Reserve_array_pyomo[i]= value[i]
list_par2.append({'name':'Reserve','type':'parameter','sets':sets,'val':Reserve_array_pyomo})
# Renewable targets:
# For each node, in MW
sets = [set_n]
list_par.append({'name':'RenewableMinimum','type':'parameter','sets':sets,'val':np.array([[ ]])})
RenewableMinimum_array_pyomo = {}
RenewableMinimum_array_pyomo[0]=0
list_par2.append({'name':'RenewableMinimum','type':'parameter','sets':sets,'val':RenewableMinimum_array_pyomo})
# Storage variables:
# get the indexes of all storage units:
idx_sto = [i for i,x in enumerate(plants_merged['type']) if x == 'HUSTO']
Nsto = len(idx_sto)

set_s = {'name':'s','type':'set','uels':sto,'index':idx_sto }

# Discharge Efficiency parameter:
# for each pumped storage unit, setting the efficiency to sqrt(0.75):
sets = [set_u]
value = np.array([idx_sto,[np.sqrt(0.75)]*Nsto])
list_par.append({'name':'StorageDischargeEfficiency','type':'parameter','sets':sets,'val':value.T})
StorageDischargeEfficiency_array_pyomo = {}
for i in idx_sto:
    StorageDischargeEfficiency_array_pyomo[i]=np.sqrt(0.75)*Nsto

list_par2.append({'name':'StorageDischargeEfficiency','type':'parameter','sets':sets,'val':StorageDischargeEfficiency_array_pyomo})
# Storage capacity:
# Assuming 4.5 hours at full load for all storage units:
value = np.array([idx_sto,[plants_merged['pmax'][x]*hours_sto for x in idx_sto]])
list_par.append({'name':'StorageCapacity','type':'parameter','sets':sets,'val':value.T})
StorageCapacity_array_pyomo = {}
for i in idx_sto:
    StorageCapacity_array_pyomo[i]=plants_merged['pmax'][i]*hours_sto

list_par2.append({'name':'StorageCapacity','type':'parameter','sets':sets,'val':StorageCapacity_array_pyomo})
# Storage inflow:
# for each hour and each storage unit
sets = [set_h,set_u]
list_par.append({'name':'StorageInflow','type':'parameter','sets':sets,'val':np.array([[ ]])})

# Storage outflow:
# for each hour and each storage unit
sets = [set_h,set_u]
list_par.append({'name':'StorageOutflow','type':'parameter','sets':sets,'val':np.array([[ ]])})
StorageInflow_array = np.zeros((set_u['uels'].__len__(),set_h['uels'].__len__()))
StorageOutflow_array = np.zeros((set_u['uels'].__len__(),set_h['uels'].__len__()))
# Initial state (assuming half-filled):
sets = [set_u]
value = np.array([idx_sto,[plants_merged['pmax'][x]*hours_sto/2 for x in idx_sto]])
list_par.append({'name':'StorageInitial','type':'parameter','sets':sets,'val':value.T})
StorageInitial_array_pyomo = {}
for i in idx_sto:
    StorageInitial_array_pyomo [i]=plants_merged['pmax'][i]*hours_sto/2
list_par2.append({'name':'StorageInitial','type':'parameter','sets':sets,'val':StorageInitial_array_pyomo})
# Minimum level:
sets = [set_u]
list_par.append({'name':'StorageMinimum','type':'parameter','sets':sets,'val':np.array([[ ]])})
StorageMinimum_array_pyomo = {}
for i in idx_sto:
    StorageMinimum_array_pyomo [i]=0
list_par2.append({'name':'StorageMinimum','type':'parameter','sets':sets,'val':StorageMinimum_array_pyomo})
# Pumping Efficiency parameter:
# for each pumped storage unit, setting the efficiency to sqrt(0.75):
sets = [set_u]
value = np.array([idx_sto,[np.sqrt(0.75)]*Nsto])
list_par.append({'name':'StorageChargingEfficiency','type':'parameter','sets':sets,'val':value.T})
StorageChargingEfficiency_array_pyomo = {}
for i in idx_sto:
    StorageChargingEfficiency_array_pyomo [i]=np.sqrt(0.75)*Nsto
list_par2.append({'name':'StorageChargingEfficiency','type':'parameter','sets':sets,'val':StorageChargingEfficiency_array_pyomo})
# Pumping Capacity parameter:
# Assuming the same power in pmping and in turbining mode:
sets = [set_u]
value = np.array([idx_sto,[plants_merged['pmax'][x] for x in idx_sto]])
list_par.append({'name':'StorageChargingCapacity','type':'parameter','sets':sets,'val':value.T})
StorageChargingCapacity_array_pyomo = {}
for i in idx_sto:
    StorageChargingCapacity_array_pyomo [i]=plants_merged['pmax'][i]

list_par2.append({'name':'StorageChargingCapacity','type':'parameter','sets':sets,'val':StorageChargingCapacity_array_pyomo})
# Technologies
# For each unit, parameter=1 if the input technology corresponds
sets = [set_u,set_t]
value = []
for i in range(Nunits):
    value.append(set_t['uels'].index(plants_merged['type'][i]))             # Find the index of the plant type in set_t:
list_par.append({'name':'Technology','type':'parameter','sets':sets,'val':np.array([range(Nunits),value,[1]*Nunits]).T})
Technology_array_pyomo={}
Technology_array_pyomo2={}
for i in range(0,Nunits,1):
    for j in range(0,Ntech,1):
        if (j==value[i]):
            Technology_array_pyomo[i,j]=1
            Technology_array_pyomo2[i]=set_t['uels'][j]
        else: Technology_array_pyomo[i,j]=0

list_par2.append({'name':'Technology','type':'parameter','sets':sets,'val':Technology_array_pyomo})

# Minimum Time up/down
# In hours for each unit (values are set arbitrarily at the moment)
# assuming same value for up and down
sets = [set_u]
value = []
for i in range(Nunits):
    if plants_merged['type'][i] in ['CCGT','CL']:
        value.append(5)
    elif plants_merged['type'][i] in ['NU']:
        value.append(24)
    else:
        value.append(0)
list_par.append({'name':'TimeDownMinimum','type':'parameter','sets':sets,'val':np.array([range(Nunits),value]).T})
list_par.append({'name':'TimeUpMinimum','type':'parameter','sets':sets,'val':np.array([range(Nunits),value]).T})
TimeDownMinimum_array_pyomo ={}
for i in range(0,Nunits,1):
    TimeDownMinimum_array_pyomo[i]=value[i]
list_par2.append({'name':'TimeDownMinimum','type':'parameter','sets':sets,'val':TimeDownMinimum_array_pyomo})
list_par2.append({'name':'TimeUpMinimum','type':'parameter','sets':sets,'val':TimeDownMinimum_array_pyomo})
# Initial state:
# Number of hours the plant has been up or down at time 0
sets = [set_u]
value = np.array([range(Nunits),[48]*Nunits])
list_par.append({'name':'TimeDownInitial','type':'parameter','sets':sets,'val':value.T})
list_par.append({'name':'TimeUpInitial','type':'parameter','sets':sets,'val':value.T})
a=[48]*Nunits
TimeDownInitial_array_pyomo ={}
for i in range(0,Nunits,1):
    TimeDownInitial_array_pyomo[i]=a[0]
list_par2.append({'name':'TimeDownInitial','type':'parameter','sets':sets,'val':TimeDownInitial_array_pyomo})
list_par2.append({'name':'TimeUpInitial','type':'parameter','sets':sets,'val':TimeDownInitial_array_pyomo})


# Initial Power
# Three possibilities to define the initial power:
# - Use a very simple merit order and dispatch to define an initial state
#   that complies with the constraints
# - Use historical data
# - Set all values to zero => in this case, the UC model decides by itself

# In this case, we assume that all CCGT and NU plants are up initially, half way between pmin and pmax:
value = []
for i in range(Nunits):
    if plants_merged['type'][i] in ['CCGT','NU']:
        value.append((plants_merged['part-load'][i]+1)/2 * plants_merged['pmax'][i] * outages_merged[i][0])
    else:
        value.append(0)
list_par.append({'name':'PowerInitial','type':'parameter','sets':sets,'val':np.array([range(Nunits),value]).T})
PowerInitial_array_pyomo ={}
CommittedInitial_array_pyomo ={}
for i in range(0,Nunits,1):
    PowerInitial_array_pyomo[i]=value[i]
    CommittedInitial_array_pyomo[i]=0
    if(value[i]>0):
        CommittedInitial_array_pyomo[i]=1
list_par2.append({'name':'PowerInitial','type':'parameter','sets':sets,'val':PowerInitial_array_pyomo})

list_par2.append({'name':'CommittedInitial','type':'parameter','sets':sets,'val':CommittedInitial_array_pyomo})



########################################################################################################################
################################################# List of set and parameters  ##########################################
########################################################################################################################


list_sets = {}



list_sets[2]=set_mk
list_sets[3]=set_b
list_sets[4]=set_n
list_sets[5]=set_c
list_sets[6]=set_g
list_sets[7]=set_p
list_sets[8]=set_l
list_sets[9]=set_f
list_sets[10]=set_s
list_sets[11]=set_t
list_sets[12]=set_tr
list_sets[13]=set_u

TimeUp = {}
TimeDown = {}

   
   
for d in range(0,Ndays-1,1):
    list_par_time = []
    AF={} # fct (u,h)
    D= {} # fct (n,mk,h)
    CV={} # fct (b,u,h)
    FP={} # fct (n,f,h)
    M= {} # fct (n,f,h)
    FM={} # fct (l,h)
    Fm={} # fct (l,h)
    OF={} # fct (u,h)
    PT={} # fct (l,h)
    SI={} # fct (u,h)
    SO={} # fct (u,h)

    TimeDownInitial = {}
    TimeUpInitial = {}
    PowerInitial={}
    CommittedInitial = {}

    for u in range(0,Nunits,1):
        for h in range(0+d*24,48+d*24,1):
            AF[u,h]= AvailabilityFactor_array[u][h]
            OF[u,h]= OutageFactor_array [u][h]
            SI[u,h]= StorageInflow_array[u][h]
            SO[u,h]= StorageOutflow_array[u][h]
    for l in range(0,set_l['uels'].__len__(),1):
        for h in range(0+d*24,48+d*24,1):
            FM[l,h]=FlowMaximum_array[l][h]
            Fm[l,h]= FlowMinimum_array[l][h]
            PT[l,h]= PriceTransmission_array[l][h]
    for n in range(0,set_n['uels'].__len__(),1):
        for mk in range(0,set_mk['uels'].__len__(),1):
            for h in range(0+d*24,48+d*24,1):
                D [n,mk,h] = Demand_array[n][mk][h]
    for b in range(0,set_b['uels'].__len__(),1):
        for u in range(0,Nunits,1):
            for h in range(0+d*24,48+d*24,1):
                CV[b,u,h] = CostVariable_array[b][u][h]
    for n in range(0,set_n['uels'].__len__(),1):
        for f in range(0,set_f['uels'].__len__(),1):
            for h in range(0+d*24,48+d*24,1):
                FP[n,f,h]=Fuel_price_array[n][f][h]
                M[n,f,h]= Markup_array[n][f][h]

    set_h_uels=range(0+d*24,48+d*24,1)
    set_h_1 = {'name':'h','type':'set','uels':[set_h_uels]}
    set_d_1 = {'name':'d','type':'set','uels':[x/24 + 1 for x in set_h_uels]}
    list_sets[0]=set_h_1
    list_sets[1]=set_d_1
    sets = [set_u, set_h]
    list_par_time.append({'name':'AvailabilityFactor','type':'parameter','sets':sets,'val':AF})
    list_par_time.append({'name':'OutageFactor','type':'parameter','sets':sets,'val':OF})
    list_par_time.append({'name':'StorageInflow','type':'parameter','sets':sets,'val':SI})
    list_par_time.append({'name':'StorageOutflow','type':'parameter','sets':sets,'val':SO})
    sets = [set_l, set_h]
    list_par_time.append({'name':'FlowMaximum','type':'parameter','sets':sets,'val':FM})
    list_par_time.append({'name':'FlowMinimum','type':'parameter','sets':sets,'val':Fm})
    list_par_time.append({'name':'PriceTransmission','type':'parameter','sets':sets,'val':PT})
    sets = [set_n,set_mk,set_h]
    list_par_time.append({'name':'Demand','type':'parameter','sets':sets,'val':D})
    sets = [set_b,set_u,set_h]
    list_par_time.append({'name':'CostVariable','type':'parameter','sets':sets,'val':CV})
    sets = [set_n,set_f,set_h]
    list_par_time.append({'name':'FuelPrice','type':'parameter','sets':sets,'val':FP})
    list_par_time.append({'name':'Markup','type':'parameter','sets':sets,'val':M})

