# -*- coding: utf-8 -*-
"""
This example shows how to run a batch of simulations by varying some of the input parameters or data
A base case is defined by the excel configuration file, then Dispa-SET functions are used to modify the input data
The input space is defined as a latin hypercube in which all the "corners" are simulated 
For each simulation, a separate simulation environment folder is created and the simualtion is run in GAMS 
Finally, all the folders with the result files are read and the results are stored in a dataframe exported to excel

TODO:
- harmonize reserves and reserve output!
- get back the solver status
- Loadshedding not activated although there is lost load
- bad hatching in dispatch plot and windows
- no horizontal line in bar plots in windows

@author: Sylvain Quoilin
"""
#%%
# Change directory to the root folder of Dispa-SET:
import os
import numpy as np
import pandas as pd
os.chdir('..')

# Import Dispa-SET
import dispaset as ds

# Define the boundaries of the inputs space:
overcapacity = [0.8,1.3]   # thermal capacity divided by peak load
share_flex = [0.1,0.9]     # share of the thermal capacity that is flexible
hours_sto = [0.1,3]        # hours of peak load that can be shifted using the storage capacity
share_sto  = [0,0.4]       # Storage Power divided by the peak load
share_wind = [0,0.3]       # Yearly wind generation divided by yearly power consumption
share_pv = [0,0.3]         # Yearly PV generation divided by yearly power consumption

# Define the folder in which all simulations should be stored:
sim_folder = 'Simulations/batch/'

# Load the configuration file
config = ds.load_config_excel('ConfigFiles/ConfigBE.xlsx')

config['SimulationType'] = 'Integer clustering'
                         # 'Integer clustering'
                         # 'No clustering'

# Build the  reference simulation environment:
SimData = ds.build_simulation(config)

# get a few important values:
load_max = SimData['parameters']['Demand']['val'].max()                     # peak load
CF_pv = SimData['parameters']['AvailabilityFactor']['val'][-2,:].mean()     # capacity factor of PV
CF_wton = SimData['parameters']['AvailabilityFactor']['val'][-1,:].mean()   # capacity factor of onshore wind

# Generate a 6-dimensional latin hypercube varying between 0 and 1:
from pyDOE import lhs
lh = lhs(6, samples=200)

# Add all the corners of the hypercube to the input space:
import itertools
lst = list(itertools.product([0, 1], repeat=6))
lh = np.concatenate((lh,np.array(lst)))

#%%
# Loop through the input space (time consuming !!):
Nruns = lh.shape[0]
for i in range(Nruns):
    print('Run ' + str(i) + ' of ' + str(Nruns))
    cap = overcapacity[0] + lh[i,0] * (overcapacity[-1] - overcapacity[0])
    flex = share_flex[0] + lh[i,1] * (share_flex[-1] - share_flex[0])
    hours = hours_sto[0] + lh[i,2] * (hours_sto[-1] - hours_sto[0])
    sto = share_sto[0] + lh[i,3] * (share_sto[-1] - share_sto[0])
    wind = share_wind[0] + lh[i,4] * (share_wind[-1] - share_wind[0])
    pv = share_pv[0] + lh[i,5] * (share_pv[-1] - share_pv[0])

    folder = sim_folder + "%.2f" % cap  + ' - ' + "%.2f" % flex + ' - ' + "%.2f" % hours + ' - ' + "%.2f" % sto + ' - ' + "%.2f" % wind + ' - ' + "%.2f" % pv

    # in the first iteration, we load the input data from the original simulation directory:
    SimData = ds.adjust_capacity(config['SimulationDirectory'],('HPHS','WAT'),value=load_max*sto)
    # then we use the dispa-set fuction to adjust the installed capacities:
    SimData = ds.adjust_capacity(SimData,('COMC','GAS'),value=load_max*cap*flex)
    SimData = ds.adjust_capacity(SimData,('STUR','NUC'),value=load_max*cap*(1-flex))
    SimData = ds.adjust_storage(SimData,('HPHS','WAT'),value=hours*load_max)
    # For wind and PV, the units should be lumped into a single unit:
    SimData = ds.adjust_capacity(SimData,('WTON','WIN'),value=load_max*cap*wind/CF_wton,singleunit=True)
    # In this last iteration, the new gdx file is written to the simulation folder:
    SimData = ds.adjust_capacity(SimData,('PHOT','SUN'),value=load_max*cap*pv/CF_pv,singleunit=True,write_gdx=True,dest_path=folder)
    
    # Finally the modified simulation environment is simulated:
    r = ds.solve_GAMS(folder, config['GAMS_folder'])

#%%
# Read all the simulation folders one by one and store key results in dataframe:
paths = os.listdir(sim_folder)
# Only take into account the ones for which a valid dispa-set result file is present
paths_ok = [x for x in paths if os.path.isfile(sim_folder + x + '/Results.gdx')]

N = len(paths_ok)
data = pd.DataFrame(index = range(N))

for i,path in enumerate(paths_ok):
    
    inputs,results = ds.get_sim_results(path=sim_folder + path,cache=True)
    FuelPower = ds.aggregate_by_fuel(results['OutputPower'], inputs, SpecifyFuels=None)
    
    # get installed capacities by fuel:
    cap = ds.plot_country_capacities(inputs,plot=False)
    
    # Capacity factors:
    CF = {}
    for f in ['GAS','NUC','WAT','WIN','SUN']:
        CF[f] = FuelPower[f].sum()/ (cap['PowerCapacity'].loc['BE',f]*8760)
    
    r = ds.get_result_analysis(inputs,results)
    
    # Compute total lost load
    LostLoad=0
    for key in ['Out_LostLoad_MaxPower','Out_LostLoad_MinPower','Out_LostLoad_Reserve2D','Out_LostLoad_Reserve2U']:
        if key in results:
            LostLoad = LostLoad = results[key].sum().sum()
    
    # Get input data from folder name (to be improved!):
    cap,flex,hours,sto,wind,pv = [float(x) for x in path.split(' - ')]
    data.loc[i,'Thermal Power [MW/MWp]'] = cap
    data.loc[i,'Flexibility [-]'] = flex
    data.loc[i,'Storage Capacity [h]'] = hours
    if sto > 0:
        data.loc[i,'Storage Hours [h]'] = hours/sto
    else:
        data.loc[i,'Storage Hours [h]'] = 0
    data.loc[i,'Storage Power [MW/MWp]'] = sto
    data.loc[i,'Wind Penetration [-]'] = wind
    data.loc[i,'PV Penetration [-]'] = pv
    data.loc[i,'cost'] = r['Cost_kwh']
    data.loc[i,'curtailment'] = r['CountryData'].loc['BE','Curtailment']
    data.loc[i,'shedding'] = r['CountryData'].loc['BE','LoadShedding']
    data.loc[i,'CF_gas'] = CF['GAS']
    data.loc[i,'CF_nuc'] = CF['NUC']
    data.loc[i,'CF_wat'] = CF['WAT']
    data.loc[i,'CF_win'] = CF['WIN']
    data.loc[i,'CF_sun'] = CF['SUN']
    data.loc[i,'LostLoad'] = LostLoad

data.fillna(0,inplace=True)
data.to_csv(sim_folder + 'dispaset_results.csv',index=False)