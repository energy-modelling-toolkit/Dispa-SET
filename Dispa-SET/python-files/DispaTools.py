# Functions used in the DispaSET pre-processing tool
__author__ = 'Sylvain Quoilin (sylvain.quoilin@ec.europa.eu)'

import os
import sys
import numpy as np
import pandas as pd 

### PYOMO TOOLS

def get_sets(instance, var):
    '''Get sets that belong to a pyomo Variable or Param

    :param instance: Pyomo Instance
    :param var: Pyomo Variable
    :return: A list with the sets that belong to this Param 
    '''
    sets = []
    var = getattr(instance, var)

    if var.dim() > 1:
        for pset in var._index.set_tuple:
            sets.append(pset.name)
    else:
        sets.append(var._index.name)
    return sets

def get_set_members(instance, sets):
    '''Get set members that belong to this set

    :param instance: Pyomo Instance
    :param set: Pyomo Set
    :return: A list with the set members  
    '''
    sm = []
    for s in sets:
        sm.append([v for v in getattr(instance, s).value])
    return sm


def pyomo_to_pandas(instance, var):
    '''
    Function converting a pyomo variable or parameter into a pandas dataframe.
    The variable must have one or two dimensions and the sets must be provided as a list of lists
    '''
    setnames = get_sets(instance, var)
    sets = get_set_members(instance, setnames)
    var = getattr(instance, var)  # Previous script used model.var instead of var
    ####
    if len(sets) != var.dim():
        sys.exit('The number of provided set lists (' + str(len(sets)) + ') does not match the dimensions of the variable (' + str(var.dim()) + ')')
    if var.dim() == 1:
        [SecondSet] = sets
        out = pd.DataFrame(columns=[var.name], index=SecondSet)
        data = var.get_values()
        for idx in data:
            out[var.name][idx] = data[idx]
        return out
        
    elif var.dim() == 2:
        [FirstSet,SecondSet] = sets
        out = pd.DataFrame(columns=FirstSet, index=SecondSet)
        data = var.get_values()
        for idx in data:
            out[idx[0]][idx[1]] = data[idx]
        return out
    else:
        print 'the pyomo_to_pandas function currently only accepts one or two-dimensional variables'
        return []


def pyomo_format(sets, param):
    '''
    Function that flattens the multidimensional dispaset input data into the pyomo format: a dicitonary with a tuple and the parameter value. 
    The tuple contains the strings of the corresponding set values
    '''    

    param_new = {}
    ndims = len(param['sets'])
    for i in range(ndims):
        if param['val'].shape[i] != len(sets[param['sets'][i]]):
            sys.exit( 'The size of the matrix and the number of set values do no match for set ' + param['sets'][i])
    if ndims == 1:
        for k in range(len(param['val'])):
            param_new[(sets[param['sets'][0]][k])] = param['val'][k]
    elif ndims == 2:
        X,Y = np.meshgrid(np.arange(param['val'].shape[1]),np.arange(param['val'].shape[0]))
        X_flat = X.flatten()
        Y_flat = Y.flatten()
        array_flat = param['val'].flatten()
        for k in range(len(array_flat)):
            i = Y_flat[k]
            j = X_flat[k]
            param_new[(sets[param['sets'][0]][i],sets[param['sets'][1]][j])] = array_flat[k]
    elif ndims == 3:
        X,Y,Z = np.meshgrid(np.arange(param['val'].shape[1]),np.arange(param['val'].shape[0]),np.arange(param['val'].shape[2]))
        X_flat = X.flatten()
        Y_flat = Y.flatten()
        Z_flat = Z.flatten()
        array_flat = param['val'].flatten()
        for k in range(len(array_flat)):
            i = Y_flat[k]
            j = X_flat[k]
            z = Z_flat[k]
            param_new[(sets[param['sets'][0]][i],sets[param['sets'][1]][j],sets[param['sets'][2]][z])] = array_flat[k]
    else:
        sys.exit('Maximum 3 dimensions')
    return param_new        


### END OF PYOMO TOOLS

### HELPER FUNCTIONS
def tuple_format(array):
    '''
    Function that flattens a n-dimensional matrix and returns a dictionary with the values and their coordinates in a tuple     
    '''
    ndims = len(array.shape)
    out = {}
    if ndims == 1:
        for k in range(len(array)):
            out[(k)] = array[k]
    elif ndims == 2:
        X,Y = np.meshgrid(np.arange(array.shape[1]),np.arange(array.shape[0]))
        X_flat = X.flatten()
        Y_flat = Y.flatten()
        array_flat = array.flatten()
        for k in range(len(array_flat)):
            out[(Y_flat[k],X_flat[k])] = array_flat[k]
    elif ndims == 3:
        X,Y,Z = np.meshgrid(np.arange(array.shape[1]),np.arange(array.shape[0]),np.arange(array.shape[2]))
        X_flat = X.flatten()
        Y_flat = Y.flatten()
        Z_flat = Z.flatten()
        array_flat = array.flatten()
        for k in range(len(array_flat)):
            out[(Y_flat[k],X_flat[k],Z_flat[k])] = array_flat[k]
    else:
        sys.exit('Maximum 3 dimensions')
    return out


def mylogspace(low,high,N):
    '''
    Self-defined logspace function in which low and high are the first and last values of the space
    '''
    # shifting all values so that low = 1
    space = np.logspace(0,np.log10(high+low+1),N)-(low+1)
    return(space)

def find_nearest(array,value):
    '''
    Self-defined function to find the index of the nearest value in a vector
    '''
    idx = (np.abs(array-value)).argmin()
    return idx

def append_to_dict(k,source,destination):
    '''
    Function to add a record to a dictionnary (e.g. a new power plant) from another dictionary
    '''
    for key in source:
        if key in destination:
            destination[key]=destination[key].append(source[key][k])
        else:
            destination[key]=[source[key][k]]

def load_var(input,string):
    '''
    Load a particular variable from the DispaSET input data structure v2.0 (Obsolete from v2.1.1 onwards)
    '''
    out = []
    for var in input:
        if var['name'] == string:
            out = var['val']
    return out

def load_set(input,string,set_name):
    '''
    Load a particular set from the DispaSET input data structure (Obsolete for v2.1.1 onwards)
    '''
    out = []
    for i in input:
        if i['name'] == string:
            for j in i['sets']:
                if j['name'] == set_name:
                    out = j['uels']
    return out


def load_csv_to_pd(path_csv,file_csv,path_pandas,file_pandas):
    '''
    Function that loads a csv sheet into a panda variable and saves it in a separate path. If the saved variable is newer
    than the sheet, do no load the sheet again.
    '''
    
    filepath_csv = os.path.join(path_csv,file_csv)
    filepath_pandas = os.path.join(path_pandas,file_pandas)
    if not os.path.isdir(path_pandas):
        os.mkdir(path_pandas)
    if not os.path.isfile(filepath_pandas):
       time_pd = 0
    else:
        time_pd = os.path.getmtime(filepath_pandas)
    if os.path.getmtime(filepath_csv) > time_pd:
        data = pd.read_csv(filepath_csv)
        data.to_pickle(filepath_pandas)
    else:
        data= pd.read_pickle(filepath_pandas)
    return data

def load_xl_to_pd(path_excel,file_excel,sheet_excel,path_pandas,file_pandas,header=0):
    '''
    Function that loads an xls sheet into a panda variable and saves it in a separate path. If the saved variable is newer
    than the sheet, do no load the sheet again.
    '''
    
    filepath_excel = os.path.join(path_excel,file_excel)
    filepath_pandas = os.path.join(path_pandas,file_pandas)
    if not os.path.isdir(path_pandas):
        os.mkdir(path_pandas)
    if not os.path.isfile(filepath_pandas):
       time_pd = 0
    else:
        time_pd = os.path.getmtime(filepath_pandas)
    if os.path.getmtime(filepath_excel) > time_pd:
        data = pd.read_excel(filepath_excel,sheet_excel,header=header)
        data.to_pickle(filepath_pandas)
    else:
        data= pd.read_pickle(filepath_pandas)
    return data
    
    
def clustering(plants,AdditionalArrays=[],Nslices=20):
    '''
    Merge excessively disaggregated power Units.
    
    :param plants: Pandas dataframe with each power plant and their characteristics (following the DispaSET format)
    :param AdditionalArrays: List of arrays to be merged. The number of rows must be equal to the number of plants. The merged values are a weighted average of the original values with respect to capacities
    :param Nslices: number slices used to fingerprint each power plant characteristics. slices in the power plant data to categorize them  (fewer slices involves that the plants will be aggregated more easily)
    :return: A list with the merged plants and a sublist of the merged additional arrays
    '''
    
    # Checking the the required columns are present in the input pandas dataframe:
    required_inputs = ['Unit','PowerCapacity','PartLoadMin','RampUpRate','RampDownRate','StartUpTime','MinUpTime','MinDownTime','NoLoadCost','StartUpCost','Efficiency']   
    for input in required_inputs:
        if input not in plants.columns:
            sys.exit("The plants dataframe requires a '" + input + "' column for clustering")
    
    # Number of units:
    Nunits = len(plants)
    plants.index = range(Nunits)
    
    # Checking the size of the additional arrays:
    for array in AdditionalArrays:
        if isinstance(array,np.ndarray):
            if array.shape[0] != Nunits:
                sys.exit("The number of rows in all the additional array must be equal to the number of units (" + str(Nunits))
        else:
            sys.exit('All AdditionalArrays objects must be numpy ndarrays')
            
    # Definition of the mapping variable, from the old power plant list the new (merged) one:
    map_old_new = np.zeros(Nunits)
    map_plant_orig = []

    # Slicing:
    bounds={'PartLoadMin':np.linspace(0,1,Nslices),'RampUpRate':np.linspace(0,1,Nslices),'RampDownRate':np.linspace(0,1,Nslices),'StartUpTime':mylogspace(0,36,Nslices),'MinUpTime':mylogspace(0,168,Nslices),'MinDownTime':mylogspace(0,168,Nslices),'NoLoadCost':np.linspace(0,50,Nslices),'StartUpCost':np.linspace(0,500,Nslices),'Efficiency':np.linspace(0,1,Nslices)}
    
    # Definition of the fingerprint value of each power plant, i.e. the pattern of the slices number in which each of
    # its characteristics falls:
    fingerprints = []
    fingerprints_merged = []
    for i in plants.index:
        fingerprints.append([find_nearest(bounds['PartLoadMin'],plants['PartLoadMin'][i]),find_nearest(bounds['RampUpRate'],plants['RampUpRate'][i]),find_nearest(bounds['RampDownRate'],plants['RampDownRate'][i]),find_nearest(bounds['StartUpTime'],plants['StartUpTime'][i]),find_nearest(bounds['MinUpTime'],plants['MinUpTime'][i]),find_nearest(bounds['MinDownTime'],plants['MinDownTime'][i]),find_nearest(bounds['NoLoadCost'],plants['NoLoadCost'][i]),find_nearest(bounds['StartUpCost'],plants['StartUpCost'][i]),find_nearest(bounds['Efficiency'],plants['Efficiency'][i])])
    
    # Definition of the merged power plants dataframe:
    plants_merged = pd.DataFrame(columns=plants.columns)
    merged_arrays = []
    
    # Find the columns containing string values (in addition to "Unit")
    string_keys = []
    for i in range(len(plants.columns)):
        if plants.columns[i] != 'Unit' and plants.dtypes[i] == np.dtype('O'):
            string_keys.append(plants.columns[i])
    
    for i in plants.index:                 # i is the plant to be added to the new list
        merged=False
        for j in plants_merged.index:    # j corresponds to the clustered plants
            same_type = all([plants[key][i] == plants_merged[key][j] for key in string_keys])
            same_fingerprint = (fingerprints[i]==fingerprints_merged[j])
            low_pmin = (plants['PartLoadMin'][i] <= 0.01 )
            low_pmax = (plants['PowerCapacity'][i] <= 30)
            highly_flexible = (plants['RampUpRate'][i] > 1/60 and (plants['RampDownRate'][i] > 1/60) and (plants['StartUpTime'][i]<1))
            if same_type and same_fingerprint and low_pmin or same_type and highly_flexible or same_type and low_pmax:     # merge the two plants in plants_merged:
                P_old = plants_merged['PowerCapacity'][j]                        # Old power in plants_merged
                P_add = plants['PowerCapacity'][i]                      # Additional power to be added
                for key in plants_merged:
                    if key in ['RampUpRate','RampDownRate','MinUpTime','MinDownTime','NoLoadCost','Efficiency','MinEfficiency']:
                        # Do a weighted average:
                        plants_merged.loc[j,key]= (plants_merged[key][j] * P_old + plants[key][i] * P_add)/(P_add + P_old)
                    elif key in ['PowerCapacity']:
                        # Do a sum:
                        plants_merged.loc[j,key] = plants_merged[key][j] + plants[key][i]
                    elif key in ['PartLoadMin','StartUpTime']:
                        # Take the minimum
                        plants_merged.loc[j,key] = np.minimum(plants_merged[key][j]*P_old, plants[key][i] * P_add)/(P_add + P_old)
                    elif key == 'RampingCost':
                        # The starting cost must be added to the ramping cost
                        Cost_to_fullload = P_add * (1 -  plants['PartLoadMin'][i]) * plants['RampingCost'][i] + plants['StartUpCost'][i]
                        plants_merged.loc[j,key] = (P_old * plants_merged[key][j] + Cost_to_fullload)/(P_old+P_add)
                map_old_new[i]=j
                map_plant_orig[j].append(i)
                merged = True
                break
        if not merged:       # Add a new plant in plants_merged:
            plants_merged = plants_merged.append(plants.loc[i],ignore_index=True)
            plants_merged = plants_merged.copy()
            map_plant_orig.append([i])
            map_old_new[i]=len(map_plant_orig)-1
            fingerprints_merged.append(fingerprints[i])
    
    Nunits_merged = len(plants_merged)
        
    for array in AdditionalArrays:
        array_merged = np.zeros([Nunits_merged,array.shape[1]])
        for j in range(Nunits_merged):
            idx_orig = map_plant_orig[j]
            value = np.zeros(array.shape[1])
            for idx in idx_orig:
                value = value + array[idx,:]*np.maximum(1e-9,plants['PowerCapacity'][idx])
            P_j = np.sum(np.maximum(1e-9,plants['PowerCapacity'].values[idx_orig]))
            array_merged[j,:] = value/P_j
        merged_arrays.append(array_merged)
        
    # Modify the Unit names with the original index number. In case of merged plants, indicate all indexes + the plant type and fuel
    for j in range(Nunits_merged):
        if len(map_plant_orig[j])==1:                                     # The plant has not been merged
            plants_merged.loc[j,'Unit'] = str(map_plant_orig[j]) + ' - ' + plants_merged['Unit'][j]
        else:
            all_stringkeys = ''
            for key in string_keys:
                all_stringkeys = all_stringkeys + ' - ' + plants_merged[key][j]  
            plants_merged.loc[j,'Unit'] = str(map_plant_orig[j]) + all_stringkeys
        
    return([plants_merged,merged_arrays])
    
    
    
    
    
    