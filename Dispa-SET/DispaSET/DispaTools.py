'''
This files gathers different functions used in the DispaSET pre-processing and post-processing tools

@author: Sylvain Quoilin (sylvain.quoilin@ec.europa.eu)
'''

from __future__ import division
import os
import sys
import numpy as np
import pandas as pd 
from warnings import warn

def clean_strings(x,exclude_digits=False,exclude_punctuation=False):
    """
    Function to convert strange unicode
    and remove characters punctuation

    :param x: any string or list of strings  

    Usage:: 
        
        df['DisplayName'].apply(clean_strings) 
        
    """    
    import unicodedata
    import string
    def clean_singlestring(x):
        if exclude_digits:            # modify the following depending on what you need to exclude 
            exclude1 = set(string.punctuation)
            # exception to the exclusion:
            exclude1.remove('_')
            exclude1.remove('-')
            exclude1.remove('[')
            exclude1.remove(']')
        else:
            exclude1 = set([])
        if exclude_punctuation:
            exclude2 = set(string.digits)
        else:
            exclude2 = set([])
        exclude = exclude1 | exclude2
    
        #http://stackoverflow.com/questions/2365411/python-convert-unicode-to-ascii-without-errors
        x = str(x).decode('utf-8') # to string byte and then unicode
        x = unicodedata.normalize('NFKD', x).encode('ascii', 'ignore') #convert utf characters and to ascii
        
        #x = x.upper() #to UPPERCASE
        x = ''.join(ch for ch in x if ch not in exclude) #remove numbers and punctuation
        return x
    if type(x) == str:
        return clean_singlestring(x)
    elif type(x) == list:
        return [clean_singlestring(xx) for xx in x]
    else:
        sys.exit('Argument type not supported')


def shrink_to_64(x,N=64):
    '''
    Function that reduces the length of the keys to be written to 64 (max admissible length for GAMS)

    :param x:   String or list of strings
    :param N:   Integer with the maximum string length (if different from 64)
    
    :returns:   Shrinked string or list of strings
    '''
    def shrink_singlestring(key,N):
        if len(key) >= N:
            return key[:20] + ' ... ' + key[-20:]
        else:
            return key
    if type(x) == str:
        return shrink_singlestring(x,N)
    elif type(x) == list:
        return [shrink_singlestring(xx,N) for xx in x]
    else:
        sys.exit('Argument type not supported')



def ds_to_df(inputs):
    '''
    Function that converts the dispaset data format into a dictionayr of dataframes

    :param list_vars: List containing the dispaset variables 
    :param format: Format version of the dispaset variables
    :param timeindex:   Time index to be applied to the dataframes. The size of the output is adapted accordingly.
    :return: dictionary of dataframes
    '''
    
    sets, parameters = inputs['sets'],inputs['parameters']
    
    #config = parameters['Config']['val']
    try:
        config = inputs['config']
        first_day = pd.datetime(config['StartDate'][0],config['StartDate'][1],config['StartDate'][2],0)
        last_day = pd.datetime(config['StopDate'][0],config['StopDate'][1],config['StopDate'][2],23)
        dates = pd.date_range(start=first_day,end=last_day,freq='1h')
    except:
        warn('Could not find the start/stop date information in the inputs. Using an integer index')
        dates = range(1,len(sets['z'])+1)
    if len(dates) > len(sets['h']):
        sys.exit('The provided index has a length of ' + str(len(dates)) + ' while the data only comprises ' +  str(len(sets['h'])) +  ' time elements')
    elif len(dates) > len(sets['z']):
        warn('The provided index has a length of ' + str(len(dates)) + ' while the simulation was designed for ' +  str(len(sets['z'])) +  ' time elements')
    elif len(dates) < len(sets['z']):
        warn('The provided index has a length of ' + str(len(dates)) + ' while the simulation was designed for ' +  str(len(sets['z'])) +  ' time elements')

    idx = range(len(dates))

    out={}
    out['sets'] = sets

    # Printing each parameter in a separate sheet and workbook:
    for p in parameters:
        var = parameters[p]
        dim = len(var['sets'])
        if var['sets'][-1] =='h' and isinstance(dates,pd.tseries.index.DatetimeIndex) and dim > 1:
            #if len(dates) != var['val'].shape[-1]:
            #    sys.exit('The date range in the Config variable (' + str(len(dates)) + ' time steps) does not match the length of the time index (' + str(var['val'].shape[-1]) + ') for variable ' + p)
            var['firstrow'] = 5
        else:
            var['firstrow']=1
        if dim == 1:
            if var['sets'][0] == 'h':
                out[p] = pd.DataFrame(var['val'][idx],columns=[p],index=dates)
            else:
                out[p] = pd.DataFrame(var['val'],columns=[p],index=sets[var['sets'][0]])
        elif dim ==2:  
            values = var['val']
            list_sets = [sets[var['sets'][0]], sets[var['sets'][1]]]
            if var['sets'][1] == 'h':
                out[p] = pd.DataFrame(values.transpose()[idx,:],index=dates,columns=list_sets[0])
            else:
                out[p] = pd.DataFrame(values.transpose(),index=list_sets[1],columns=list_sets[0])
        elif dim ==3:    
            list_sets = [sets[var['sets'][0]], sets[var['sets'][1]], sets[var['sets'][2]]]
            values = var['val']
            values2 = np.zeros([len(list_sets[0])*len(list_sets[1]),len(list_sets[2])])
            cols = np.zeros([2,len(list_sets[0])*len(list_sets[1])])
            for i in range(len(list_sets[0])):
                values2[i*len(list_sets[1]):(i+1)*len(list_sets[1]),:] = values[i,:,:]
                cols[0,i*len(list_sets[1]):(i+1)*len(list_sets[1])] = i
                cols[1,i*len(list_sets[1]):(i+1)*len(list_sets[1])] = range(len(list_sets[1]))
                
            columns = pd.MultiIndex([list_sets[0],list_sets[1]],cols)
            if var['sets'][2] == 'h':
                out[p] = pd.DataFrame(values2.transpose()[idx,:],index=dates,columns=columns)
            else:
                out[p] = pd.DataFrame(values2.transpose(),index=list_sets[2],columns=columns)
        else:
            sys.exit('Only three dimensions currently supported. Parameter ' + p + ' has ' + str(dim) + ' dimensions.' )
    return out


def InputDir(msg=None):
    '''
    Function that requires the user to input a directory
    '''
    if msg == None:
        msg = 'Please enter a directory: '
    while True:
        x = raw_input(msg)
        if os.path.isdir(x):
            break
        print x + 'is not a not a valid Directory.  Try again...'
    return x
    
def InputFile(msg=None):
    '''
    Function that requires the user to input a file path
    '''
    if msg == None:
        msg = 'Please enter a directory: '
    while True:
        x = raw_input(msg)
        if os.path.isfile(x):
            break
        print x + 'is not a file.  Try again...'
    return x


def ParamDefinition(sets_in,sets,value=0):
    '''
    Function to define a Dispa-SET parameter and fill it with a constant value

    :param sets_in:     List with the labels of the sets corresponding to the parameter
    :param sets:        Dictionnary containing the definition of all the sets (must comprise those referenced in sets_in)
    :param value:       Default value to attribute to the parameter
    '''
    if value == 'bool':
        values = np.zeros([len(sets[setx]) for setx in sets_in],dtype='bool')
    elif value == 0:
        values = np.zeros([len(sets[setx]) for setx in sets_in])
    elif value == 1:    
        values = np.ones([len(sets[setx]) for setx in sets_in])
    else:
        values = np.ones([len(sets[setx]) for setx in sets_in])*value
    return {'sets':sets_in,'val':values}



def incidence_matrix(sets,set_used,parameters,param_used):
    ''' 
    This function generates the incidence matrix of the lines within the nodes
    A particular case is considered for the node "Rest Of the World", which is no explicitely defined in Dispa-SET
    '''

    for i in range(len(sets[set_used])):
        if 'RoW' not in sets[set_used][i]:
            first_country = sets[set_used][i][0:2]    
            second_country = sets[set_used][i][6:8]
        elif 'RoW' == sets[set_used][i][0:3]:
            first_country = sets[set_used][i][0:3]
            second_country = sets[set_used][i][7:9]
        elif 'RoW' == sets[set_used][i][6:9]:
            first_country = sets[set_used][i][0:2]
            second_country = sets[set_used][i][6:9]
        else:
            sys.exit('The format of the interconnection is not admitted.')
            
        for j in range(len(sets['n'])):
            if first_country == sets['n'][j]:
                parameters[param_used]['val'][i,j] = -1
            elif second_country == sets['n'][j]:
                parameters[param_used]['val'][i,j] = 1
                
    return parameters[param_used]




def interconnections(Simulation_list,list_all_countries,NTC_inter,Historical_flows):
    '''
    Function that checks for the possible interconnections of the countries included
    in the simulation. If the interconnections occurs between two of the countries
    defined by the user to perform the simulation with, it extracts the NTC between 
    those two countries. If the interconnection occurs between one of the countries
    selected by the user and one country outside the simulation, it extracts the 
    physical flows; it does so for each pair (country inside-country outside) and 
    sums them together creating the interconnection of this country with the RoW.

    :param Simulation_list:     List of simulated countries
    :param list_all_countries:  List of all countries
    :param NTC:                 Day-ahead net transfer capacities (pd dataframe)
    :param Historical_flows:    Historical flows (pd dataframe)
    '''
    if len(NTC_inter.index) != len(Historical_flows.index) or NTC_inter.index[0] != Historical_flows.index[0]:
        sys.exit('The two input dataframes must have the same index')
    else: 
        index= NTC_inter.index
    all_connections = []
    simulation_connections = []
    for i in Simulation_list:
        for j in list_all_countries:
            if i != j:
                all_connections.append(i + ' -> ' + j)
                all_connections.append(j + ' -> ' + i)
        for k in Simulation_list:
            if i != k:
                simulation_connections.append(i + ' -> ' + k)       
    
    df_countries_simulated = pd.DataFrame(index=index)
    for interconnection in simulation_connections:
        if interconnection in NTC_inter.columns:
            df_countries_simulated[interconnection] = NTC_inter[interconnection]
    interconnections1 = df_countries_simulated.columns        
        
    df_RoW_temp = pd.DataFrame(index=index)
    connNames = []
    for interconnection in all_connections:
        if interconnection in Historical_flows.columns and interconnection not in simulation_connections:
            df_RoW_temp[interconnection] = Historical_flows[interconnection]
            connNames.append(interconnection)
            
    compare_set = set()
    for k in connNames:
        if not k[0:2] in compare_set and k[0:2] in Simulation_list:
            compare_set.add(k[0:2])
    
    df_countries_RoW = pd.DataFrame(index=index)
    while compare_set:
        nameToCompare = compare_set.pop()
        exports = []
        imports = []
        for name in connNames:
            if nameToCompare[0:2] in name[0:2]:
                exports.append(connNames.index(name))
            elif nameToCompare[0:2] in name[6:8]:
                imports.append(connNames.index(name))
                    
        flows_out = pd.concat(df_RoW_temp[connNames[exports[i]]] for i in range(len(exports)))
        flows_out = flows_out.groupby(flows_out.index).sum()
        flows_out.name = nameToCompare + ' -> RoW'
        df_countries_RoW[nameToCompare + ' -> RoW'] = flows_out
        flows_in = pd.concat(df_RoW_temp[connNames[imports[j]]] for j in range(len(imports)))
        flows_in = flows_in.groupby(flows_in.index).sum()
        flows_in.name = 'RoW -> ' + nameToCompare
        df_countries_RoW['RoW -> ' + nameToCompare] = flows_in
    interconnections2 = df_countries_RoW.columns   
    inter = list(interconnections1) + list(interconnections2)
    return(df_countries_simulated,df_countries_RoW,inter)




def get_gams_path():
    ''' 
    Function that attemps to search for the GAMS installation path (required to write the GDX or run gams)

    It returns the path if it has been found, or an empty string otherwise.

    Currently works for Windows and Linux. More searching rules and patterns should be added in the future    
    '''
    import sys
    import subprocess
    import os
    out = ''
    if sys.platform == 'linux2':
        try:
            proc = subprocess.Popen(['locate', '-i', 'libgamscall64.so'], stdout=subprocess.PIPE)
            tmp = proc.stdout.read()
        except:
            tmp = ''
        lines = tmp.split('\n')
        for line in lines:
            path = line.strip('libgamscall64.so')
            if os.path.exists(path):
                out = path
            break
    elif sys.platform == 'win32':
        paths = ['C:\\GAMS','C:\\Program Files\\GAMS','C:\\Program Files (x86)\\GAMS']
        lines = []
        for path in paths:
            if os.path.exists(path):
                paths1 = [path + os.sep + tmp for tmp in os.listdir(path) if tmp.startswith('win') and os.path.exists(path + os.sep + tmp)]
                for path1 in paths1:
                    lines = lines + [path1 + os.sep + tmp for tmp in os.listdir(path1) if tmp.startswith('24') and os.path.isfile(path1 + os.sep + tmp +  os.sep + 'gams.exe')]
        for line in lines:
            if os.path.exists(line):
                out = line
            break     
    elif sys.platform == 'darwin':
        paths = ['/Applications/']
        lines = []
        for path in paths:
            if os.path.exists(path):
                paths1 = [path + os.sep + tmp for tmp in os.listdir(path) if tmp.startswith('GAMS') and os.path.exists(path + os.sep + tmp)]
                for path1 in paths1:
                    lines = lines + [path1 + os.sep + tmp for tmp in os.listdir(path1) if tmp.startswith('sysdir') and os.path.isfile(path1 + os.sep + tmp +  os.sep + 'gams')]
        if len(lines)==0:
            proc = subprocess.Popen(['mdfind', '-name', 'libgamscall64.dylib'], stdout=subprocess.PIPE)
            tmp = proc.stdout.read()
            lines = [x.strip('libgamscall64.dylib') for x in tmp.split('\n')]
        for line in lines:
            if os.path.exists(line):
                out = line
            break                    
            
    if out != '':
        print 'Detected ' + out + ' as GAMS path on this computer'
    else:
        tmp = input('Specify the path to GAMS within quotes (e.g. "C:\\\\GAMS\\\\win64\\\\24.3"): ')
        if os.path.isdir(tmp):
            if sys.platform == 'win32':
                if os.path.isfile(tmp + os.sep + 'gams.exe'):
                    out = tmp
                else:
                    sys.exit('The provided path is not a valid windows gams folder')
            elif sys.platform == 'linux2':
                if os.path.isfile(tmp + os.sep + 'gamslib'):
                    out = tmp
                else:
                    sys.exit('The provided path is not a valid linux gams folder')
            else:
                if os.path.isdir(tmp):
                    out = tmp

    return out


def invert_dic_df(dic):
    '''
    Function that takes as input a dictionnary of dataframes, and inverts the key of
    the dictionnary with the columns headers of the dataframes

    :param dic: Dictionnary of dataframes, with the same columns headers and the same index
    :returns: Dictionnary of dataframes, with swapped headers
    '''
    cols = [aa for aa in dic[dic.keys()[0]].columns]
    index = dic[dic.keys()[0]].index
    cols_out = dic.keys()
    dic_out = {}
    for col in cols:
        dic_out[col] = pd.DataFrame(index=index,columns=cols_out)
    for key in dic:
        if [aa for aa in dic[key].columns] != cols:
            sys.exit('The columns of the dataframes within the dictionnary are not all equal')
        if len(dic[key].index) != len(index):
            sys.exit('The indexes of the dataframes within the dictionnary are not all equal')
        for col in cols:
            dic_out[col][key] = dic[key][col]
    return dic_out

def get_sets(instance, varname):
    '''Get sets that belong to a pyomo Variable or Param

    :param instance: Pyomo Instance
    :param varname: Name of the Pyomo Variable (string) 
    :return: A list with the sets that belong to this Param 
    '''
    var = getattr(instance, varname)

    if var.dim() > 1:
        sets = [pset.cname() for pset in var._index.set_tuple]
    else:
        sets = [var._index.name]
    return sets

def get_set_members(instance, sets):
    '''Get set members relative to a list of sets

    :param instance: Pyomo Instance
    :param sets: List of strings with the set names
    :return: A list with the set members  
    '''
    sm = []
    for s in sets:
        sm.append([v for v in getattr(instance, s).value])
    return sm


def pyomo_to_pandas(instance, varname):
    '''
    Function converting a pyomo variable or parameter into a pandas dataframe.
    The variable must have one or two dimensions and the sets must be provided as a list of lists
    
    :param instance: Pyomo Instance
    :param varname: Name of the Pyomo Variable (string)     
    '''
    setnames = get_sets(instance, varname)
    sets = get_set_members(instance, setnames)
    var = getattr(instance, varname)  # Previous script used model.var instead of var
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
    Function that flattens the multidimensional dispaset input data into the pyomo format: a dictionnary with a tuple and the parameter value. 
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

    
def load_csv(filename,TempPath='.pickle',header=0,skiprows=[],skip_footer=0,index_col=None,parse_dates=False):
    '''
    Function that loads an xls sheet into a dataframe and saves a temporary pickle version of it. 
    If the pickle is newer than the sheet, do no load the sheet again.

    :param file_excel: path to the excel file
    :param TempPath: path to store the temporary data files
    '''
    import os
    import pandas as pd
    
    filepath_pandas = TempPath + '/' + filename.replace('/','-') + '-' + '.p'
    if not os.path.isdir(TempPath):
        os.mkdir(TempPath)
    if not os.path.isfile(filepath_pandas):
       time_pd = 0
    else:
        time_pd = os.path.getmtime(filepath_pandas)
    if os.path.getmtime(filename) > time_pd:
        data = pd.read_csv(filename,header=header,skiprows=skiprows,skip_footer=skip_footer,index_col=index_col,parse_dates=parse_dates)
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
    
 
def clustering(plants,method='Standard',Nslices=20,PartLoadMax=0.1,Pmax=30):
    '''
    Merge excessively disaggregated power Units.

    :param plants:          Pandas dataframe with each power plant and their characteristics (following the DispaSET format)
    :param method:          Select clustering method ('Standard'/'LP'/None)
    :param Nslices:         Number of slices used to fingerprint each power plant characteristics. slices in the power plant data to categorize them  (fewer slices involves that the plants will be aggregated more easily)
    :param PartLoadMax:     Maximum part-load capability for the unit to be clustered
    :param Pmax:            Maximum power for the unit to be clustered
    :return:                A list with the merged plants and the mapping between the original and merged units
    '''
    from DispaSET_io_data import shrink_to_64
    
    if method == 'Standard':
        cluster=True
        LP=False
    elif method =='LP':
        cluster=True
        LP=True
    elif method==None:
        cluster=False
        LP=False
    else:
        sys.exit('Method argument not recognized in the clustering function')
    
    # Checking the the required columns are present in the input pandas dataframe:
    required_inputs = ['Unit','PowerCapacity','PartLoadMin','RampUpRate','RampDownRate','StartUpTime','MinUpTime','MinDownTime','NoLoadCost','StartUpCost','Efficiency']   
    for input in required_inputs:
        if input not in plants.columns:
            sys.exit("The plants dataframe requires a '" + input + "' column for clustering")
    
    # Number of units:
    Nunits = len(plants)
    plants.index = range(Nunits)
    
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
    
    # Find the columns containing string values (in addition to "Unit")
#    string_keys = []
#    for i in range(len(plants.columns)):
#        if plants.columns[i] != 'Unit' and plants.dtypes[i] == np.dtype('O'):
#            string_keys.append(plants.columns[i])
    string_keys = ['Zone','Technology','Fuel']
    
    for i in plants.index:                 # i is the plant to be added to the new list
        merged=False
        for j in plants_merged.index:    # j corresponds to the clustered plants
            same_type = all([plants[key][i] == plants_merged[key][j] for key in string_keys]) and cluster       # if clustering is off, all plants will be considered as different and will therefore not be merged
            same_fingerprint = (fingerprints[i]==fingerprints_merged[j])
            low_pmin = (plants['PartLoadMin'][i] <= PartLoadMax)
            low_pmax = (plants['PowerCapacity'][i] <= Pmax)
            highly_flexible = plants['RampUpRate'][i] > 1/60 and (plants['RampDownRate'][i] > 1/60) and (plants['StartUpTime'][i]<1) and (plants['MinDownTime'][i]<=1) and (plants['MinUpTime'][i]<=1)
            if (same_type and same_fingerprint and low_pmin) or (same_type and highly_flexible) or (same_type and low_pmax) or (same_type and LP):     # merge the two plants in plants_merged:
                P_old = plants_merged['PowerCapacity'][j]               # Old power in plants_merged
                P_add = plants['PowerCapacity'][i]                      # Additional power to be added
                for key in plants_merged:
                    if key in ['RampUpRate','RampDownRate','MinUpTime','MinDownTime','NoLoadCost','Efficiency','MinEfficiency','STOChargingEfficiency','CO2Intensity','STOSelfDischarge']:
                        # Do a weighted average:
                        plants_merged.loc[j,key]= (plants_merged[key][j] * P_old + plants[key][i] * P_add)/(P_add + P_old)
                    elif key in ['PowerCapacity','STOCapacity','STOMaxChargingPower']:
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
    mapping = {'NewIndex':{},'FormerIndexes':{} }
#    mapping['NewIdx'] = map_plant_orig
#    mapping['OldIdx'] = map_old_new
    # Modify the Unit names with the original index number. In case of merged plants, indicate all indexes + the plant type and fuel
    for j in range(Nunits_merged):
        if len(map_plant_orig[j])==1:                                     # The plant has not been merged
            NewName = str(map_plant_orig[j]) + ' - ' + plants_merged['Unit'][j]
            NewName = shrink_to_64(clean_strings(NewName))
            plants_merged.loc[j,'Unit'] = NewName
            mapping['FormerIndexes'][NewName] = [map_plant_orig[j][0]]
            mapping['NewIndex'][map_plant_orig[j][0]] = NewName
        else:
            all_stringkeys = ''
            for key in string_keys:
                all_stringkeys = all_stringkeys + ' - ' + plants_merged[key][j]  
            NewName = str(map_plant_orig[j]) + all_stringkeys
            NewName = shrink_to_64(clean_strings(NewName))
            plants_merged.loc[j,'Unit'] = NewName
            list_oldplants = [x for x in map_plant_orig[j]]
            mapping['FormerIndexes'][NewName] = list_oldplants
            for oldplant in list_oldplants:
                mapping['NewIndex'][oldplant] = NewName           
            
    if LP:
        for i in range(Nunits_merged):
            if plants_merged['RampingCost'][i] == 0:
                Power = plants_merged['PowerCapacity'][i]
                Start_up = plants_merged['StartUpCost'][i]
                plants_merged.loc[i,'RampingCost'] = Start_up/Power
    
    # Updating the index of the merged plants dataframe with the new unit names, after some cleaning:
    plants_merged.index = plants_merged['Unit']
    
    if Nunits != len(plants_merged):
        print('Clustered ' + str(Nunits) + ' original units into ' + str(len(plants_merged)) + ' new units')
    else:
        print('Did not cluster any unit')
    return plants_merged,mapping    
    



def MergeSeries(plants,data,mapping,method='WeightedAverage'):
    '''
    Function that merges the times series corresponding to the merged units (e.g. outages, inflows, etc.)

    :param plants:      Pandas dataframe with the information relative to the original units
    :param data:        Pandas dataframe with the time series and the original unit names as column header
    :param mapping:     Mapping between the merged units and the original units. Output of the clustering function
    :param method:      Select the merging method ('WeightedAverage'/'Sum') 
    :return merged:     Pandas dataframe with the merged time series when necessary
    '''

    Nunits = len(plants)
    plants.index = range(Nunits)
    merged = pd.DataFrame(index = data.index)
    unitnames = [plants['Unit'][x] for x in mapping['NewIndex']]
    for key in data:
        if key in unitnames:
            i = unitnames.index(key)
            newunit = mapping['NewIndex'][i]
            if newunit not in merged:    # if the columns name is in the mapping and the new unit has not been processed yet
                oldindexes = mapping['FormerIndexes'][newunit]
                oldnames = [plants['Unit'][x] for x in oldindexes]
                if all([name in data for name in oldnames]):
                    subunits = data[oldnames]
                else:
                    for name in oldnames:
                        if name not in data:
                            sys.exit('The column "' + name + '" is required for the aggregation of unit "' + key + '", but it has not been found in the input data')
                value = np.zeros(len(data))
                #Renaming the subunits df headers with the old plant indexes instead of the unit names:
                subunits.columns = mapping['FormerIndexes'][newunit]
                if method=='WeightedAverage':
                    for idx in oldindexes:
                        name = plants['Unit'][idx]
                        value = value + subunits[idx]*np.maximum(1e-9,plants['PowerCapacity'][idx])
                    P_j = np.sum(np.maximum(1e-9,plants['PowerCapacity'][oldindexes]))
                    merged[newunit] = value/P_j
                elif method=='Sum':
                    merged[newunit] = subunits.sum(axis=1)
                else:
                    sys.exit('Method "'+ str(method) + '" unknown in function MergeSeries')        
        else:
            if not type(key) == tuple:  # if the columns header is a tuple, it does not come from the data and has been added by Dispa-SET
                warn('Column ' + str(key) + ' not found in the mapping between original and clustered units. Skipping')
    return merged
    
    