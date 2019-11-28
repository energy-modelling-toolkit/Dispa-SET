import logging
import os
import sys

import numpy as np
import pandas as pd

from six.moves import reload_module
from ..common import commons 
try:
    from future.builtins import int
except ImportError:
    pass


def NodeBasedTable(varname,config,default=None):
    '''
    This function loads the tabular data stored in csv files relative to each
    zone of the simulation.

    :param varname:             Variable name (as defined in config)
    :param idx:                 Pandas datetime index to be used for the output
    :param zones:               List with the zone codes to be considered
    :param fallback:            List with the order of data source.
    :param default:             Default value to be applied if no data is found

    :return:           Dataframe with the time series for each unit
    '''
    
    path = config[varname]
    zones=config['zones']       
    paths = {}
    if os.path.isfile(path):
        paths['all'] = path
        SingleFile=True
    elif '##' in path:
        for z in zones:
            path_c = path.replace('##', str(z))
            if os.path.isfile(path_c):
                paths[str(z)] = path_c
            else:
                logging.error('No data file found for the table ' + varname + ' and zone ' + z + '. File ' + path_c + ' does not exist')
                sys.exit(1)
        SingleFile=False
    data = pd.DataFrame(index=config['idx_long'])
    if len(paths) == 0:
        logging.info('No data file found for the table ' + varname + '. Using default value ' + str(default))
        if default is None:
            pass
        elif isinstance(default,(float,int)):
            data = pd.DataFrame(default,index=config['idx_long'],columns=zones)
        else:
            logging.error('Default value provided for table ' + varname + ' is not valid')
            sys.exit(1)
    elif SingleFile:
        # If it is only one file, there is a header with the zone code
        tmp = load_time_series(config,paths['all'])
           
        if len(tmp.columns) == 1:    # if there is only one column, assign its value to all the zones, whatever the header
            try:    # if the column header is numerical, there was probably no header. Load the file again.
                float(tmp.columns[0])   # this will fail if the header is not numerical
                tmp = load_time_series(config,paths['all'],header=None)
            except:
                pass
            for key in zones:
                data[key] = tmp.iloc[:,0]
        else:
            for key in zones:
                if key in tmp:
                    data[key] = tmp[key]
                else:
                    logging.error('Zone ' + key + ' could not be found in the file ' + path + '. Using default value ' + str(default))
                    if default is None:
                        pass
                    elif isinstance(default,(float,int)):
                        data[key] = default
                    else:
                        logging.error('Default value provided for table ' + varname + ' is not valid')
                        sys.exit(1)
    else: # assembling the files in a single dataframe:
        for z in paths:
            # In case of separated files for each zone, there is no header
            tmp = load_time_series(config,paths[z])  
            data[z] = tmp.iloc[:,0]

    return data



def UnitBasedTable(plants,varname,config,fallbacks=['Unit'],default=None,RestrictWarning=None):
    '''
    This function loads the tabular data stored in csv files and assigns the
    proper values to each unit of the plants dataframe. If the unit-specific 
    value is not found in the data, the script can fallback on more generic
    data (e.g. fuel-based, technology-based, zone-based) or to the default value.
    The order in which the data should be loaded is specified in the fallback
    list. For example, ['Unit','Technology'] means that the script will first
    try to find a perfect match for the unit name in the data table. If not found,
    a column with the unit technology as header is search. If not found, the
    default value is assigned.

    :param plants:              Dataframe with the units for which data is required
    :param varname:             Variable name (as defined in config)
    :param idx:                 Pandas datetime index to be used for the output
    :param zones:           List with the zone codes to be considered
    :param fallback:            List with the order of data source. 
    :param default:             Default value to be applied if no data is found
    :param RestrictWarning:     Only display the warnings if the unit belongs to the list of technologies provided in this parameter
    
    :return:           Dataframe with the time series for each unit
    '''
    path = config[varname]   
    zones = config['zones']    
    paths = {}
    if os.path.isfile(path):
        paths['all'] = path
        SingleFile=True
    elif '##' in path:
        for z in zones:
            path_c = path.replace('##', str(z))
            if os.path.isfile(path_c):
                paths[str(z)] = path_c
            else:
                logging.critical('No data file found for the table ' + varname + ' and zone ' + z + '. File ' + path_c + ' does not exist')
#                sys.exit(1)
        SingleFile=False
    data = pd.DataFrame(index=config['idx_long'])
    if len(paths) == 0:
        logging.info('No data file found for the table ' + varname + '. Using default value ' + str(default))
        if default is None:
            out = pd.DataFrame(index=config['idx_long'])
        elif isinstance(default,(float,int)):
            out = pd.DataFrame(default,index=config['idx_long'],columns=plants['Unit'])
        else:
            logging.error('Default value provided for table ' + varname + ' is not valid')
            sys.exit(1)
    else: # assembling the files in a single dataframe:
        columns = []
        for z in paths:
            tmp = load_time_series(config,paths[z])
            if SingleFile:
                for key in tmp:
                    data[key] = tmp[key]                
            else:    # use the multi-index header with the zone
                for key in tmp:
                    columns.append((z,key))
                    data[z+','+key] = tmp[key]
        if not SingleFile:
            data.columns = pd.MultiIndex.from_tuples(columns, names=['Zone', 'Data'])
        # For each plant and each fallback key, try to find the corresponding column in the data
        out = pd.DataFrame(index=config['idx_long'])
        for j in plants.index:
            warning = True
            if not RestrictWarning is None:
                warning = False
                if plants.loc[j,'Technology'] in RestrictWarning:
                    warning=True
            u = plants.loc[j,'Unit']
            found = False
            for i,key in enumerate(fallbacks):
                if SingleFile:
                    header = plants.loc[j,key]
                else:
                    header = (plants.loc[j,'Zone'],plants.loc[j,key])
                if header in data:
                    out[u] = data[header]
                    found = True
                    if i > 0 and warning:
                        logging.warning('No specific information was found for unit ' + u + ' in table ' + varname + '. The generic information for ' + str(header) + ' has been used')
                    break
            if not found:
                if warning:
                    logging.info('No specific information was found for unit ' + u + ' in table ' + varname + '. Using default value ' + str(default))
                if not default is None:
                    out[u] = default
    if not out.columns.is_unique:
        logging.error('The column headers of table "' + varname + '" are not unique!. The following headers are duplicated: ' + str(out.columns.get_duplicates()))
        sys.exit(1)
    return out



def merge_series(plants, data, mapping, method='WeightedAverage', tablename=''):
    """
    Function that merges the times series corresponding to the merged units (e.g. outages, inflows, etc.)

    :param plants:      Pandas dataframe with the information relative to the original units
    :param data:        Pandas dataframe with the time series and the original unit names as column header
    :param mapping:     Mapping between the merged units and the original units. Output of the clustering function
    :param method:      Select the merging method ('WeightedAverage'/'Sum')
    :param tablename:   Name of the table being processed (e.g. 'Outages'), used in the warnings
    :return merged:     Pandas dataframe with the merged time series when necessary
    """
    # backward compatibility:
    if not "Nunits" in plants:
        plants['Nunits'] = 1

    plants.index = range(len(plants))
    merged = pd.DataFrame(index=data.index)
    unitnames = plants.Unit.values.tolist()
    # First check the data:
    if not isinstance(data,pd.DataFrame):
        logging.error('The input "' + tablename + '" to the merge_series function must be a dataframe')
        sys.exit(1)
    for key in data:
        if str(data[key].dtype) not in ['bool','int','float','float16', 'float32', 'float64', 'float128','int8', 'int16', 'int32', 'int64']:
            logging.critical('The column "' + str(key) + '" of table + "' + tablename + '" is not numeric!')
    for key in data:
        if key in unitnames:
            i = unitnames.index(key)
            newunit = mapping['NewIndex'][i]
            if newunit not in merged:  # if the columns name is in the mapping and the new unit has not been processed yet
                oldindexes = mapping['FormerIndexes'][newunit]
                oldnames = [plants['Unit'][x] for x in oldindexes]
                if all([name in data for name in oldnames]):
                    subunits = data[oldnames]
                else:
                    for name in oldnames:
                        if name not in data:
                            logging.critical('The column "' + name + '" is required for the aggregation of unit "' + key +
                                             '", but it has not been found in the input data')
                            sys.exit(1)
                value = np.zeros(len(data))
                # Renaming the subunits df headers with the old plant indexes instead of the unit names:
                subunits.columns = mapping['FormerIndexes'][newunit]
                if method == 'WeightedAverage':
                    for idx in oldindexes:
                        name = plants['Unit'][idx]
                        value = value + subunits[idx] * np.maximum(1e-9, plants['PowerCapacity'][idx]*plants['Nunits'][idx])
                    P_j = np.sum(np.maximum(1e-9, plants['PowerCapacity'][oldindexes]*plants['Nunits'][oldindexes]))
                    merged[newunit] = value / P_j
                elif method == 'Sum':
                    merged[newunit] = subunits.sum(axis=1)
                else:
                    logging.critical('Method "' + str(method) + '" unknown in function MergeSeries')
                    sys.exit(1)
        elif key in plants['Unit']:
            if not isinstance(key, tuple):  # if the columns header is a tuple, it does not come from the data and has been added by Dispa-SET
                logging.warning('Column ' + str(key) + ' present in the table "' + tablename + '" not found in the mapping between original and clustered units. Skipping')
        else:
            if not isinstance(key, tuple):  # if the columns header is a tuple, it does not come from the data and has been added by Dispa-SET
                logging.warning('Column ' + str(key) + ' present in the table "' + tablename + '" not found in the table of power plants. Skipping')
    return merged


def define_parameter(sets_in, sets, value=0):
    """
    Function to define a DispaSET parameter and fill it with a constant value

    :param sets_in:     List with the labels of the sets corresponding to the parameter
    :param sets:        dictionary containing the definition of all the sets (must comprise those referenced in sets_in)
    :param value:       Default value to attribute to the parameter
    """
    if value == 'bool':
        values = np.zeros([len(sets[setx]) for setx in sets_in], dtype='bool')
    elif value == 0:
        values = np.zeros([len(sets[setx]) for setx in sets_in])
    elif value == 1:
        values = np.ones([len(sets[setx]) for setx in sets_in])
    else:
        values = np.ones([len(sets[setx]) for setx in sets_in]) * value
    return {'sets': sets_in, 'val': values}


def invert_dic_df(dic,tablename=''):
    """
    Function that takes as input a dictionary of dataframes, and inverts the key of
    the dictionary with the columns headers of the dataframes

    :param dic: dictionary of dataframes, with the same columns headers and the same index
    :param tablename: string with the name of the table being processed (for the error msg)
    :returns: dictionary of dataframes, with swapped headers
    """
    # keys are defined as the keys of the original dictionary, cols are the columns of the original dataframe
    # items are the keys of the output dictionary, i.e. the columns of the original dataframe
    dic_out = {}
    # First, check that all indexes have the same length:
    index = dic[dic.keys()[0]].index
    for key in dic:
        if len(dic[key].index) != len(index):
            logging.error('The indexes of the data tables "' + tablename + '" are not equal in all the files')
            sys.exit(1)
    # Then put the data in a panda Panel with minor orientation:
    panel = pd.Panel.fromDict(dic, orient='minor')
    # Display a warning if some items are missing in the original data:
    for item in panel.items:
        for key in dic.keys():
            if item not in dic[key].columns:
                logging.warning('The column "' + item + '" is not present in "' + key + '" for the "' + tablename + '" data. Zero will be assumed')
        dic_out[item] = panel[item].fillna(0)
    return dic_out



def load_time_series(config,path,header='infer'):
    """
    Function that loads time series data, checks the compatibility of the indexes
    and guesses when no exact match between the required index and the data is 
    present
    """

    data = pd.read_csv(path, index_col=0, parse_dates=True, header=header)
    
    if not data.index.is_unique:
        logging.error('The index of data file ' + path + ' is not unique. Please check the data')
        sys.exit(1)  
        
    if not data.index.is_monotonic_increasing:
        logging.error('The index of data file ' + path + ' is not monotoneously increasing. Please check that you have used the proper american date format (yyyy-mm-dd hh:mm:ss)')
        sys.exit(1)  
        
    # First convert numerical indexes into datetimeindex:
    if data.index.is_numeric():
        if len(data) == len(config['idx']):  # The data has the same length as the provided index range
            logging.info('A numerical index has been found for file ' + path + 
                         '. It has the same length as the target simulation and is therefore considered valid')
            data.index=config['idx']
        elif len(data) == 8760:
            logging.info('A numerical index has been found for file ' + path + 
                         '. Since it contains 8760 elements, it is assumed that it corresponds to a whole year')
            data.index = pd.DatetimeIndex(start=pd.datetime(*(config['idx'][0].year,1,1,0,0)),
                                                        end=pd.datetime(*(config['idx'][0].year,12,31,23,59,59)),
                                                        freq=commons['TimeStep'])
        else:
            logging.error('A numerical index has been found for file ' + path + 
                         '. However, its length does not allow guessing its timestamps. Please use a 8760 elements time series')
            sys.exit(1)

    if data.index.is_all_dates:   
        data.index = data.index.tz_localize(None)   # removing locational data
        # Checking if the required index entries are in the data:
        common = data.index.tz_localize(None).intersection(config['idx'])
        if len(common) == 0:
            # try to see if it is just a year-mismatch
            index2 = data.index.shift(8760 * (config['idx'][0].year - data.index[0].year),freq=commons['TimeStep'])
            common2 = index2.intersection(config['idx'])
            if len(common2) == len(config['idx']):
                logging.warn('File ' + path + ': data for year '+ str(data.index[0].year) + ' is used instead of year ' + str(config['idx'][0].year))
                data.index=index2
        elif len(common) == len(config['idx'])-1:  # there is only one data point missing. This is deemed acceptable
            logging.warn('File ' + path + ': there is one data point missing in the time series. It will be filled with the nearest data')
        elif len(common) < len(config['idx'])-1:
            logging.error('File ' + path + ': the index does not contain the necessary time range (from ' + str(config['idx'][0]) + ' to ' + str(config['idx'][-1]) + ')')
            sys.exit(1)
    else:
        logging.error('Index for file ' + path + ' is not valid')
        sys.exit(1)
        
    # re-indexing with the longer index (including look-ahead) and filling possibly missing data at the beginning and at the end::
    return data.reindex(config['idx_long'], method='nearest').fillna(method='bfill')



def load_config_excel(ConfigFile,AbsPath=True):
    """
    Function that loads the DispaSET excel config file and returns a dictionary
    with the values

    :param ConfigFile: String with (relative) path to the DispaSET excel configuration file
    :param AbsPath:    If true, relative paths are automatically changed into absolute paths (recommended)
    """
    import xlrd
    wb = xlrd.open_workbook(filename=ConfigFile)  # Option for csv to be added later
    sheet = wb.sheet_by_name('main')

    config = {}
    config['SimulationDirectory'] = sheet.cell_value(17, 2)
    config['WriteExcel'] = sheet.cell_value(18, 2)
    config['WriteGDX'] = sheet.cell_value(19, 2)
    config['WritePickle'] = sheet.cell_value(20, 2)
    config['GAMS_folder'] = sheet.cell_value(21, 2)
    config['cplex_path'] = sheet.cell_value(22, 2)

    config['StartDate'] = xlrd.xldate_as_tuple(sheet.cell_value(30, 2), wb.datemode)
    config['StopDate'] = xlrd.xldate_as_tuple(sheet.cell_value(31, 2), wb.datemode)
    config['HorizonLength'] = int(sheet.cell_value(32, 2))
    config['LookAhead'] = int(sheet.cell_value(33, 2))

    config['SimulationType'] = sheet.cell_value(46, 2)
    config['ReserveCalculation'] = sheet.cell_value(47, 2)
    config['AllowCurtailment'] = sheet.cell_value(48, 2)

    # List of parameters for which an external file path must be specified:
    params = ['Demand', 'Outages', 'PowerPlantData', 'RenewablesAF', 'LoadShedding', 'NTC', 'Interconnections',
              'ReservoirScaledInflows', 'PriceOfNuclear', 'PriceOfBlackCoal', 'PriceOfGas', 'PriceOfFuelOil',
              'PriceOfBiomass', 'PriceOfCO2', 'ReservoirLevels', 'PriceOfLignite', 'PriceOfPeat','HeatDemand',
              'CostHeatSlack','CostLoadShedding']
    for i, param in enumerate(params):
        config[param] = sheet.cell_value(61 + i, 2)

    if AbsPath:
    # Changing all relative paths to absolute paths. Relative paths must be defined 
    # relative to the parent folder of the config file.
        abspath = os.path.abspath(ConfigFile)
        basefolder = os.path.abspath(os.path.join(os.path.dirname(abspath),os.pardir))
        if not os.path.isabs(config['SimulationDirectory']):
            config['SimulationDirectory'] = os.path.join(basefolder,config['SimulationDirectory'])
        for param in params:
            if not os.path.isabs(config[param]):
                config[param] = os.path.join(basefolder,config[param])

    config['default'] = {}
    config['default']['PriceOfNuclear'] = sheet.cell_value(69, 5)
    config['default']['PriceOfBlackCoal'] = sheet.cell_value(70, 5)
    config['default']['PriceOfGas'] = sheet.cell_value(71, 5)
    config['default']['PriceOfFuelOil'] = sheet.cell_value(72, 5)
    config['default']['PriceOfBiomass'] = sheet.cell_value(73, 5)
    config['default']['PriceOfCO2'] = sheet.cell_value(74, 5)
    config['default']['PriceOfLignite'] = sheet.cell_value(76, 5)
    config['default']['PriceOfPeat'] = sheet.cell_value(77, 5)
    config['default']['LoadShedding'] = sheet.cell_value(65, 5)
    config['default']['CostHeatSlack'] = sheet.cell_value(79, 5)
    config['default']['CostLoadShedding'] = sheet.cell_value(80, 5)

    # read the list of zones to consider:
    def read_truefalse(sheet, rowstart, colstart, rowstop, colstop):
        """
        Function that reads a two column format with a list of strings in the first
        columns and a list of true false in the second column
        The list of strings associated with a True value is returned
        """
        out = []
        for i in range(rowstart, rowstop):
            if sheet.cell_value(i, colstart + 1) == 1:
                out.append(sheet.cell_value(i, colstart))
        return out

    config['zones'] = read_truefalse(sheet, 86, 1, 109, 3)
    config['zones'] = config['zones'] + read_truefalse(sheet, 86, 4, 109, 6)

    config['modifiers'] = {}
    config['modifiers']['Demand'] = sheet.cell_value(111, 2)
    config['modifiers']['Wind'] = sheet.cell_value(112, 2)
    config['modifiers']['Solar'] = sheet.cell_value(113, 2)
    config['modifiers']['Storage'] = sheet.cell_value(114, 2)

    # Read the technologies participating to reserve markets:
    config['ReserveParticipation'] = read_truefalse(sheet, 131, 1, 145, 3)

    logging.info("Using config file " + ConfigFile + " to build the simulation environment")
    logging.info("Using " + config['SimulationDirectory'] + " as simulation folder")

    return config

def load_config_yaml(filename,AbsPath=True):
    """ Loads YAML file to dictionary"""
    import yaml
    with open(filename, 'r') as f:
        try:
            config = yaml.load(f)
        except yaml.YAMLError as exc:
            logging.error('Cannot parse config file: {}'.format(filename))
            raise exc

    # List of parameters for which an external file path must be specified:
    params = ['Demand', 'Outages', 'PowerPlantData', 'RenewablesAF', 'LoadShedding', 'NTC', 'Interconnections',
              'ReservoirScaledInflows', 'PriceOfNuclear', 'PriceOfBlackCoal', 'PriceOfGas', 'PriceOfFuelOil',
              'PriceOfBiomass', 'PriceOfCO2', 'ReservoirLevels', 'PriceOfLignite', 'PriceOfPeat','HeatDemand',
              'CostHeatSlack','CostLoadShedding']

    if AbsPath:
    # Changing all relative paths to absolute paths. Relative paths must be defined 
    # relative to the parent folder of the config file.
        abspath = os.path.abspath(filename)
        basefolder = os.path.abspath(os.path.join(os.path.dirname(abspath),os.pardir))
        if not os.path.isabs(config['SimulationDirectory']):
            config['SimulationDirectory'] = os.path.join(basefolder,config['SimulationDirectory'])
        for param in params:
            if not os.path.isabs(config[param]):
                config[param] = os.path.join(basefolder,config[param])

    return config

def export_yaml_config(ExcelFile,YAMLFile):
    """
    Function that loads the DispaSET excel config file and dumps it as a yaml file.

    :param ExcelFile:   Path to the Excel config file
    :param YAMLFile:    Path to the YAML config file to be written
    """
    import yaml
    config = load_config_excel(ExcelFile,AbsPath=False)
    with open(YAMLFile, 'w') as outfile:
        yaml.dump(config, outfile, default_flow_style=False)
    return True
