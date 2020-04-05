import datetime as dt
import logging
import os
import sys

import numpy as np
import pandas as pd

from ..common import commons
try:
    from future.builtins import int
except ImportError:
    pass

DEFAULTS = {'ReservoirLevelInitial':0.5,'ReservoirLevelFinal':0.5,'ValueOfLostLoad':1E5,
                'PriceOfSpillage':1,'WaterValue':100,'ShareOfQuickStartUnits':0.5,
                'PriceOfNuclear':0,'PriceOfBlackCoal':0,'PriceOfGas':0,'PriceOfFuelOil':0,'PriceOfBiomass':0,
                'PriceOfCO2':0,'PriceOfLignite':0,'PriceOfPeat':0,'LoadShedding':0,'CostHeatSlack':0,
                'CostLoadShedding':100,'ShareOfFlexibleDemand':0,'DemandFlexibility':0,'PriceTransmission':0,
                'CostH2Slack':75}

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
                logging.critical('No data file found for the table ' + varname + ' and zone ' + z + '. File ' + path_c + ' does not exist')
                sys.exit(1)
        SingleFile=False
    elif path != '':
        logging.critical('A path has been specified for table ' + varname + ' (' + path + ') but no file has been found')
        sys.exit(1)
    data = pd.DataFrame(index=config['idx_long'])
    if len(paths) == 0:
        logging.info('No data file specified for the table ' + varname + '. Using default value ' + str(default))
        if default is None:
            pass
        elif isinstance(default,(float,int)):
            data = pd.DataFrame(default,index=config['idx_long'],columns=zones)
        else:
            logging.critical('Default value provided for table ' + varname + ' is not valid')
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
                        logging.critical('Default value provided for table ' + varname + ' is not valid')
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
                logging.error('No data file found for the table ' + varname + ' and zone ' + z + '. File ' + path_c + ' does not exist')
#                sys.exit(1)
        SingleFile=False
    elif path != '':
        logging.critical('A path has been specified for table ' + varname + ' (' + path + ') but no file has been found')
        sys.exit(1)

    data = pd.DataFrame(index=config['idx_long'])
    if len(paths) == 0:
        logging.info('No data file specified for the table ' + varname + '. Using default value ' + str(default))
        if default is None:
            out = pd.DataFrame(index=config['idx_long'])
        elif isinstance(default,(float,int)):
            out = pd.DataFrame(default,index=config['idx_long'],columns=plants['Unit'])
        else:
            logging.critical('Default value provided for table ' + varname + ' is not valid')
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
        logging.critical('The column headers of table "' + varname + '" are not unique!. The following headers are duplicated: ' + str(out.columns.get_duplicates()))
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

    plants.index = range(len(plants))
    merged = pd.DataFrame(index=data.index)
    unitnames = plants.Unit.values.tolist()
    # First check the data:
    if not isinstance(data,pd.DataFrame):
        logging.critical('The input "' + tablename + '" to the merge_series function must be a dataframe')
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


def load_time_series(config,path,header='infer'):
    """
    Function that loads time series data, checks the compatibility of the indexes
    and guesses when no exact match between the required index and the data is 
    present
    """

    data = pd.read_csv(path, index_col=0, parse_dates=True, header=header)
    
    if not data.index.is_unique:
        logging.critical('The index of data file ' + path + ' is not unique. Please check the data')
        sys.exit(1)

    if not data.index.is_monotonic_increasing:
        logging.error('The index of data file ' + path + ' is not monotoneously increasing. Trying to check if it can be parsed with a "day first" format ')
        data = pd.read_csv(path, index_col=0, parse_dates=True, header=header, dayfirst=True)
        if not data.index.is_monotonic_increasing:
            logging.critical('Could not parse index of ' + path + '. To avoid problems make sure that you use the proper american date format (yyyy-mm-dd hh:mm:ss)')
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
            data.index = pd.date_range(start=dt.datetime(*(config['idx'][0].year,1,1,0,0)),
                                                        end=dt.datetime(*(config['idx'][0].year,12,31,23,59,59)),
                                                        freq=commons['TimeStep'])
        else:
            logging.critical('A numerical index has been found for file ' + path + 
                         '. However, its length does not allow guessing its timestamps. Please use a 8760 elements time series')
            sys.exit(1)

    if data.index.is_all_dates:   
        data.index = data.index.tz_localize(None)   # removing locational data
        # Checking if the required index entries are in the data:
        common = data.index.intersection(config['idx'])
        if len(common) == 0:
            # check if original year is leap year and destination year is not (remove leap date)
            if (data.index[0].is_leap_year is True) and (config['idx'][0].is_leap_year is False):
                data = data[~((data.index.month == 2) & (data.index.day == 29))]
                logging.warning('File ' + path + ': data for year ' + str(data.index[0].year) +
                                ' is used instead of year ' + str(config['idx'][0].year))
                data.index = data.index.map(lambda t: t.replace(year=config['idx'][0].year))
            # check if both years are either leap or non leap
            elif (data.index[0].is_leap_year is True) and (config['idx'][0].is_leap_year is True) or \
                    (data.index[0].is_leap_year is False) and (config['idx'][0].is_leap_year is False):
                logging.warning('File ' + path + ': data for year ' + str(data.index[0].year) +
                                ' is used instead of year ' + str(config['idx'][0].year) +
                                '. Leap year date is removed from the original DataFrame.')
                data.index = data.index.map(lambda t: t.replace(year=config['idx'][0].year))
            # check if original year is not a leap year and destination year is a leap year (add leap date and take average hourly values between 28.02. and 1.3.
            elif (data.index[0].is_leap_year is False) and (config['idx'][0].is_leap_year is True):
                logging.warning('File ' + path + ': data for year ' + str(data.index[0].year) +
                                ' is used instead of year ' + str(config['idx'][0].year) +
                                '. Leap year date is interpolated between the two neighbouring days.')
                data.index = data.index.map(lambda t: t.replace(year=config['idx'][0].year))
                mask = data.loc[str(config['idx'][0].year)+'-2-28': str(config['idx'][0].year)+'-3-1']
                mask = mask.groupby(mask.index.hour).mean()
                time = pd.date_range(str(config['idx'][0].year)+'-2-29', periods=24, freq='H')
                mask = mask.set_index(time)
                data = data.reindex(config['idx'])
                data.update(mask)
        # recompute common index entries, and check again:
        common = data.index.intersection(config['idx'])
        if len(common) < len(config['idx'])-1:
            logging.critical('File ' + path + ': the index does not contain the necessary time range (from ' + str(config['idx'][0]) + ' to ' + str(config['idx'][-1]) + ')')
            sys.exit(1)
        elif len(common) == len(config['idx'])-1:  # there is only one data point missing. This is deemed acceptable
            logging.warning('File ' + path + ': there is one data point missing in the time series. It will be filled with the nearest data')
        else:
            pass              # the defined simulation index is found within the data. No action required
    else:
        logging.critical('Index for file ' + path + ' is not valid')
        sys.exit(1)
        
    # re-indexing with the longer index (including look-ahead) and filling possibly missing data at the beginning and at the end::
    return data.reindex(config['idx_long'], method='nearest').fillna(method='bfill')


def load_config(ConfigFile,AbsPath=True):
    """
    Wrapper function around load_config_excel and load_config_yaml
    """
    if ConfigFile.endswith(('.xlsx','.xls')):
        config = load_config_excel(ConfigFile,AbsPath=True)
    elif ConfigFile.endswith(('.yml','.yaml')):
        config = load_config_yaml(ConfigFile,AbsPath=True)
    else:
        logging.critical('The extension of the config file should be .xlsx or .yml')
        sys.exit(1)
    return config

def read_truefalse(sheet, rowstart, colstart, rowstop, colstop, colapart=1):
    """
    Function that reads a two column format with a list of strings in the first
    columns and a list of true false in the second column
    The list of strings associated with a True value is returned
    """
    out = []
    for i in range(rowstart, rowstop):
        if sheet.cell_value(i, colstart + colapart) == 1:
            out.append(sheet.cell_value(i, colstart))
    return out

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
    
    if sheet.cell_value(0,0) == 'Dispa-SET Configuration File (v20.01)':
        config['Description'] = sheet.cell_value(5, 1)
        config['StartDate'] = xlrd.xldate_as_tuple(sheet.cell_value(56, 2), wb.datemode)
        config['StopDate'] = xlrd.xldate_as_tuple(sheet.cell_value(57, 2), wb.datemode)
        config['HorizonLength'] = int(sheet.cell_value(58, 2))
        config['LookAhead'] = int(sheet.cell_value(59, 2))
        
        # Defning the input locations in the config file:
        StdParameters={'SimulationDirectory':33,'WriteGDX':34,'WritePickle':35,'GAMS_folder':36,
                          'cplex_path':37,'DataTimeStep':60,'SimulationTimeStep':61,
                          'SimulationType':76,'ReserveCalculation':77,'AllowCurtailment':78,
                          'HydroScheduling':98,'HydroSchedulingHorizon':99,'InitialFinalReservoirLevel':100}
        PathParameters={'Demand':124, 'Outages':126, 'PowerPlantData':127, 'RenewablesAF':128, 
                          'LoadShedding':129, 'NTC':130, 'Interconnections':131, 'ReservoirScaledInflows':132, 
                          'PriceOfNuclear':180, 'PriceOfBlackCoal':181, 'PriceOfGas':182, 
                          'PriceOfFuelOil':183,'PriceOfBiomass':184, 'PriceOfCO2':166, 
                          'ReservoirLevels':133, 'PriceOfLignite':185, 'PriceOfPeat':186,
                          'HeatDemand':134,'CostHeatSlack':165,'CostLoadShedding':168,'ShareOfFlexibleDemand':125,
                          'Temperatures':135,'PriceTransmission':169,'Reserve2U':160,'Reserve2D':161,
                          'H2Demand':136,'CostH2Slack':170}
        modifiers= {'Demand':274,'Wind':275,'Solar':276,'Storage':277}
        default = {'ReservoirLevelInitial':101,'ReservoirLevelFinal':102,'PriceOfNuclear':180,'PriceOfBlackCoal':181,
                    'PriceOfGas':182,'PriceOfFuelOil':183,'PriceOfBiomass':184,'PriceOfCO2':166,'PriceOfLignite':185,
                    'PriceOfPeat':186,'LoadShedding':129,'CostHeatSlack':167,'CostLoadShedding':168,'ValueOfLostLoad':204,
                    'PriceOfSpillage':205,'WaterValue':206,'ShareOfQuickStartUnits':163,'ShareOfFlexibleDemand':125,
                    'DemandFlexibility':162,'PriceTransmission':169,'CostH2Slack':170}
        for p in StdParameters:
            config[p] = sheet.cell_value(StdParameters[p], 2)
        for p in PathParameters:
            config[p] = sheet.cell_value(PathParameters[p], 2)
        config['modifiers'] = {}
        for p in modifiers:
            config['modifiers'][p] = sheet.cell_value(modifiers[p], 2)
        config['default'] = {}
        for p in default:
            config['default'][p] = sheet.cell_value(default[p], 5)
            
        #True/Falst values:
        config['zones'] = read_truefalse(sheet, 225, 1, 246, 3)
        config['zones'] = config['zones'] + read_truefalse(sheet, 225, 4, 246, 6)
        config['mts_zones'] = read_truefalse(sheet, 225, 1, 246, 3, 2)
        config['mts_zones'] = config['mts_zones'] + read_truefalse(sheet, 225, 4, 246, 6, 2)
        config['ReserveParticipation'] = read_truefalse(sheet, 305, 1, 319, 3)

        # Set default values (for backward compatibility):
        for param in DEFAULTS:
            if config['default'][param]=='':
                config['default'][param]=DEFAULTS[param]
                logging.warning('No value was provided in config file for {}. Will use {}'.format(param, DEFAULTS[param]))
                config['default'][param] = DEFAULTS[param]

        if AbsPath:
        # Changing all relative paths to absolute paths. Relative paths must be defined 
        # relative to the parent folder of the config file.
            abspath = os.path.abspath(ConfigFile)
            basefolder = os.path.abspath(os.path.join(os.path.dirname(abspath),os.pardir))
            if not os.path.isabs(config['SimulationDirectory']):
                config['SimulationDirectory'] = os.path.join(basefolder,config['SimulationDirectory'])
            for param in PathParameters:
                if config[param] == '' or config[param].isspace():
                    config[param] = ''
                elif not os.path.isabs(config[param]):
                    config[param] = os.path.join(basefolder,config[param])

        logging.info("Using config file (v20.01) " + ConfigFile + " to build the simulation environment")
        logging.info("Using " + config['SimulationDirectory'] + " as simulation folder")
        logging.info("Description of the simulation: "+ config['Description'])
        
        return config        
        

    elif sheet.cell_value(0,0) == 'Dispa-SET Configuration File':
        config['Description'] = sheet.cell_value(5, 1)
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
        config['DataTimeStep'] = sheet.cell_value(34, 2)
        config['SimulationTimeStep'] = sheet.cell_value(35, 2)
    
        config['SimulationType'] = sheet.cell_value(46, 2)
        config['ReserveCalculation'] = sheet.cell_value(47, 2)
        config['AllowCurtailment'] = sheet.cell_value(48, 2)
    
        config['HydroScheduling'] = sheet.cell_value(53, 2)
        config['HydroSchedulingHorizon'] = sheet.cell_value(54, 2)
        config['InitialFinalReservoirLevel'] = sheet.cell_value(55, 2)
    
        # Set default values (for backward compatibility):
        NonEmptyarameters = {'DataTimeStep':1,'SimulationTimeStep':1,'HydroScheduling':'Off','HydroSchedulingHorizon':'Annual','InitialFinalReservoirLevel':True}
        for param in NonEmptyarameters:
            if config[param]=='':
                config[param]=NonEmptyarameters[param]   
    
        # List of parameters for which an external file path must be specified:
        PARAMS = ['Demand', 'Outages', 'PowerPlantData', 'RenewablesAF', 'LoadShedding', 'NTC', 'Interconnections',
          'ReservoirScaledInflows', 'PriceOfNuclear', 'PriceOfBlackCoal', 'PriceOfGas', 'PriceOfFuelOil',
          'PriceOfBiomass', 'PriceOfCO2', 'ReservoirLevels', 'PriceOfLignite', 'PriceOfPeat','HeatDemand',
          'CostHeatSlack','CostLoadShedding','ShareOfFlexibleDemand']
        for i, param in enumerate(PARAMS):
            config[param] = sheet.cell_value(61 + i, 2)
    
        # List of new parameters for which an external file path must be specified:
        params2 = ['Temperatures','PriceTransmission','Reserve2D','Reserve2U','H2Demand','CostH2Slack']
        if sheet.nrows>150:                 # for backward compatibility (old excel sheets had less than 150 rows)
            for i, param in enumerate(params2):
                config[param] = sheet.cell_value(156 + i, 2)
        else:
            for param in params2:
                config[param] = ''
    
        if AbsPath:
        # Changing all relative paths to absolute paths. Relative paths must be defined 
        # relative to the parent folder of the config file.
            abspath = os.path.abspath(ConfigFile)
            basefolder = os.path.abspath(os.path.join(os.path.dirname(abspath),os.pardir))
            if not os.path.isabs(config['SimulationDirectory']):
                config['SimulationDirectory'] = os.path.join(basefolder,config['SimulationDirectory'])
            for param in PARAMS+params2:
                if config[param] == '' or config[param].isspace():
                    config[param] = ''
                elif not os.path.isabs(config[param]):
                    config[param] = os.path.join(basefolder,config[param])
    
        config['default'] = {}
        config['default']['ReservoirLevelInitial'] = sheet.cell_value(56, 5)
        config['default']['ReservoirLevelFinal'] = sheet.cell_value(57, 5)
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
        config['default']['ValueOfLostLoad'] = sheet.cell_value(81, 5)
        config['default']['PriceOfSpillage'] = sheet.cell_value(82, 5)
        config['default']['WaterValue'] = sheet.cell_value(83, 5)
        config['default']['ShareOfQuickStartUnits'] = 0.5          # to be added to xlsx file
        
        # Set default values (for backward compatibility):
        for param in DEFAULTS:
            if config['default'].get(param,'')=='':
                config['default'][param]=DEFAULTS[param]
    
        config['zones'] = read_truefalse(sheet, 86, 1, 109, 3)
        config['zones'] = config['zones'] + read_truefalse(sheet, 86, 4, 109, 6)
    
        config['mts_zones'] = read_truefalse(sheet, 86, 1, 109, 3, 2)
        config['mts_zones'] = config['mts_zones'] + read_truefalse(sheet, 86, 4, 109, 6, 2)
    
        config['modifiers'] = {}
        config['modifiers']['Demand'] = sheet.cell_value(111, 2)
        config['modifiers']['Wind'] = sheet.cell_value(112, 2)
        config['modifiers']['Solar'] = sheet.cell_value(113, 2)
        config['modifiers']['Storage'] = sheet.cell_value(114, 2)
    
        # Read the technologies participating to reserve markets:
        config['ReserveParticipation'] = read_truefalse(sheet, 131, 1, 145, 3)
    
        logging.info("Using config file " + ConfigFile + " to build the simulation environment")
        logging.info("Using " + config['SimulationDirectory'] + " as simulation folder")
        logging.info("Description of the simulation: "+ config['Description'])
        
        return config
    
    else:
        logging.critical('The format of the excel config file (defined by its main title) is not recognized')
        sys.exit(1)

def load_config_yaml(filename, AbsPath=True):
    """ Loads YAML file to dictionary"""
    import yaml
    with open(filename, 'r') as f:
        try:
            config = yaml.full_load(f)
        except yaml.YAMLError as exc:
            logging.error('Cannot parse config file: {}'.format(filename))
            raise exc
            
    # List of parameters to be added with a default value if not present (for backward compatibility):
    
    params_to_be_added = {'Temperatures':'','DataTimeStep':1,'SimulationTimeStep':1,'HydroScheduling':'Off','HydroSchedulingHorizon':'Annual','InitialFinalReservoirLevel':True}
    for param in params_to_be_added:
        if param not in config:
            config[param] = params_to_be_added[param]
                        
    # Set default values (for backward compatibility):
    NonEmptyDefaultss = {'ReservoirLevelInitial':0.5,'ReservoirLevelFinal':0.5,'ValueOfLostLoad':1E5,'PriceOfSpillage':1,'WaterValue':100,'ShareOfQuickStartUnits':0.5}
    for param in NonEmptyDefaultss:
        if param not in config['default']:
            config['default'][param]=NonEmptyDefaultss[param]


    # Define missing parameters if they were not provided in the config file
    PARAMS = ['Demand', 'Outages', 'PowerPlantData', 'RenewablesAF', 'LoadShedding', 'NTC', 'Interconnections',
          'ReservoirScaledInflows', 'PriceOfNuclear', 'PriceOfBlackCoal', 'PriceOfGas', 'PriceOfFuelOil',
          'PriceOfBiomass', 'PriceOfCO2', 'ReservoirLevels', 'PriceOfLignite', 'PriceOfPeat','HeatDemand',
          'CostHeatSlack','CostLoadShedding','ShareOfFlexibleDemand','Temperatures','PriceTransmission',
          'Reserve2D','Reserve2U','H2Demand','CostH2Slack']
    for param in PARAMS:
        if param not in config:
            config[param] = ''    
    global DEFAULTS
    for key in DEFAULTS:
        if key not in config['default']:
            config['default'][key]=DEFAULTS[key]

    if AbsPath:
    # Changing all relative paths to absolute paths. Relative paths must be defined 
    # relative to the parent folder of the config file.
        abspath = os.path.abspath(filename)
        basefolder = os.path.abspath(os.path.join(os.path.dirname(abspath),os.pardir))
        if not os.path.isabs(config['SimulationDirectory']):
            config['SimulationDirectory'] = os.path.join(basefolder,config['SimulationDirectory'])
        for param in PARAMS:
            if not os.path.isabs(config[param]):
                if config[param] == '' or config[param].isspace():
                    config[param] = ''
                elif not os.path.isabs(config[param]):
                    config[param] = os.path.join(basefolder,config[param])
    return config

def export_yaml_config(ExcelFile, YAMLFile):
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
