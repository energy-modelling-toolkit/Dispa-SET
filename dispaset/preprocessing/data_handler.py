import logging
import os
import sys

import numpy as np
import pandas as pd

from six.moves import reload_module
try:
    from future.builtins import int
except ImportError:
    pass


def NodeBasedTable(path,idx,countries,tablename='',default=None):
    '''
    This function loads the tabular data stored in csv files relative to each
    zone (a.k.a node, country) of the simulation.

    :param path:                Path to the data to be loaded
    :param idx:                 Pandas datetime index to be used for the output
    :param countries:           List with the country codes to be considered
    :param fallback:            List with the order of data source.
    :param tablename:           String with the name of the table being processed
    :param default:             Default value to be applied if no data is found

    :return:           Dataframe with the time series for each unit
    '''
              
    paths = {}
    if os.path.isfile(path):
        paths['all'] = path
        SingleFile=True
    elif '##' in path:
        for c in countries:
            path_c = path.replace('##', str(c))
            if os.path.isfile(path_c):
                paths[str(c)] = path_c
            else:
                logging.error('No data file found for the table ' + tablename + ' and country ' + c + '. File ' + path_c + ' does not exist')
                sys.exit(1)
        SingleFile=False
    data = pd.DataFrame(index=idx)
    if len(paths) == 0:
        logging.info('No data file found for the table ' + tablename + '. Using default value ' + str(default))
        if default is None:
            pass
        elif isinstance(default,(float,int)):
            data = pd.DataFrame(default,index=idx,columns=countries)
        else:
            logging.error('Default value provided for table ' + tablename + ' is not valid')
            sys.exit(1)
    elif SingleFile:
        # If it is only one file, there is a header with the country code
        tmp = load_csv(paths['all'], index_col=0, parse_dates=True)
        if not tmp.index.is_unique:
            logging.error('The index of data file ' + paths['all'] + ' is not unique. Please check the data')
            sys.exit(1)
        for key in countries:
            if key in tmp:
                data[key] = tmp[key]
            elif len(tmp.columns) == 1:    # if the country code is not in the header, it can also be because it is a single country simulation and no header is needed:
                data[key] = tmp.iloc[:,0]
            else:
                logging.error('Country ' + key + ' could not be found in the file ' + path + '. Using default value ' + str(default))
                if default is None:
                    pass
                elif isinstance(default,(float,int,long)):
                    data[key] = default
                else:
                    logging.error('Default value provided for table ' + tablename + ' is not valid')
                    sys.exit(1)
    else: # assembling the files in a single dataframe:
        for c in paths:
            path = paths[c]
            # In case of separated files for each country, there is no header
            tmp = load_csv(path, index_col=0, parse_dates=True)
            # check that the loaded file is ok:
            if not tmp.index.is_unique:
                logging.error('The index of data file ' + paths['all'] + ' is not unique. Please check the data')
                sys.exit(1)
            data[c] = tmp.iloc[:,0]
     
    return data



def UnitBasedTable(plants,path,idx,countries,fallbacks=['Unit'],tablename='',default=None,RestrictWarning=None):
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
    :param path:                Path to the data to be loaded
    :param idx:                 Pandas datetime index to be used for the output
    :param countries:           List with the country codes to be considered
    :param fallback:            List with the order of data source. 
    :param tablename:           String with the name of the table being processed
    :param default:             Default value to be applied if no data is found
    :param RestrictWarning:     Only display the warnings if the unit belongs to the list of technologies provided in this parameter
    
    :return:           Dataframe with the time series for each unit
    '''
              
    paths = {}
    if os.path.isfile(path):
        paths['all'] = path
        SingleFile=True
    elif '##' in path:
        for c in countries:
            path_c = path.replace('##', str(c))
            if os.path.isfile(path_c):
                paths[str(c)] = path_c
            else:
                logging.critical('No data file found for the table ' + tablename + ' and country ' + c + '. File ' + path_c + ' does not exist')
#                sys.exit(1)
        SingleFile=False
    data = pd.DataFrame(index=idx)
    if len(paths) == 0:
        logging.info('No data file found for the table ' + tablename + '. Using default value ' + str(default))
        if default is None:
            out = pd.DataFrame(index=idx)
        elif isinstance(default,(float,int)):
            out = pd.DataFrame(default,index=idx,columns=plants['Unit'])
        else:
            logging.error('Default value provided for table ' + tablename + ' is not valid')
            sys.exit(1)
    else: # assembling the files in a single dataframe:
        columns = []
        for c in paths:
            path = paths[c]
            tmp = load_csv(path, index_col=0, parse_dates=True)
            # check that the loaded file is ok:
            if not tmp.index.is_unique:
                logging.error('The index of data file ' + path + ' is not unique. Please check the data')
                sys.exit(1)
            if SingleFile:
                for key in tmp:
                    data[key] = tmp[key]                
            else:    # use the multi-index header with the country
                for key in tmp:
                    columns.append((c,key))
                    data[c+','+key] = tmp[key]
        if not SingleFile:
            data.columns = pd.MultiIndex.from_tuples(columns, names=['Country', 'Data'])
        # For each plant and each fallback key, try to find the corresponding column in the data
        out = pd.DataFrame(index=idx)
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
                        logging.warning('No specific information was found for unit ' + u + ' in table ' + tablename + '. The generic information for ' + str(header) + ' has been used')
                    break
            if not found:
                if warning:
                    logging.info('No specific information was found for unit ' + u + ' in table ' + tablename + '. Using default value ' + str(default))
                if not default is None:
                    out[u] = default
    if not out.columns.is_unique:
        logging.error('The column headers of table "' + tablename + '" are not unique!. The following headers are duplicated: ' + str(out.columns.get_duplicates()))
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


def write_to_excel(xls_out, list_vars):
    """
    Function that reads all the variables (in list_vars) and inserts them one by one to excel

    :param xls_out: The path of the folder where the excel files are to be written
    :param list_vars: List containing the dispaset variables
    :returns: Binary variable (True)
    """


    reload_module(sys)
    try: # Hack needed in 2.7
        sys.setdefaultencoding("utf-8")
    except:
        pass

    if not os.path.exists(xls_out):
        os.mkdir(xls_out)


        # Printing all sets in one sheet:
    writer = pd.ExcelWriter(os.path.join(xls_out, 'InputDispa-SET - Sets.xlsx'), engine='xlsxwriter')

    [sets, parameters] = list_vars

    try:
        config = parameters['Config']['val']
        first_day = pd.datetime(config[0, 0], config[0, 1], config[0, 2], 0)
        last_day = pd.datetime(config[1, 0], config[1, 1], config[1, 2], 23)
        dates = pd.date_range(start=first_day, end=last_day, freq='1h')
    except:
        dates = []

    i = 0
    for s in sets:
        df = pd.DataFrame(sets[s], columns=[s])
        df.to_excel(writer, sheet_name='Sets', startrow=1, startcol=i, header=True, index=False)
        i += 1
    writer.save()
    logging.info('All sets successfully written to excel')

    # Printing each parameter in a separate sheet and workbook:
    for p in parameters:
        var = parameters[p]
        dim = len(var['sets'])
        if var['sets'][-1] == 'h' and isinstance(dates, pd.DatetimeIndex) and dim > 1:
            if len(dates) != var['val'].shape[-1]:
                logging.critical('The date range in the Config variable (' + str(
                    len(dates)) + ' time steps) does not match the length of the time index (' + str(
                    var['val'].shape[-1]) + ') for variable ' + p)
                sys.exit(1)
            var['firstrow'] = 5
        else:
            var['firstrow'] = 1
        writer = pd.ExcelWriter(os.path.join(xls_out, 'InputDispa-SET - ' + p + '.xlsx'), engine='xlsxwriter')
        if dim == 1:
            df = pd.DataFrame(var['val'], columns=[p], index=sets[var['sets'][0]])
            df.to_excel(writer, sheet_name=p, startrow=var['firstrow'], startcol=0, header=True, index=True)
            worksheet = writer.sheets[p]
            worksheet.write_string(0, 0, p + '(' + var['sets'][0] + ')')
            worksheet.set_column(0, 0, 30)
        elif dim == 2:
            list_sets = [sets[var['sets'][0]], sets[var['sets'][1]]]
            values = var['val']
            df = pd.DataFrame(values, columns=list_sets[1], index=list_sets[0])
            df.to_excel(writer, sheet_name=p, startrow=var['firstrow'], startcol=0, header=True, index=True)
            worksheet = writer.sheets[p]
            if var['firstrow'] == 5:
                worksheet.write_row(1, 1, dates.year)
                worksheet.write_row(2, 1, dates.month)
                worksheet.write_row(3, 1, dates.day)
                worksheet.write_row(4, 1, dates.hour + 1)
            worksheet.write_string(0, 0, p + '(' + var['sets'][0] + ',' + var['sets'][1] + ')')
            worksheet.freeze_panes(var['firstrow'] + 1, 1)
            worksheet.set_column(0, 0, 30)
        elif dim == 3:
            list_sets = [sets[var['sets'][0]], sets[var['sets'][1]], sets[var['sets'][2]]]
            values = var['val']
            for i in range(len(list_sets[0])):
                key = list_sets[0][i]
                Nrows = len(list_sets[1])
                df = pd.DataFrame(values[i, :, :], columns=list_sets[2], index=list_sets[1])
                df.to_excel(writer, sheet_name=p, startrow=var['firstrow'] + 1 + i * Nrows, startcol=1, header=False,
                            index=True)
                df2 = pd.DataFrame(np.array([key]).repeat(Nrows))
                df2.to_excel(writer, sheet_name=p, startrow=var['firstrow'] + 1 + i * Nrows, startcol=0, header=False,
                             index=False)
            worksheet = writer.sheets[p]
            if var['firstrow'] == 5:
                worksheet.write_row(1, 2, dates.year)
                worksheet.write_row(2, 2, dates.month)
                worksheet.write_row(3, 2, dates.day)
                worksheet.write_row(4, 2, dates.hour + 1)
            worksheet.write_string(0, 0, p + '(' + var['sets'][0] + ',' + var['sets'][1] + ',' + var['sets'][2] + ')')
            worksheet.write_string(var['firstrow'] - 1, 0, var['sets'][0])
            worksheet.write_string(var['firstrow'] - 1, 1, var['sets'][1])
            worksheet.freeze_panes(var['firstrow'], 2)
            worksheet.set_column(0, 1, 30)
            df = pd.DataFrame(columns=list_sets[2])
            df.to_excel(writer, sheet_name=p, startrow=var['firstrow'], startcol=2, header=True, index=False)
        else:
            logging.error('Only three dimensions currently supported. Parameter ' + p + ' has ' + str(dim) + ' dimensions.')
        writer.save()
        logging.info('Parameter ' + p + ' successfully written to excel')


    # Writing a gams file to process the excel sheets:
    gmsfile = open(os.path.join(xls_out, 'make_gdx.gms'), 'w')
    i = 0

    for s in sets:
        gmsfile.write('\n')
        gmsfile.write('$CALL GDXXRW "InputDispa-SET - Sets.xlsx" Set=' + s + ' rng=' + chr(
            i + ord('A')) + '3 Rdim=1  O=' + s + '.gdx \n')
        gmsfile.write('$GDXIN ' + s + '.gdx \n')
        gmsfile.write('Set ' + s + '; \n')
        gmsfile.write('$LOAD ' + s + '\n')
        gmsfile.write('$GDXIN \n')
        i = i + 1

    for p in parameters:
        var = parameters[p]
        dim = len(var['sets'])
        gmsfile.write('\n')
        if dim == 1:
            gmsfile.write('$CALL GDXXRW "InputDispa-SET - ' + p + '.xlsx" par=' + p + ' rng=A' + str(
                var['firstrow'] + 1) + ' Rdim=1 \n')
        elif dim == 2:
            gmsfile.write('$CALL GDXXRW "InputDispa-SET - ' + p + '.xlsx" par=' + p + ' rng=A' + str(
                var['firstrow'] + 1) + ' Rdim=1 Cdim=1 \n')
        elif dim == 3:
            gmsfile.write('$CALL GDXXRW "InputDispa-SET - ' + p + '.xlsx" par=' + p + ' rng=A' + str(
                var['firstrow'] + 1) + ' Rdim=2 Cdim=1 \n')
        gmsfile.write('$GDXIN "InputDispa-SET - ' + p + '.gdx" \n')
        gmsfile.write('Parameter ' + p + '; \n')
        gmsfile.write('$LOAD ' + p + '\n')
        gmsfile.write('$GDXIN \n')

    gmsfile.write('\n')
    gmsfile.write('Execute_Unload "Inputs.gdx"')
    gmsfile.close()

    logging.info('Data Successfully written to the ' + xls_out + ' directory.')


def load_csv(filename, TempPath='.pickle', header=0, skiprows=None, skipfooter=0, index_col=None, parse_dates=False):
    """
    Function that loads acsv sheet into a dataframe and saves a temporary pickle version of it.
    If the pickle is newer than the sheet, do no load the sheet again.

    :param filename: path to csv file
    :param TempPath: path to store the temporary data files
    """

    import hashlib
    m = hashlib.new('md5', filename.encode('utf-8'))
    resultfile_hash = m.hexdigest()
    filepath_pandas = TempPath + os.sep + resultfile_hash + '.p'
    
    if not os.path.isdir(TempPath):
        os.mkdir(TempPath)
    if not os.path.isfile(filepath_pandas):
        time_pd = 0
    else:
        time_pd = os.path.getmtime(filepath_pandas)
    if os.path.getmtime(filename) > time_pd:
        data = pd.read_csv(filename, header=header, skiprows=skiprows, skipfooter=skipfooter, index_col=index_col,
                           parse_dates=parse_dates)
        data.to_pickle(filepath_pandas)
    else:
        data = pd.read_pickle(filepath_pandas)
    return data


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

    try: 
        config["CEP"] = sheet.cell_value(151, 2)
    except: # config sheet without CEP cell leads to out of bounds error
        config["CEP"] = None


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

    # read the list of countries to consider:
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

    config['countries'] = read_truefalse(sheet, 86, 1, 101, 3)
    config['countries'] = config['countries'] + read_truefalse(sheet, 86, 4, 102, 6)

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