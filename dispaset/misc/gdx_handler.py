"""
Collection of functions to write Dispa-SET input data to a gdx file and/or to a simulation directory
with one excel file per parameter.

Example:
    read gdx file::

        data = GdxToList(gams_dir,'Results.gdx',varname='all',verbose=True)

    write it to a dictionary of dataframes::

        dataframes = GdxToDataframe(data,fixindex=True,verbose=True)

@author: Sylvain Quoilin (sylvain.quoilin@ec.europa.eu)

"""
import copy
import logging
import os
import platform
import sys
import time as tm

import numpy as np
import pandas as pd
from datetime import datetime, timedelta

from .str_handler import shrink_to_64, force_str


def package_exists(package_name):
    """
    Function that checks if a package is installed
    """
    try:
        __import__(package_name)
        return True
    except ImportError:
        return False


def import_local_lib(lib):
    """
    Try to import the GAMS api and gdxcc to write gdx files
    """
    # First define the path to the 'Externals' folder.
    # This path must be defined relatively to the current script location
    path_script = os.path.dirname(__file__)
    path_ext = os.path.join(path_script, '../../Externals')

    if sys.platform == 'win32' and platform.architecture()[0] == '64bit' and sys.version[:3] == '3.7':
        sys.path.append(os.path.join(path_ext, 'gams_api/win64/'))
    else:
        logging.error(
            'Pre-compiled GAMS libraries are only available for python 3.7 64 bits under windows. '
            'You are using platform ' + sys.platform + ' and architecture ' + platform.architecture()[0] +
            'Please install the gams API using: "pip install gamsxcc gdxcc optcc"')

    if lib == 'gams':
        try:
            import gams
            return True
        except ImportError:
            logging.error(
                'Could not load the gams high-level api. '
                'The gams library is required to run the GAMS versions of DispaSET.'
                'Please install the gams API using: "python setup.py install" in the gams api folder')
            sys.exit(1)
    elif lib == 'lowlevel':
        try:
            import gdxcc, gamsxcc, optcc
            return True
        except ImportError:
            logging.error(
                'Could not load the gams low-level api. '
                'The gams library is required to run the GAMS versions of DispaSET.'
                'Please install the gams API using: "pip install gamsxcc gdxcc optcc"')
            sys.exit(1)
    elif lib == 'gdxcc':
        try:
            import gdxcc
            return True
        except ImportError:
            logging.critical("gdxcc module could not be imported from Externals. GDX cannot be produced or read"
                             'Please install the gams API using: "pip install gamsxcc gdxcc optcc"')
            sys.exit(1)
    else:
        logging.error('Only "gams" and "gdxcc" are present')


if package_exists('gdxcc'):
    import gdxcc
else:
    logging.warning('Could not import gdxcc. Trying to use pre-compiled libraries')
    try:
        if sys.platform == 'win32' and platform.architecture()[0] == '64bit' and sys.version[:3] == '3.7':
            sys.path.append(os.path.join(path_ext, 'gams_api/win64/'))
        import gdxcc
    except ImportError:
        logging.critical('Importing gdxcc from the new gams api')
        try:
            import gams.core.gdx as gdxcc
        except ImportError:
            logging.critical("gdxcc module could not be imported from Externals. GDX cannot be produced or read"
                             'Please install the gams API using: "pip install gamsxcc gdxcc optcc"')
            sys.exit(1)    


#####################


def _insert_symbols(gdxHandle, sets, parameters):
    """
    Function that writes all sets and parameters to the gdxHandle

    :param sets: dictionary with all the sets
    :param parameters: dictionary with all the parameters
    """

    # It is essential to write the sets first, otherwise h might be written in the wrong order
    for s in sets:
        gdxSymbolType = gdxcc.GMS_DT_SET
        dims = 1

        gdxcc.gdxDataWriteStrStart(gdxHandle, s, "", dims, gdxSymbolType, 0)
        gdxValues = gdxcc.doubleArray(5)
        gdxValues[gdxcc.GMS_VAL_LEVEL] = 0.0  # 0.0 == Y (explanatory text of set in gdx)

        Nrows = len(sets[s])

        for row in range(Nrows):
            gdxKeys = [str(ss) for ss in shrink_to_64([sets[s][row]])]  # Reduce the size if bigger than 64 characters
            try:
                success = gdxcc.gdxDataWriteStr(gdxHandle, gdxKeys, gdxValues)
            except:
                success = False
            if not success:
                logging.error('Key ' + gdxKeys[0] + ' of set ' + s + ' could not be written')

        gdxcc.gdxDataWriteDone(gdxHandle)

    # Check array sizes for parameters:
    for p in parameters:
        variable = parameters[p]

        # Check that the required fields are present:
        dims = len(variable['sets'])
        shape = variable['val'].shape
        Nrows = variable['val'].shape[0]
        gdxSymbolType = gdxcc.GMS_DT_PAR

        gdxcc.gdxDataWriteStrStart(gdxHandle, p, "", dims, gdxSymbolType, 0)
        gdxValues = gdxcc.doubleArray(5)
        gdxValues[gdxcc.GMS_VAL_LEVEL] = 0.0  # 0.0 == Y (explanatory text of set in gdx)

        if len(shape) != dims:
            logging.error('Variable ' + p + ': The \'val\' data matrix has ' + str(
                len(shape)) + ' dimensions and should have ' + str(dims))
            sys.exit(1)
        for i in range(dims):
            if shape[i] != len(sets[variable['sets'][i]]):
                logging.error(
                    'Variable ' + p + ': The \'val\' data matrix has ' + str(shape[i]) + ' elements for dimention ' +
                    str(variable['sets'][i]) + ' while there are ' + str(
                        len(variable['sets'])) + ' set values')
                sys.exit(1)

        for index, value in np.ndenumerate(variable['val']):
            # Write line by line if value is non null
            if value != 0 and not pd.isnull(value):
                gdxKeys = []  # All the set values for this line
                for i in range(dims):
                    key = sets[variable['sets'][i]][
                        index[i]]  # Get the string value of the set by using the indice in the val matrix
                    gdxKeys.append(str(key))
                gdxKeys = shrink_to_64(gdxKeys)  # Reduce the size if bigger than 64 characters
                gdxValues[gdxcc.GMS_VAL_LEVEL] = float(value)
                try:
                    success = gdxcc.gdxDataWriteStr(gdxHandle, gdxKeys, gdxValues)
                except:
                    logging.error("Didn't work")
                    success = False
                if not success:
                    logging.error('Key ' + gdxKeys[0] + ' of parameter ' + p + ' could not be written')
        gdxcc.gdxDataWriteDone(gdxHandle)
        logging.debug('Parameter ' + p + ' successfully written')

    logging.debug('Set ' + s + ' successfully written')


def write_variables(config, gdx_out, list_vars):
    """
    This function performs the following:
    * Use the gdxcc library to create a gdxHandle instance
    * Check that the gams path is well defined
    * Call the 'insert_symbols' function to write all sets and parameters to gdxHandle

    :param config:          Main config dictionary
    :param gdx_out:         (Relative) path to the gdx file to be written
    :param list_vars:       List with the sets and parameters to be written
    """
    gams_dir = get_gams_path(config.get('GAMS_folder'))
    if not gams_dir:  # couldn't locate
        logging.critical('GDXCC: Could not find a valid GAMS installation. Please check your configuration and environment variables.')
        sys.exit(1)
    gams_dir = force_str(gams_dir)
    config['GAMS_folder'] = gams_dir  # updating the config dictionary
    gdx_out = force_str(gdx_out)

    gdxHandle = gdxcc.new_gdxHandle_tp()
    gdxcc.gdxCreateD(gdxHandle, gams_dir, gdxcc.GMS_SSSIZE)  # it accepts only str type
    gdxcc.gdxOpenWrite(gdxHandle, gdx_out, "")

    [sets, parameters] = list_vars
    _insert_symbols(gdxHandle, sets, parameters)

    gdxcc.gdxClose(gdxHandle)

    logging.info('Data Successfully written to ' + gdx_out)


def gdx_to_list(gams_dir, filename, varname='all', verbose=False):
    """original
    This function loads the gdx with the results of the simulation
    All results are stored in an unordered list

    :param gams_dir:    Gams working directory
    :param filename:    Path to the gdx file to be read
    :param varname:     In case online one variable is needed, specify it name (otherwise specify 'all')
    :returns:        Dictionary with all the collected values (within lists)
    """

    try:
        from gdxcc import gdxSymbolInfo, gdxCreate, gdxCreateD, gdxOpenRead, GMS_SSSIZE, gdxDataReadDone, \
            new_gdxHandle_tp, gdxDataReadStr, gdxFindSymbol, gdxErrorStr, gdxDataReadStrStart, gdxGetLastError
    except ImportError: 
        try:
            from gams.core.gdx import gdxSymbolInfo, gdxCreate, gdxCreateD, gdxOpenRead, GMS_SSSIZE, gdxDataReadDone, \
                new_gdxHandle_tp, gdxDataReadStr, gdxFindSymbol, gdxErrorStr, gdxDataReadStrStart, gdxGetLastError
        except ImportError:
            logging.critical("gdxcc module could not be imported. GDX cannot read"
                             'Please install the gams API"')
            sys.exit(1)


    out = {}
    tgdx = tm.time()
    gdxHandle = new_gdxHandle_tp()
    if gams_dir == None:
        gdxCreate(gdxHandle, GMS_SSSIZE)
    else:
        gdxCreateD(gdxHandle, force_str(gams_dir), GMS_SSSIZE)

    # make sure the file path is properly formatted:
    filename = filename.replace('/', os.path.sep).replace('\\\\', os.path.sep).replace('\\', os.path.sep)
    filename = str(filename)  # removing possible unicode formatting

    if not os.path.isfile(filename):
        logging.critical('Gdx file "' + filename + '" does not exist')
        sys.exit(1)

    gdxOpenRead(gdxHandle, filename)

    if varname == 'all':
        # go through all the symbols one by one and add their data to the dict
        symNr = 0
        SymbolInfo = gdxSymbolInfo(gdxHandle, 0)
        while SymbolInfo[0] > 0:
            ret, nrRecs = gdxDataReadStrStart(gdxHandle, symNr)
            assert ret, "Error in gdx data string" + gdxErrorStr(gdxHandle, gdxGetLastError(gdxHandle))[1]

            res = []
            for i in range(nrRecs):
                ret, elements, values, afdim = gdxDataReadStr(gdxHandle)
                res.append(elements + [values[0]])
            out[SymbolInfo[1]] = res
            symNr += 1
            SymbolInfo = gdxSymbolInfo(gdxHandle, symNr)
    else:
        # find the number of the required symbol:
        ret, symNr = gdxFindSymbol(gdxHandle, varname)
        assert ret, "Symbol not found"

        ret, nrRecs = gdxDataReadStrStart(gdxHandle, symNr)
        assert ret, "Error in gdx data string" + gdxErrorStr(gdxHandle, gdxGetLastError(gdxHandle))[1]

        res = []
        for i in range(nrRecs):
            ret, elements, values, afdim = gdxDataReadStr(gdxHandle)
            res.append(elements + [values[0]])
        out[varname] = res

    gdxDataReadDone(gdxHandle)
    if verbose:
        logging.info("Loading gdx file " + filename + " took {}s".format(tm.time() - tgdx))
    return out


def gdx_to_dataframe(data, fixindex=False, verbose=False, inputs=False):
    """
    This function structures the raw data extracted from a gdx file (using the function GdxToList)
    and outputs it as a dictionary of pandas dataframes (or series)

    :param data:        Dictionary with all the collected values (within lists), from GdxToList function
    :param fixindex:    This flag allows converting string index into integers and sort the data
    :returns:        dictionary of dataframes
    """
    out = {}
    tc = tm.time()
    for symbol in data:
        if len(data[symbol]) > 0:
            dim = len(data[symbol][0])
            if dim == 3:
                vars1 = set()
                for element in data[symbol]:
                    if not element[0] in vars1:
                        vars1.add(element[0])
                vals = {}
                while vars1:
                    vars2 = {}
                    var1 = vars1.pop()
                    for element in data[symbol]:
                        if var1 == element[0]:
                            vars2[element[1]] = element[2]
                    vals[var1] = vars2
                out[symbol] = pd.DataFrame(vals)
                logging.debug('Successfully loaded variable ' + symbol)
            elif dim == 2:
                vals = {}
                for element in data[symbol]:
                    vals[element[0]] = element[1]
                out[symbol] = pd.Series(vals)
                logging.debug('Successfully loaded variable ' + symbol)
            elif dim == 1:
                logging.warning('Variable ' + symbol + ' has dimension 0, which should not occur. Skipping')
            elif dim > 4:
                logging.warning('Variable ' + symbol + ' has more than 2 dimensions, which is very tiring. Skipping')
            elif dim == 4:
                vars1 = set()
                vars2 = set()
                for element in data[symbol]:
                    if not element[0] in vars1:
                        vars1.add(element[0])
                    if not element[1] in vars2:
                        vars2.add(element[1])
                vars1 = list(vars1)
                vars2 = list(vars2)
                pd_index = pd.MultiIndex.from_product([vars1, vars2], names=('Zones', 'Emissions'))
                if inputs:
                    out[symbol] = pd.DataFrame(columns=pd_index)
                else:
                    out[symbol] = pd.DataFrame(columns=pd_index, index=out['OutputPower'].index)
                for element1 in data[symbol]:
                    element = copy.deepcopy(element1)
                    for var1 in vars1:
                        for var2 in vars2:
                            if (var2 == element[1]) and (var1 == element[0]):
                                out[symbol].loc[element[2], (var1, var2)] = element[3]
                out[symbol].index = out[symbol].index.astype(int)
                out[symbol].sort_index(inplace=True)
                out[symbol] = out[symbol].fillna(0)

                logging.warning('Successfully loaded variable ' + symbol)
        else:
            logging.debug('Variable ' + symbol + ' is empty. Skipping')
    for symbol in out:
        try:
            out[symbol].fillna(value=0, inplace=True)
        except:
            logging.error('Error while trying to remove nan')
            pass
    if fixindex:
        for symbol in out:
            try:
                index_int = [int(idx) for idx in out[symbol].index]
                out[symbol].index = index_int
                out[symbol].sort_index(inplace=True)
            except:
                pass
    if verbose:
        logging.info("Time to convert to dataframes: {}s".format(tm.time() - tc))
    return out


def get_gdx(gams_dir, resultfile):
    """
    Short wrapper of the two gdx reading functions (GdxToDataframe and GdxToList)

    :param gams_dir:    Gams working directory
    :param resultfile:  Path to the gdx file to be read
    :returns:           dictionary of dataframes
    """
    return gdx_to_dataframe(gdx_to_list(gams_dir, resultfile,
                                        varname='all', verbose=True),
                            fixindex=True, verbose=True)


def get_gams_path(gams_folder=None):
    """
    Get the GAMS installation path using the following priority:
    1. User provided path (if valid)
    2. GAMSPATH environment variable
    3. GAMSDIR environment variable
    4. Try to locate gams executable in common paths
    5. On Linux, try to locate using the 'locate' command
    6. Prompt user for path

    :param gams_folder: Optional user-provided GAMS path
    :return: Path to GAMS installation or None if not found
    """
    # If a path is provided, try to use it first
    if gams_folder:
        if isinstance(gams_folder, bytes):
            gams_folder = gams_folder.decode()
        if os.path.exists(gams_folder):
            if platform.system() == 'Windows':
                if os.path.exists(os.path.join(gams_folder, 'gams.exe')):
                    return gams_folder
            else:
                if os.path.exists(os.path.join(gams_folder, 'gams')):
                    return gams_folder
        logging.warning(f'The provided GAMS path ({gams_folder}) is not valid. Trying to locate automatically...')

    # Try environment variables
    gams_path = os.environ.get('GAMSPATH')
    if gams_path and os.path.exists(gams_path):
        return gams_path

    gams_path = os.environ.get('GAMSDIR')
    if gams_path and os.path.exists(gams_path):
        return gams_path

    # Try to locate gams executable
    if platform.system() == 'Windows':
        # Common Windows paths
        common_paths = [
            r'C:\GAMS',
            r'C:\Program Files\GAMS',
            r'C:\Program Files (x86)\GAMS',
            os.path.expanduser('~\GAMS'),
        ]
        
        for base_path in common_paths:
            if os.path.exists(base_path):
                # Look for the most recent version
                versions = [d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d))]
                if versions:
                    latest_version = sorted(versions)[-1]
                    gams_path = os.path.join(base_path, latest_version)
                    if os.path.exists(os.path.join(gams_path, 'gams.exe')):
                        return gams_path
    else:
        # Linux/Mac paths
        common_paths = [
            '/opt/gams',
            '/usr/local/gams',
            os.path.expanduser('~/gams'),
        ]
        
        for base_path in common_paths:
            if os.path.exists(base_path):
                # Look for the most recent version
                versions = [d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d))]
                if versions:
                    latest_version = sorted(versions)[-1]
                    gams_path = os.path.join(base_path, latest_version)
                    if os.path.exists(os.path.join(gams_path, 'gams')):
                        return gams_path

        # On Linux, try to locate using the 'locate' command
        if platform.system() == 'Linux':
            try:
                from subprocess import check_output
                tmp = check_output(['locate', '-i', 'libgamscall64.so']).decode()
                lines = tmp.split('\n')
                for line in lines:
                    path = line.strip('libgamscall64.so')
                    if os.path.exists(path):
                        # Try to find the gams executable in the parent directory
                        parent_dir = os.path.dirname(path)
                        if os.path.exists(os.path.join(parent_dir, 'gams')):
                            return parent_dir
            except Exception as e:
                logging.debug(f"Could not use 'locate' command to find GAMS: {str(e)}")

    # If not found, prompt user
    print("\nGAMS installation not found. Please provide the path to your GAMS installation.")
    print("This should be the directory containing the GAMS executable (gams.exe on Windows, gams on Linux/Mac)")
    while True:
        gams_path = input("GAMS path: ").strip()
        if os.path.exists(gams_path):
            if platform.system() == 'Windows':
                if os.path.exists(os.path.join(gams_path, 'gams.exe')):
                    return gams_path
            else:
                if os.path.exists(os.path.join(gams_path, 'gams')):
                    return gams_path
        print("Invalid path. Please try again.")
