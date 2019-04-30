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
import platform
import os
import sys
import time as tm
import numpy as np
import pandas as pd
import logging

from .str_handler import shrink_to_64

def package_exists(package):
    # http://stackoverflow.com/questions/14050281/how-to-check-if-a-python-module-exists-without-importing-it
    import pkgutil
    package_loader = pkgutil.find_loader(package)
    return package_loader is not None


def import_local_lib(lib):
    '''
    Try to import the GAMS api and gdxcc to write gdx files
    '''
    # First define the path to the 'Externals' folder. This path must be defined relatively to the current script location
    path_script = os.path.dirname(__file__)
    path_ext = os.path.join(path_script,'../../Externals')

    if lib == 'gams':
        if sys.platform == 'linux2' and platform.architecture()[0] == '64bit':
            sys.path.append(os.path.join(path_ext,'gams_api/linux64/'))
        elif sys.platform == 'win32' and platform.architecture()[0] == '64bit':
            sys.path.append(os.path.join(path_ext,'gams_api/win64/'))
        #elif sys.platform == 'darwin':
        #    sys.path.append(os.path.join(path_ext, 'gams_api/linux64/'))
        try:
            import gams
            return True
        except ImportError:
            logging.error('Could not find gams. The gams library is required to run the GAMS versions of DispaSET.'
                          'Please install it from the /apifiles/Python/api/ folder in the GAMS directory')
            sys.exit(1)
    elif lib == 'gdxcc':
        if sys.platform == 'linux2' and platform.architecture()[0] == '32bit':
            sys.path.append(os.path.join(path_ext,'gdxcc/linux32/'))
        elif sys.platform == 'linux2' and platform.architecture()[0] == '64bit':
            sys.path.append(os.path.join(path_ext,'gams_api/linux64/'))
        elif sys.platform == 'win32' and platform.architecture()[0] == '32bit':
            sys.path.append(os.path.join(path_ext,'gdxcc/win32/'))
        elif sys.platform == 'win32' and platform.architecture()[0] == '64bit':
            sys.path.append(os.path.join(path_ext,'gams_api/win64/'))
        elif sys.platform == 'darwin':
            sys.path.append(os.path.join(path_ext,'gdxcc/osx64/'))
        try:
            import gdxcc
            return True
        except ImportError:
            print [x for x in sys.path]
            logging.critical("gdxcc module could not be found. GDX cannot be produced or read")
            sys.exit(1)
    else:
        logging.error('Only "gams" and "gdxcc" are present')

if package_exists('gdxcc'):
    import gdxcc
else:
    logging.warning('Could not import gdxcc. Trying to use pre-compiled libraries')
    if import_local_lib('gdxcc'):
        import gdxcc

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


def write_variables(gams_dir, gdx_out, list_vars):
    # type: (object, object, object) -> object
    """
    This function performs the following:
    * Use the gdxcc library to create a gdxHandle instance
    * Check that the gams path is well defined
    * Call the 'insert_symbols' function to write all sets and parameters to gdxHandle

    :rtype: object
    :param gams_dir:        (Relative) path to the gams directory
    :param gdx_out:         (Relative) path to the gdx file to be written
    :param list_vars:       List with the sets and parameters to be written
    """
    if not os.path.isdir(gams_dir):
        gams_dir = get_gams_path()
    if not os.path.isdir(gams_dir):
        logging.critical('GDXCC: Could not find the specified gams directory: ' + gams_dir)
        sys.exit(1)
    gams_dir = gams_dir.encode()

    gdxHandle = gdxcc.new_gdxHandle_tp()
    gdxcc.gdxCreateD(gdxHandle, gams_dir, gdxcc.GMS_SSSIZE)
    gdxcc.gdxOpenWrite(gdxHandle, gdx_out, "")

    [sets, parameters] = list_vars
    _insert_symbols(gdxHandle, sets, parameters)

    gdxcc.gdxClose(gdxHandle)

    logging.info('Data Successfully written to ' + gdx_out)



def gdx_to_list(gams_dir, filename, varname='all', verbose=False):
    """
    This function loads the gdx with the results of the simulation
    All results are stored in an unordered list

    :param gams_dir:    Gams working directory
    :param filename:    Path to the gdx file to be read
    :param varname:     In case online one variable is needed, specify it name (otherwise specify 'all')
    :returns:        Dictionary with all the collected values (within lists)
    """



    from gdxcc import gdxSymbolInfo, gdxCreateD, gdxOpenRead, GMS_SSSIZE, gdxDataReadDone, new_gdxHandle_tp, \
        gdxDataReadStr, gdxFindSymbol, gdxErrorStr, gdxDataReadStrStart, gdxGetLastError
    out = {}
    tgdx = tm.time()
    gdxHandle = new_gdxHandle_tp()
    gdxCreateD(gdxHandle, gams_dir, GMS_SSSIZE)

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


def gdx_to_dataframe(data, fixindex=False, verbose=False):
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
            #if dim != 4: continue
            if dim == 4:
                vals = {}
                for i, j, n, m in data[symbol]:
                    if (i,j) not in vals:
                        vals[(i,j)] = {}
                    vals[(i,j)][n] = m
                out[symbol] = pd.DataFrame(vals)
                logging.debug('Successfully loaded variable ' + symbol)
                """
                s1_unique, idx1 = np.unique(col1, return_inverse=True)
                s2_unique, idx2 = np.unique(col2, return_inverse=True)
                M = s2_unique.shape[0]
                idx_combined = idx1 * M + idx2
                print idx1.shape, idx2.shape, idx_combined.shape
                series = pd.DataFrame(data=[col4], index=col3, columns=idx_combined )
                """

                """
                vars1 = set(col1)
                vars2 = set(col2)
                N, M = len(vars1), len(vars2)
                series1, series2 = np.meshgrid(vars1, vars2)
                s1_unique, idx1 = np.unique(series1, return_inverse=True)
                s2_unique, idx2 = np.unique(series2, return_inverse=True)
                idx_combwined = idx1 * M + idx2
                """

                """
                for element in data[symbol]:
                    if not element[0] in vars1:
                        vars1.add(element[0])
                vars2 = set()
                for element in data[symbol]:
                    if not element[1] in vars2:
                        vars2.add(element[1])
                vars1 = list(vars1)
                series_1 = list(vars1 * len(vars2))
                series_2 = []
                for var2 in vars2:
                    series_2 = series_2 + [var2] * len(vars1)
                combined_array = [series_1, series_2]
                out[symbol] = pd.DataFrame(columns=combined_array)
                idx = pd.IndexSlice
                #series=pd.DataFrame(columns=combined_array)
                idx1, idx2, values = [], [], []
                for i, j, n, m in data[symbol]:
                    if i in vars1 and j in vars2:
                        idx1.append( idx[n] )
                        idx2.append(dx[i,j] i )
                        values.append( m )
                        #series.loc[idx[n], idx[i, j]] = m
                idx1, idx2, values = np.asarray(idx1), np.asarray(idx2), np.asarray(values)
                #series.loc[idx1, idx2] =  values
                print values.shape, idx1.shape, idx2.shape
                series = pd.DataFrame(data=[values], index=idx1, columns=idx2)
                """

                """
                for element in data[symbol]:
                    for i in vars1:
                        for j in vars2:
                            if i == element[0] and j== element[1]:
                                series.loc[idx[element[2]], idx[i, j]] = element[3]
                """
            elif dim == 3:
                vals = {}
                for i, n, m in data[symbol]:
                    if i not in vals:
                        vals[i] = {}
                    vals[i][n] = m
                out[symbol] = pd.DataFrame(vals)
                logging.debug('Successfully loaded variable ' + symbol)
                """
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
                """
            elif dim == 2:
                vals = {}
                for element in data[symbol]:
                    vals[element[0]] = element[1]
                out[symbol] = pd.Series(vals)
                logging.debug('Successfully loaded variable ' + symbol)
            elif dim == 1:
                logging.warn('Variable ' + symbol + ' has dimension 0, which should not occur. Skipping')
            elif dim > 4:
                logging.warn('Variable ' + symbol + ' has more than 3 dimensions, which is very tiring. Skipping')
        else:
             logging.debug('Variable ' + symbol + ' is empty. Skipping')
    for symbol in out:
        try:
            out[symbol].fillna(value=0, inplace=True)
        except:
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


def get_gams_path():
    """
    Function that attempts to search for the GAMS installation path (required to write the GDX or run gams)

    It returns the path if it has been found, or an empty string otherwise.

    Currently works for Windows, Linux and OSX. More searching rules and patterns should be added in the future
    """
    import subprocess
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
        paths = ['C:\\GAMS', 'C:\\Program Files\\GAMS', 'C:\\Program Files (x86)\\GAMS']
        lines_32 = []
        lines_64 = []
        for path in paths:
            if os.path.exists(path):
                paths_32 = [path + os.sep + tmp for tmp in os.listdir(path) if
                          tmp.startswith('win32') and os.path.exists(path + os.sep + tmp)]
                paths_64 = [path + os.sep + tmp for tmp in os.listdir(path) if
                          tmp.startswith('win64') and os.path.exists(path + os.sep + tmp)]
                for path1 in paths_32:
                    lines_32 = lines_32 + [path1 + os.sep + tmp for tmp in os.listdir(path1) if
                                     tmp.startswith('24') and os.path.isfile(
                                         path1 + os.sep + tmp + os.sep + 'gams.exe')]
                for path1 in paths_64:
                    lines_64 = lines_64 + [path1 + os.sep + tmp for tmp in os.listdir(path1) if
                                     tmp.startswith('24') and os.path.isfile(
                                         path1 + os.sep + tmp + os.sep + 'gams.exe')]
        for line in lines_64:
            if os.path.exists(line):
                out = line
            break
        if out == '':    # The 32-bit version of gams should never be preferred
            for line in lines_32:
                if os.path.exists(line):
                    out = line
                    logging.critical('It seems that the installed version of gams is 32-bit, which might cause consol crashing and compatibility problems. Please consider using GAMS 64-bit')
                break
    elif sys.platform == 'darwin':
        paths = ['/Applications/']
        lines = []
        for path in paths:
            if os.path.exists(path):
                paths1 = [path + os.sep + tmp for tmp in os.listdir(path) if
                          tmp.startswith('GAMS') and os.path.exists(path + os.sep + tmp)]
                for path1 in paths1:
                    lines = lines + [path1 + os.sep + tmp for tmp in os.listdir(path1) if
                                     tmp.startswith('sysdir') and os.path.isfile(
                                         path1 + os.sep + tmp + os.sep + 'gams')]
        if len(lines) == 0:
            proc = subprocess.Popen(['mdfind', '-name', 'libgamscall64.dylib'], stdout=subprocess.PIPE)
            tmp = proc.stdout.read()
            lines = [x.strip('libgamscall64.dylib') for x in tmp.split('\n')]
        for line in lines:
            if os.path.exists(line):
                out = line
            break

    if out != '':
        logging.info('Detected ' + out + ' as GAMS path on this computer')
    else:
        tmp = input('Specify the path to GAMS within quotes (e.g. "C:\\\\GAMS\\\\win64\\\\24.3"): ')
        if os.path.isdir(tmp):
            if sys.platform == 'win32':
                if os.path.isfile(tmp + os.sep + 'gams.exe'):
                    out = tmp
                else:
                    logging.critical('The provided path is not a valid windows gams folder')
                    sys.exit(1)
            elif sys.platform == 'linux2':
                if os.path.isfile(tmp + os.sep + 'gamslib'):  # does not always work... gamslib_ml
                    out = tmp
                else:
                    logging.critical('The provided path is not a valid linux gams folder')
                    sys.exit(1)
            else:
                if os.path.isdir(tmp):
                    out = tmp

    return out
