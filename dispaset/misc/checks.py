# -*- coding: utf-8 -*-
"""
Set of functions used to test the availability of different libraries required by Dispa-SET

@author: Sylvain Quoilin
"""
import platform
import sys
import os
import subprocess


def get_gams_path():
    """
    Function that attempts to search for the GAMS installation path (required to write the GDX or run gams)

    It returns the path if it has been found, or an empty string otherwise.

    Currently works for Windows, Linux and OSX. More searching rules and patterns should be added in the future
    """
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
                    print('It seems that the installed version of gams is 32-bit, which might cause consol crashing and compatibility problems. Please consider using GAMS 64-bit')
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
    return out



def cmd_exists(cmd):
    return subprocess.call("type " + cmd, shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0




# ## First check logging:

print('CHECK LOGGING')
import logging
import logging.config
# select ColorStreamHandler based on platform

_LOGCONFIG = {
     "version": 1,
     "disable_existing_loggers": False,
     'formatters': {
        'standard': {
            'format': '%(asctime)s [%(levelname)-8s] (%(funcName)s): %(message)s',
            'datefmt': '%y/%m/%d %H:%M:%S'
        },
        'notime': {
            'format': '[%(levelname)-8s] (%(funcName)s): %(message)s',
            'datefmt': '%y/%m/%d %H:%M:%S'
        },
     },
     "handlers": {
         "console": {
            "class": "logging.StreamHandler",
             "stream": "ext://sys.stderr",
#             "stream": "sys.stdout",
             "level": "INFO",
             'formatter': 'notime',
         },

         "error_file": {
             "class": "logging.FileHandler",
             "level": "WARNING",
             'formatter': 'standard',
             'filename': 'warn.log',
             'encoding': 'utf8'

         }
     },

     "root": {
         "level": "INFO",
         "handlers": ["console", "error_file"],
     }
}

try:
    logging.config.dictConfig(_LOGCONFIG)
        # Display a few messages:
    logging.warn('This is a warning')
    logging.info('This is information')
    logging.critical('This is critical')
    logging.error('This is an error')
    print('Logging seems to be fine')
except Exception as e:
    print('Logging cannot be handled. The following error msg was generated:')
    print(e)






print('\n \nCHECK THE GDXCC LIBRARY')
success_gdxcc = True
try:
    import gdxcc
except ImportError as e:
    if str(e) == 'No module named gdxcc':
        print('Could not find the gdxcc library in the standard python PATH. Trying to import the pre-compiled libraries')
        path_script = os.path.dirname(__file__)
        path_ext = os.path.join(path_script,'../../Externals')
        if sys.platform == 'linux2' and platform.architecture()[0] == '32bit':
            path = os.path.join(path_ext,'gdxcc/linux32/')
            print('Adding the following folder to the system PATH: ' + path)
            sys.path.append(path)
        elif sys.platform == 'linux2' and platform.architecture()[0] == '64bit':
            path = os.path.join(path_ext,'gams_api/linux64/')
            print('Adding the following folder to the system PATH: ' + path)
            sys.path.append(path)
        elif sys.platform == 'win32' and platform.architecture()[0] == '32bit':
            path = os.path.join(path_ext,'gdxcc/win32/')
            print('Adding the following folder to the system PATH: ' + path)
            sys.path.append(path)
        elif sys.platform == 'win32' and platform.architecture()[0] == '64bit':
            path = os.path.join(path_ext,'gams_api/win64/')
            print('Adding the following folder to the system PATH: ' + path)
            sys.path.append(path)
        elif sys.platform == 'darwin':
            path = os.path.join(path_ext,'gdxcc/osx64/')
            print('Adding the following folder to the system PATH: ' + path)
            sys.path.append(path)
        try:
            import gdxcc
        except ImportError as ee:
            print('ERROR: Could not load the precompiled gdxcc library. The following error was issued: ' + str(ee))
            success_gdxcc = False
    else:
        print('ERROR: ' + str(e))
        success_gdxcc = False
if success_gdxcc:
    print('gdxcc library successfully loaded at the following location: ' + gdxcc.__file__)



# ## Check the gams library:

print('\n \nCHECK GAMS library')
success_lib = True
success_path = True
success_sim = True
try:
    import gams
except Exception as e:
    if str(e) == 'No module named gams':
        print('Could not find the gams library in the standard python PATH. Trying to import the pre-compiled libraries')
        path_script = os.path.dirname(__file__)
        path_ext = os.path.join(path_script,'../../Externals')
    
        if sys.platform == 'linux2' and platform.architecture()[0] == '64bit':
            path = os.path.join(path_ext,'gams_api/linux64/')
            print('Adding the following folder to the system PATH: ' + path)
            sys.path.append(path)
        elif sys.platform == 'win32' and platform.architecture()[0] == '64bit':
            path = os.path.join(path_ext,'gams_api/win64/')
            print('Adding the following folder to the system PATH: ' + path)
            sys.path.append(path)
        try:
            import gams
        except ImportError as ee:
            print('ERROR: Could not load the precompiled gams library. The following error was issued: ' + str(ee))
            success_lib = False
    else:
        print('ERROR: ' + str(e))
        success_lib = False
if success_lib:
    print('GAMS library successfully loaded at the following location: ' + gams.__file__)


print('\n \nCHECK GAMS INSTALLATION FOLDER')

gamspath = get_gams_path()
if gamspath == '':
    print('ERROR: The GAMS installation folder could not be found')
    success_path = False
else:
    print('GAMS folder found: ' + gamspath)





if success_path and success_lib:
    print('\n \nTRY TO RUN A SIMPLE GAMS MODEL:')
    try:
        ws = gams.GamsWorkspace(system_directory=gamspath,debug=1)
        ws.gamslib("trnsport")
        t1 = ws.add_job_from_file("trnsport.gms")
        t1.run()
        for rec in t1.out_db["x"]:
            print ("x(" + rec.keys[0] + "," + rec.keys[1] + "): level=" + str(rec.level) + " marginal=" + str(rec.marginal))
        print('The optimization seems to have run properly')
    except Exception as e:
        print('ERROR while trying to run the optimization: ' + str(e))
        success_sim = False







if success_gdxcc and success_path:
    print('\n \nTRY TO GENERATE GDX FILE')
    try:
        gdxHandle = gdxcc.new_gdxHandle_tp()
        gdxcc.gdxCreateD(gdxHandle, gamspath, gdxcc.GMS_SSSIZE)
        gdxcc.gdxOpenWrite(gdxHandle, 'test.gdx', "")
    
        # Write a set:
        dims = 1
        setname = 'set_test'
        keys = ['aa']
        gdxcc.gdxDataWriteStrStart(gdxHandle, setname, "", dims, gdxcc.GMS_DT_SET, 0)
        gdxValues = gdxcc.doubleArray(5)
        gdxValues[gdxcc.GMS_VAL_LEVEL] = 0.0  # 0.0 == Y (explanatory text of set in gdx)
    
        try:
            success = gdxcc.gdxDataWriteStr(gdxHandle, keys, gdxValues)
        except Exception as e:
            success = False
            print('ERROR: the set could not be written to the gdx file. Error msg: ' + str(e))
        gdxcc.gdxDataWriteDone(gdxHandle)
        gdxcc.gdxClose(gdxHandle)
    except Exception as ee:
        print('ERROR: the gdxfile could not be created. Error msg: ' + str(ee))
        success = False
    if success and os.path.isfile('test.gdx'):
        print('GDX successfully written. Cleaning up')
        os.remove('test.gdx')




print('\n \nCHECK PYOMO')
try:
    import pyomo
    print('Pyomo is available')
except ImportError as e:
    print('ERROR: Pyomo could not be load. Erros msg: ' + str(e))


print('\n \nCHECK CPLEX')
if cmd_exists('cplex'):
    print('cplex is available from the command prompt!')
else:
    print('ERROR: the cplex command is not available. It should be in a location referenced in the PATH environment variable')

