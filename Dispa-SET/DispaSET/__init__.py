import sys
import platform


# Try to import the GAMS api and gdxcc to write gdx files:
try:
    import gams
    import gdxcc
    gams_ok = True
    gdxcc_ok = True
except:
    if sys.platform == 'linux2' and platform.architecture()[0] == '64bit':
        sys.path.append('../Externals/gams_api/linux64/')
    elif sys.platform == 'linux2' and platform.architecture()[0] == '32bit':
        sys.path.append('../Externals/gdxcc/linux32/')
    elif sys.platform == 'win32' and platform.architecture()[0] == '64bit':
        sys.path.append('../Externals/gams_api/win64/')
    elif sys.platform == 'win32' and platform.architecture()[0] == '32bit':
        sys.path.append('../Externals/gdxcc/win32/')
    elif sys.platform == 'darwin':
        sys.path.append('../Externals/gdxcc/osx64/')        
    try:
        import gams
        gams_ok = True
    except:
        gams_ok = False     
    try:
        import gdxcc
        gdxcc_ok = True
    except:
        gdxcc_ok = False
     
       
#from .DispaCheck import *
from .DispaSET_io_data import *
#from .DispaSolve import *      # importing dispasolve causes errors if pyomo is not installed
from .DispaTools import *
from .GenerateInputs import *
from .PostProcessing import *

