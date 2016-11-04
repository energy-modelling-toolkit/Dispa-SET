'''
Collection of functions to write Dispa-SET input data to a gdx file and/or to a simulation directory
with one excel file per parameter.

Example: 
    read gdx file::
    
        data = GdxToList(gams_dir,'Results.gdx',varname='all',verbose=True)

    write it to a dictionary of dataframes::
    
        dataframes = GdxToDataframe(data,fixindex=True,verbose=True)

@author: Sylvain Quoilin (sylvain.quoilin@ec.europa.eu)

'''

from DispaTools import shrink_to_64
import numpy as np
import os
import sys
import pandas as pd
import time as tm
try:
    import gdxcc
    gdxcc_available = True
except ImportError:
    gdxcc_available = False



def LoadConfig(ConfigFile):
    '''
    Function that loads the DispaSET excel config file and returns a dictionnary
    with the values
    
    :param ConfigFile: String with (relative) path to the DispaSET excel configuration file
    '''
    import xlrd
    wb = xlrd.open_workbook(filename=ConfigFile)   # Option for csv to be added later
    sheet = wb.sheet_by_name('main')
    
    config = {}
    config['SimulationDirectory'] = sheet.cell_value(17,2)
    config['WriteExcel'] = sheet.cell_value(18,2)
    config['WriteGDX'] = sheet.cell_value(19,2)
    config['WritePickle'] = sheet.cell_value(20,2)
    config['GAMS_folder']= sheet.cell_value(21,2)
    
    config['StartDate'] = xlrd.xldate_as_tuple(sheet.cell_value(30,2),wb.datemode)
    config['StopDate'] = xlrd.xldate_as_tuple(sheet.cell_value(31,2),wb.datemode)
    config['HorizonLength'] = int(sheet.cell_value(32,2))
    config['LookAhead'] = int(sheet.cell_value(33,2))
    
    config['Clustering'] = sheet.cell_value(45,2)
    config['SimulationType'] = sheet.cell_value(46,2)
    config['ReserveCalculation'] = sheet.cell_value(47,2)
    config['AllowCurtailment'] = sheet.cell_value(48,2)
    
    params = ['Demand', 'Outages', 'PowerPlantData','RenewablesAF','LoadShedding','NTC','Interconnections','ReservoirScaledInflows','PriceOfNuclear','PriceOfBlackCoal','PriceOfGas','PriceOfFuelOil','PriceOfBiomass','PriceOfCO2','ReservoirLevels']
    for i,param in enumerate(params):
        config[param] = sheet.cell_value(61+i,2)
        
    config['default'] ={}
    config['default']['PriceOfNuclear'] = sheet.cell_value(69,5)
    config['default']['PriceOfBlackCoal'] = sheet.cell_value(70,5)
    config['default']['PriceOfGas'] = sheet.cell_value(71,5)
    config['default']['PriceOfFuelOil'] = sheet.cell_value(72,5)
    config['default']['PriceOfBiomass'] = sheet.cell_value(73,5)
    config['default']['PriceOfCO2'] = sheet.cell_value(74,5)
    config['default']['LoadShedding'] = sheet.cell_value(65,5)
    
    # read the list of countries to consider:
    def read_truefalse(sheet,rowstart,colstart,rowstop,colstop):
        ''' 
        Function that reads a two column format with a list of strings in the first
        columns and a list of true false in the second column
        The list of strings associated with a True value is returned
        '''
        out = []
        for i in range(rowstart,rowstop):
            if sheet.cell_value(i,colstart+1) == 1:
                out.append(sheet.cell_value(i,colstart))
        return out
                
    config['countries'] = read_truefalse(sheet,86,1,101,3)
    config['countries'] = config['countries'] + read_truefalse(sheet,86,4,101,6)
    
    config['modifiers']= {}
    config['modifiers']['Demand'] = sheet.cell_value(111,2)
    config['modifiers']['Wind'] = sheet.cell_value(112,2)
    config['modifiers']['Solar'] = sheet.cell_value(113,2)
    config['modifiers']['Storage'] = sheet.cell_value(114,2)
    
    # Read the technologies participating to reserve markets:
    config['ReserveParticipation'] = read_truefalse(sheet,131,1,145,3)
    
    return config


def insert_symbols(gdxHandle,sets,parameters):

    '''
    Function that writes all sets and parameters to the gdxHandle
    
    :param sets: Dictionnary with all the sets
    :param parameters: Dictionnary with all the parameters
    '''

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
        gdxValues[gdxcc.GMS_VAL_LEVEL] = 0.0 #0.0 == Y (explanatory text of set in gdx)        
        
        if len(shape) != dims:
            sys.exit('Variable ' + p + ': The \'val\' data matrix has ' + str(len(shape)) + ' dimensions and should have ' + str(dims))
        for i in range(dims):
            if shape[i] != len(sets[variable['sets'][i]]):
                sys.exit('Variable ' + p + ': The \'val\' data matrix has ' + str(shape[i]) + ' elements for dimention ' + sets[variable['sets'][i]] + ' while there are ' + str(len(variable['sets'][i]['uels'])) + ' set values')

        for index,value in np.ndenumerate(variable['val']):
            # Write line by line if value is non null
            if value != 0 and not np.isnan(value):
                gdxKeys = []                                                        # All the set values for this line
                for i in range(dims):
                    key=sets[variable['sets'][i]][index[i]]                       # Get the string value of the set by using the indice in the val matrix
                    gdxKeys.append(str(key))
                gdxKeys = shrink_to_64(gdxKeys)                                     # Reduce the size if bigger than 64 characters
                gdxValues[gdxcc.GMS_VAL_LEVEL] = float(value)
                try:
                    success = gdxcc.gdxDataWriteStr(gdxHandle, gdxKeys, gdxValues)
                except:
                    print "Didn't work"
                    success= False
                if not success:
                    print('Key ' + gdxKeys[0] + ' of parameter ' + p + ' could not be written')
        gdxcc.gdxDataWriteDone(gdxHandle)
        print('Parameter ' + p + ' successfully written')


    for s in sets:
        gdxSymbolType = gdxcc.GMS_DT_SET
        dims = 1
        
        gdxcc.gdxDataWriteStrStart(gdxHandle, s, "", dims, gdxSymbolType, 0)
        gdxValues = gdxcc.doubleArray(5)
        gdxValues[gdxcc.GMS_VAL_LEVEL] = 0.0 #0.0 == Y (explanatory text of set in gdx)     
        
        Nrows = len(sets[s])

        for row in range(Nrows):
            gdxKeys = [str(ss) for ss in shrink_to_64([sets[s][row]])]                            # Reduce the size if bigger than 64 characters
            try:
                success = gdxcc.gdxDataWriteStr(gdxHandle, gdxKeys, gdxValues)
            except:
                success = False
            if not success:
                print('Key ' + gdxKeys[0] + ' of set ' + s + ' could not be written')

        gdxcc.gdxDataWriteDone(gdxHandle)
        print('Set ' + s + ' successfully written')

 
def write_variables(gams_dir,gdx_out,list_vars):
    '''
    This function performs the following:
    * Use the gdxcc library to create a gdxHandle instance
    * Check that the gams path is well defined
    * Call the 'insert_symbols' function to write all sets and parameters to gdxHandle
    
    :param gams_dir:        (Relative) path to the gams directory
    :param gdx_out:         (Relative) path to the gdx file to be written
    :param list_vars:       List with the sets and parameters to be written    
    '''
    if gdxcc_available:
        #Check gams_dir:
        if not os.path.isdir(gams_dir):
            sys.exit('Could not find the specified gams directory: ' + gams_dir)

        gdxHandle = gdxcc.new_gdxHandle_tp()
        gdxcc.gdxCreateD(gdxHandle, gams_dir, gdxcc.GMS_SSSIZE)
        gdxcc.gdxOpenWrite(gdxHandle, gdx_out, "")
            

        [sets, parameters] = list_vars
        insert_symbols(gdxHandle,sets,parameters)

        gdxcc.gdxClose(gdxHandle)

        print 'Data Successfully written to ' + gdx_out
    else: 
        sys.exit("gdxcc module not available. GDX cannot be produced")



def write_toexcel(xls_out,list_vars):
    '''
    Function that reads all the variables (in list_vars) and inserts them one by one to excel
    
    :param xls_out: The path of the folder where the excel files are to be written
    :param list_vars: List containing the dispaset variables 
    :returns: Binary variable (True)
    '''
    
    import pandas as pd
    
    #import sys
    reload(sys)
    sys.setdefaultencoding("utf-8")
    
    if not os.path.exists(xls_out):
        os.mkdir(xls_out)    
     
    
    # Printing all sets in one sheet:    
    writer = pd.ExcelWriter(os.path.join(xls_out,'InputDispa-SET - Sets.xlsx'), engine='xlsxwriter')

    [sets, parameters] = list_vars
    
    try:
        config = parameters['Config']['val']
        first_day = pd.datetime(config[0,0],config[0,1],config[0,2],0)
        last_day = pd.datetime(config[1,0],config[1,1],config[1,2],23)
        dates = pd.date_range(start=first_day,end=last_day,freq='1h')
    except:
        dates=[]
    
    i = 0
    for s in sets:
        df = pd.DataFrame(sets[s],columns=[s])
        df.to_excel(writer, sheet_name='Sets',startrow=1,startcol=i,header=True,index=False)
        i += 1
    writer.save()    
    print 'All sets successfully written to excel'
    
    # Printing each parameter in a separate sheet and workbook:
    for p in parameters:
        var = parameters[p]
        dim = len(var['sets'])
        if var['sets'][-1] =='h' and isinstance(dates,pd.tseries.index.DatetimeIndex) and dim > 1:
            if len(dates) != var['val'].shape[-1]:
                sys.exit('The date range in the Config variable (' + str(len(dates)) + ' time steps) does not match the length of the time index (' + str(var['val'].shape[-1]) + ') for variable ' + p)
            var['firstrow'] = 5
        else:
            var['firstrow']=1
        writer = pd.ExcelWriter(os.path.join(xls_out,'InputDispa-SET - ' + p + '.xlsx'), engine='xlsxwriter')
        if dim == 1:
            df = pd.DataFrame(var['val'],columns=[p],index=sets[var['sets'][0]])
            df.to_excel(writer, sheet_name=p,startrow=var['firstrow'],startcol=0,header=True,index=True)
            worksheet = writer.sheets[p]
            worksheet.write_string(0, 0, p+'(' + var['sets'][0] + ')' )
            worksheet.set_column(0,0,30)
        elif dim ==2:    
            list_sets = [sets[var['sets'][0]], sets[var['sets'][1]]]
            values = var['val']
            df = pd.DataFrame(values,columns=list_sets[1],index=list_sets[0])
            df.to_excel(writer, sheet_name=p,startrow=var['firstrow'],startcol=0,header=True,index=True)
            worksheet = writer.sheets[p]
            if var['firstrow'] == 5:
                worksheet.write_row(1,1,dates.year)
                worksheet.write_row(2,1,dates.month)
                worksheet.write_row(3,1,dates.day)
                worksheet.write_row(4,1,dates.hour+1)
            worksheet.write_string(0, 0, p+'(' + var['sets'][0] + ',' + var['sets'][1] + ')'  )
            worksheet.freeze_panes(var['firstrow']+1,1)
            worksheet.set_column(0,0,30)
        elif dim ==3:    
            list_sets = [sets[var['sets'][0]], sets[var['sets'][1]], sets[var['sets'][2]]]
            values = var['val']
            for i in range(len(list_sets[0])):
                key = list_sets[0][i]
                Nrows = len(list_sets[1])
                df = pd.DataFrame(values[i,:,:],columns=list_sets[2],index=list_sets[1])
                df.to_excel(writer, sheet_name=p,startrow=var['firstrow']+1+i*Nrows,startcol=1,header=False,index=True)
                df2 = pd.DataFrame(np.array([key]).repeat(Nrows))
                df2.to_excel(writer, sheet_name=p,startrow=var['firstrow']+1+i*Nrows,startcol=0,header=False,index=False)
            worksheet = writer.sheets[p]
            if var['firstrow'] == 5:
                worksheet.write_row(1,2,dates.year)
                worksheet.write_row(2,2,dates.month)
                worksheet.write_row(3,2,dates.day)
                worksheet.write_row(4,2,dates.hour+1)                
            worksheet.write_string(0, 0, p+'(' + var['sets'][0] + ',' + var['sets'][1] + ',' + var['sets'][2] + ')'  )
            worksheet.write_string(var['firstrow']-1, 0, var['sets'][0])
            worksheet.write_string(var['firstrow']-1, 1, var['sets'][1])
            worksheet.freeze_panes(var['firstrow'],2)
            worksheet.set_column(0,1,30)
            df = pd.DataFrame(columns=list_sets[2])
            df.to_excel(writer, sheet_name=p,startrow=var['firstrow'],startcol=2,header=True,index=False)
        else:
            print 'Only three dimensions currently supported. Parameter ' + p + ' has ' + str(dim) + ' dimensions.'
        writer.save() 
        print 'Parameter ' + p + ' successfully written to excel'        


    # Writing a gams file to process the excel sheets:
    gmsfile = open(os.path.join(xls_out,'make_gdx.gms'),'w')
    i= 0
        
    for s in sets:   
        gmsfile.write('\n')
        gmsfile.write('$CALL GDXXRW "InputDispa-SET - Sets.xlsx" Set=' + s + ' rng=' + chr(i + ord('A')) + '3 Rdim=1  O=' + s + '.gdx \n')
        gmsfile.write('$GDXIN ' + s + '.gdx \n')
        gmsfile.write('Set ' + s + '; \n')
        gmsfile.write('$LOAD ' + s + '\n') 
        gmsfile.write('$GDXIN \n')
        i = i+1
            
    for p in parameters:
        var = parameters[p]
        dim = len(var['sets'])
        gmsfile.write('\n')
        if dim ==1:
            gmsfile.write('$CALL GDXXRW "InputDispa-SET - ' + p +'.xlsx" par=' + p + ' rng=A' + str(var['firstrow']+1) + ' Rdim=1 \n')
        elif dim ==2:
            gmsfile.write('$CALL GDXXRW "InputDispa-SET - ' + p +'.xlsx" par=' + p + ' rng=A' + str(var['firstrow']+1) + ' Rdim=1 Cdim=1 \n')
        elif dim ==3:
            gmsfile.write('$CALL GDXXRW "InputDispa-SET - ' + p +'.xlsx" par=' + p + ' rng=A' + str(var['firstrow']+1) + ' Rdim=2 Cdim=1 \n')
        gmsfile.write('$GDXIN "InputDispa-SET - ' + p + '.gdx" \n')
        gmsfile.write('Parameter ' + p + '; \n')
        gmsfile.write('$LOAD ' + p + '\n') 
        gmsfile.write('$GDXIN \n')

    gmsfile.write('\n')
    gmsfile.write('Execute_Unload "Inputs.gdx"')             
    gmsfile.close()

    print 'Data Successfully written to the ' + xls_out + ' directory.'


def GdxToList(gams_dir,filename,varname='all',verbose=False):
    '''
    This function loads the gdx with the results of the simulation 
    All results are stored in an unordered list

    :param gams_dir:    Gams working directory
    :param filename:    Path to the gdx file to be read
    :param varname:     In case online one variable is needed, specify it name (otherwise specify 'all')
    :returns:        Dictionary with all the collected values (within lists)
    '''
    
    if not gdxcc_available:
        sys.exit("gdxcc module not available. GDX cannot be read")    
    
    from gdxcc import gdxSymbolInfo,gdxCreateD,gdxOpenRead,GMS_SSSIZE,gdxDataReadDone,new_gdxHandle_tp,gdxDataReadStr,gdxFindSymbol,gdxErrorStr,gdxDataReadStrStart,gdxGetLastError
    out = {}
    tgdx = tm.time()
    gdxHandle = new_gdxHandle_tp()
    gdxCreateD(gdxHandle, gams_dir, GMS_SSSIZE)

    # make sure the file path is properly formatted:
    filename = filename.replace('/',os.path.sep).replace('\\\\',os.path.sep).replace('\\',os.path.sep)
    filename = str(filename)            # removing possible unicode formatting    
    
    if not os.path.isfile(filename):
        sys.exit('Gdx file "' + filename + '" does not exist')
    
    gdxOpenRead(gdxHandle, filename)

    if varname == 'all':
        # go through all the symbols one by one and add their data to the dict
        symNr = 0
        SymbolInfo = gdxSymbolInfo(gdxHandle,0)
        while SymbolInfo[0] > 0:
            ret, nrRecs = gdxDataReadStrStart(gdxHandle,symNr)
            assert ret,"Error in gdx data string" + gdxErrorStr(gdxHandle,gdxGetLastError(gdxHandle))[1]
            
            res = []
            for i in range(nrRecs):
                ret,elements,values,afdim = gdxDataReadStr(gdxHandle)
                res.append(elements + [values[0]])  
            out[SymbolInfo[1]] = res
            symNr += 1
            SymbolInfo = gdxSymbolInfo(gdxHandle,symNr)
    else:
        # find the number of the required symbol:
        ret, symNr = gdxFindSymbol(gdxHandle,varname)
        assert ret,"Symbol not found"
        
        ret, nrRecs = gdxDataReadStrStart(gdxHandle,symNr)
        assert ret,"Error in gdx data string" + gdxErrorStr(gdxHandle,gdxGetLastError(gdxHandle))[1]
        
        res = []
        for i in range(nrRecs):
            ret,elements,values,afdim = gdxDataReadStr(gdxHandle)
            res.append(elements + [values[0]])
        out[varname] = res
        
    gdxDataReadDone(gdxHandle)
    if verbose:
        print("Loading gdx file " + filename + " took {}s".format(tm.time()-tgdx))  
    return out

def GdxToDataframe(data,fixindex=False,verbose=False):
    '''
    This function structures the raw data extracted from a gdx file (using the function GdxToList)
    and outputs it as a dictionnary of pandas dataframes (or series)

    :param data:        Dictionary with all the collected values (within lists), from GdxToList function
    :param fixindex:    This flag allows converting string index into integers and sort the data
    :returns:        Dictionnary of dataframes
    '''
    out = {}
    tc = tm.time()
    for symbol in data:
        if len(data[symbol])>0:
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
                if verbose:
                    print 'Successfully loaded variable ' + symbol
            elif dim==2:
                vals = {}
                for element in data[symbol]:
                    vals[element[0]] = element[1]          
                out[symbol] = pd.Series(vals)
                if verbose:
                    print 'Successfully loaded variable ' + symbol
            elif dim==1:
                print 'Variable ' + symbol + ' has dimension 0, which should not occur. Skipping'
            elif dim>3:
                print 'Variable ' + symbol + ' has more than 2 dimensions, which is very tiring. Skipping'
        else:
            if verbose: 
                print 'Variable ' + symbol + ' is empty. Skipping'
    for symbol in out:
        try:
            out[symbol].fillna(value=0,inplace=True)
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
        print("Time to convert to dataframes: {}s".format(tm.time()-tc)) 
    return out

def GetGdx(gams_dir,resultfile):
    ''' 
    Short wrapper of the two gdx reading functions (GdxToDataframe and GdxToList)
    
    :param gams_dir:    Gams working directory
    :param resultfile:  Path to the gdx file to be read
    :returns:           Dictionnary of dataframes
    '''
    return GdxToDataframe(GdxToList(gams_dir,resultfile,varname='all',verbose=True),fixindex=True,verbose=True)


