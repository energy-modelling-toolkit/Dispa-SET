# Collection of functions to write Dispa-SET input data to a gdx file and/or to a simulation directory
# with one excel file per parameter.

__author__ = 'Sylvain Quoilin (sylvain.quoilin@ec.europa.eu)'

import numpy as np
import os
import sys
 

def shrink_to_64(list_keys):
    '''
    Function that reduces the length of the keys to be written to 64 (max admissible length for GAMS)
    '''
    for i in range(len(list_keys)):
        key = list_keys[i]
        if len(key) > 63:
            list_keys[i] = key[:20] + ' ... ' + key[-20:]
    return list_keys


def insert_symbol(gdxHandle,variable):
    '''
    Function that writes a variable (either set or parameter) to the gdxHandle
    The format of the dispaset variable must be the version 2.0
    '''
    from gdxcc import *  

    domains = []                # Domains contains the names of the input sets for the parameter (empty for a set)
    description = None          # Describes the variable being written (to be added)

    # Check array sizes for parameters:
    if variable['type']=='parameter':
        # Check that the required fields are present:

        dims = len(variable['sets'])
        Nrows = variable['val'].shape[0]
        gdxSymbolType = GMS_DT_PAR
        if len(variable['val'].shape) != 2:
            sys.exit('Variable ' + variable['name'] + ': The \'val\' data matrix should be a 2-dimensional numpy array (NB: it can also be empty)')
        if variable['val'].shape[1] == 0:
            empty_val = True
            Nrows = 0
            max_uels = [-1 for x in range(dims)]
        else:
            empty_val = False
            max_uels = np.max(variable['val'],axis=0)          # Maximum value of the indices referring to the uels of each set
        # check that there are at least two rows in 'val'
        if variable['val'].shape[1] <= 1 and not empty_val:
            sys.exit('Variable ' + variable['name'] + ': The \'val\' data matrix should contain at least two rows (one for the uels, on for the parameter values')
        if dims != variable['val'].shape[1] -1 and not empty_val:
            sys.exit('The number of uels values (' +str(len(variable['sets'])) + ') does not match the array size for variable ' + variable['name'])
        for i in range(dims):
            set = variable['sets'][i]
            if len(variable['sets'][i]['uels']) < max_uels[i]+1:
                sys.exit('Variable ' + variable['name'] + ': The number of rows in the uels is smaller than some indices in the data table')
        for i in range(dims):
            domains.append(variable['sets'][i]['name'])

    elif variable['type']=='set':
        gdxSymbolType = GMS_DT_SET
        dims = 1
        Nrows = len(variable['uels'])
        variable['val'] = np.zeros([Nrows,1])
        variable['val'][:,0] = range(Nrows)

    else:
        sys.exit('Variable ' + variable['name'] + ': Type (' + variable['type'] + ') not supported')

    gdxDataWriteStrStart(gdxHandle, variable['name'], "", dims, gdxSymbolType, 0)
    gdxValues = doubleArray(5)
    gdxValues[GMS_VAL_LEVEL] = 0.0 #0.0 == Y (explanatory text of set in gdx)

    if variable['type']=='parameter':                                           # Write a parameter
        for row in range(Nrows):                                                # Write line by line
            gdxKeys = []                                                        # All the set values for this line
            for j in range(dims):
                key=variable['sets'][j]['uels'][int(variable['val'][row,j])]    # Get the string value of the uels from the set and using the indice in the val matrix
                gdxKeys.append(key)
            gdxKeys = shrink_to_64(gdxKeys)                                     # Reduce the size if bigger than 64 characters
            gdxValues[GMS_VAL_LEVEL] = variable['val'][row,dims]
            success = gdxDataWriteStr(gdxHandle, gdxKeys, gdxValues)
            if not success:
                print('Key ' + gdxKeys[0] + ' of ' + variable['type'] + ' ' + variable['name'] + ' could not be written')

    elif variable['type']=='set':                                               # Do the same if it a set
        for row in range(Nrows):
            gdxKeys = []
            gdxKeys.append(variable['uels'][row])
            #gdxAddSetText(gdxHandle, 'description')
            gdxKeys = shrink_to_64(gdxKeys)                                 # Reduce the size if bigger than 64 characters
            success = gdxDataWriteStr(gdxHandle, gdxKeys, gdxValues)
            if not success:
                print('Key ' + gdxKeys[0] + ' of ' + variable['type'] + ' ' + variable['name'] + ' could not be written')

    gdxDataWriteDone(gdxHandle)
    print(variable['type'] + ' ' + variable['name'] + ' successfully written')


def insert_symbol_21(gdxHandle,variable):
    '''
    Function that writes a variable (either set or parameter) to the gdxHandle
    The format of the dispaset variable must be the version 2.1
    '''
    from gdxcc import *  

    domains = []                # Domains contains the names of the input sets for the parameter (empty for a set)
    description = None          # Describes the variable being written (to be added)

    # Check array sizes for parameters:
    if variable['type']=='parameter':
        # Check that the required fields are present:

        dims = len(variable['sets'])
        shape = variable['val'].shape
        Nrows = variable['val'].shape[0]
        gdxSymbolType = GMS_DT_PAR
        if len(shape) != dims:
            sys.exit('Variable ' + variable['name'] + ': The \'val\' data matrix has ' + str(len(shape)) + ' dimensions and should have ' + str(dims))
        for i in range(dims):
            if shape[i] != len(variable['sets'][i]['uels']):
                sys.exit('Variable ' + variable['name'] + ': The \'val\' data matrix has ' + str(shape[i]) + ' elements for dimention ' + variable['sets'][i]['name'] + ' while there are ' + str(len(variable['sets'][i]['uels'])) + ' set values')
        else:
            empty_val = False

        for i in range(dims):
            domains.append(variable['sets'][i]['name'])

    elif variable['type']=='set':
        gdxSymbolType = GMS_DT_SET
        dims = 1
        Nrows = len(variable['uels'])
        variable['val'] = np.zeros([Nrows,1])
        variable['val'][:,0] = range(Nrows)

    else:
        print('Variable ' + variable['name'] + ': Type (' + variable['type'] + ') not supported')
        return()

    gdxDataWriteStrStart(gdxHandle, variable['name'], "", dims, gdxSymbolType, 0)
    gdxValues = doubleArray(5)
    gdxValues[GMS_VAL_LEVEL] = 0.0 #0.0 == Y (explanatory text of set in gdx)

    if variable['type']=='parameter':    
        for index,value in np.ndenumerate(variable['val']):
            # Write line by line if value is non null
            if value != 0:
                gdxKeys = []                                                        # All the set values for this line
                for j in range(dims):
                    key=variable['sets'][j]['uels'][index[j]]    # Get the string value of the uels from the set and using the indice in the val matrix
                    gdxKeys.append(key)
                gdxKeys = shrink_to_64(gdxKeys)                                     # Reduce the size if bigger than 64 characters
                gdxValues[GMS_VAL_LEVEL] = float(value)
                success = gdxDataWriteStr(gdxHandle, gdxKeys, gdxValues)
                if not success:
                    print('Key ' + gdxKeys[0] + ' of ' + variable['type'] + ' ' + variable['name'] + ' could not be written')

    elif variable['type']=='set':                                               # Do the same if it a set
        for row in range(Nrows):
            gdxKeys = []
            gdxKeys.append(variable['uels'][row])
            #gdxAddSetText(gdxHandle, 'description')
            gdxKeys = shrink_to_64(gdxKeys)                                 # Reduce the size if bigger than 64 characters
            success = gdxDataWriteStr(gdxHandle, gdxKeys, gdxValues)
            if not success:
                print('Key ' + gdxKeys[0] + ' of ' + variable['type'] + ' ' + variable['name'] + ' could not be written')

    gdxDataWriteDone(gdxHandle)
    print(variable['type'] + ' ' + variable['name'] + ' successfully written')

def insert_symbols_211(gdxHandle,sets,parameters):
    '''
    Function that writes a variable (either set or parameter) to the gdxHandle
    The format of the dispaset variable must be the version 2.1.1
    '''
    from gdxcc import *  
    
    # Check array sizes for parameters:
    for p in parameters:
        variable = parameters[p]
        
        # Check that the required fields are present:
        dims = len(variable['sets'])
        shape = variable['val'].shape
        Nrows = variable['val'].shape[0]
        gdxSymbolType = GMS_DT_PAR
        
        gdxDataWriteStrStart(gdxHandle, p, "", dims, gdxSymbolType, 0)
        gdxValues = doubleArray(5)
        gdxValues[GMS_VAL_LEVEL] = 0.0 #0.0 == Y (explanatory text of set in gdx)        
        
        if len(shape) != dims:
            sys.exit('Variable ' + variable['name'] + ': The \'val\' data matrix has ' + str(len(shape)) + ' dimensions and should have ' + str(dims))
        for i in range(dims):
            if shape[i] != len(sets[variable['sets'][i]]):
                sys.exit('Variable ' + variable['name'] + ': The \'val\' data matrix has ' + str(shape[i]) + ' elements for dimention ' + sets[variable['sets'][i]] + ' while there are ' + str(len(variable['sets'][i]['uels'])) + ' set values')

        for index,value in np.ndenumerate(variable['val']):
            # Write line by line if value is non null
            if value != 0:
                gdxKeys = []                                                        # All the set values for this line
                for i in range(dims):
                    key=sets[variable['sets'][i]][index[i]]                       # Get the string value of the set by using the indice in the val matrix
                    gdxKeys.append(key)
                gdxKeys = shrink_to_64(gdxKeys)                                     # Reduce the size if bigger than 64 characters
                gdxValues[GMS_VAL_LEVEL] = float(value)
                success = gdxDataWriteStr(gdxHandle, gdxKeys, gdxValues)
                if not success:
                    print('Key ' + gdxKeys[0] + ' of parameter ' + p + ' could not be written')
        gdxDataWriteDone(gdxHandle)
        print('Parameter ' + p + ' successfully written')


    for s in sets:
        gdxSymbolType = GMS_DT_SET
        dims = 1
        
        gdxDataWriteStrStart(gdxHandle, s, "", dims, gdxSymbolType, 0)
        gdxValues = doubleArray(5)
        gdxValues[GMS_VAL_LEVEL] = 0.0 #0.0 == Y (explanatory text of set in gdx)     
        
        Nrows = len(sets[s])

        for row in range(Nrows):
            gdxKeys = shrink_to_64([sets[s][row]])                                 # Reduce the size if bigger than 64 characters
            success = gdxDataWriteStr(gdxHandle, gdxKeys, gdxValues)
            if not success:
                print('Key ' + gdxKeys[0] + ' of set ' + s + ' could not be written')

        gdxDataWriteDone(gdxHandle)
        print('Set ' + s + ' successfully written')



# 
def write_variables(gams_dir,gdx_out,list_vars,format='2.0'):
    '''
    Function that reads all the variables (in list_vars) and inserts them one by one into gdx_out
    Currently accepts Dispaset data formats 2.0, 2.1 and 2.1.1
    '''
    from gdxcc import *  
    #Check gams_dir:
    if not os.path.isdir(gams_dir):
        sys.exit('Could not find the specific gams directory' + gams_dir)

    gdxHandle = new_gdxHandle_tp()
    gdxCreateD(gdxHandle, gams_dir, GMS_SSSIZE)
    gdxOpenWrite(gdxHandle, gdx_out, "")
    

    if format == '2.0':
        for var in list_vars:
            insert_symbol(gdxHandle,var)
    elif format == '2.1':
        for var in list_vars:
            insert_symbol_21(gdxHandle,var)
    elif format == '2.1.1':
        [sets, parameters] = list_vars
        insert_symbols_211(gdxHandle,sets,parameters)
    else:
        print 'Format ' + format + ' not supported'
        return()

    gdxClose(gdxHandle)

    #pickle.dump(list_vars, open( gdx_out[:-4]+'.p', "wb" ) )

    print 'Data Successfully written to ' + gdx_out


def write_toexcel(xls_out,list_vars,format='2.1',firstrow=7):
    '''
    Function that reads all the variables (in list_vars) and inserts them one by one into gdx_out
    :param xls_out: The path of the folder where the excel files are to be written
    :param list_vars: List containing the dispaset variables 
    :param format: Format version of the dispaset variables. Currently supports 2.1 and 2.1.1
    :param firstrow: Excel row number containing the first line of the DATA (not the headers)
    :return: Binary variable (True)
    '''
    
    import pandas as pd
    
    #import sys
    reload(sys)
    sys.setdefaultencoding("utf-8")
    
    if format == '2.1':
        print 'Using data format version 2.1'
        from DispaTools import load_var
        config = load_var(list_vars,'Config')
    elif format == '2.1.1':
        print 'Using data format version 2.1.1'       
        config = list_vars[1]['Config']['val']
    else:
        sys.exit('Format ' + format + ' not supported')
    
    if not os.path.exists(xls_out):
        os.mkdir(xls_out)    
    
    
    if config != []:      # if there is a config variable
        first_day = pd.datetime(config[0,0],config[0,1],config[0,2],0)
        last_day = pd.datetime(config[1,0],config[1,1],config[1,2],23)
        dates = pd.date_range(start=first_day,end=last_day,freq='1h')
    else:
        dates = []
        
    
    # Printing all sets in one sheet:    
    writer = pd.ExcelWriter(os.path.join(xls_out,'InputDispa-SET - Sets.xlsx'), engine='xlsxwriter')

    if format == '2.1':
        i = 0
        headers = []
        for var in list_vars:
            if var['type'] == 'set':
                df = pd.DataFrame(var['uels'],columns=[var['name']])
                df.to_excel(writer, sheet_name='Sets',startrow=1,startcol=i,header=True,index=False)
                headers.append(var['description'])
                i = i+1
        # write description headers:
        df = pd.DataFrame(columns = headers)
        df.to_excel(writer, sheet_name='Sets',startrow=0,startcol=0,header=True,index=False)
        writer.save()    
        print 'All sets successfully written to excel'
        
        # Printing each parameter in a separate sheet and workbook:
        for var in list_vars:
            if var['type'] == 'parameter':    
                dim = len(var['sets'])
                sets = [var['sets'][i]['name'] for i in range(dim)][:]
                if sets[-1] =='h' and isinstance(dates,pd.tseries.index.DatetimeIndex) and dim > 1:
                    if len(dates) != var['val'].shape[-1]:
                        sys.exit('The date range in the Config variable (' + str(len(dates)) + ' time steps) does not match the length of the time index (' + str(values.shape[2]) + ') for variable ' + var['name'])
                    var['firstrow'] = 5
                else:
                    var['firstrow']=1
                writer = pd.ExcelWriter(os.path.join(xls_out,'InputDispa-SET - ' + var['name'] + '.xlsx'), engine='xlsxwriter')
                if dim == 1:
                    df = pd.DataFrame(var['val'],columns=[var['name']],index=var['sets'][0]['uels'])
                    df.to_excel(writer, sheet_name=var['name'],startrow=var['firstrow'],startcol=0,header=True,index=True)
                    worksheet = writer.sheets[var['name']]
                    worksheet.write_string(0, 0, var['name']+'(' + var['sets'][0]['name'] + ')' )
                    worksheet.set_column(0,0,30)
                elif dim ==2:    
                    values = var['val']
                    df = pd.DataFrame(values,columns=var['sets'][1]['uels'],index=var['sets'][0]['uels'])
                    df.to_excel(writer, sheet_name=var['name'],startrow=var['firstrow'],startcol=0,header=True,index=True)
                    worksheet = writer.sheets[var['name']]
                    if var['firstrow'] == 5:
                        worksheet.write_row(1,1,dates.year)
                        worksheet.write_row(2,1,dates.month)
                        worksheet.write_row(3,1,dates.day)
                        worksheet.write_row(4,1,dates.hour+1)
                    worksheet.write_string(0, 0, var['name']+'(' + sets[0] + ',' + sets[1] + ')'  )
                    worksheet.freeze_panes(var['firstrow']+1,1)
                    worksheet.set_column(0,0,30)
                elif dim ==3:    
                    list_sets = var['sets'][:]
                    values = var['val']
                    for i in range(len(list_sets[0]['uels'])):
                        key = list_sets[0]['uels'][i]
                        Nrows = len(list_sets[1]['uels'])
                        df = pd.DataFrame(values[i,:,:],columns=list_sets[2]['uels'],index=list_sets[1]['uels'])
                        df.to_excel(writer, sheet_name=var['name'],startrow=var['firstrow']+1+i*Nrows,startcol=1,header=False,index=True)
                        df2 = pd.DataFrame(np.array([key]).repeat(Nrows))
                        df2.to_excel(writer, sheet_name=var['name'],startrow=var['firstrow']+1+i*Nrows,startcol=0,header=False,index=False)
                    worksheet = writer.sheets[var['name']]
                    if var['firstrow'] == 5:
                        worksheet.write_row(1,2,dates.year)
                        worksheet.write_row(2,2,dates.month)
                        worksheet.write_row(3,2,dates.day)
                        worksheet.write_row(4,2,dates.hour+1)                
                    worksheet.write_string(0, 0, var['name']+'(' + list_sets[0]['name'] + ',' + list_sets[1]['name'] + ',' + list_sets[2]['name']+ ')'  )
                    worksheet.write_string(var['firstrow']-1, 0, list_sets[0]['name'])
                    worksheet.write_string(var['firstrow']-1, 1, list_sets[1]['name'])
                    worksheet.freeze_panes(var['firstrow'],2)
                    worksheet.set_column(0,1,30)
                    df = pd.DataFrame(columns=list_sets[2]['uels'])
                    df.to_excel(writer, sheet_name=var['name'],startrow=var['firstrow'],startcol=2,header=True,index=False)
                else:
                    print 'Only three dimensions currently supported. Parameter ' + var['name'] + ' has ' + str(dim) + ' dimensions.'
                writer.save() 
                print 'Parameter ' + var['name'] + ' successfully written to excel'
    elif format == '2.1.1':
        [sets, parameters] = list_vars
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
                    sys.exit('The date range in the Config variable (' + str(len(dates)) + ' time steps) does not match the length of the time index (' + str(values.shape[2]) + ') for variable ' + p)
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

    if format == '2.1':
        for var in list_vars:   
            if var['type'] == 'set':
                gmsfile.write('\n')
                gmsfile.write('$CALL GDXXRW "InputDispa-SET - Sets.xlsx" Set=' + var['name'] + ' rng=' + chr(i + ord('A')) + '3 Rdim=1  O=' + var['name'] + '.gdx \n')
                gmsfile.write('$GDXIN ' + var['name'] + '.gdx \n')
                gmsfile.write('Set ' + var['name'] + '; \n')
                gmsfile.write('$LOAD ' + var['name']+ '\n') 
                gmsfile.write('$GDXIN \n')
                i = i+1
                
        for var in list_vars:
            if var['type'] == 'parameter':
                dim = len(var['sets'])
                gmsfile.write('\n')
                if dim ==1:
                    gmsfile.write('$CALL GDXXRW "InputDispa-SET - ' + var['name'] +'.xlsx" par=' + var['name'] + ' rng=A' + str(var['firstrow']+1) + ' Rdim=1 \n')
                elif dim ==2:
                    gmsfile.write('$CALL GDXXRW "InputDispa-SET - ' + var['name'] +'.xlsx" par=' + var['name'] + ' rng=A' + str(var['firstrow']+1) + ' Rdim=1 Cdim=1 \n')
                elif dim ==3:
                    gmsfile.write('$CALL GDXXRW "InputDispa-SET - ' + var['name'] +'.xlsx" par=' + var['name'] + ' rng=A' + str(var['firstrow']+1) + ' Rdim=2 Cdim=1 \n')
                gmsfile.write('$GDXIN "InputDispa-SET - ' + var['name'] + '.gdx" \n')
                gmsfile.write('Parameter ' + var['name'] + '; \n')
                gmsfile.write('$LOAD ' + var['name']+ '\n') 
                gmsfile.write('$GDXIN \n')
    
        gmsfile.write('\n')
        gmsfile.write('Execute_Unload "Inputs.gdx"')             
        gmsfile.close()
        
    elif format == '2.1.1':
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






