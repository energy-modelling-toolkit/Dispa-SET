# -*- coding: utf-8 -*-
"""
This file defines the global variables to be used in Dispa-SET such as fluids, technologies, etc.

@author: 'Sylvain Quoilin'
"""
import itertools

def commonvars():
    '''
    Function providing a dictionary with the common variable definitions
    to be used in Dispa-SET
    '''
    commons={}
    # Timestep
    commons['TimeStep'] = '1h'

    # DispaSET technologies:
    commons['Technologies'] = ['COMC', 'GTUR', 'HDAM', 'HROR', 'HPHS', 'ICEN', 'PHOT', 'STUR', 'WTOF', 'WTON', 'CAES', 'BATS',
                    'BEVS', 'THMS', 'P2GS']
    # List of renewable technologies:
    commons['tech_renewables'] = ['WTON', 'WTOF', 'PHOT', 'HROR']
    # List of storage technologies:
    commons['tech_storage'] = ['HDAM', 'HPHS', 'BATS', 'BEVS', 'CAES', 'THMS']
    # List of CHP types:
    commons['types_CHP'] = ['extraction','back-pressure', 'p2h']
    # DispaSET fuels:
    commons['Fuels'] = ['BIO', 'GAS', 'HRD', 'LIG', 'NUC', 'OIL', 'PEA', 'SUN', 'WAT', 'WIN', 'WST', 'OTH', 'GEO']
    # Ordered list of fuels for plotting (the first ones are negative):
    commons['MeritOrder'] = ['Storage','FlowOut','GEO','NUC', 'LIG', 'HRD', 'BIO', 'GAS', 'OIL', 'PEA', 'WST', 'OTH', 'SUN', 'WIN', 'FlowIn', 'WAT']
    # Colors associated with each fuel:
    commons['colors'] = {'NUC': 'orange', 'LIG': 'brown', 'HRD': 'grey', 'BIO': 'darkgreen', 'GAS': 'lightcoral', 
           'OIL': 'chocolate','PEA':'green', 'WST': 'dodgerblue', 'OTH':'grey', 'SUN': 'yellow', 'WIN': 'red', 'FlowIn': 'green', 'WAT': 'blue', 
           'Storage': 'blue','FlowOut': 'green','GEO':'grey'}
    # Hatches associated with each fuel (random):
    hatches = itertools.cycle(['x', '//', '\\', '/'])
    commons['hatches'] = {}
    for x in commons['colors']:
        commons['hatches'][x] = next(hatches)

    return commons