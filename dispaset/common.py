# -*- coding: utf-8 -*-
"""
This file defines a dictionary with global variables to be used in Dispa-SET such as fluids, technologies, etc.
"""
import datetime

commons={}
# Timestep
commons['TimeStep'] = '1h'

# DispaSET technologies:
commons['Technologies'] = ['COMC', 'GTUR', 'HDAM', 'HROR', 'HPHS', 'ICEN', 'PHOT', 'STUR', 'WTOF', 'WTON', 'CAES',
                           'BATS', 'BEVS', 'THMS', 'P2GS','P2HT','SCSP', 'COMC_CCS']
# List of VRES technologies:
commons['tech_renewables'] = ['WTON', 'WTOF', 'PHOT', 'HROR']
# List of storage technologies:
commons['tech_storage'] = ['HDAM', 'HPHS', 'BATS', 'BEVS', 'CAES', 'THMS','P2GS','SCSP']
# List of CHP types:
commons['types_CHP'] = ['extraction','back-pressure', 'p2h']
# DispaSET fuels:
commons['Fuels'] = ['BIO', 'GAS', 'HRD', 'LIG', 'NUC', 'OIL', 'PEA', 'SUN', 'WAT', 'WIN', 'WST', 'OTH', 'GEO','HYD']
# Ordered list of fuels for plotting (the first ones are negative):
# commons['MeritOrder'] = ['Storage','FlowOut','GEO','NUC', 'LIG', 'HRD', 'BIO', 'GAS', 'OIL', 'PEA', 'WST', 'OTH', 'SUN', 'WIN', 'FlowIn', 'WAT']
commons['MeritOrder'] = ['THMS', 'SCSP', 'BATS', 'BEVS', 'HDAM', 'HPHS', 'FlowOut', 'GEO', 'NUC', 'LIG',
                         'HRD', 'BIO', 'GAS', 'OIL', 'PEA', 'WST', 'OTH', 'SUN', 'WIN', 'FlowIn', 'WAT']

# Colors associated with each fuel:
# commons['colors'] = {'LIG': '#af4b9180', 'PEA': '#af4b9199', 'HRD': '#af4b91b2', 'OIL': '#af4b91ff',
# #                      'GAS': '#d7642dff',
# #                      'NUC': '#466eb4ff',
# #                      'SUN': '#e6a532ff',
# #                      'WIN': '#41afaaff',
# #                      'WAT': '#00a0e1ff',
# #                      'BIO': '#7daf4bff', 'GEO': '#7daf4bbf',
# #                      'Storage': '#b93c46ff', 'FlowIn': '#b93c46b2', 'FlowOut': '#b93c4666',
# #                      'OTH': '#b9c33799', 'WST': '#b9c337ff',
# #                      'HDAM': '#00a0e1ff',
# #                      'HPHS': '#3090C7ff',
# #                      'THMS': '#C04000ff',
# #                      'BATS': '#41A317ff',
# #                      'BEVS': '#CC80FFff'}
commons['colors'] = {'LIG': '#af4b9180', 'PEA': '#af4b9199', 'HRD': 'darkviolet', 'OIL': 'magenta',
                     'GAS': '#d7642dff',
                     'NUC': '#466eb4ff',
                     'SUN': '#e6a532ff',
                     'WIN': '#41afaaff',
                     'WAT': '#00a0e1ff',
                     'BIO': '#7daf4bff', 'GEO': '#7daf4bbf',
                     'Storage': '#b93c46ff', 'FlowIn': '#b93c46b2', 'FlowOut': '#b93c4666',
                     'OTH': '#b9c33799', 'WST': '#b9c337ff',
                     'HDAM': '#00a0e1ff',
                     'HPHS': 'blue',
                     'THMS': '#C04000ff',
                     'BATS': '#41A317ff',
                     'BEVS': '#b9c33799',
                     'SCSP': '#e6a532ff'}

commons['colors']['curtailment'] = 'red'
# Hatches associated with each fuel:
commons['hatches'] = {'LIG': '', 'PEA': '', 'HRD': '', 'OIL': '', 'GAS': '', 'NUC': '', 'SUN': '', 'WIN': '', 'WAT': '',
                      'BIO': '', 'GEO': '', 'Storage': '', 'WST': '', 'OTH': '',
                      'FlowIn': '/', 'FlowOut': '\\', 'HDAM': '/','HPHS': '/','SCSP': '/', 'THMS': '','BATS': '/','BEVS': '/'
                      }

commons['logfile'] = str(datetime.datetime.now()).replace(':','-').replace(' ','_') + '.dispa.log'