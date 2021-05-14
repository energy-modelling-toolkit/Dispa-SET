# -*- coding: utf-8 -*-
"""
This file defines a dictionary with global variables to be used in Dispa-SET such as fluids, technologies, etc.
"""
import datetime

commons = {}
# Timestep
commons['TimeStep'] = '1h'

# DispaSET technologies:
# 1) Hydro and renewables
# 2) Thermal
# 3) Storage
# 4) P2HT
# 5) Fuel cells
# 6) Electrolyzers
commons['Technologies'] = ['HDAM', 'HROR', 'HPHS', 'PHOT', 'WAVE', 'WHEN', 'WTOF', 'WTON',
                           'COMC', 'GTUR', 'ICEN', 'SCSP', 'STUR',
                           'BATS', 'BEVS', 'CAES', 'THMS',
                           'ABHP', 'ASHP', 'GETH', 'GSHP', 'HOBO', 'HYHP', 'P2HT', 'REHE', 'SOTH', 'WSHP',
                           'PEFC', 'DMFC', 'ALFC', 'PAFC', 'MCFC', 'SOFC', 'REFC',
                           'P2GS', 'ALKE', 'PEME', 'SOXE']
# List of VRES technologies:
commons['tech_renewables'] = ['HROR', 'PHOT', 'WAVE', 'WTOF', 'WTON', 'SOTH']
# List of storage technologies:
commons['tech_storage'] = ['HDAM', 'HPHS', 'BATS', 'BEVS', 'CAES', 'P2GS', 'SCSP']
# List of power to heat technologies:
commons['tech_p2ht'] = ['P2HT', 'ABHP', 'ASHP', 'GSHP', 'HYHP', 'WSHP', 'REHE']
# List of heat only technologies:
commons['tech_heat'] = ['GETH', 'HOBO', 'SOTH']
# List of thermal storage technologies:
commons['tech_thermal_storage'] = ['THMS']
# List of CHP types:
commons['types_CHP'] = ['extraction', 'back-pressure', 'p2h']
# DispaSET fuels:
commons['Fuels'] = ['AIR', 'BIO', 'GAS', 'HRD', 'LIG', 'NUC', 'OIL', 'PEA', 'SUN', 'WAT', 'WIN', 'WST', 'OTH', 'GEO',
                    'HYD', 'WHT']
# Ordered list of fuels for plotting (the first ones are negative):
commons['MeritOrder'] = ['SCSP', 'BATS', 'BEVS', 'HDAM', 'HPHS', 'P2GS', 'FlowOut', 'GEO', 'NUC', 'LIG',
                         'HRD', 'BIO', 'GAS', 'OIL', 'PEA', 'WST', 'OTH', 'SUN', 'WIN', 'FlowIn', 'WAT',
                         'HYD', 'AIR', 'WHT']
commons['MeritOrderHeat'] = ['GEO', 'NUC', 'LIG', 'HRD', 'BIO', 'GAS', 'OIL', 'PEA', 'WST', 'OTH', 'SUN', 'WIN', 'WAT',
                             'HYD', 'AIR', 'WHT', 'HeatSlack', 'THMS']

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
commons['colors'] = {'LIG': '#af4b9180',
                     'PEA': '#af4b9199',
                     'HRD': 'darkviolet',
                     'OIL': 'magenta',
                     'GAS': '#d7642dff',
                     'NUC': '#466eb4ff',
                     'SUN': '#e6a532ff',
                     'WIN': '#41afaaff',
                     'WAT': '#00a0e1ff',
                     'HYD': '#A0522D',
                     'BIO': '#7daf4bff',
                     'GEO': '#7daf4bbf',
                     'Storage': '#b93c46ff',
                     'FlowIn': '#b93c46b2',
                     'FlowOut': '#b93c4666',
                     'OTH': '#57D53B',
                     'WST': '#b9c337ff',
                     'HDAM': '#00a0e1ff',
                     'HPHS': 'blue',
                     'THMS': '#C04000ff',
                     'BATS': '#41A317ff',
                     'BEVS': '#b9c33799',
                     'SCSP': '#e6a532ff',
                     'P2GS': '#A0522D',
                     'ShedLoad': '#ffffffff',
                     'AIR': '#aed6f1ff',
                     'WHT': '#a93226ff',
                     'HeatSlack': '#943126ff'}

commons['colors']['curtailment'] = 'red'
# Hatches associated with each fuel:
commons['hatches'] = {'LIG': '', 'PEA': '', 'HRD': '', 'OIL': '', 'GAS': '', 'NUC': '', 'SUN': '', 'WIN': '', 'WAT': '',
                      'BIO': '', 'GEO': '', 'Storage': '', 'WST': '', 'OTH': '', 'HYD': '',
                      'FlowIn': '/', 'FlowOut': '\\', 'HDAM': '/', 'HPHS': '/', 'SCSP': '/', 'THMS': '', 'BATS': '/',
                      'BEVS': '/', 'P2GS': '/', 'AIR': '', 'WHT': '', 'HeatSlack': '/'
                      }

commons['logfile'] = str(datetime.datetime.now()).replace(':', '-').replace(' ', '_') + '.dispa.log'

# Specifying the na value is required to avoid configusion with the 'NA' (Namibia) country code
commons['na_values']=['', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', '-NaN', '-nan',
                                        '1.#IND', '1.#QNAN', 'N/A', 'NULL', 'NaN', 'nan']