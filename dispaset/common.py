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
                           # BS Version
                           'HDAMC', 'HRORC', 'HPHSC', 'COMCX', 'GTURX', 'ICENX', 'STURX', 'HOBOX',
                           # Storage
                           'BATS', 'BEVS', 'CAES', 'THMS', 'H2ST',
                           'ABHP', 'ASHP', 'GETH', 'GSHP', 'HOBO', 'HYHP', 'P2HT', 'REHE', 'SOTH', 'WSHP',
                           'PEFC', 'DMFC', 'ALFC', 'PAFC', 'MCFC', 'SOFC', 'REFC',
                           'P2GS', 'ALKE', 'PEME', 'SOXE',
                           'P2BS', 'BSPG',
                           'HDLZ',  # Bio-Hydrolysis
                           'HBBS'  # Haber-Bosch
                           ]
# List of VRES technologies:
commons['tech_renewables'] = ['HROR', 'PHOT', 'WAVE', 'WTOF', 'WTON', 'SOTH']
# List of Conventional technologies:
#commons['tech_conventional'] = ['HDAM', 'HROR', 'COMC', 'GTUR', 'STUR']
# List of storage technologies:
commons['tech_storage'] = ['HDAM', 'HPHS', 'BATS', 'BEVS', 'CAES', 'SCSP']
# commons['tech_storage'] = ['HDAM', 'HPHS', 'BATS', 'BEVS', 'CAES', 'SCSP', 'HDAMC']
# List of power to heat technologies:
commons['tech_p2ht'] = ['P2HT', 'ASHP', 'GSHP', 'HYHP', 'WSHP', 'REHE']
# List of power to boundary sector technologies
commons['tech_p2bs'] = []#'P2GS', 'ALKE', 'PEME', 'SOXE', 'P2BS', 'PEFC',
                        #'DMFC', 'ALFC', 'PAFC', 'MCFC', 'SOFC', 'REFC', 'HDAMC', 'HRORC', 'HDLZ',
                        #'COMCX', 'GTURX', 'ICENX', 'STURX'
# # List of heat only technologies:
# commons['tech_heat'] = ['GETH', 'HOBO', 'SOTH', 'ABHP']
# List of noundary sector only technologies:
commons['tech_boundary_sector'] = []#'BSPG', 'GETH', 'HOBO', 'SOTH', 'ABHP', 'HOBOX'
# List of thermal storage technologies:
commons['tech_thermal_storage'] = ['THMS']
# List of CHP types:
commons['types_CHP'] = ['extraction', 'back-pressure', 'p2h']
# DispaSET fuels:
commons['Fuels'] = ['AIR', 'AMO', 'BIO', 'GAS', 'HRD', 'LIG', 'NUC', 'OIL', 'PEA', 'SUN', 'WAT', 'WIN', 'WST', 'OTH',
                    'GEO', 'HYD', 'WHT', 'ELE', 'THE']
# Ordered list of fuels for plotting (the first ones are negative):
commons['MeritOrder'] = ['SCSP', 'BATS', 'BEVS', 'HDAM', 'HPHS', 'P2GS', 'FlowOut', 'GEO', 'NUC', 'LIG',
                         'HRD', 'BIO', 'AMO', 'GAS', 'OIL', 'PEA', 'WST', 'OTH', 'SUN', 'WIN', 'FlowIn', 'WAT',
                         'HYD', 'AIR', 'WHT', 'ELE']
commons['MeritOrderHeat'] = ['GEO', 'NUC', 'LIG', 'HRD', 'BIO', 'AMO', 'GAS', 'OIL', 'PEA', 'WST', 'OTH', 'SUN', 'WIN',
                             'WAT', 'HYD', 'AIR', 'WHT', 'ELE', 'THE', 'HeatSlack', 'THMS']

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
                     'AMO': '#ffff00ff',
                     'GEO': '#7daf4bbf',
                     'Storage': '#b93c46ff',
                     'FlowIn': 'red',
                     'FlowOut': 'green',
                     # 'FlowIn': '#b93c46b2',
                     # 'FlowOut': '#b93c4666',
                     'OTH': '#57D53B',
                     'WST': '#b9c337ff',
                     'HDAM': '#00a0e1ff',
                     'HDAMC': '#00a0e1ff',
                     'HPHS': 'blue',
                     'THMS': '#C04000ff',
                     'BATS': '#41A317ff',
                     'BEVS': '#b9c33799',
                     'SCSP': '#e6a532ff',
                     'P2GS': '#A0522D',
                     'ShedLoad': '#ffffffff',
                     'AIR': '#aed6f1ff',
                     'WHT': '#a93226ff',
                     'ELE': '#2C75FFff',
                     'THE': '#c70509ff',
                     'HeatSlack': '#943126ff',
                     'curtailment': 'red'}

commons['colors']['curtailment'] = 'red'
commons['colors']['shed load'] = 'white'
commons['colors']['reserves'] = 'black'
# Hatches associated with each fuel:
commons['hatches'] = {'LIG': '', 'PEA': '', 'HRD': '', 'OIL': '', 'GAS': '', 'NUC': '', 'SUN': '', 'WIN': '', 'WAT': '',
                      'BIO': '', 'AMO': '', 'GEO': '', 'Storage': '', 'WST': '', 'OTH': '', 'HYD': '',
                      'FlowIn': '//', 'FlowOut': '//', 'HDAM': '/', 'HPHS': '/', 'SCSP': '/', 'THMS': '', 'BATS': '/',
                      'BEVS': '/', 'P2GS': '/', 'AIR': '', 'WHT': '', 'HeatSlack': '/', 'ELE': '', 'THE': ''
                      }

commons['logfile'] = str(datetime.datetime.now()).replace(':', '-').replace(' ', '_') + '.dispa.log'

# Specifying the na value is required to avoid configusion with the 'NA' (Namibia) country code
commons['na_values']=['', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', '-NaN', '-nan',
                                        '1.#IND', '1.#QNAN', 'N/A', 'NULL', 'NaN', 'nan']

commons['StdParameters'] = {
    # Scenario options
    'SimulationDirectory': 33, 'WriteGDX': 34, 'WritePickle': 35, 'GAMS_folder': 36,
    'cplex_path': 37,
    # Horizon Settings
    'DataTimeStep': 60, 'SimulationTimeStep': 61,
    # Simulation Options
    'SimulationType': 76, 'ReserveCalculation': 77, 'AllowCurtailment': 78, 'TransmissionGridType': 79, 'FrequencyServices': 80,
    # Mid-term scheduling related
    'HydroScheduling': 98, 'HydroSchedulingHorizon': 99, 'InitialFinalReservoirLevel': 100
}

commons['PathParameters'] = {
    # Power system data
    'Demand': 124, 'ShareOfFlexibleDemand': 125, 'Outages': 126, 'PowerPlantData': 127,
    'RenewablesAF': 128, 'LoadShedding': 129,
    # Interconnection data
    'NTC': 130, 'Interconnections': 131,
    # Hydro data
    'ReservoirScaledInflows': 132, 'ReservoirLevels': 133,
    # Heat data
    'HeatDemand': 134, 'Temperatures': 135,
    # Geo data
    'GeoData': 136,
    # DC-Power Flow data
    'GridData': 145,
    # 'PTDFMatrix': 145,
    # Inertia Limit data
    'InertiaLimit': 146,
    # # Hydrogen data
    # 'H2RigidDemand': 137, 'H2FlexibleDemand': 138, 'H2FlexibleCapacity': 139,
    # Reserves input data
    'PrimaryResponse': 159, 'Reserve2U': 160, 'Reserve2D': 161, 
    # Other costs related data
    'PriceOfCO2': 166, 'CostHeatSlack': 167, 'CostLoadShedding': 168, 'PriceTransmission': 169, 'CostCurtailment': 171,
    # Fuel price related data
    'PriceOfNuclear': 180, 'PriceOfBlackCoal': 181, 'PriceOfGas': 182, 'PriceOfFuelOil': 183,
    'PriceOfBiomass': 184, 'PriceOfLignite': 185, 'PriceOfPeat': 186, 'PriceOfAmmonia': 187 

}

commons['modifiers'] = {'Demand': 274, 'Wind': 275, 'Solar': 276, 'Storage': 277}

commons['default'] = {
    # Hydro scheduling defaults
    'ReservoirLevelInitial': 101, 'ReservoirLevelFinal': 102,
    # Fuel price defaults
    'PriceOfNuclear': 180, 'PriceOfBlackCoal': 181, 'PriceOfGas': 182, 'PriceOfFuelOil': 183,
    'PriceOfBiomass': 184, 'PriceOfLignite': 185, 'PriceOfPeat': 186, 'PriceOfAmmonia': 187,
    # Other price defaults
    'PriceOfCO2': 166, 'CostHeatSlack': 167, 'CostLoadShedding': 168, 'PriceTransmission': 169, 'CostCurtailment': 171,
    # Optimization and infeasibility cost data
    'ShareOfFlexibleDemand': 125, 'LoadShedding': 129,
    'DemandFlexibility': 162, 'ShareOfQuickStartUnits': 163,
    'ValueOfLostLoad': 204, 'PriceOfSpillage': 205, 'WaterValue': 206,
    # Inertia requirement default
    'InertiaLimit': 146
}
