# -*- coding: utf-8 -*-
"""
This file defines a dictionary with global variables to be used in Dispa-SET such as fluids, technologies, etc.
"""
import datetime


class DispaSETValidationError(Exception):
    """Raised when Dispa-SET input data or configuration fails validation.

    Replaces the previous pattern of ``logging.critical(...); sys.exit(1)``
    so that the error can be caught programmatically (e.g. in Jupyter or
    in tests) rather than killing the interpreter.
    """

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
# NOTE (2026-08-20): SOTH removed -- per Docs/data.rst it is "Solar thermal district heating"
# (Heat only category, no electrical output), not a power-producing renewable; it was
# misclassified here (it was also already correctly listed in tech_boundary_sector below,
# so this was an inconsistent double-classification). SCSP added -- Docs/data.rst flags it
# VRES=Y ("Concentrated Solar Power").
commons['tech_renewables'] = ['HROR', 'PHOT', 'WAVE', 'WTOF', 'WTON', 'SCSP']
# List of Conventional technologies (real synchronous/turbine-based generation):
# NEW (2026-08-20): WHEN ("Waste heat engine", Power-only per Docs/data.rst), SCSP
# (Concentrated Solar Power with storage -- steam turbine), CAES (traditional CAES uses a
# real synchronous generator on discharge), and the CHP/boundary-sector variants of hydro/
# thermal units (same electrical-side machine as their non-CHP counterparts -- these were
# previously excluded from reserve eligibility purely because of this classification gap).
commons['tech_conventional'] = ['HDAM', 'HROR', 'HPHS',
                                'COMC', 'GTUR', 'ICEN', 'STUR',
                                'WHEN', 'SCSP', 'CAES',
                                'HDAMC', 'HRORC', 'HPHSC',
                                'COMCX', 'GTURX', 'ICENX', 'STURX',
                                ]#MARCO
# List of Batteries technologies:
# PENDING (not added, 2026-08-20): BEVS was considered for inclusion here (same underlying
# battery/power-electronics technology as BATS), but deliberately left out. The real Bolivia
# scenarios (S0-S5) have no BEVS units at all, so this has no effect on current results.
# Reasons to hold off until a BEVS-containing study is actually built: (1) `ba(au)` (and every
# reserve equation gated by it) is derived directly from this list, so adding BEVS here would
# immediately credit it with full battery-like reserve capability; (2) that is only realistic
# if the unit's AvailabilityFactor reflects actual vehicle plug-in patterns, not a flat 1 --
# there is real precedent for this in this repo (Database/AvailabilityFactors/NL/1h/
# 2015_withBEVS.csv, a hourly profile ranging 0.34-1.0), so the mechanism exists, it just
# needs real data prepared alongside enabling it; (3) V2G participation in ancillary-service
# markets is still regulatorily immature/pilot-stage in most jurisdictions, a scope question
# for the specific study, not something code can decide. Revisit together, not separately, if
# BEVS units are ever added to a study.
commons['tech_batteries'] = ['BATS']#MARCO
# List of storage technologies:
commons['tech_storage'] = ['HDAM', 'HPHS', 'BATS', 'BEVS', 'CAES', 'SCSP']
# commons['tech_storage'] = ['HDAM', 'HPHS', 'BATS', 'BEVS', 'CAES', 'SCSP', 'HDAMC']

# List of power to boundary sector technologies
# NEW (2026-08-20): HPHSC added -- it was previously unclassified anywhere in this file (not
# even here), unlike its hydro siblings HDAMC/HRORC.
commons['tech_p2bs'] = ['P2GS', 'ALKE', 'PEME', 'SOXE', 'P2BS', 'PEFC',
                        'DMFC', 'ALFC', 'PAFC', 'MCFC', 'SOFC', 'REFC', 'HDAMC', 'HRORC', 'HPHSC', 'HDLZ',
                        'COMCX', 'GTURX', 'ICENX', 'STURX',
                        'P2HT', 'ASHP', 'GSHP', 'HYHP', 'WSHP', 'REHE']
# PENDING (not resolved, 2026-08-20): fuel cells (PEFC/DMFC/ALFC/PAFC/MCFC/SOFC/REFC),
# electrolyzers (P2GS/ALKE/PEME/SOXE), bio-hydrolysis (HDLZ), and Haber-Bosch (HBBS, still
# unclassified anywhere in this file) are custom additions to this project, not part of the
# officially-documented Dispa-SET technology set (Docs/data.rst does not cover them), so no
# literature backing specific to *this model's* representation of them was found this
# session. Whether/how they should participate in frequency reserves (e.g. fuel cells as slow
# X2P generation with limited FCR/mFRR capability; electrolyzers as demand-response-style
# downward reserve) is an open question requiring the owner's domain input -- not decided.
# List of boundary sector to power technologies
commons['tech_bs2p'] = ['BSPG']
# List of boundary sector only technologies:
commons['tech_boundary_sector'] = ['BSPG', 'GETH', 'HOBO', 'SOTH', 'ABHP', 'HOBOX', 'P2BS']#MARCO , 'HDAMC', 'HRORC'
# List of CHP types:
commons['types_CHP'] = ['extraction', 'back-pressure', 'p2h']

# Physical constants used to size droop-based reserve-participation factors. Centralized here
# instead of a local dict in build.py -- same precedent as ReserveDuration/FullActivationTime/
# default_StartUpTime below. SystemFrequency is region-specific (Bolivia = 50 Hz); move to
# config instead if this is ever reused for a 60 Hz study. RoCoF_max is kept for reference/
# possible future use -- it is not currently used since VIRU participation is binary, not
# RoCoF-derated.
commons['FrequencyResponseConstants'] = {
    'SystemFrequency': 50,        # Hz
    'DeltaFrequencyMax': 0.8,     # Hz -- reference frequency deviation for droop-based sizing
    'RoCoF_max': 0.5,             # Hz/s -- reference RoCoF (see note above)
}

# Which technologies may provide each reserve product. Derived from the technology-class
# lists above (single source of truth -- do not hardcode a separate flat technology list
# here). VIRU/FFRU/FFRD: batteries only. Real rotating-mass inertia from synchronous machines
# is handled separately via EQ_Inertia_Delivery_Limit/cu(au); VIRU specifically represents
# *synthetic* inertia from non-synchronous resources. Wind (WTOF/
# WTON) was considered and deliberately NOT added to VIRU/FFRU this round: grid-following wind
# turbines have zero inertia by default (their power electronics decouple the rotor from grid
# frequency), synthetic-inertia control is a non-standard add-on not universally installed,
# and even where present it causes a "recovery dip" (temporary output reduction to re-
# accelerate the rotor) that this model has no way to represent -- crediting it by default
# would repeat the same class of overstatement this session's audit was about fixing.
commons['ReserveEligibleTechnologies'] = {
    'VIRU':  commons['tech_batteries'],
    'FFRU':  commons['tech_batteries'],
    'FFRD':  commons['tech_batteries'],
    'FCRU':  commons['tech_batteries'] + commons['tech_conventional'] + commons['tech_renewables'],
    'FCRD':  commons['tech_batteries'] + commons['tech_conventional'] + commons['tech_renewables'],
    'aFRRU': commons['tech_batteries'] + commons['tech_conventional'] + commons['tech_renewables'],
    'aFRRD': commons['tech_batteries'] + commons['tech_conventional'] + commons['tech_renewables'],
    'mFRRU': commons['tech_batteries'] + commons['tech_conventional'] + commons['tech_renewables'],
}

# Fraction of zonal demand sheddable per reserve category for emergency UFLS (Under-Frequency
# Load Shedding) / OFDM (downward emergency action). Centralized here instead of a hardcoded
# loop in build.py. Values unchanged from before; still proposed/order-of-magnitude, not tied
# to a specific relay-staging scheme -- these should also vary by zone eventually.
commons['UFLS_Participation'] = {'FFRU': 0.1, 'FCRU': 0.2, 'aFRRU': 0.2, 'mFRRU': 0.2}
commons['OFDM_Participation'] = {'FFRD': 0.1, 'FCRD': 0.2, 'aFRRD': 0.2}

# Default StartUpTime [h], used ONLY when a unit's PowerPlantData value is missing (NaN) --
# an explicitly-entered value (including 0) is real data and is never overridden. Proposed,
# order-of-magnitude, technology-representative placeholders (not a substitute for real
# per-unit data).
# NOTE: deriving a missing value from RampUpRate/PartLoadMin instead was considered and
# rejected -- that ramp-implied time is a lower bound under an unrealistically optimistic
# assumption (no separate thermal/mechanical start-up delay), not an estimate.
commons['default_StartUpTime'] = {
    # Batteries / power-electronics-coupled storage: near-instantaneous
    'BATS': 0, 'BEVS': 0, 'CAES': 0,
    # Variable renewables: no thermal start-up process
    'PHOT': 0, 'WTON': 0, 'WTOF': 0, 'WAVE': 0, 'HROR': 0,
    # Hydro: mechanically fast to start
    'HDAM': 0.083, 'HPHS': 0.083,
    # Fast-start thermal
    'GTUR': 0.17, 'ICEN': 0.083,
    # Slower thermal (representative hot/warm start; cold start is materially longer)
    'COMC': 2, 'STUR': 4,
    # Boundary-sector prime movers: electrically driven, fast-responding
    'ASHP': 0.083, 'GSHP': 0.083, 'HYHP': 0.083, 'WSHP': 0.083, 'ABHP': 0.083,
    'REHE': 0.083, 'P2HT': 0.083,
}
# Conservative fallback for any technology not listed above (always paired with a logged
# warning naming the unlisted technology, so the assumption is never silent).
commons['default_StartUpTime_fallback'] = 0.17

# Reserve-product duration and full-activation-time requirements [h], derived from
# commons['ReserveTiming'] below (ramp/3600 = FullActivationTime; (deact-ramp)/3600 = ReserveDuration).
commons['ReserveDuration'] = {
    'VIRU': 0.0000833, 'FFRU': 0.008333, 'FCRU': 0.25, 'aFRRU': 0.125, 'mFRRU': 1.0,
    'FFRD': 0.008333, 'FCRD': 0.25, 'aFRRD': 0.125,
}
commons['FullActivationTime'] = {
    'FFRU': 0.000361, 'FCRU': 0.008333, 'aFRRU': 0.083333, 'mFRRU': 0.208333,
    'FFRD': 0.000361, 'FCRD': 0.008333, 'aFRRD': 0.083333,
}

# Reserve-product activation timeline [s] from contingency onset (t=0): prep/ramp/delivery/deact,
# per ENTSO-E FCR/aFRR/mFRR (box-car, delivery=deact); VIRU decays until FFRU's ramp instead, for
# lack of a published duration. Override per-case via config['ReserveTiming'].
commons['ReserveTiming'] = {
    'VIRU':  dict(prep=0.50, ramp=1.00, delivery=1.00, deact=1.30),
    'FFRU':  dict(prep=0.90, ramp=1.30, delivery=31.3, deact=31.3),
    'FCRU':  dict(prep=1.10, ramp=30.00, delivery=930, deact=930),
    'aFRRU': dict(prep=31.00, ramp=300.00, delivery=750, deact=750),
    'mFRRU': dict(prep=300.00, ramp=750.00),
}
# DispaSET fuels:
commons['Fuels'] = ['AIR', 'AMO', 'BIO', 'GAS', 'HRD', 'LIG', 'NUC', 'OIL', 'PEA', 'SUN', 'WAT', 'WIN', 'WST', 'OTH',
                    'GEO', 'HYD', 'WHT', 'ELE', 'THE']
# Ordered list of fuels for plotting (the first ones are negative):
commons['MeritOrder'] = ['SCSP', 'BATS', 'BEVS', 'HDAM', 'HPHS', 'P2X', 'FlowOut', 'GEO', 'NUC', 'LIG',
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
                     'HPHS': '#00a0e1ff',
                     'THMS': '#C04000ff',
                     'BATS': '#41A317ff',
                     'BEVS': '#b9c33799',
                     'SCSP': '#e6a532ff',
                     'P2X': '#A0522D',
                     'X2P': '#A0522D',
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
                      'FlowIn': '', 'FlowOut': '', 'HDAM': '', 'HPHS': '', 'SCSP': '', 'THMS': '', 'BATS': '',
                      'BEVS': '', 'P2X': '', 'X2P': '','AIR': '', 'WHT': '', 'ELE': '', 'THE': ''
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
    'SimulationType': 76, 'ReserveCalculation': 77, 'AllowCurtailment': 78, 'TransmissionGridType': 79, 'FrequencyStability': 80,
    'SectorCoupling': 81,
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
    'GridData': 154,
    # 'PTDFMatrix': 145,
    # Inertia Limit data
    'InertiaLimit': 155,
    # Hydrogen data
    'H2RigidDemand': 137, 'H2FlexibleDemand': 138, 'H2FlexibleCapacity': 139,
    # Storage data
    'StorageAlertLevels': 140, 'StorageFloodControl': 141, 
    # Reserves input data
    'FFRLimit': 158, 'PrimaryReserveLimit': 159, 'Reserve2U': 160, 'Reserve2D': 161, 
    # Other costs related data
    'PriceOfCO2': 166, 'CostHeatSlack': 167, 'CostLoadShedding': 168, 'PriceTransmission': 169, 'CostH2Slack': 170, 
    'CostCurtailment': 171, 'CostNotServed': 172,
    # Fuel price related data
    'PriceOfNuclear': 180, 'PriceOfBlackCoal': 181, 'PriceOfGas': 182, 'PriceOfFuelOil': 183,
    'PriceOfBiomass': 184, 'PriceOfLignite': 185, 'PriceOfPeat': 186, 'PriceOfAmmonia': 187,
    'CostOfSpillage': 205

}

commons['modifiers'] = {'Demand': 274, 'Wind': 275, 'Solar': 276, 'Storage': 277}

commons['default'] = {
    # Hydro scheduling defaults
    'ReservoirLevelInitial': 101, 'ReservoirLevelFinal': 102,
    # Fuel price defaults
    'PriceOfNuclear': 180, 'PriceOfBlackCoal': 181, 'PriceOfGas': 182, 'PriceOfFuelOil': 183,
    'PriceOfBiomass': 184, 'PriceOfLignite': 185, 'PriceOfPeat': 186, 'PriceOfAmmonia': 187,
    # Other price defaults
    'PriceOfCO2': 166, 'CostHeatSlack': 167, 'CostLoadShedding': 168, 'PriceTransmission': 169, 'CostH2Slack': 170,
    'CostCurtailment': 171, 'CostNotServed': 172,
    # Optimization and infeasibility cost data
    'ShareOfFlexibleDemand': 125, 'LoadShedding': 129,
    'DemandFlexibility': 162, 'ShareOfQuickStartUnits': 163,
    'ValueOfLostLoad': 204, 'CostOfSpillage': 205, 'WaterValue': 206,
    # Inertia requirement default
    'InertiaLimit': 155,
    # Inertia requirement default
    'FFRLimit': 158, 'PrimaryReserveLimit': 159
}
