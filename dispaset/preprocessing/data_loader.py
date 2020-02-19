
import datetime as dt
import logging
import os
import sys

import numpy as np
import pandas as pd
try:
    from future.builtins import int
except ImportError:
    logging.warning("Couldn't import future package. Numeric operations may differ among different versions due to incompatible variable types")
    pass

from .data_check import check_units, check_chp, check_sto, check_heat_demand, check_df, isStorage, check_MinMaxFlows,check_AvailabilityFactors, check_clustering
from .utils import clustering, interconnections, incidence_matrix
from .data_handler import UnitBasedTable,NodeBasedTable,merge_series, load_csv

from ..common import commons  # Load fuel types, technologies, timestep, etc:



def get_indices(config):
    # Indexes of the simulation:
    idx_std = pd.DatetimeIndex(pd.date_range(start=pd.datetime(*config['StartDate']),
                                            end=pd.datetime(*config['StopDate']),
                                            freq=commons['TimeStep'])
                            ) #todo check brackets on master

    idx_utc_noloc = idx_std - dt.timedelta(hours=1)
    idx_utc = idx_utc_noloc.tz_localize('UTC')
    # Indexes for the whole year considered in StartDate
    idx_utc_year_noloc = pd.DatetimeIndex(pd.date_range(start=pd.datetime(*(config['StartDate'][0],1,1,0,0)),
                                                        end=pd.datetime(*(config['StartDate'][0],12,31,23,59,59)),
                                                        freq=commons['TimeStep'])
                                        )
    return idx_utc, idx_utc_noloc, idx_utc_year_noloc

def load_loads(config, idx_utc_noloc, idx_utc_year_noloc):

    # Load :
    Load = NodeBasedTable(config['Demand'],idx_utc_noloc,config['countries'],tablename='Demand')
    # For the peak load, the whole year is considered:
    PeakLoad = NodeBasedTable(config['Demand'],idx_utc_year_noloc,config['countries'],tablename='PeakLoad').max()
    if config['modifiers']['Demand'] != 1:
        logging.info('Scaling load curve by a factor ' + str(config['modifiers']['Demand']))
        Load = Load * config['modifiers']['Demand']
        PeakLoad = PeakLoad * config['modifiers']['Demand']
    return Load, PeakLoad

def load_interconnections(config, idx_utc_noloc):

    # Interconnections:
    if os.path.isfile(config['Interconnections']):
        flows = load_csv(config['Interconnections'], index_col=0, parse_dates=True).fillna(0)
    else:
        logging.warning('No historical flows will be considered (no valid file provided)')
        flows = pd.DataFrame(index=idx_utc_noloc)
    if os.path.isfile(config['NTC']):
        NTC = load_csv(config['NTC'], index_col=0, parse_dates=True).fillna(0)
    else:
        logging.warning('No NTC values will be considered (no valid file provided)')
        NTC = pd.DataFrame(index=idx_utc_noloc)

    # Interconnections:
    [Interconnections_sim, Interconnections_RoW, Interconnections] = interconnections(config['countries'], NTC, flows)

    if len(Interconnections_sim.columns) > 0:
        NTCs = Interconnections_sim.reindex(idx_utc_noloc)
    else:
        NTCs = pd.DataFrame(index=idx_utc_noloc)
    Inter_RoW = Interconnections_RoW.reindex(idx_utc_noloc)
    return flows, NTC, Interconnections, NTCs, Inter_RoW

def load_load_shedding(config, idx_utc_noloc):
    # Load Shedding:
    LoadShedding = NodeBasedTable(config['LoadShedding'],idx_utc_noloc,config['countries'],tablename='LoadShedding',default=config['default']['LoadShedding'])
    CostLoadShedding = NodeBasedTable(config['CostLoadShedding'],idx_utc_noloc,config['countries'],tablename='CostLoadShedding',default=config['default']['CostLoadShedding'])
    return LoadShedding, CostLoadShedding

def load_fuel_prices(config, idx_utc_noloc):

    # Fuel prices:
    fuels = ['PriceOfNuclear', 'PriceOfBlackCoal', 'PriceOfGas', 'PriceOfFuelOil', 'PriceOfBiomass', 'PriceOfCO2', 'PriceOfLignite', 'PriceOfPeat']
    FuelPrices = pd.DataFrame(columns=fuels, index=idx_utc_noloc)
    for fuel in fuels:
        if os.path.isfile(config[fuel]):
            tmp = load_csv(config[fuel], header=None, index_col=0, parse_dates=True)
            FuelPrices[fuel] = tmp[1][idx_utc_noloc].values
        elif isinstance(config['default'][fuel], (int, float, complex)):
            logging.warning('No data file found for "' + fuel + '. Using default value ' + str(config['default'][fuel]) + ' EUR')
            FuelPrices[fuel] = pd.Series(config['default'][fuel], index=idx_utc_noloc)
        # Special case for lignite and peat, for backward compatibility
        elif fuel == 'PriceOfLignite':
            logging.warning('No price data found for "' + fuel + '. Using the same value as for Black Coal')
            FuelPrices[fuel] = FuelPrices['PriceOfBlackCoal']
        elif fuel == 'PriceOfPeat':
            logging.warning('No price data found for "' + fuel + '. Using the same value as for biomass')
            FuelPrices[fuel] = FuelPrices['PriceOfBiomass']
        else:
            logging.warning('No data file or default value found for "' + fuel + '. Assuming zero marginal price!')
            FuelPrices[fuel] = pd.Series(0, index=idx_utc_noloc)
    return FuelPrices


#todo NOT ACTIVE
def load_cep_parameters_cart_prod(Plants_merged, countries, expandable_units=['HRD-STUR', 'LIG-STUR', 'NUC-STUR','OIL-STUR', 'GAS-GTUR']):
    logging.info("Capacity Expansion used!")
    all_cost = load_csv('Database/CapacityExpansion/TechsCost.csv') #TODO
    plant_new = load_csv('Database/CapacityExpansion/techs_cap.csv')
    plant_new = plant_new[plant_new.Unit.isin(expandable_units)]

    # create variables (cartesian product of tech x country)
    n_countries = len(countries)
    n_technologies = plant_new.shape[0]
    plant_new = pd.concat([plant_new] * n_countries) # for each zone create new uc
    plant_new['Zone'] = np.repeat(countries, n_technologies) # create zone column
    plant_new['Unit'] = plant_new.apply(lambda x:  x['Zone'] + "-" + x['Unit'], axis=1) # naming
    
    plant_new = plant_new.set_index('Unit', drop=False)

    ## Cost of new technologies
    index = plant_new[['Fuel', 'Technology']].reset_index().set_index('Unit', drop=False)
    Plants_merged = Plants_merged.merge(index[['Fuel', 'Technology']],  how='outer', on = ['Fuel', 'Technology'],
        indicator=True).query('_merge == "left_only"')
    del Plants_merged['_merge']
    Plants_merged = Plants_merged.set_index('Unit', drop=False)
    Plants_merged = Plants_merged.append(plant_new)
    plant_new_cost = all_cost[all_cost.Unit.isin(expandable_units)]
    df_expanded = pd.merge(index, plant_new_cost, on=['Fuel', 'Technology'], how='left')
    return df_expanded


def _fill_missing_cols_by_mean(df, df_mean, on_cols, merge_cols):
    df.loc[:,"_missing"] = df.apply(lambda x: x[["Investment", "FixedCost", "EconomicLifetime"]].isna().any(), axis=1)
    
    df_merged = df.merge(df_mean[on_cols + merge_cols],  how='outer', 
                         on = on_cols,  suffixes = ["", "_mean"], indicator=True)
    df_merged = df_merged[(df_merged["_merge"] == "left_only") | (df_merged["_merge"] == "both")]
    df_only_left = df_merged[(df_merged["_merge"] == "left_only") & (df_merged["_missing"] == True)]
    
    if df_only_left.shape[0] > 0:
        print("No values found for filling missing capacity values of: %s" % ", ".join(df_only_left.Unit.values))
    
    for col in merge_cols:
        merge_col_name = col + "_mean"
        df_merged.loc[:, col] = df_merged[col].fillna(df_merged[merge_col_name])
        del df_merged[merge_col_name]
    
    del df_merged["_missing"]
    del df_merged["_merge"]
    return df_merged


def load_cep_parameters(config, Plants_merged):
    cap_cols = ["Investment", "EconomicLifetime", "FixedCost"]
    for col in cap_cols:
        if col not in Plants_merged:
            Plants_merged[col] = np.nan
    df_cap = Plants_merged[Plants_merged['Extendable'].astype(str) == "x"]
    
    if df_cap.shape[0] > 0: #any extendable power plant technology
        logging.info("Capacity Expansion used!")
        


        if df_cap[cap_cols].isnull().values.any():
            logging.warning("Merging missing values in the columns: %s" % ", ".join(["Investment", "EconomicLifetime", "FixedCost"]))
            
            ## Cost of new technologies
            all_cost = load_csv('Database/CapacityExpansion/techs_cost.csv') #basic cost data
            df_cap = _fill_missing_cols_by_mean(df_cap, all_cost, ["Technology", "Fuel"], ["Investment", "EconomicLifetime", "FixedCost"])
        plant_new = df_cap[:]
        plant_new = plant_new.set_index('Unit', drop=False)
    else:
        plant_new = pd.DataFrame() # empty dataframe -> empty set in optimization model
        
    return plant_new

def cluster_plants(config, plants):

    # Clustering of the plants:
    Plants_merged, mapping = clustering(plants, method=config['SimulationType'])
    # Check clustering:
    check_clustering(plants, Plants_merged)
    return Plants_merged, mapping

def get_unit_based_tables(idx_utc_noloc, config, plants, plants_sto, plants_chp):

    Outages = UnitBasedTable(plants,config['Outages'],idx_utc_noloc,config['countries'],fallbacks=['Unit','Technology'],tablename='Outages')
    AF = UnitBasedTable(plants,config['RenewablesAF'],idx_utc_noloc,config['countries'],fallbacks=['Unit','Technology'],tablename='AvailabilityFactors',default=1,RestrictWarning=commons['tech_renewables'])
    ReservoirLevels = UnitBasedTable(plants_sto,config['ReservoirLevels'],idx_utc_noloc,config['countries'],fallbacks=['Unit','Technology','Zone'],tablename='ReservoirLevels',default=0)
    ReservoirScaledInflows = UnitBasedTable(plants_sto,config['ReservoirScaledInflows'],idx_utc_noloc,config['countries'],fallbacks=['Unit','Technology','Zone'],tablename='ReservoirScaledInflows',default=0)
    HeatDemand = UnitBasedTable(plants_chp,config['HeatDemand'],idx_utc_noloc,config['countries'],fallbacks=['Unit'],tablename='HeatDemand',default=0)
    CostHeatSlack = UnitBasedTable(plants_chp,config['CostHeatSlack'],idx_utc_noloc,config['countries'],fallbacks=['Unit','Zone'],tablename='CostHeatSlack',default=config['default']['CostHeatSlack'])

    # data checks:
    check_AvailabilityFactors(plants, AF)
    check_heat_demand(plants, HeatDemand)
    return Outages, AF, ReservoirLevels, ReservoirScaledInflows, HeatDemand, CostHeatSlack


def load_plants(config):
    # Power plants:
    plants = pd.DataFrame()
    if os.path.isfile(config['PowerPlantData']):
        plants = load_csv(config['PowerPlantData'])
    elif '##' in config['PowerPlantData']:
        for c in config['countries']:
            path = config['PowerPlantData'].replace('##', str(c))
            tmp = load_csv(path)
            plants = plants.append(tmp, ignore_index=True)
    plants = plants[plants['Technology'] != 'Other']
    plants = plants[pd.notnull(plants['PowerCapacity'])]
    plants.index = range(len(plants))

    # Some columns can be in two format (absolute or per unit). If not specified, they are set to zero:
    for key in ['StartUpCost', 'NoLoadCost']:
        if key in plants:
            pass
        elif key+'_pu' in plants:
            plants[key] = plants[key+'_pu'] * plants['PowerCapacity']
        else:
            plants[key] = 0
    # check plant list:
    check_units(config, plants)
    # If not present, add the non-compulsory fields to the units table:
    for key in ['CHPPowerLossFactor','CHPPowerToHeat','CHPType','STOCapacity','STOSelfDischarge','STOMaxChargingPower','STOChargingEfficiency', 'CHPMaxHeat']:
        if key not in plants.columns:
            plants[key] = np.nan


    # Defining the hydro storages:
    plants_sto = plants[[u in commons['tech_storage'] for u in plants['Technology']]]
    # check storage plants:
    check_sto(config, plants_sto)
    plants_chp = plants[[str(x).lower() in commons['types_CHP'] for x in plants['CHPType']]]

    return plants, plants_sto, plants_chp


class DataLoader(object):

    def __rename_plant_columns(self):
            # Renaming the columns to ease the production of parameters:
        self.Plants_merged.rename(columns={ 'StartUpCost': 'CostStartUp',
                                            'RampUpMax': 'RampUpMaximum',
                                            'RampDownMax': 'RampDownMaximum',
                                            'MinUpTime': 'TimeUpMinimum',
                                            'MinDownTime': 'TimeDownMinimum',
                                            'RampingCost': 'CostRampUp',
                                            'STOCapacity': 'StorageCapacity',
                                            'STOMaxChargingPower': 'StorageChargingCapacity',
                                            'STOChargingEfficiency': 'StorageChargingEfficiency',
                                            'STOSelfDischarge': 'StorageSelfDischarge',
                                            'CO2Intensity': 'EmissionRate'}, inplace=True)

    def __init__(self, config):

        # loading/assigning basic data
        self.config = config
        self.idx_utc, self.idx_utc_noloc, self.idx_utc_year_noloc = get_indices(self.config) #todo do i really need you all?
        self.Load, self.PeakLoad = load_loads(self.config, self.idx_utc_noloc, self.idx_utc_year_noloc)
        self.flows, self.NTC, self.Interconnections, self.NTCs, self.Inter_RoW = load_interconnections(config, self.idx_utc_noloc)
        self.LoadShedding = NodeBasedTable(config['LoadShedding'],self.idx_utc_noloc,config['countries'],tablename='LoadShedding',default=config['default']['LoadShedding'])
        self.CostLoadShedding = NodeBasedTable(config['CostLoadShedding'],self.idx_utc_noloc,config['countries'],tablename='CostLoadShedding',default=config['default']['CostLoadShedding'])
        #self.LoadShedding, self.CostLoadShedding = load_load_shedding(self.config, self.idx_utc_noloc)
        self.plants, self.plants_sto, self.plants_chp = load_plants(self.config)
        self.FuelPrices = load_fuel_prices(self.config, self.idx_utc_noloc)
        self.Plants_merged, self.mapping = cluster_plants(self.config, self.plants)
        self.Outages, self.AF, self.ReservoirLevels, self.ReservoirScaledInflows, self.HeatDemand, self.CostHeatSlack = get_unit_based_tables(self.idx_utc_noloc, self.config, self.Plants_merged, self.plants_sto, self.plants_chp)

        # preprocess data
        self.__rename_plant_columns()
        self.__merge_time_series()
        self.__check_plants()
        self.__check_dfs()
        self.__extend_data_with_lookahead()
        self.__prepare_plant_data()

        
        # adding cep
        self.plants_expanded = load_cep_parameters(self.config, self.Plants_merged) #todo
        self.extend_cep_cart_prod = False

        # currently not used:
        # if self.extend_cep_cart_prod:
        #     self.plants_expanded = load_cep_parameters_cart_prod(self.plants_expanded, self.config["countries"])

    def __check_plants(self):
        # data checks:
        check_AvailabilityFactors(self.plants, self.AF)
        check_heat_demand(self.plants, self.HeatDemand)

    def __merge_time_series(self):

        plants = self.plants
        mapping = self.mapping
        # Merging the time series relative to the clustered power plants:
        self.ReservoirScaledInflows = merge_series(plants, self.ReservoirScaledInflows, mapping, method='WeightedAverage', tablename='ScaledInflows')
        self.ReservoirLevels = merge_series(plants, self.ReservoirLevels, mapping, tablename='ReservoirLevels')
        self.Outages = merge_series(plants, self.Outages, mapping, tablename='Outages')
        self.HeatDemand = merge_series(plants, self.HeatDemand, mapping, tablename='HeatDemand',method='Sum')
        self.AF = merge_series(plants, self.AF, mapping, tablename='AvailabilityFactors')
        self.CostHeatSlack = merge_series(plants, self.CostHeatSlack, mapping, tablename='CostHeatSlack')
        #return ReservoirScaledInflows, self.ReservoirLevels, Outages, HeatDemand, AF, CostHeatSlack


    def __check_dfs(self):
        # checking data
        idx_utc_noloc = self.idx_utc_noloc

        check_df(self.Load, StartDate=idx_utc_noloc[0], StopDate=idx_utc_noloc[-1], name='Load')
        check_df(self.AF, StartDate=idx_utc_noloc[0], StopDate=idx_utc_noloc[-1], name='AF')
        check_df(self.Outages, StartDate=idx_utc_noloc[0], StopDate=idx_utc_noloc[-1], name='Outages')
        check_df(self.Inter_RoW, StartDate=idx_utc_noloc[0], StopDate=idx_utc_noloc[-1], name='Inter_RoW')
        check_df(self.FuelPrices, StartDate=idx_utc_noloc[0], StopDate=idx_utc_noloc[-1], name='FuelPrices')
        check_df(self.NTCs, StartDate=idx_utc_noloc[0], StopDate=idx_utc_noloc[-1], name='NTCs')
        check_df(self.ReservoirLevels, StartDate=idx_utc_noloc[0], StopDate=idx_utc_noloc[-1], name='ReservoirLevels')
        check_df(self.ReservoirScaledInflows, StartDate=idx_utc_noloc[0], StopDate=idx_utc_noloc[-1], name='ReservoirScaledInflows')
        check_df(self.HeatDemand, StartDate=idx_utc_noloc[0], StopDate=idx_utc_noloc[-1], name='HeatDemand')
        check_df(self.CostHeatSlack, StartDate=idx_utc_noloc[0], StopDate=idx_utc_noloc[-1], name='CostHeatSlack')
        check_df(self.LoadShedding, StartDate=idx_utc_noloc[0], StopDate=idx_utc_noloc[-1], name='LoadShedding')
        check_df(self.CostLoadShedding, StartDate=idx_utc_noloc[0], StopDate=idx_utc_noloc[-1], name='CostLoadShedding')


    def __extend_data_with_lookahead(self):

        # Extending the data to include the look-ahead period (with constant values assumed)
        enddate_long = self.idx_utc_noloc[-1] + dt.timedelta(days=self.config['LookAhead'])
        idx_long = pd.DatetimeIndex(pd.date_range(start=self.idx_utc_noloc[0], end=enddate_long, freq=commons['TimeStep']))
        Nhours_long = len(idx_long)

        # re-indexing with the longer index and filling possibly missing data at the beginning and at the end::
        self.Load = self.Load.reindex(idx_long, method='nearest').fillna(method='bfill')
        self.AF = self.AF.reindex(idx_long, method='nearest').fillna(method='bfill')
        self.Inter_RoW = self.Inter_RoW.reindex(idx_long, method='nearest').fillna(method='bfill')
        self.NTCs = self.NTCs.reindex(idx_long, method='nearest').fillna(method='bfill')
        self.FuelPrices = self.FuelPrices.reindex(idx_long, method='nearest').fillna(method='bfill')
        self.Load = self.Load.reindex(idx_long, method='nearest').fillna(method='bfill')
        self.Outages = self.Outages.reindex(idx_long, method='nearest').fillna(method='bfill')
        self.ReservoirLevels = self.ReservoirLevels.reindex(idx_long, method='nearest').fillna(method='bfill')
        self.ReservoirScaledInflows = self.ReservoirScaledInflows.reindex(idx_long, method='nearest').fillna(method='bfill')
        self.LoadShedding = self.LoadShedding.reindex(idx_long, method='nearest').fillna(method='bfill')
        self.CostLoadShedding = self.CostLoadShedding.reindex(idx_long, method='nearest').fillna(method='bfill')
    #    for tr in Renewables:
    #        Renewables[tr] = Renewables[tr].reindex(idx_long, method='nearest').fillna(method='bfill')

    def __prepare_plant_data(self):

        # using references for prettier coding
        config = self.config
        plants = self.plants
        Plants_merged = self.Plants_merged

        for key in ['TimeUpMinimum','TimeDownMinimum']:
            if any([not x.is_integer() for x in Plants_merged[key].fillna(0).values.astype('float')]):
                logging.warning(key + ' in the power plant data has been rounded to the nearest integer value')
                Plants_merged.loc[:,key] = Plants_merged[key].fillna(0).values.astype('int32')

            if not len(Plants_merged.index.unique()) == len(Plants_merged):
                # Very unlikely case:
                logging.error('plant indexes not unique!')
                sys.exit(1)

            # Apply scaling factors:
            if config['modifiers']['Solar'] != 1:
                logging.info('Scaling Solar Capacity by a factor ' + str(config['modifiers']['Solar']))
                for u in Plants_merged.index:
                    if Plants_merged.Technology[u] == 'PHOT':
                        Plants_merged.loc[u, 'PowerCapacity'] = Plants_merged.loc[u, 'PowerCapacity'] * config['modifiers']['Solar']
            if config['modifiers']['Wind'] != 1:
                logging.info('Scaling Wind Capacity by a factor ' + str(config['modifiers']['Wind']))
                for u in Plants_merged.index:
                    if Plants_merged.Technology[u] == 'WTON' or Plants_merged.Technology[u] == 'WTOF':
                        Plants_merged.loc[u, 'PowerCapacity'] = Plants_merged.loc[u, 'PowerCapacity'] * config['modifiers']['Wind']
            if config['modifiers']['Storage'] != 1:
                logging.info('Scaling Storage Power and Capacity by a factor ' + str(config['modifiers']['Storage']))
                for u in Plants_merged.index:
                    if isStorage(Plants_merged.Technology[u]):
                        Plants_merged.loc[u, 'PowerCapacity'] = Plants_merged.loc[u, 'PowerCapacity'] * config['modifiers']['Storage']
                        Plants_merged.loc[u, 'StorageCapacity'] = Plants_merged.loc[u, 'StorageCapacity'] * config['modifiers']['Storage']
                        Plants_merged.loc[u, 'StorageChargingCapacity'] = Plants_merged.loc[u, 'StorageChargingCapacity'] * config['modifiers']['Storage']

            # Defining the hydro storages:
            self.plants_sto = Plants_merged[[u in commons['tech_storage'] for u in Plants_merged['Technology']]]
            # check storage plants:
            check_sto(config, self.plants_sto, raw_data=False)
            # Defining the CHPs:
            self.plants_chp = Plants_merged[[x.lower() in commons['types_CHP'] for x in Plants_merged['CHPType']]].copy()
            # check chp plants:
            check_chp(config, self.plants_chp)
            # For all the chp plants correct the PowerCapacity, which is defined in cogeneration mode in the inputs and in power generation model in the optimization model
            for u in self.plants_chp.index:
                PowerCapacity = self.plants_chp.loc[u, 'PowerCapacity']

                if self.plants_chp.loc[u,'CHPType'].lower() == 'p2h':
                    PurePowerCapacity = PowerCapacity
                else:
                    if pd.isnull(self.plants_chp.loc[u,'CHPMaxHeat']):  # If maximum heat is not defined, then it is defined as the intersection between two lines
                        MaxHeat = PowerCapacity / self.plants_chp.loc[u,'CHPPowerToHeat']
                        self.plants_chp.loc[u, 'CHPMaxHeat'] = 'inf'
                    else:
                        MaxHeat = self.plants_chp.loc[u, 'CHPMaxHeat']
                    PurePowerCapacity = PowerCapacity + self.plants_chp.loc[u,'CHPPowerLossFactor'] * MaxHeat
                Plants_merged.loc[u, 'PartLoadMin'] = Plants_merged.loc[u,'PartLoadMin'] * PowerCapacity / PurePowerCapacity  # FIXME: Is this correct?
                Plants_merged.loc[u, 'PowerCapacity'] = PurePowerCapacity


            # Get the hydro time series corresponding to the original plant list: #FIXME Unused variable ?
            #StorageFormerIndexes = [s for s in plants.index if
            #                        plants['Technology'][s] in commons['tech_storage']]


            # Same with the CHPs:
            # Get the heat demand time series corresponding to the original plant list:
            CHPFormerIndexes = [s for s in plants.index if
                                    plants['CHPType'][s] in commons['types_CHP']]
            for s in CHPFormerIndexes:  # for all the old plant indexes
                # get the old plant name corresponding to s:
                oldname = plants['Unit'][s]
                # newname = mapping['NewIndex'][s] #FIXME Unused variable ?
                if oldname not in self.HeatDemand:
                    logging.warning('No heat demand profile found for CHP plant "' + str(oldname) + '". Assuming zero')
                    self.HeatDemand[oldname] = 0
                if oldname not in self.CostHeatSlack:
                    logging.warning('No heat cost profile found for CHP plant "' + str(oldname) + '". Assuming zero')
                    self.CostHeatSlack[oldname] = 0

            # # merge the outages:
            # for i in plants.index:  # for all the old plant indexes
            #     # get the old plant name corresponding to s:
            #     print(i)
            #     print(self.mapping['NewIndex'])
            #     oldname = plants['Unit'][i]
            #     newname = self.mapping['NewIndex'][i]
