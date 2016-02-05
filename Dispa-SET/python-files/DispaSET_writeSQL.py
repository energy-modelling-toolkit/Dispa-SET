# -*- coding: utf-8 -*-
__author__ = 'Sylvain Quoilin <sylvain.quoilin@ec.europa.eu>'

# This scripts is used to read raw data from excel files and write it in the Dispaset data template format in a mysql database


import DispaSET_io_data
from DispaSET_tools import *
import numpy as np
from sqlalchemy import create_engine

path_data_csv = '/home/sylvain/Dropbox/power_system_modelling-master2/Data/BE/data_2012-2013_csv'
path_data_pandas = '/home/sylvain/Dropbox/power_system_modelling-master2/Data/BE/data_2012-2013_pandas'


# Connect to database: 
engine = create_engine('mysql://dispaset:aaa@localhost/dispaset')


# Load the timestamps in string format:
dates_pd = load_csv_to_pd(path_data_csv,'daysofyear.csv',path_data_pandas,'daysofyear')
load_data_pd = load_csv_to_pd(path_data_csv,'load_data.csv',path_data_pandas,'load_data')
Nrows = len(dates_pd)



firstday = pd.datetime(load_data_pd.year[0],load_data_pd.month[0],load_data_pd.day[0],load_data_pd.hour[0],load_data_pd['min'][0])
dates= dates_pd.Date.tolist()
dates_index = [pd.datetime(load_data_pd.year[i],load_data_pd.month[i],load_data_pd.day[i],load_data_pd.hour[i],load_data_pd['min'][i]) for i in range(Nrows)]

# The dates must be register in the UTC timezone to avoid duplicates in the sql database. This should be improved
dates_index = pd.DatetimeIndex(start=firstday,periods=Nrows,freq='15Min',tz='UTC')
dates_index_hour = pd.DatetimeIndex(start=firstday,periods=Nrows/4,freq='1h',tz='UTC')

# Load plant data:
plants = load_csv_to_pd(path_data_csv,'Production_park_data_v2.1.csv',path_data_pandas,'Production_park_data_v2.1').set_index('Index')

# Load demand data:
load_data_pd.index = dates_index
load_data= np.array(load_data_pd['Vertical load [kW]'])

# Outage data:
outages_pd = load_csv_to_pd(path_data_csv,'outages.csv',path_data_pandas,'outages')
outages= np.array(outages_pd)[:,1:]
outages_pd2 =pd.DataFrame(np.repeat(outages.T,24,axis=0),index=dates_index_hour)

# Load load_forecast
load_forecast_pd = load_csv_to_pd(path_data_csv,'load_forecast.csv',path_data_pandas,'load_forecast')
load_forecast= np.array(load_forecast_pd)
load_forecast_pd.index = dates_index

# Load interconnections
inter_fr_pd = load_csv_to_pd(path_data_csv,'interconnections_fr.csv',path_data_pandas,'interconnections_fr')
inter_nl_pd = load_csv_to_pd(path_data_csv,'interconnections_nl.csv',path_data_pandas,'interconnections_nl')
interconnections = {'fr':np.array(inter_fr_pd['Power exchange']),'nl':np.array(inter_nl_pd['Power exchange'])}
inter_fr_pd.index = dates_index
inter_nl_pd.index = dates_index

# Load solar data:
solar_data_pd = load_csv_to_pd(path_data_csv,'solar_data.csv',path_data_pandas,'solar_data')
solar_data= np.array(solar_data_pd)
solar_data_pd.index =dates_index

# Load wind data:
wind_data_pd = load_csv_to_pd(path_data_csv,'WindData.csv',path_data_pandas,'WindData')
wind_data= np.array(wind_data_pd)
wind_data_pd.index = dates_index

# Fuel price data from IEA, monthly, in eur/MWhth
# From october 2012 to september 2013
price = {'coal_month':np.array([12.2430352051,12.2430352051,12.2430352051,10.2058479078,10.2058479078,10.2058479078,10.5107415473,10.5107415473,10.5107415473,9.8076194807,9.8076194807,9.8076194807])}
price['oil_month'] = np.array([73.7754656768,73.7754656768,73.7754656768,69.6732200066,69.6732200066,69.6732200066,67.2718323463,67.2718323463,67.2718323463,79.5365975377,79.5365975377,79.5365975377])
price['gas_month'] = np.array([28.85,28.85,28.85,30.05,30.05,30.05,29.9,29.9,29.9,30.2,30.2,30.2])
# Transforming the fuel price (monthly data) into hourly data:
x = np.array(range(1,13))
xvals = np.linspace(1,12,8760)
price_h = {}
price_h['Oil'] = np.interp(xvals,x,price['oil_month'])
price_h['Gas'] = np.interp(xvals,x,price['gas_month'])
price_h['Coal'] = np.interp(xvals,x,price['coal_month'])



plants.to_sql('plants_be_2012',engine,if_exists='replace')

load = pd.DataFrame({'vertical': load_data_pd['Vertical load [kW]'], 'forecast':load_forecast_pd['Vertical load [kW]']},index = dates_index)
load.to_sql('load_be',engine,if_exists='replace')

inter = pd.DataFrame({'fr-be':inter_fr_pd['Power exchange'], 'nl-be':inter_nl_pd['Power exchange']},index = dates_index)
inter.to_sql('interconnections_be',engine,if_exists='replace')

solar = pd.DataFrame({'Day-Ahead forecast':solar_data_pd['Day-Ahead forecast [MW]'], 'Intraday forecast':solar_data_pd['Intraday forecast [MW]'], 'Measurement':solar_data_pd['Corrected Upscaled Measurement [MW]'],'Capacity':solar_data_pd['Monitored Capacity [MWp]']},index = dates_index)
solar['Measurement'][0:4451] = solar_data_pd['Real-time Upscaled Measurement [MW]'][0:4451]    # missing data => using the non corrected values
solar.to_sql('solar_be',engine,if_exists='replace')

wind = pd.DataFrame({'Day-Ahead forecast':wind_data_pd['Total Forecast'],'Measurement':wind_data_pd['Total Measured & Upscaled'],'Capacity':wind_data_pd['Monitored Capacity']},index = dates_index)
wind.to_sql('wind_be',engine,if_exists='replace')

outages_pd2.to_sql('outages_be',engine,if_exists='replace')

prices = pd.DataFrame(price_h,index=dates_index_hour)
prices.to_sql('fuelprice_be',engine,if_exists='replace')

