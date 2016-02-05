# analysis of dispaset results files

__author__ = 'Sylvain Quoilin (sylvain.quoilin@ec.europa.eu)'

import matplotlib.pyplot as plt
import numpy as np
import pickle
from load_dispaset_vars import *
from scipy import signal

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning,
                        module="pandas", lineno=276)


result_id = 'sim1'


#path_results = '/home/sylvain/Dropbox/power_system_modelling-master2/Data/BE'
path_results = '/home/sylvain/JRC/DispaSET/simulations papier EEM/' + result_id + '/Data/BE'

#load summary file of the results:
try:
    all_results = pickle.load(open('summary.p','rb'))
except:
    print 'Could not load summary file, creating a new one'
    all_results = {}

# Selecting the data range to be considered:
# 0: Winter day
# 1: Winter night
# 2: winter peak
# 3: spring day
# 4: srping night
# 5: spring peak
# 6: summer day
# 7: summer night
# 8: summer peak
# 9: fall day
# 10: fall night
# 11: fall peak
pos = -1




# Load inputs:
demand_pd = load_xl_to_pd(path_results,'ResultTest.xlsx','Demand',path_results,'Demand.p')
techno_pd = load_xl_to_pd(path_results,'ResultTest.xlsx','Technology',path_results,'Technology.p',header=None)
powercapacity_pd = load_xl_to_pd(path_results,'ResultTest.xlsx','PowerCapacity',path_results,'PowerCapacity.p',header=None)
AF_pd = load_xl_to_pd(path_results,'ResultTest.xlsx','AvailabilityFactor',path_results,'AvailabilityFactor.p').transpose()
OF_pd = load_xl_to_pd(path_results,'ResultTest.xlsx','OutageFactor',path_results,'OutageFactor.p').transpose()
powerminstable_pd = load_xl_to_pd(path_results,'ResultTest.xlsx','PartLoadMin',path_results,'PartLoadMin.p',header=None)
storagecapacity_pd = load_xl_to_pd(path_results,'ResultTest.xlsx','StorageCapacity',path_results,'StorageCapacity.p',header=None)



load = demand_pd['DA']



techno = techno_pd[1].values
PowerCapacity = powercapacity_pd[1].values
PowerMinStable = powerminstable_pd[1].values
Nhours = len(demand_pd['DA'])
load = demand_pd['DA'].values
Nunits = len(powercapacity_pd)

# Load outputs:
results_power_pd = load_xl_to_pd(path_results,'ResultTest.xlsx','Power',path_results,'results_power.p')
results_power = np.array(results_power_pd)
results_pumping_pd = load_xl_to_pd(path_results,'ResultTest.xlsx','StorageInput',path_results,'results_pumping.p')
results_pumping = np.array(results_pumping_pd['[45, 46, 49] - HUSTO - WA'])
results_curt_pd = load_xl_to_pd(path_results,'ResultTest.xlsx','CurtailedPower',path_results,'results_curt.p')
results_curt = np.zeros(Nhours)
results_curt[np.array(results_curt_pd.index) -1] = results_curt_pd['BE']

# Setting the season boundaries (see above):
# In Times, starting with winter
bounds_seas = [range(1945,4105), range(4105, 7032), range(7032, 8520), range(8520,8760)+range(1945)]
lengths = []
length_inday = [11, 12, 1]*4    # Day, night, peak
for i in range(12):
    s = i/4
    lengths.append(len(bounds_seas[s])/24*length_inday[i])

# for each day, find the maximum load and define de timeslices:
ranges = []
for bounds in bounds_seas:
    Nrows = len(bounds)
    Ndays = Nrows/24
    day_all = []
    max_all = []
    night_all = []
    for d in range(Ndays):
        day_range = range(d*24,(d+1)*24)
        night_range = range(d*24,d*24+9)+range(d*24+21,(d+1)*24)
        pos_max = load[range(d*24+9,d*24+21)].argmax() + day_range[0]
        day_range11 = range(d*24+9,pos_max) + range(pos_max+1,d*24+21)
        day_all = day_all + [bounds[x] for x in day_range11]
        max_all = max_all + [bounds[pos_max]]
        night_all = night_all + [bounds[x] for x in night_range]
    ranges.append(day_all)
    ranges.append(night_all)
    ranges.append(max_all)

# Find VRE entries to get the capacity factors (if several lines, only the first one):
idx_wind = 44               # !! to be changed if plant definition changes!
idx_sol = 45                # !! to be changed if plant definition changes!
idx_sto = 0
AF_wind = AF_pd['[109, 110, 111, 112, ...  118, 119] - WT - WI'].values
AF_sol = AF_pd['[120] - AllSol'].values
OF_wind = OF_pd['[109, 110, 111, 112, ...  118, 119] - WT - WI'].values
OF_sol = OF_pd['[120] - AllSol'].values

Pcap_sol = PowerCapacity[idx_sol]
Pcap_wind = PowerCapacity[idx_wind]
P_sto = PowerCapacity[idx_sto]

cap_sto = np.sum(storagecapacity_pd[1])

Psol = Pcap_sol * AF_sol * (1 - OF_sol)
Pwind = Pcap_wind * AF_wind * (1 - OF_wind)

# Aggregation (careful because GAMS changes the order):
t = ['HUSTO','CCGT','CL','D','HU','IS','NU','TJ','WKK','WT','PV','GT']
P={}
CF={}
P_hist = {}
for tech in t:
    P[tech] = 0
    P_hist[tech] = 0
    for j in range(Nunits):
        if techno[j]==tech:
             P[tech]=P[tech]+PowerCapacity[j]
             P_hist[tech] = P_hist[tech] + results_power[j,:]
P_tot = P['CCGT'] + P['D'] + P['GT']+P['IS'] + P['WKK'] + P['CL'] + P['NU']
P_dispatch = P['CCGT'] + P['D'] + P['GT']
P_baseload = P['IS'] + P['WKK'] + P['CL'] + P['NU']
P_hist_tot = P_hist['CCGT'] + P_hist['D'] + P_hist['GT'] + P_hist['IS'] + P_hist['WKK'] + P_hist['CL'] + P_hist['NU']  + P_hist['HU']  + P_hist['TJ'] + P['HUSTO']

P_hist_tot = np.sum(results_power[0:-2,:],axis=0)

residu = load - P_hist_tot - Psol - Pwind + results_pumping + results_curt
#plt.plot(residu)

Pvre = Psol + Pwind
Ppump_sol = results_pumping * Psol/Pvre
Ppump_wind = results_pumping * Pwind/Pvre
Pcurt = results_curt
Pbase = Pvre - Pcurt - Ppump_sol - Ppump_wind

if pos >=0 and pos <12:
    rng = ranges[pos]
else:
    rng = range(Nhours)

result_id2 = result_id + ' - ' + str(pos)
if not all_results.has_key(result_id2):
    all_results[result_id2] = {}

all_results[result_id2]['E_vre'] = np.sum(Pvre[rng])
all_results[result_id2]['E_pp'] = np.sum(results_pumping[rng])
all_results[result_id2]['E_pp_sol'] = np.sum(Ppump_sol[rng])
all_results[result_id2]['E_pp_wind'] = np.sum(Ppump_wind[rng])
all_results[result_id2]['E_curt'] = np.sum(results_curt[rng])
all_results[result_id2]['E_load'] = np.sum(load[rng])
all_results[result_id2]['load_min'] = np.min(load[rng])
all_results[result_id2]['load_max'] = np.max(load[rng])
all_results[result_id2]['Pcap_sol'] = Pcap_sol
all_results[result_id2]['Pcap_wind'] = Pcap_wind
all_results[result_id2]['P_dispatch'] = P_dispatch
all_results[result_id2]['P_baseload'] = P_baseload
all_results[result_id2]['P_sto'] = P_sto
all_results[result_id2]['cap_sto'] = cap_sto
all_results[result_id2]['Nhours_wind'] = len(np.nonzero(Pwind)[0])
all_results[result_id2]['Nhours_sol'] = len(np.nonzero(Psol)[0])
all_results[result_id2]['timeslice'] = pos



####################################### Plotting: ######################################################

idx = [rng[i] for i in np.argsort(-Pvre[rng])]

idx2 = idx
x = range(len(idx2))
P1 = Pbase[idx2]
P2 = Pbase[idx2]+Pcurt[idx2]
P3 = Pbase[idx2]+Pcurt[idx2]+Ppump_wind[idx2]
P4 = Pvre[idx2]
demand = np.mean(load[idx2])
D = np.ones(len(idx2)) * demand

print np.sum(Pcurt)
print np.sum(Pcurt[25:])


'''
plt.plot(Pvre[idx])
plt.fill_between(x,P3,P4,linewidth=0,color='#000066')
plt.fill_between(x,0,P1,linewidth=0,color='grey',alpha=0.5)
plt.fill_between(x,P2,P3,linewidth=0,color='#001400')
plt.fill_between(x,P1,P2,linewidth=0,color='red',alpha=0.7)
plt.plot(D)
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,15000))
plt.title('Non filtered: base, curtailment, PV pumping, wind pumping')
plt.show()
'''


b, a = signal.butter(4, 0.2, 'low')
P1_smooth = P1[0] + signal.lfilter(b, a, P1-P1[0])
P2_smooth = P2[0] + signal.lfilter(b, a, P2-P2[0])
P3_smooth = P3[0] + signal.lfilter(b, a, P3-P3[0])

'''
plt.fill_between(x,P3_smooth,P4,linewidth=0,color='blue')
plt.fill_between(x,0,P1_smooth,linewidth=0,color='grey',alpha=0.5)
plt.fill_between(x,P2_smooth,P3_smooth,linewidth=0,color='green')
plt.fill_between(x,P1_smooth,P2_smooth,linewidth=0,color='red',alpha=0.7)
plt.plot(D)
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,15000))
plt.title('Filtered: base, curtailment, PV pumping, wind pumping')
plt.show()
'''

# Same thing, but only with Psol:
Pvre = Psol
Ppump_sol = results_pumping
Pcurt = results_curt
Pbase = Pvre - Pcurt - Ppump_sol

idx = [rng[i] for i in np.argsort(-Pvre[rng])]

idx2 = idx
x = range(len(idx2))
P1 = Pbase[idx2]
P2 = Pbase[idx2]+Pcurt[idx2]
P3 = P2
P4 = Pvre[idx2]
demand = np.mean(load[idx2]-Pwind[idx2])
D = np.ones(len(idx2)) * demand

P1_smooth = P1[0] + signal.lfilter(b, a, P1-P1[0])
P2_smooth = P2[0] + signal.lfilter(b, a, P2-P2[0])
P3_smooth = P3[0] + signal.lfilter(b, a, P3-P3[0])

'''
plt.fill_between(x,P3_smooth,P4,linewidth=0,color='blue')
plt.fill_between(x,0,P1_smooth,linewidth=0,color='grey',alpha=0.5)
plt.fill_between(x,P2_smooth,P3_smooth,linewidth=0,color='green')
plt.fill_between(x,P1_smooth,P2_smooth,linewidth=0,color='red',alpha=0.7)
plt.plot(D)
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,0,8000))
plt.title('Filtered: base PV, curtailment, Pumping')
plt.show()
'''

# Same thing, residual load:
Pres = load - Psol - Pwind

idx = [rng[i] for i in np.argsort(Pres[rng])]

idx2 = idx
x = range(len(idx2))
P1 = Pres[idx2]
P2 = Pres[idx2] + results_curt[idx2]
P3 = Pres[idx2] + results_curt[idx2]  + results_pumping[idx2]
P4 = load[idx2]
demand = np.mean(load[idx2])
D = np.ones(len(idx2)) * demand

P1_smooth = P1[0] + signal.lfilter(b, a, P1-P1[0])
P2_smooth = P2[0] + signal.lfilter(b, a, P2-P2[0])
P3_smooth = P3[0] + signal.lfilter(b, a, P3-P3[0])
P4_smooth = P4[0] + signal.lfilter(b, a, P4-P4[0])

plt.fill_between(x,P3_smooth,P4_smooth,linewidth=0,color='blue')
plt.fill_between(x,0,P1_smooth,linewidth=0,color='grey',alpha=0.5)
plt.fill_between(x,P2_smooth,P3_smooth,linewidth=0,color='green')
plt.fill_between(x,P1_smooth,P2_smooth,linewidth=0,color='red',alpha=0.7)
plt.plot(D,color='black')
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,np.minimum(x2,8760),-6500,11500))
#plt.title('Filtered: Residual load, curtailment, pumping, VRE')
# Creating the legend items for the fill between:
p1 = plt.Rectangle((0, 0), 1, 1, fc="red",alpha=0.7)
p2 = plt.Rectangle((0, 0), 1, 1, fc="green")
p3 = plt.Rectangle((0, 0), 1, 1, fc="blue")
plt.legend([p1, p2, p3], ['Curtailment', 'Storage','VRE'],'lower right')
plt.ylabel('Power (MW)')
plt.xlabel('Time (h)')
plt.show()

pickle.dump(all_results,open( "summary.p", "wb" ))
