.. _model:

Model Description
=================

The model is expressed as a MILP problem. Continuous variables include the individual unit dispatched power, the shedded load and the curtailed power generation. The binary variables are the commitment status of each unit. The main model features can be summarized as follows:


Variables
^^^^^^^^^

Sets
----

.. table:: 

	======= =================================================================================
	Name	Description
	======= =================================================================================
	f	Fuel types
	h	Hours
	i	Time step in the solving loop
	l	Transmission lines between nodes
	mk	{DA: Day-Ahead, 2U: Reserve up, 2D: Reserve Down}
	n	Zones within each country (currently one zone, or node, per country)
	p	Pollutants
	t	Power generation technologies
	tr	Renewable power generation technologies
	u	Units
	s	Storage units (including hydro reservoirs)
	i	Subset of all simulated hours
	======= =================================================================================

Parameters
----------

.. table:: 

	======================================= ======= =============================================================
	Name					Units	Description
	======================================= ======= =============================================================
	AvailabilityFactor(u,i)			%	Percentage of nominal capacity available
	CommittedInitial(u)			n.a.	Initial commitment status
	CostFixed(u)		 		€/h	Fixed costs
	CostLoadShedding(n,h)			€/MWh	Shedding costs
	CostRampDown(u)				€/MW	Ramp-down costs
	CostRampUp(u)				€/MW	Ramp-up costs
	CostShutDown(u)				€/h	Shut-down costs
	CostStartUp(u)				€/h	Start-up costs
	CostVariableH(u,i)			€/MWh	Variable costs
	Curtailment(n)				n.a.	Curtailment {binary: 1 allowed}
	Demand(mk,n,i)				MW	Hourly demand in each zone
	Efficiency(u)				%	Power plant efficiency
	EmissionMaximum(n,p)			€/tP	Emission limit per zone for pollutant p
	EmissionRate(u,p)			tP/MW	Emission rate of pollutant p from unit u
	FlexibilityDown(u)			MW/h	Available fast shut-down ramping capacity
	FlexibilityUp(u)			MW/h	Available fast start-up ramping capacity
	Fuel(u,f)				n.a.	Fuel type used by unit u {binary: 1 u uses f}
	LineNode(l,n)				n.a.	Line-zone incidence matrix {-1,+1}
	LoadMaximum(u,h)			%	Maximum load for each unit
	LoadShedding(n,h)			MW	Load that may be shed per zone in 1 hour
	Location(u,n)				n.a.	Location {binary: 1 u located in n}
	OutageFactor(u,h)			%	Outage factor (100 % = full outage) per hour
	PartLoadMin(u)				%	Percentage of minimum nominal capacity
	PermitPrice(p)				€/tP	Permit price for pollutant p
	PowerCapacity(u)			MW	Installed capacity
	PowerInitial(u)				MW	Power output before initial period
	PowerMinStable(u)			MW	Minimum power for stable generation
	PowerMustRun(u)				MW	Minimum power output
	PriceTransmission(l,h)			€/MWh	Price of transmission between zones
	RampDownMaximum(u)			MW/h	Ramp down limit
	RampShutDownMaximum(u)			MW/h	Shut-down ramp limit
	RampStartUpMaximum(u)			MW/h	Start-up ramp limit
	RampUpMaximum(u)			MW/h	Ramp up limit
	Reserve(t)				n.a.	Reserve provider {binary}
	StorageCapacity(s)			MWh 	Storage capacity (reservoirs)
	StorageChargingCapacity(s)		MW	Maximum charging capacity
	StorageChargingEfficiency(s)		%	Charging efficiency
	StorageDischargeEfficiency(s)		%	Discharge efficiency
	StorageInflow(s,h)			MWh 	Storage inflows
	StorageInitial(s)			MWh 	Storage level before initial period
	StorageMinimum(s)			MWh 	Minimum storage level
	StorageOutflow(s,h)			MWh	Storage outflows (spills) 
	Technology(u,t)				n.a.	Technology type {binary: 1: u belongs to t}
	TimeDownInitial(u)			h	Hours down before initial period
	TimeDownLeftInitial(u)			h	Time down remaining at initial time
	TimeDownLeftJustStopped(u,i)		h	Time down remaining if started at time i
	TimeDownMinimum(u)			h	Minimum down time
	TimeDown(u,h)				h	Number of hours down
	TimeUpInitial(u)			h	Number of hours up before initial period
	TimeUpLeftInitial(u)			h	Time up remaining at initial time
	TimeUpLeftJustStarted(u,i)		h	Time up remaining if started at time i
	TimeUpMinimum(u)			h	Minimum up time
	TimeUp(u,h)				h	Number of hours up
	VOLL)					€/MWh	Value of lost load
        ======================================= ======= =============================================================

Optimization Variables
----------------------

.. table:: 

	======================= ======= =============================================================
	Name			Units	Description
	======================= ======= =============================================================
	Committed(u,h)		n.a.	Unit committed at hour h {1,0}
	CostStartUpH(u,h)	[EUR]	Cost of starting up
	CostShutDownH(u,h)	[EUR]	cost of shutting down
	CostRampUpH(u,h)	[EUR]	Ramping cost
	CostRampDownH(u,h)	[EUR]	Ramping cost
	CurtailedPower(n,h)	[MW]	Curtailed power at node n
	Flow(l,h)		[MW]	Flow through lines
	MaxRamp2U(u,h)		[MW/h]	Maximum 15-min Ramp-up capbility
	MaxRamp2D(u,h)		[MW/h]	Maximum 15-min Ramp-down capbility
	Power(u,h)		[MW]	Power output
	PowerMaximum(u,h)	[MW]	Power output
	PowerMinimum(u,h)	[MW]	Power output
	ShedLoad(n,h)		[MW]	Shed load
	StorageInput(s,h)	[MWh]	Charging input for storage units
	StorageLevel(s,h)	[MWh]	Storage level of charge
	SystemCostD		[EUR]	Total system cost  for one optimization period
	LostLoadMaxPower(n,h)	[MW]	Deficit in terms of maximum power
	LostLoadRampUp(u,h)	[MW]	Deficit in terms of ramping up for each plant
	LostLoadRampDown(u,h)	[MW]	Deficit in terms of ramping down
	LostLoadMinPower(n,h)	[MW]	Power exceeding the demand
	LostLoadReserve2U(n,h)	[MW]	Deficit in reserve up
	======================= ======= =============================================================




Equations
^^^^^^^^^

The aim of this model is to represent with a high level of detail the short-term operation of large-scale power systems solving the so-called unit commitment problem. To that aim we consider that the system is managed by a central operator with full information on the technical and economic data of the generation units, the demands in each node, and the transmission network.

The unit commitment problem considered in this report is a simplified instance of the problem faced by the operator in charge of clearing the competitive bids of the participants into a wholesale day-ahead power market. In the present formulation the demand side is an aggregated input for each node, while the transmission network is modelled as a transport problem between the nodes (that is, the problem is network-constrained but the model does not include the calculation of the optimal power flows).

The unit commitment problem consists of two parts: i) scheduling the start-up, operation, and shut down of the available generation units, and ii) allocating (for each period of the simulation horizon of the model) the total power demand among the available generation units in such a way that the overall power system costs is minimized. The first part of the problem, the unit scheduling during several periods of time, requires the use of binary variables in order to represent the start-up and shut down decisions, as well as the consideration of constraints linking the commitment status of the units in different periods. The second part of the problem is the so-called economic dispatch problem, which determines the continuous output of each and every generation unit in the system. Therefore, given all the features of the problem mentioned above, it can be naturally formulated as a mixed-integer linear program (MILP). The formulation of the model presented in this report is based upon publicly available modelling approaches [ CITATION Arr00 {\textbackslash}l 1033  {\textbackslash}m Car06 {\textbackslash}m Mor13]. Since our goal is to model a large European interconnected power system, we have implemented a so-called tight and compact formulation, in order to simultaneously reduce the region where the solver searches for the solution and increase the speed at which the solver carries out that search. Tightness refers to the distance between the relaxed and integer solutions of the MILP and therefore defines the search space to be explored by the solver, while compactness is related to the amount of data to be processed by the solver and thus determines the speed at which the solver searches for the optimum. Usually tightness is increased by adding new constraints, but that also increases the size of the problem (decreases compactness), so both goals contradict each other and a trade-off must be found.

Objective function
------------------

The goal of the unit commitment problem is to minimize the total power system costs (expressed in {\texteuro} in equation ), which are defined as the sum of different cost items, namely: start-up and shut-down, fixed, variable, ramping, transmission-related and load shedding (voluntary and involuntary) costs.

.. math::
	 i=1:

	 \mathit{CostStartUp}_{u,i} \geq \mathit{CostStartUp}_u \cdot \left(\mathit{Committed}_{u,i}-\mathit{CommittedInitial}_u\right)

	 \mathit{CostShutDown}_{u,i} \geq \mathit{CostShutDown}_u \cdot (\mathit{CommittedInitial}_u-\mathit{Committed}_{u,i})

	 i>1:

	 \mathit{CostStartUp}_{u,i} \geq \mathit{CostStartUp}_u \cdot \left(\mathit{Committed}_{u,i}-\mathit{Committed}_{u,i-1}\right)

	 \mathit{CostShutDown}_{u,i} \geq \mathit{CostShutDown}_u \cdot (\mathit{Committed}_{u,i-1}-\mathit{Committed}_{u,i})
	 

In the previous equation, as in some of the following, a distinction is made between the equation for the first and subsequent periods. The equation for the first period takes into account the commitment status of the unit before the beginning of the simulation, which is part of the information fed into the model.

Ramping costs are computed in the same manner:

.. math:: 
	 i=1:

	 \mathit{CostRampUp}_{u,i} \geq \mathit{CostRampUp}_u \cdot \left(\mathit{Power}_{u,i}-\mathit{PowerInitial}_u\right)

	 \mathit{CostRampDown}_{u,i} \geq \mathit{CostRampDown}_u \cdot (\mathit{PowerInitial}_u-\mathit{Power}_{u,i})

	 i>1:

	 \mathit{CostRampUp}_{u,i} \geq \mathit{CostRampUp}_u \cdot \left(\mathit{Power}_{u,i}-\mathit{Power}_{u,i-1}\right)

	 \mathit{CostRampDown}_{u,i} \geq \mathit{CostRampDown}_u \cdot (\mathit{Power}_{u,i-1}-\mathit{Power}_{u,i})


It should be noted that in case of start-up and shut-down, the ramping costs are added to the objective function. Using start-up, shut-down and ramping costs at the same time should therefore be performed with care.

In the current formulation all other costs (fixed and variable) are considered as exogenous parameters. The variable production costs (in {\texteuro}/MW), are determined by fuel and emission prices corrected by the efficiency (which is considered to be constant for all levels of output in this version of the model) and the emission rate of the unit (equation ):

.. math::
	 \mathit{CostVariable}_{u,h}=

	 \mathit{Markup}_{u,h} + \sum _{n,f}\left(\frac{\mathit{Fuel}_{u,f} \cdot \mathit{FuelPrice}_{n,f,h} \cdot \mathit{Location}_{u,n}}{\mathit{Efficiency}_u}\right)

	  + \sum _p\left(\mathit{EmissionRate}_{u,p} \cdot \mathit{PermitPrice}_p\right)

The previous equation includes an additional mark-up parameter that is used for calibration and validation purposes.

Transmission costs are also considered to be exogenous, and they result from multiplying the energy flows through the network by the corresponding transmission price (exogenous).

As regards load shedding, the model considers the possibility of voluntary load shedding resulting from contractual arrangements between generators and consumers. Additionally, in order to facilitate tracking and debugging of errors, the model also considers some variables representing the capacity the system is not able to provide when the minimum/maximum power, reserve, or ramping constraints are reached. These lost loads are a very expensive last resort of the system used when there is no other choice available. The different lost loads are assigned very high values (with respect to any other costs). This allows running the simulation without infeasibilities, thus helping to detect the origin of the loss of load. In a normal run of the model, without errors, all these variables are expected to be equal to zero.

Demand-related constraints
--------------------------

The main constraint to be met is the supply-demand balance, for each period and each zone, in the day-ahead market (equation ). According to this restriction, the sum of all the power produced by all the units present in the node (including the power generated by the storage units), the power injected from neighbouring nodes, and the curtailed power from intermittent sources is equal to the load in that node, plus the power consumed for energy storage, minus the load interrupted and the load shed.

.. math::
	 \sum _u\left(\mathit{Power}_{u,i} \cdot \mathit{Location}_{u,n}\right)

	  + \sum _l\left(\mathit{Flow}_{l,i} \cdot \mathit{LineNode}_{l,n}\right)

	 -\mathit{Curtailment}_n \cdot \mathit{CurtailedPower}_{n,i}

	 =\mathit{Demand}_{\mathit{DA},n,h} + \sum _r\left(\mathit{StorageInput}_{s,h} \cdot \mathit{Location}_{r,n}\right)
	
	  -\mathit{ShedLoad}_{n,i} - \mathit{LostLoadMaxPower}_{n,i} 
	  
	  + \mathit{LostLoadMinPower}_{n,i}

Besides that balance, the reserve requirements (upwards and downwards) in each node must be met as well. The reserve requirements considered in this model are an aggregation of secondary and tertiary reserves, which are typically brought online in periods shorter than an hour, the time step of this model. Therefore, additional equations and constraints must be defined for representing the up/down ramping requirements, by computing the ability of each unit to adapt its power output in periods below 60 minutes.

For each power plant, the ability to increase its power is the ramp-up capability if it is already committed or the nominal power if it is stopped and its starting time is lower than M minutes (equation ). This is to take into account that fast starting units could provide reserve (hydro units for secondary reserve, gas turbine for tertiary reserve).

.. math::

	\mathit{MaxRamp}2U_{u,i} 

	\leq \mathit{RampUpMaximum}_u  \cdot  \mathit{Committed}_{u,i} + \mathit{FlexibilityUp}_u 

	 \cdot  \left(1-\mathit{Committed}_{u,i} \right)

The parameter FlexibilityUpu is the maximum flexibility (in terms of ramping rate) that can be provided by the unit in case of cold start:

.. math::

	 \mathit{If}\mathit{RampStartUpMaximum}_u \geq \mathit{PowerMinStable}_u  \cdot  60

	 \mathit{Then}\mathit{FlexibilityUp}_u = \mathit{RampStartUpMaximum}_u

	 \mathit{Else}\mathit{FlexibilityUp}_u = 0

The maximum ramping rate is also limited by the available capacity margin between current and maximum power output (equation ).

.. math::

 	\mathit{MaxRamp2U}_{u,i} \leq (\mathit{PowerCapacit}y_u \cdot \mathit{AvailabilityFactor}_{u,i}

	 \cdot  (1-\mathit{OutageFactor}_{u,i})-\mathit{Power}_{u,i}) \cdot 60

The same applies to ramping down capabilities within periods below 60 minutes.

.. math::

	\mathit{MaxRamp}2D_{i,u}
	
	 \leq \mathit{max}\left(\mathit{RampDownMaximu}m_u,\mathit{Flexibility}\mathit{Down}_u\right) \cdot \mathit{Committed}_{u,i}

The parameter FlexibilityDownu is defined as the maximum ramp down rate at which the unit can shut down in M minutes.

In case the unit cannot be shut-down in M minutes (and only in this case) the maximum ramping down capability is limited by the capacity margin between actual and minimum power:

.. math::

	 If \mathit{RampShutDownMaximu}m_u<\mathit{PowerMinStabl}e_u \cdot 60 :

	 \mathit{Then}\mathit{MaxRamp}2D_{u,i} \leq \left(\mathit{Power}_{u,i}-\mathit{PowerMinStable}_u \cdot \mathit{Committed}_{u,i}\right) \cdot 60

	 Else :

	\mathit{MaxRamp}2D_{u,i} \leq \mathit{Power}_{u,i} \cdot 60 

The reserve requirements are defined by the users. In case no input is provided a default formula is used to evaluate the needs for secondary reserves as a function of the maximum expected load for each day. The default formula is described by equation , as defined by ENTSO-E [CITATION Mil13 {\textbackslash}m Eur04 {\textbackslash}l 1033 ]:

.. math::

	\mathit{Demand}_{2U,n,i}=\sqrt{10 \cdot \underset h{\mathit{max}}\left(\mathit{Demand}_{\mathit{DA},n,h}\right) + 150^2}-150

Downward reserves are defined as 50\% of the upward margin:

.. math::

	\mathit{Demand}_{2D,n,h}=0.5 \cdot \mathit{Demand}_{2U,n,h}

The reserve demand should be fulfilled at all times by all the plants allowed to participate in the reserve market:

.. math::

	\mathit{Demand}_{2U,n,h}
	
	 \leq \sum _{u,t}\left(\mathit{MaxRamp}2U_{u,i} \cdot \mathit{Technology}_{u,t} \cdot \mathit{Reserv}e_t \cdot \mathit{Locatio}n_{u,n}\right)

	+ \mathit{LostLoadReserve2UH}_{n,i}

The same equation applies to downward reserve requirements (2D).


Power output bounds
-------------------

The minimum power output is determined by the must-run or stable generation level of the unit if it is committed:

.. math::

	\mathit{Power}\mathit{MustRun}_{u,i} \cdot \mathit{Committed}_{u,i}

	 \leq \mathit{Power}_{u,i}

On the other hand, the output is limited by the available capacity, if the unit is committed:

.. math::

	\mathit{Power}_{u,i}

	 \leq \mathit{PowerCapacity}_u \cdot \mathit{AvailabilityFactor}_{u,i}

	 \cdot (1-\mathit{OutageFactor}_{u,i}) \cdot \mathit{Committed}_{u,i}

The power output in a given period also depends on the output levels in the previous and the following periods and on the ramping capabilities of the unit. If the unit was down, the ramping capability is given by the maximum start up ramp, while if the unit was online the limit is defined by the maximum ramp up rate. Those bounds are given by equation :

.. math::

	 i=1:

	 \mathit{Power}_{u,i} \leq 

	 \mathit{PowerInitial}_u

	  + \mathit{CommittedInitial}_u \cdot \mathit{RampUpMaximum}_u

	  + \left(1-\mathit{CommittedInitial}_u\right) \cdot \mathit{RampStartUpMaximum}_u

	  + \mathit{LostLoadRampUp}_{u,i}

	 i>1:

	 \mathit{Power}_{u,i} \leq 

	 \mathit{Power}_{u,i-1}

	  + \mathit{Committed}_{u,i-1} \cdot \mathit{RampUpMaximum}_u

	  + \left(1-\mathit{Committed}_{u,i-1}\right) \cdot \mathit{RampStartUpMaximum}_u

	  + \mathit{LostLoadRampUp}_{u,i}

And by equation :

.. math::

	 i=1:

	 \mathit{Power}_{u,i} \leq 

	 \mathit{PowerCapacity}_u \cdot \mathit{LoadMaximum}_{u,i} \cdot \mathit{Committed}_{u,i}

	  + \left(1-\mathit{Committed}_{u,i}\right) \cdot \mathit{RampShutDownMaximum}_u

	  + \mathit{LostLoadRampDown}_{u,i}

	 i<N:

	 \mathit{Power}_{u,i} \leq 

	 \mathit{PowerCapacity}_u \cdot \mathit{LoadMaximum}_{u,i} \cdot \mathit{Committed}_{u,i + 1}

	  + \left(1-\mathit{Committed}_{u,i + 1}\right) \cdot \mathit{RampShutDownMaximum}_u

	  + \mathit{LostLoadRampDown}_{u,i} 

Where the :math:`LoadMaximum_{u,i}` parameter is calculated taking into account the availability factor and the outage factor:

.. math::

	\mathit{LoadMaximum}_{u,h}=\mathit{AvailabilityFactor}_{u,h} \cdot (1-\mathit{OutageFactor}_{u,g})

Similarly, the ramp down capability is limited by the maximum ramp down or the maximum shut down ramp rate:

.. math::

	 i=1:

	 \mathit{PowerInitial}_u-\mathit{Power}_{u,i} \leq 

	 \mathit{Committed}_{u,i} \cdot \mathit{RampDownMaximum}_u

	  + \left(1-\mathit{Committed}_{u,i}\right) \cdot \mathit{RampShutDownMaximum}_u

	  + \mathit{LostLoadRampDown}_{u,i}

	 i>1:

	 \mathit{Power}_{u,i-1}-\mathit{Power}_{u,i} \leq 

	 \mathit{Committed}_{u,i} \cdot \mathit{RampDownMaximum}_u

	  + \left(\mathit{Committed}_{u,i-1}-\mathit{Committed}_{u,i}\right) \cdot \mathit{RampShutDownMaximum}_u

	  + \mathit{LostLoadRampDown}_{u,i}

While the ramp up limitation is defined by:

.. math::

	 i=1:

	 \mathit{PowerInitial}_u-\mathit{Power}_{u,i} \leq 

	 \mathit{CommittedInitial}_u \cdot \mathit{RampUpMaximum}_u

	  + \left(\mathit{Committed}_{u,i}-\mathit{CommittedInitial}_u\right) \cdot \mathit{RampStartUpMaximum}_u

	  + \mathit{LostLoadRampUp}_{u,i}

	 i>1:

	 \mathit{Power}_{u,i}-\mathit{Power}_{u,i-1} \leq 

	 \mathit{Committed}_{u,i-1} \cdot \mathit{RampUpMaxim}\mathit{um}_u

	  + \left(\mathit{Committed}_{u,i}-\mathit{Committed}_{u,i-1}\right) \cdot \mathit{RampStartUpMaximum}_u

	  + \mathit{LostLoadRampUp}_{u,i}

Minimum up and down times
-------------------------

The operation of the generation units is limited as well by the amount of time the unit has been running or stopped. Due to the physical characteristics of the generators, once a unit is started up it cannot be shut down immediately, while if the unit is shut down it may not be started immediately. 

That is, the value of the time counter with respect to the minimum up time and down times determines the commitment status of the unit. In order to model theses constraints linearly, it is necessary to keep track of the number of hours the unit must be online at the beginning of the simulation for having been online less than the minimum up time:

.. math::

	\mathit{TimeUpLeftInitial}_u =

	\mathit{min}\left\{N,\left(\mathit{TimeUpMinimum}_u - \mathit{TimeUpInitial}_u\right) \cdot \mathit{CommittedInitial}_u\right\}

If the unit is initially started up, it has to remain committed until reaching the minimum up time:

.. math::

	\sum _{i=1}^{\mathit{TimeUpLeftInitial}_u}\left(1-\mathit{Committed}_{u,i}\right)=0

If the unit is started during the considered horizon, the time it has to remain online is TimeUpMinimum, but cannot exceed the time remaining in the simulated period. This is expressed in equation  and is pre-calculated for each time step of the period.

.. math::

	\mathit{TimeUpLeftJustStarted}_{u,i}=

	\mathit{min}\left\{N -i + 1,\mathit{TimeUpMinimum}_u\right\}

The equation imposing the unit to remain committed is written:

.. math::

	 i=1:

	 \sum _{\mathit{ii}=i}^{i + \mathit{TimeUpLeftJustStarted}_{u,i}-1}\mathit{Committed}_{u,\mathit{ii}} \geq 

	 \mathit{TimeUpLeftJustStarted}_{u,i} \cdot \left(\mathit{Committed}_{u,i}-\mathit{CommittedInitial}_u\right)

	 i>1:

	 \sum _{\mathit{ii}=i}^{i + \mathit{TimeUpLeftJustStarted}_u-1}\mathit{Committed}_{u,\mathit{ii}} \geq 

	\mathit{TimeUpLeftJustStarted}_{u,i} \cdot \left(\mathit{Committed}_{u,i}-\mathit{Committed}_{u,i-1}\right)

The same method can be applied to the minimum down time constraint:

.. math::

	 \mathit{TimeDownLeft}_u = 

	 \mathit{min}\{N,(\mathit{TimeDownMinimum}_u-\mathit{TimeDownInitial}_u) 

	 \cdot (1-\mathit{CommittedInitial}_u)\}

Related to the initial status of the unit:

.. math::

	\sum _{i=1}^{\mathit{TimeDownLeft}_u}\mathit{Committed}_{u,i}=0

The TimeDownLeftJustStopped parameter is computed by:

.. math::

	\mathit{TimeDownLeftJustStopped}_{u,i} = 

	\mathit{min}\left\{N - i + 1,\mathit{TimeDownMinimum}_u\right\}

Finally, the equation imposing the time the unit has to remain de-committed is defined as:

.. math:: 

	 i=1:

	 \sum _{\mathit{ii}=i}^{i + \mathit{TimeDownLeftJustStopped}_{i,u}-1}\left(1-\mathit{Committed}_{u,i}\right) \geq 

	 \mathit{TimeDownLeftJustStopped}_{u,i} \cdot \left(\mathit{CommittedInitial}_u-\mathit{Committed}_{u,i}\right)

	 i>1:

	 \sum _{\mathit{ii}=i}^{i + \mathit{TimeDownLeftJustStopped}_u-1}\left(1-\mathit{Committed}_{u,i}\right) \geq 

	 \mathit{TimeDownLeftJustStopped}_{u,i} \cdot \left(\mathit{Committed}_{u,i-1}-\mathit{Committed}_{u,i}\right)

Storage-related constraints
---------------------------

Generation units with energy storage capabilities (mostly large hydro reservoirs and pumped hydro storage units) must meet additional restrictions related to the amount of energy stored. Storage units are considered to be subject to the same constraints as non-storage power plants. In addition to those constraints, storage-specific restrictions are added for the set of storage units (i.e. a subset of all units). These restrictions include the storage capacity, inflow, outflow, charging, charging capacity, charge/discharge efficiencies, etc. Discharging is considered as the standard operation mode and is therefore linked to the Power variable, common to all units.

The first constrain imposes that the energy stored by a given unit is bounded by a minimum value:

.. math::

 	\mathit{StorageMinimum}_s \leq \mathit{StorageLevel}_{s,i}

And the storage capacity:

.. math::

	\mathit{StorageLevel}_{s,i} \leq \mathit{StorageCapacity}_s

The energy added to the storage unit is limited by the charging capacity. Charging is allowed only if the unit is not producing (discharging) at the same time (i.e. if Committed, corresponding to the {\textquotedbl}normal{\textquotedbl} mode, is equal to 0).

.. math::

	\mathit{StorageInput}_{s,i} \leq \mathit{StorageChargingCapacity}_s \cdot (1-\mathit{Committed}_{s,i})

Discharge is limited by the level of charge of the storage unit:

.. math::

	\frac{\mathit{Power}_{i,s}}{\mathit{StorageDischargeEfficienc}y_s} + \mathit{StorageOutflow}_{s,i}-\mathit{StorageInflow}_{s,i} 

	\leq \mathit{StorageLevel}_{s,i}


Charge is limited by the level of charge of the storage unit:

.. math::

	\mathit{StorageInput}_{s,i} \cdot \mathit{StorageChargingEfficienc}y_s-\mathit{StorageOutflow}_{s,i} 
	
	+ \mathit{Storag}\mathit{eInflow}_{s,i} \leq \mathit{StorageCapacity}_s-\mathit{StorageLevel}_{s,i}

Besides, the energy stored in a given period is given by the energy stored in the previous period, net of charges and discharges:

.. math::

	i=1:

	\mathit{StorageLevelInitial}_s + \mathit{StorageInflow}_{s,i} 

	+ \mathit{StorageInput}_{s,i} \cdot \mathit{StorageChargingEfficiency}_s

	= \mathit{StorageLevel}_{s,i} + \mathit{StorageOutflow}_{s,i} + \frac{\mathit{Power}_{s,i}}{\mathit{StorageDischargeEfficienc}y_s}

	i>1:
	
	\mathit{StorageLevel}_{s,i-1} + \mathit{StorageInflow}_{s,i} 

	+ \mathit{StorageInput}_{s,i} \cdot \mathit{StorageChargingEfficiency}_s

	= \mathit{StorageLevel}_{s,i} + \mathit{StorageOutflow}_{s,i} + \frac{\mathit{Power}_{s,i}}{\mathit{StorageDischargeEfficienc}y_s}

Emission limits
---------------

The operating schedule also needs to take into account any cap on the emissions (not only CO2) from the generation units existing in each node:

.. math::

	\sum _u\left(\mathit{Power}_{u,i} \cdot \mathit{EmisionRate}_{u,p} \cdot \mathit{Location}_{u,n}\right)

	\leq \mathit{EmisionMaximum}_{n,p}

Network-related constraints
---------------------------

The flow of power between nodes is limited by the capacities of the transmission lines:

.. math::

	\mathit{FlowMinimum}_{l,i} \leq \mathit{Flow}_{l,i}

	\mathit{Flow}_{l,i} \leq \mathit{FlowMaximum}_{l,i}

In this model a simple transport-problem approach is followed.

Curtailment
-----------

If curtailment of intermittent generation sources is allowed in one node, the amount of curtailed power is bounded by the output of the renewable (tr) units present in that node:

.. math::

	\mathit{CurtailedPower}_{n,i}

	\leq \sum _{u,\mathit{tr}}\left(\mathit{Power}_{u,i} \cdot \mathit{Technology}_{u,\mathit{tr}} \cdot \mathit{Location}_{u,n}\right) \cdot \mathit{Curtailment}_n


Load shedding
-------------

If load shedding is allowed in a node, the amount of shed load is limited by the shedding capacity contracted on that particular node (e.g. through interruptible industrial contracts

.. math::

	\mathit{ShedLoad}_{n,i} \leq \mathit{LoadShedding}_n	



Interface
^^^^^^^^^

This section describes the different simulation files, templates and scripts required to run the DispaSET model. For each simulation, these files are included into a single directory corresponding to a self-sufficient simulation environment.
The typical step-by-step procedure to parametrize and run a DispaSET simulation is the following:

1. Fill the data template (the different files InputDispa-SET – ParameterName.xlsx) with properly formatted data (time series, power plant data, etc.)
2. Configure the simulation parameter (rolling horizon, data slicing) in the InputDispa-SET - Config.xlsx file.
3. Generate the Inputs.gdx file using the bat script makeGDX.bat 
4. Open the GAMS simulation files (project: UCM.gpr and model: UCM_h.gms). Run the model.
5. Read and display results in Results.xlsx

A more comprehensive description of these different files is provided hereunder.

UCM_h.gms and UCM.gpr
---------------------

UCM_h.gms is the main GAMS model described in Chapter 1. A copy of this file is included in each simulation environment, allowing keeping track of the exact version of the model used for the simulation. The model must be run in GAMS and requires a proper input file (Inputs.gdx). 

.. table:: 

	=============== =============================== =====================================
	Requires: 	Inputs.gdx			Input file for the simulation.
	Generates:	Results.gdx			Simulation results in gdx format	
	. 		Results.xlsx			Simulation results in xlsx format.
	=============== =============================== =====================================

UCM.gpr is the GAMS project file which should be opened before UCM_h.gms.

make_gdx.gms
------------

GAMS file that reads the different template excel files and generates the Inputs.gdx file. This file should be opened in GAMS.

.. table:: 

	=============== =============================== =====================================
	Requires: 	InputDispa-SET – xxx.xlsx	DispaSET template files
	Generates:	Inputs.gdx			Input file for the simulation	
	=============== =============================== =====================================
			

makeGDX.bat
-----------

Batch script that generates the input file from the template without requiring opening GAMS. The first time it is executed, the path of the GAMS folder must be provided.

.. table:: 

	=============== =============================== =====================================
	Requires: 	InputDispa-SET – xxx.xlsx	DispaSET template files
	.		make_gdx.gms			GAMS file to generate Inputs.gdx
	Generates:	Inputs.gdx			Input file for the simulation	
	=============== =============================== =====================================


writeresults.gms
----------------

GAMS file to generate the excel Results.xlsx file from the Results.gdx generated by GAMS (in case the write_excel function was deactivated in GAMS. 

.. table:: 

	=============== =============================== =====================================
	Requires: 	Results.gdx			Simulation results in gdx format
	Generates:	Results.xlsx			Simulation results in xlsx format	
	=============== =============================== =====================================
			

Inputs.gdx
----------

All the inputs of the model must be stored in the Inputs.gdx file since it is the only file read by the main GAMS model. This file is generated from the DispaSET template.

.. table:: 

	=============== =============================== =====================================
	Requires: 	InputDispa-SET – xxx.xlsx	DispaSET template files
	Generates:					 
	=============== =============================== =====================================


InputDispa-SET -  [ParameterName].xlsx
--------------------------------------

Series of 42 excel files, each corresponding to a parameter of the DispaSET model (see Chapter 1). The files must be formatted according to section 2.2.

InputDispa-SET -  Sets.xlsx
---------------------------

Single excel file that contains all the sets used in the model in a column format. 

InputDispa-SET -  Config.xlsx
-----------------------------

Single excel file that contains simulation metadata in the form of a Table. This metadata allows setting the rolling horizon parameter and slicing the input data to simulate a subset only.

.. table:: Config

	=============================== ======= ======= ======= =================================================
					Year	Month	Day	Description
	=============================== ======= ======= ======= =================================================
	FirstDay			2012	10	1	First day of the simulation in the template data
	LastDay				2013	9	30	Last day of the simulation in the template data
	DayStart			2013	5	29	First day to be simulated
	DayStop				2013	7	2	Last day to be simulated
	RollingHorizon Length		0	0	3	Length of the rolling horizons 
	RollingHorizon LookAhead	0	0	1	Overlap period of the rolling horizon 
	=============================== ======= ======= ======= =================================================



Structure of the Excel template
-------------------------------

The name of the input files are "Input Dispa-SET – [Parameter name].xlsx". These files contain the data to be read by the model, after conversion into a GDX file. 

The structure of all input files follows the following rules: 

1. There is one file per model parameter 
2. Each file contains only one sheet 
3. The first row is left blank for non-time series data (i.e. data starts at A2)
4. For time series data, the rows are organized as follows:
	a. The first row is left blank
	b. Rows 2 to 5 contains the year, month, day and hour of each data
	c. Row 6 contains the time index of the data, which will be used in DispaSET
	d. The data therefore starts at A6
5. If one of the input sets of the data is u (the unit name), it is always defined as the first column of the data (column A)
6. If one of the input sets of the data is h (the time index), it is always defined as the only horizontal input in row 6

In the case of the file "Input Dispa-SET – Sets.xlsx", all the required sets are written in columns with the set name in row 2.








References
^^^^^^^^^^

.. [1] Hidalgo González, I., Quoilin, S., & Zucker, A. (2014). Dispa-SET 2.0: unit commitment and power dispatch model (EUR 27015 EN). Petten, Netherlands: European Commission. 
.. [2] Quoilin, S., Nijs, W., Hidalgo, I., & Thiel, C. (2015). Evaluation of simplified flexibility evaluation tools using a unit commitment model. IEEE Digital Library. 
.. [3] Quoilin, S., Gonzalez Vazquez, I., Zucker, A., & Thiel, C. (2014). Available technical flexibility for balancing variable renewable energy sources: case study in Belgium. Proceedings of the 9th Conference on Sustainable Development of Energy, Water and Environment Systems. 
