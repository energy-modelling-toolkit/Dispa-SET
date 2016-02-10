
$CALL GDXXRW "InputDispa-SET - Sets.xlsx" Set=d rng=A3 Rdim=1  O=d.gdx 
$GDXIN d.gdx 
Set d; 
$LOAD d
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Sets.xlsx" Set=y_config rng=B3 Rdim=1  O=y_config.gdx 
$GDXIN y_config.gdx 
Set y_config; 
$LOAD y_config
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Sets.xlsx" Set=f rng=C3 Rdim=1  O=f.gdx 
$GDXIN f.gdx 
Set f; 
$LOAD f
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Sets.xlsx" Set=h rng=D3 Rdim=1  O=h.gdx 
$GDXIN h.gdx 
Set h; 
$LOAD h
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Sets.xlsx" Set=tr rng=E3 Rdim=1  O=tr.gdx 
$GDXIN tr.gdx 
Set tr; 
$LOAD tr
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Sets.xlsx" Set=l rng=F3 Rdim=1  O=l.gdx 
$GDXIN l.gdx 
Set l; 
$LOAD l
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Sets.xlsx" Set=mk rng=G3 Rdim=1  O=mk.gdx 
$GDXIN mk.gdx 
Set mk; 
$LOAD mk
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Sets.xlsx" Set=n rng=H3 Rdim=1  O=n.gdx 
$GDXIN n.gdx 
Set n; 
$LOAD n
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Sets.xlsx" Set=p rng=I3 Rdim=1  O=p.gdx 
$GDXIN p.gdx 
Set p; 
$LOAD p
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Sets.xlsx" Set=s rng=J3 Rdim=1  O=s.gdx 
$GDXIN s.gdx 
Set s; 
$LOAD s
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Sets.xlsx" Set=u rng=K3 Rdim=1  O=u.gdx 
$GDXIN u.gdx 
Set u; 
$LOAD u
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Sets.xlsx" Set=t rng=L3 Rdim=1  O=t.gdx 
$GDXIN t.gdx 
Set t; 
$LOAD t
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Sets.xlsx" Set=x_config rng=M3 Rdim=1  O=x_config.gdx 
$GDXIN x_config.gdx 
Set x_config; 
$LOAD x_config
$GDXIN 

$CALL GDXXRW "InputDispa-SET - TimeDownInitial.xlsx" par=TimeDownInitial rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - TimeDownInitial.gdx" 
Parameter TimeDownInitial; 
$LOAD TimeDownInitial
$GDXIN 

$CALL GDXXRW "InputDispa-SET - RampUpMaximum.xlsx" par=RampUpMaximum rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - RampUpMaximum.gdx" 
Parameter RampUpMaximum; 
$LOAD RampUpMaximum
$GDXIN 

$CALL GDXXRW "InputDispa-SET - CostFixed.xlsx" par=CostFixed rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - CostFixed.gdx" 
Parameter CostFixed; 
$LOAD CostFixed
$GDXIN 

$CALL GDXXRW "InputDispa-SET - EmissionRate.xlsx" par=EmissionRate rng=A2 Rdim=1 Cdim=1 
$GDXIN "InputDispa-SET - EmissionRate.gdx" 
Parameter EmissionRate; 
$LOAD EmissionRate
$GDXIN 

$CALL GDXXRW "InputDispa-SET - RampStartUpMaximum.xlsx" par=RampStartUpMaximum rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - RampStartUpMaximum.gdx" 
Parameter RampStartUpMaximum; 
$LOAD RampStartUpMaximum
$GDXIN 

$CALL GDXXRW "InputDispa-SET - PermitPrice.xlsx" par=PermitPrice rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - PermitPrice.gdx" 
Parameter PermitPrice; 
$LOAD PermitPrice
$GDXIN 

$CALL GDXXRW "InputDispa-SET - PowerInitial.xlsx" par=PowerInitial rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - PowerInitial.gdx" 
Parameter PowerInitial; 
$LOAD PowerInitial
$GDXIN 

$CALL GDXXRW "InputDispa-SET - CostStartUp.xlsx" par=CostStartUp rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - CostStartUp.gdx" 
Parameter CostStartUp; 
$LOAD CostStartUp
$GDXIN 

$CALL GDXXRW "InputDispa-SET - TimeDownMinimum.xlsx" par=TimeDownMinimum rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - TimeDownMinimum.gdx" 
Parameter TimeDownMinimum; 
$LOAD TimeDownMinimum
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Fuel.xlsx" par=Fuel rng=A2 Rdim=1 Cdim=1 
$GDXIN "InputDispa-SET - Fuel.gdx" 
Parameter Fuel; 
$LOAD Fuel
$GDXIN 

$CALL GDXXRW "InputDispa-SET - LineNode.xlsx" par=LineNode rng=A2 Rdim=1 Cdim=1 
$GDXIN "InputDispa-SET - LineNode.gdx" 
Parameter LineNode; 
$LOAD LineNode
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Config.xlsx" par=Config rng=A2 Rdim=1 Cdim=1 
$GDXIN "InputDispa-SET - Config.gdx" 
Parameter Config; 
$LOAD Config
$GDXIN 

$CALL GDXXRW "InputDispa-SET - PartLoadMin.xlsx" par=PartLoadMin rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - PartLoadMin.gdx" 
Parameter PartLoadMin; 
$LOAD PartLoadMin
$GDXIN 

$CALL GDXXRW "InputDispa-SET - FuelPrice.xlsx" par=FuelPrice rng=A6 Rdim=2 Cdim=1 
$GDXIN "InputDispa-SET - FuelPrice.gdx" 
Parameter FuelPrice; 
$LOAD FuelPrice
$GDXIN 

$CALL GDXXRW "InputDispa-SET - OutageFactor.xlsx" par=OutageFactor rng=A6 Rdim=1 Cdim=1 
$GDXIN "InputDispa-SET - OutageFactor.gdx" 
Parameter OutageFactor; 
$LOAD OutageFactor
$GDXIN 

$CALL GDXXRW "InputDispa-SET - RampShutDownMaximum.xlsx" par=RampShutDownMaximum rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - RampShutDownMaximum.gdx" 
Parameter RampShutDownMaximum; 
$LOAD RampShutDownMaximum
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Efficiency.xlsx" par=Efficiency rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - Efficiency.gdx" 
Parameter Efficiency; 
$LOAD Efficiency
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Curtailment.xlsx" par=Curtailment rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - Curtailment.gdx" 
Parameter Curtailment; 
$LOAD Curtailment
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Demand.xlsx" par=Demand rng=A6 Rdim=2 Cdim=1 
$GDXIN "InputDispa-SET - Demand.gdx" 
Parameter Demand; 
$LOAD Demand
$GDXIN 

$CALL GDXXRW "InputDispa-SET - StorageOutflow.xlsx" par=StorageOutflow rng=A6 Rdim=1 Cdim=1 
$GDXIN "InputDispa-SET - StorageOutflow.gdx" 
Parameter StorageOutflow; 
$LOAD StorageOutflow
$GDXIN 

$CALL GDXXRW "InputDispa-SET - CostShutDown.xlsx" par=CostShutDown rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - CostShutDown.gdx" 
Parameter CostShutDown; 
$LOAD CostShutDown
$GDXIN 

$CALL GDXXRW "InputDispa-SET - FlowMinimum.xlsx" par=FlowMinimum rng=A6 Rdim=1 Cdim=1 
$GDXIN "InputDispa-SET - FlowMinimum.gdx" 
Parameter FlowMinimum; 
$LOAD FlowMinimum
$GDXIN 

$CALL GDXXRW "InputDispa-SET - StorageChargingEfficiency.xlsx" par=StorageChargingEfficiency rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - StorageChargingEfficiency.gdx" 
Parameter StorageChargingEfficiency; 
$LOAD StorageChargingEfficiency
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Markup.xlsx" par=Markup rng=A6 Rdim=1 Cdim=1 
$GDXIN "InputDispa-SET - Markup.gdx" 
Parameter Markup; 
$LOAD Markup
$GDXIN 

$CALL GDXXRW "InputDispa-SET - CostVariable.xlsx" par=CostVariable rng=A6 Rdim=1 Cdim=1 
$GDXIN "InputDispa-SET - CostVariable.gdx" 
Parameter CostVariable; 
$LOAD CostVariable
$GDXIN 

$CALL GDXXRW "InputDispa-SET - EmissionMaximum.xlsx" par=EmissionMaximum rng=A2 Rdim=1 Cdim=1 
$GDXIN "InputDispa-SET - EmissionMaximum.gdx" 
Parameter EmissionMaximum; 
$LOAD EmissionMaximum
$GDXIN 

$CALL GDXXRW "InputDispa-SET - PriceTransmission.xlsx" par=PriceTransmission rng=A6 Rdim=1 Cdim=1 
$GDXIN "InputDispa-SET - PriceTransmission.gdx" 
Parameter PriceTransmission; 
$LOAD PriceTransmission
$GDXIN 

$CALL GDXXRW "InputDispa-SET - StorageCapacity.xlsx" par=StorageCapacity rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - StorageCapacity.gdx" 
Parameter StorageCapacity; 
$LOAD StorageCapacity
$GDXIN 

$CALL GDXXRW "InputDispa-SET - StorageChargingCapacity.xlsx" par=StorageChargingCapacity rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - StorageChargingCapacity.gdx" 
Parameter StorageChargingCapacity; 
$LOAD StorageChargingCapacity
$GDXIN 

$CALL GDXXRW "InputDispa-SET - PowerCapacity.xlsx" par=PowerCapacity rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - PowerCapacity.gdx" 
Parameter PowerCapacity; 
$LOAD PowerCapacity
$GDXIN 

$CALL GDXXRW "InputDispa-SET - StorageMinimum.xlsx" par=StorageMinimum rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - StorageMinimum.gdx" 
Parameter StorageMinimum; 
$LOAD StorageMinimum
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Reserve.xlsx" par=Reserve rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - Reserve.gdx" 
Parameter Reserve; 
$LOAD Reserve
$GDXIN 

$CALL GDXXRW "InputDispa-SET - RampDownMaximum.xlsx" par=RampDownMaximum rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - RampDownMaximum.gdx" 
Parameter RampDownMaximum; 
$LOAD RampDownMaximum
$GDXIN 

$CALL GDXXRW "InputDispa-SET - TimeUpMinimum.xlsx" par=TimeUpMinimum rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - TimeUpMinimum.gdx" 
Parameter TimeUpMinimum; 
$LOAD TimeUpMinimum
$GDXIN 

$CALL GDXXRW "InputDispa-SET - StorageInitial.xlsx" par=StorageInitial rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - StorageInitial.gdx" 
Parameter StorageInitial; 
$LOAD StorageInitial
$GDXIN 

$CALL GDXXRW "InputDispa-SET - StorageDischargeEfficiency.xlsx" par=StorageDischargeEfficiency rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - StorageDischargeEfficiency.gdx" 
Parameter StorageDischargeEfficiency; 
$LOAD StorageDischargeEfficiency
$GDXIN 

$CALL GDXXRW "InputDispa-SET - AvailabilityFactor.xlsx" par=AvailabilityFactor rng=A6 Rdim=1 Cdim=1 
$GDXIN "InputDispa-SET - AvailabilityFactor.gdx" 
Parameter AvailabilityFactor; 
$LOAD AvailabilityFactor
$GDXIN 

$CALL GDXXRW "InputDispa-SET - StorageInflow.xlsx" par=StorageInflow rng=A6 Rdim=1 Cdim=1 
$GDXIN "InputDispa-SET - StorageInflow.gdx" 
Parameter StorageInflow; 
$LOAD StorageInflow
$GDXIN 

$CALL GDXXRW "InputDispa-SET - FlowMaximum.xlsx" par=FlowMaximum rng=A6 Rdim=1 Cdim=1 
$GDXIN "InputDispa-SET - FlowMaximum.gdx" 
Parameter FlowMaximum; 
$LOAD FlowMaximum
$GDXIN 

$CALL GDXXRW "InputDispa-SET - LoadShedding.xlsx" par=LoadShedding rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - LoadShedding.gdx" 
Parameter LoadShedding; 
$LOAD LoadShedding
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Technology.xlsx" par=Technology rng=A2 Rdim=1 Cdim=1 
$GDXIN "InputDispa-SET - Technology.gdx" 
Parameter Technology; 
$LOAD Technology
$GDXIN 

$CALL GDXXRW "InputDispa-SET - TimeUpInitial.xlsx" par=TimeUpInitial rng=A2 Rdim=1 
$GDXIN "InputDispa-SET - TimeUpInitial.gdx" 
Parameter TimeUpInitial; 
$LOAD TimeUpInitial
$GDXIN 

$CALL GDXXRW "InputDispa-SET - Location.xlsx" par=Location rng=A2 Rdim=1 Cdim=1 
$GDXIN "InputDispa-SET - Location.gdx" 
Parameter Location; 
$LOAD Location
$GDXIN 

Execute_Unload "Inputs.gdx"