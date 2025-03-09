.. _sector_coupling:

Sector Coupling
==============

The sector coupling model extends the power system optimization to include interactions with other energy sectors through storage and conversion processes. This formulation allows modeling sector coupling with various energy carriers (e.g., heat, hydrogen, etc.) while maintaining the temporal and operational constraints of the power system.

Model Structure
--------------

Sets
~~~~

Additional sets for the sector coupling model include:

.. table::

    ======= =================================================================================
    Name    Description
    ======= =================================================================================
    nx      Boundary sector nodes
    nx_CC   Boundary sector nodes without the EndCascade node
    lx      Boundary sector lines
    slx     Boundary sector spillage lines
    xu(au)  Boundary sector only units
    ======= =================================================================================

Parameters
~~~~~~~~~~

Key parameters specific to the sector coupling model:

.. table::

    ======================================= ======= =============================================================
    Name                                    Units   Description
    ======================================= ======= =============================================================
    SectorXDemand(nx,h)                     MWh     Demand profile in boundary sectors
    X2PowerConversionMultiplier(nx,au,h)    %       Discharge efficiency for boundary sector
    Power2XConversionMultiplier(nx,au,h)    %       Charging efficiency for boundary sector
    SectorXStorageCapacity(nx)              MWh     Storage capacity of the boundary sector
    SectorXStorageSelfDischarge(nx)         %       Boundary sector storage self discharge rate
    SectorXStorageMinimum(nx)               MWh     Boundary sector storage minimum level
    SectorXStorageProfile(nx,h)             %       Required storage level profile
    SectorXAlertLevel(nx,h)                 MWh     Storage alert level - violated only to avoid rationing
    SectorXFloodControl(nx,h)               MWh     Storage flood control level
    ======================================= ======= =============================================================

Sector Coupling Units
--------------------

Sector coupling units are specific units that enable the conversion between the power sector and boundary sectors. These units can operate in two modes:

1. Power-to-X: Converting electrical power to another energy carrier (e.g., power-to-heat, power-to-hydrogen)
2. X-to-Power: Converting energy from a boundary sector back to electrical power (e.g., fuel cells)

The conversion process is governed by efficiency parameters:

.. math::
    \begin{split}
    \mathit{SectorXStorageInput}_{nx,i} = & \sum_{au} \mathit{Power2XConversionMultiplier}_{nx,au,i} \cdot \mathit{PowerConsumption}_{au,i} \cdot \mathit{LocationX}_{au,nx} \\
    & - \sum_{au} \mathit{X2PowerConversionMultiplier}_{nx,au,i} \cdot \mathit{Power}_{au,i} \cdot \mathit{LocationX}_{au,nx}
    \end{split}

This equation represents the net storage input as the difference between:
- Charging: Power consumption multiplied by the Power-to-X conversion efficiency
- Discharging: Power generation multiplied by the X-to-Power conversion efficiency

Boundary Sector Storage
----------------------

Each boundary sector can include storage capabilities, subject to the following constraints:

1. Storage Level Minimum:

.. math::
    \mathit{SectorXStorageMinimum}_{nx} \leq \mathit{SectorXStorageLevel}_{nx,i}

2. Storage Level Maximum:

.. math::
    \mathit{SectorXStorageLevel}_{nx,i} \leq \mathit{SectorXStorageCapacity}_{nx}

3. Storage Balance:

.. math::
    \begin{split}
    \mathit{SectorXStorageLevel}_{nx,i} = & \mathit{SectorXStorageInitial}_{nx} \cdot \delta_{i=1} + \mathit{SectorXStorageLevel}_{nx,i-1} \cdot (1-\delta_{i=1}) \\
    & + \mathit{SectorXStorageInput}_{nx,i} \cdot \mathit{TimeStep} \\
    & - \mathit{SectorXStorageSelfDischarge}_{nx} \cdot \mathit{SectorXStorageLevel}_{nx,i} \cdot \mathit{TimeStep}
    \end{split}

where :math:`\delta_{i=1}` is 1 for the first time step and 0 otherwise.

4. Storage Final Level:

.. math::
    \mathit{SectorXStorageFinalMin}_{nx} \leq \mathit{SectorXStorageLevel}_{nx,i} + \mathit{SectorXStorageLevelViolation}_{nx}

for the final time step i.

Boundary Sector Demand
---------------------

The demand in boundary sectors can be categorized into two types:

Non-Flexible Demand
~~~~~~~~~~~~~~~~~~

Non-flexible demand (:math:`\mathit{SectorXDemand}_{nx,h}`) represents the fixed energy requirements that must be met at each time step. This demand is typically defined through input time series and must be satisfied either through:
- Direct conversion from power using sector coupling units
- Discharge from boundary sector storage
- Alternative supply options (at a high cost penalty)

Flexible Demand
~~~~~~~~~~~~~~

Flexible demand (:math:`\mathit{SectorXFlexDemand}_{nx,h}`) represents demand that can be shifted in time, subject to constraints:

.. math::
    \sum_{h} \mathit{SectorXFlexDemand}_{nx,h} = \mathit{SectorXFlexDemandInput}_{nx}

This ensures that the total flexible demand over the optimization horizon matches the required input while allowing temporal flexibility in when it is met.

Additional constraints limit the maximum flexible demand at each time step:

.. math::
    \mathit{SectorXFlexDemand}_{nx,h} \leq \mathit{SectorXFlexMaxCapacity}_{nx}

Cost Terms
---------

The sector coupling model introduces additional cost terms to the objective function:

.. math::
    \begin{split}
    & + \sum_{nx,i} \mathit{CostXSpillage}_{slx,i} \cdot \mathit{SectorXSpillage}_{slx,i} \cdot \mathit{TimeStep} \\
    & + \sum_{nx,i} \mathit{CostXStorageAlert}_{nx,i} \cdot \mathit{SectorXStorageLevelViolation\_H}_{nx,i} \cdot \mathit{TimeStep} \\
    & + \sum_{nx,i} \mathit{CostXFloodControl}_{nx,i} \cdot \mathit{SectorXStorageLevelViolation}_{nx} \cdot \mathit{TimeStep} \\
    & + \sum_{nx,i} \mathit{CostXNotServed}_{nx,i} \cdot (\mathit{LL\_SectorXFlexDemand}_{nx} + \mathit{LL\_SectorXFlexSupply}_{nx}) \cdot \mathit{TimeStep}
    \end{split}

These terms account for:
- Spillage between boundary sector nodes
- Violations of storage alert levels and flood control levels
- Deficits in meeting flexible demand or supply requirements
- Energy not served in boundary sectors 