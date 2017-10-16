.. Dispa-SET documentation master file, created by
   sphinx-quickstart on Mon Feb  8 16:23:20 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The Dispa-SET model
===================

The Dispa-SET model is an open-source unit commitment and optimal dispatch model focused on the balancing and flexibility problems in European grids. Its pre and post-processing tools are written in Python and the main solver can be called via GAMS or via PYOMO. The selected Mixed-Integer Linear Programming (MILP) solver is CPLEX.
 
Dispa-SET is mainly developed within the Joint Research Centre of the EU Commission, in close collaboration with the University of Liège (Belgium).
 
Model description and philosophy
--------------------------------
A comprehensive description of the model is available in the 2017 JRC technical report: `Modelling Future EU Power Systems Under High Shares of Renewables`_.

.. image:: figures/report2.jpg
 
Downloading Dispa-SET
---------------------
The public version of Dispa-SET can be downloaded in the :ref:`releases` section or from its github repository (using the Clone or Download button on the right side of the screen):
https://github.com/energy-modelling-toolkit/Dispa-SET
 
Documentation
-------------
The model documentation is available by running sphinx in the Docs folder of the project or by consulting the online documentation. This documentation corresponds to the latest available public version of Dispa-SET: 
http://www.dispaset.eu/latest/index.html
 
Main contributors:

* Sylvain Quoilin (Main Developper) 

* Konstantinos Kavvadias  

* Andreas Zucker 


Contents
--------

.. toctree::
   :maxdepth: 1

   overview
   workflow
   model
   implementation
   data
   developers
   DispaSET
   releases


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. _Modelling Future EU Power Systems Under High Shares of Renewables: https://ec.europa.eu/jrc/en/publication/eur-scientific-and-technical-research-reports/modelling-future-eu-power-systems-under-high-shares-renewables-dispa-set-21-open-source
