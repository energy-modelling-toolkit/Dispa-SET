.. Dispa-SET documentation master file, created by
   sphinx-quickstart on Mon Feb  8 16:23:20 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Dispa-SET's documentation!
=====================================
:Organization:  `Joint Research Centre`_, 
		`Institute of Energy and Transport`_,
                Energy Technology Policy Outlook Unit
:Version: |version|
:Date: |today|

The Dispa-SET model is a unit commitment and dispatch model developed within the “Joint Research Centre” and focused on the balancing and flexibility problems in European grids [1]_. It is written in GAMS an Python (Pyomo) and Excel for input/output data handling and visualization. The selected Mixed-Integer Linear Programming (MILP) solver is CPLEX.



Libraries used
--------------

* `pyomo`_ Optimization object library, interface to LP solver (e.g. CPLEX)
* `pandas`_ for input and result data handling 
* `matplotlib`_ for plotting

  
Contents
--------

.. toctree::
   :maxdepth: 1

   overview
   model
   workflow
   data
   API
   developers


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

References
----------
.. [1] Hidalgo González, I., Quoilin, S., & Zucker, A. (2014). Dispa-SET 2.0: unit commitment and power dispatch model (EUR 27015 EN). Petten, Netherlands: European Commission. 

.. _matplotlib: http://matplotlib.org
.. _pandas: http://pandas.pydata.org
.. _pyomo: http://www.pyomo.org/
.. _Institute of Energy and Transport: https://ec.europa.eu/jrc/en/institutes/iet
.. _Joint Research Centre: https://ec.europa.eu/jrc/en




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

