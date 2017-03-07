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

The Dispa-SET model is a unit commitment and dispatch model developed within the “Joint Research Centre” and focused on the balancing and flexibility problems in European grids [1]_. It is written in GAMS an Python (Pyomo) and uses csv files for input data handling. The optimisation is defined as a Linear Programming (LP) or Mixed-Integer Linear Programming (MILP) problem, depending on the desired level of accuracy and complexity. 



Libraries used
--------------

* `pyomo`_ Optimization object library, interface to LP solver (e.g. CPLEX)
* `pandas`_ for input and result data handling 
* `matplotlib`_ for plotting
* `GAMS_api`_ for the communication with GAMS

  
Contents
--------

.. toctree::
   :maxdepth: 3

   overview
   workflow
   model
   implementation
   data
   developers
   DispaSET


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

References
----------
.. [1] Quoilin, S., Hidalgo Gonzalez, I., & Zucker, A. (2017). Modelling Future EU Power Systems Under High Shares of Renewables: The Dispa-SET 2.1 open-source model. Publications Office of the European Union.  

.. _matplotlib: http://matplotlib.org
.. _pandas: http://pandas.pydata.org
.. _pyomo: http://www.pyomo.org/
.. _GAMS_api: http://www.gams.com/help/index.jsp?topic=%2Fgams.doc%2Fapis%2Findex.html
.. _Institute of Energy and Transport: https://ec.europa.eu/jrc/en/institutes/iet
.. _Joint Research Centre: https://ec.europa.eu/jrc/en




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

