.. _developers:

Developers
===========
Developers' section.


Folders organization
--------------------

- GAMS code and scripts are included within the "GAMS-files" folder. The code should only be modified in that folder!
- Python code and scripts are included within the "python-files" folder. 
- Pre-processing scripts and raw data are included within the "Pre-processing" folder.

A sample Simulation Environment Folder is available at SimulationReferenceTestCase. The files in this folders are generated automatically by the pre-processing scripts. They should not be modified. The folder should be replaced and updated at every major Dispa-SET release.

By default, the pre-processing scripts are set to generate the simulation environment within the "Simulation" folder. This folder is excluded in .gitignore and can therefore be used for testing purposes.


Math equations in the Docs
--------------------------

- To use online mathjax (default), there is nothing to do but displaying the equation requires internet connection
- To use pngmath (for Linux)::
	
	sudo apt-get install dvipng
	
	In conf.py, add 'sphinx.ext.pngmath' in the extensions

	in Makefile: $(SPHINXBUILD) -D pngmath_latex=latex -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html


Clone git repository to svn
---------------------------
* Install git svn ::

	sudo apt-get install git-svn
* Create svn branch ::

	git branch svn
* Connect to svn repository ::

	git svn init -s --prefix=svn/ --username <user> https://joinup.ec.europa.eu/svn/dispaset
* Checkout svn banch ::

	git checkout svn
* Fetch remote content ::

	git svn fetch
* Reset repository ::

	git reset --hard remotes/svn/trunk
* Merge master into svn ::

	git merge master
* Commit to the remote repository :: 

	git svn dcommit

The folder can finally optionally be deleted to avoid any confusion.
