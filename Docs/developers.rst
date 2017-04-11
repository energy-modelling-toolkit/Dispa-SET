.. _developers:

Developers' section
===================


Folders organization
--------------------

* GAMS code and scripts are included within the "GAMS-files" folder. The code should only be modified in that folder!
* Python code and scripts are included within the "DispaSET" folder. 
* DispaSET configuration files are included within the "ConfigFiles" folder (one file per simulation).
* Input files for each country are stored in the "Database" folder
* Simulation directories can be written into the "Simulations" folder, which is not tracked by git
 

A sample Simulation Environment Folder is available at SimulationReferenceTestCase. The files in this folders are generated automatically by the pre-processing scripts. They should not be modified. The folder should be replaced and updated at every major Dispa-SET release.

By default, the pre-processing scripts are set to generate the simulation environment within the "Simulations" folder. This folder is excluded in .gitignore and can therefore be used for testing purposes.


Math equations in the Docs
--------------------------

- To use online mathjax (default), there is nothing to do but displaying the equation requires an internet connection. In linux, the Makefile call to build with mathjax is written:: 

	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html

	In conf.py, the only math equation should be: 'sphinx.ext.mathjax'	
- To use pngmath (for Linux)::
	
	sudo apt-get install dvipng
	
	In conf.py, add 'sphinx.ext.pngmath' in the extensions

	in Makefile: $(SPHINXBUILD) -D pngmath_latex=latex -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
- To use mathjax in offline mode, download the latest release from the mathjax github, copy it to Docs/_static/, and include it in the conf.py::

	mathjax_path = "Mathjax/MathJax.js"  

Using Autodoc
-------------
The "API" section of the Docs uses uses the sphinx autodoc extension to scan the source code of Dispa-SET and display the relevant functions together with their description, parameters and outputs.
In the Sphinx "conf.py", the path to the source file must be added::

	sys.path.insert(0, os.path.abspath('../DispaSET'))

If the API documentation is generated with sphinx-apidoc, from the Docs folder, use::

	sphinx-apidoc -o . ../DispaSET/

Add a link to "DispaSET" in the table of content of index.rst and include the root folder in conf.py::

	sys.path.insert(0, os.path.abspath('../'))


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


Issue with the compilation of the GAMS API
------------------------------------------
First, check that the installed version of GAMS is the 64 bit. 32 bit versions tend to generated compatibility issues.

When the pre-compiled libraries do not work, they must be re-compiled from the GAMS apifile folder. In Windows, this generally raises the issue of the missing vcvarsall.bat file. If the issue persists after installing the Microsoft C++ compiler for Python 2.7, try the following:

1. Enter MSVC for Python command prompt
2. SET DISTUTILS_USE_SDK=1
3. SET MSSdk=1
4. python.exe gdxsetup.py install


Public version of Dispa-SET
---------------------------
Because some input files are subject to intellectual property and copyrights, some folders available on the private repository cannot by uploaded to the public GitHub repository. The script DispaSync.sh in the root folder has been written to synchronize the subset of publicly available folders and files with an external folder. This folder can then be committed and pushed to the public repository. 

By default, all file and folders are synchronized. In order to add a private path (to a file or to a folder), edit the DispaSync.sh file and add an entry to the "--exclude" argument.

The "rsync" software is required for the synchronization and must be installed on the local machine. The script can be run in any UNIX terminal (i.e. it cannot be run in Windows).
