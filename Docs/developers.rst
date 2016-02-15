.. _developers:

Developers
===========
Developers' section.


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
