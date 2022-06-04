.. transformato documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to transformato's documentation!
=========================================================

:code:`transformato` is a lightweight implementation of the Boresch-Braunsfeld Serial-Atom-Insertion Common-Core (SAI-CC) approach [#f1]_ to relative binding/solvation free energy calculations [#f2]_. 
It does so by connecting the physical endstates of two molecules via two separate alchemical paths. The common core is defined as the maximum topology mapping between the two molecular graphs.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Installation
   System_Setup
   Running_Simulation
   Analysis
   Troubleshooting
   api
   







.. rubric:: References

.. [#f1] Wieder, M., Fleck, M., Braunsfeld, B., and Boresch, S. (2022). *Alchemical free energy simulations without speed limits. A generic framework to calculate free energy differences independent of the underlying molecular dynamics program.* J. Comput. Chem. 43, 1151–1160
.. [#f2] Karwounopoulos, J., Wieder, M., and Boresch, S. (2022). *Relative binding free energy calculations with Transformato: a molecular dynamics engine-independent tool.* frontiers, *submitted*.

.. rubric:: Maintainers

+ Markus Wieder, Department of Pharmaceutical Sciences, University of Vienna
+ Stefan Boresch, Department of Computational Biological Chemistry, University of Vienna
+ Johannes Karwounopoulos, Department of Computational Biological Chemistry, University of Vienna

* :ref:`genindex`
