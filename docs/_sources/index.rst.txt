.. transformato documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to transformato's documentation!
=========================================================

:code:`transformato` is a lightweight implementation of the Boresch-Bruckner Serial-Atom-Insertion\ [1]_ Common-Core (SAI-CC) approach\ [2]_ to relative binding/solvation free energy calculations\ [3]_. 
It does so by connecting the physical endstates of two molecules via two separate alchemical paths. 


It requires very little set-up time and is designed to work directly with output from `CHARMM - GUI <http://www.charmm-gui.org>`_\ , 
making it a quick and simple way to run accurate RBFE/RSFE calculations in either CHARMM or openMM.


If you'd like to take a look at the code, head over to our `github repository <https://github.com/wiederm/transformato>`_\ . 


.. toctree::
   :maxdepth: 2
   :caption: Contents:
   

   whatistransformato
   Installation
   System_Setup
   Running_Simulation
   Analysis
   Troubleshooting
   api
   







.. rubric:: References


.. [1] Boresch, S., Bruckner, S. (2011). *Avoiding the van der Waals endpoint problem using serial atomic insertion* J.Comput. Chem. 2011,32, 2449–2458, `DOI ⤶ <https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.21829>`_

.. [2] Wieder, M., Fleck, M., Braunsfeld, B., and Boresch, S. (2022). *Alchemical free energy simulations without speed limits. A generic framework to calculate free energy differences independent of the underlying molecular dynamics program.* J. Comput. Chem. 43, 1151–1160, `DOI ⤶ <https://doi.org/10.1002/jcc.26877>`_

.. [3] Karwounopoulos, J., Wieder, M., and Boresch, S. (2022). *Relative binding free energy calculations with Transformato: a molecular dynamics engine-independent tool.* frontiers, *submitted*.

.. rubric:: Maintainers

+ Markus Wieder, Department of Pharmaceutical Sciences, University of Vienna
+ Stefan Boresch, Department of Computational Biological Chemistry, University of Vienna
+ Johannes Karwounopoulos, Department of Computational Biological Chemistry, University of Vienna

* :ref:`genindex`
