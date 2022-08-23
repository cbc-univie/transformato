.. transformato documentation master file, created by
  sphinx-quickstart on Thu Mar 15 13:55:56 2018.
  You can adapt this file completely to your liking, but it should at least
  contain the root `toctree` directive.
 
Welcome to the documentation for |trafo|!
=========================================================
 
The python package |trafo| is an implementation of the Common Core / Serial-Atom-Insertion
(CC-SAI) approach\ [1]_ for calculating free energy differences\.
It does so by connecting the physical endstates of two molecules via alchemical pathways.
 
It requires very little set-up time and is designed to work directly with output from
`CHARMM - GUI <http://www.charmm-gui.org>`_\. Currently either **relative solvation free energies** (RSFE) [2]_
or **relative binding free energies** (RBFE) [3]_ can be calculated. For the production runs either
`CHARMM <https://academiccharmm.org/>`_\ or `openMM <https://openmm.org/>`_\ are required.
 
If you'd like to take a look at the code,
head over to our `github repository <https://github.com/wiederm/transformato>`_\ .
 
 
.. toctree::
  :maxdepth: 2
  :caption: Contents:
  :glob:
 
  whatistransformato
  Installation
  System_Setup
  Additional_Settings
  Running_Simulation
  Analysis
  Troubleshooting
  api
 
 
 
 
 
 
 
 
 
.. rubric:: References
.. [1] Boresch, S., Bruckner, S. (2011). *Avoiding the van der Waals endpoint problem using serial atomic insertion* J.Comput. Chem. 2011,32, 2449–2458, `DOI ⤶ <https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.21829>`_
 
.. [2] Wieder, M., Fleck, M., Braunsfeld, B., and Boresch, S. (2022). *Alchemical free energy simulations without speed limits. A generic framework to calculate free energy differences independent of the underlying molecular dynamics program.* J. Comput. Chem. 43, 1151–1160, `DOI ⤶ <https://doi.org/10.1002/jcc.26877>`_
 
.. [3] Karwounopoulos, J., Wieder, M., and Boresch, S. (2022). *Relative binding free energy calculations with Transformato: a molecular dynamics engine-independent tool.* Front. Mol. Biosic., *submitted*.
 
.. rubric:: Maintainers
 
+ `Markus Wieder <marcus.wieder@univie.ac.at>`_\ , Department of Pharmaceutical Sciences, University of Vienna
+ `Johannes Karwounopoulos <johannes.karwounopoulos@univie.ac.at>`_\, `Institute of Computational Biological Chemistry <https://www.mdy.univie.ac.at/index.html>`_, University of Vienna
+ `Stefan Boresch <stefan@mdy.univie.ac.at>`_\ , `Institute of Computational Biological Chemistry <https://www.mdy.univie.ac.at/index.html>`_, University of Vienna
 
 
.. |trafo| replace:: :math:`\texttt{TRANSFORMATO}`
 
 
* :ref:`genindex`
 
 
 

