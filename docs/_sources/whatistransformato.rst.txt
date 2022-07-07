Theoretical background
==========================
 
:math:`\texttt{TRANSFORMATO}` [#fspeedlimits]_ is a python package designed to calculate
either **relative binding free energies** (RBFE) or **relative solvation free energies** (RSFE) between two similar
ligands. The two ligands are the two physical endstates which are mutated to a common core structure.
 
For that reason, |trafo| first searches for the common core (Common Core approach). All atoms which are not part of this
CC and thus differ between the two ligands are then turned off step-wise. First,
electrostatic contributions of the non-CC atoms are turned off, followed by turning off the Lennard-Jones (LJ) interactions
of the non-CC atoms. For the non-hydrogen atoms this is done on an atom-by-atom basis (Serial-Atom-Insertion approach).
This is done in both environments, for RSFE; ligand in vacuum and ligand solvated in water. For RBFE; ligand solvated in
water and ligand solvated in water and bound to the receptor (protein).
 
|trafo| will create intermediate state folders to connect each ligand with its common core structure.
Each folder will contain a plain MD simulation, which can be run individually from all other simulations. Finally, for
each ligand the free energy according to the common core structure can be calculated. Since both common core structures
are the same, the resulting free energy difference :math:`\Delta\Delta G^{bind}_{L1\rightarrow L2}` can  be obtained.
 
Interplay with CHARMM-GUI
###########################
 
|trafo| works closely together with CHARMM-GUI. For starting RBFE simulations, one needs output from the CHARMM-GUI
solution builder.
 
.. figure:: assets/images/tf_overview_v3.png
   :scale: 8%
   :alt: scheme of transformato
      
   Figure: Interplay between |trafo| and CHARMM-GUI.
 
The common-core approach
###########################
 
With standard alchemical FES, the two endpoints are generally transformed directly into each other,
with nonbonded forces (electrostatics) scaled according to the coupling parameter :math:`\lambda`,
and bonded LJ interactions being turned off on an atom by atom basis. [#fshirts]_\ .
 
With |trafo| however, the structures are not transformed into each other directly.
Instead, each structure is transformed into a 'common core structure', a common topology of both systems
(ideally - and usually - the *maximum* common topology between systems), with the free energy
being the sum of free energies necessary to reach the common core:
 
.. math::
  
   \Delta\Delta G^{bind}_{L1\rightarrow L2} = \Delta\Delta G^{bind}_{L1\rightarrow |DL_1| - R_{CC}} - \Delta\Delta G^{bind}_{L2\rightarrow\ |DL_2| - R_{CC}}
 
:math:`|DL_1|` and :math:`|DL_2|` refer to the non-CC (dummy region) of Endstate 1 and 2, respectively.
:math:`R_{CC}` refers to the common core.
 
.. figure:: assets/images/partA.png
   :alt: alchemical pathway
      
   Figure: Alchemical pathway as implemented in the common core approach.
   Free energies are calculated relative to a common core structure in the
   two environments (for RBFE: ligand solvated in water and ligand bound to a protein and solvated in water).
   This common core structure of the two ligands differs only in the number of dummy atoms,
   thus, the contribution cancels out when calculating relative free energies.
 
Serial-Atom-Insertion approach
################################
 
Lennard-Jones interactions are turned off on an atom-by-atom basis. This means that LJ parameters are either fully
interacting (1) or non-interacting (0). For the heavy atoms (non hydrogen atoms) the interaction is turned off
for each atom one by one.
 
 
Citations, License and Contact
##################################
 
Transformato is released under the MIT License. For further information, visit the `repository license page <https://github.com/wiederm/transformato/blob/master/LICENSE>`_\ .
 
Acknowledgements
####################
 
Apart from the default python packages, |trafo| relies in large parts on the following packages:
 
+ MDAnalysis\ [#fMDAnalysis1]_ [#fMDAnalysis2]_ : Used for general purpose protein and ligand analysis, atom selection and especially trajectory analysis.
 
+ Numpy\ [#fNumpy1]_ : Well, it is numpy. We use it for math.
 
+ parmed\ [#fparmed1]_ : Protein structure manipulations.
 
+ pymbar\ [#fpymbar]_ : mBAR analysis.
 
+ RDKit: Visualisation of molecular structures.
 
.. rubric:: References
 
 
.. [#fspeedlimits] Wieder, M., Fleck, M., Braunsfeld, B., and Boresch, S. (2022). *Alchemical free energy simulations without speed limits. A generic framework to calculate free energy differences independent of the underlying molecular dynamics program.* J. Comput. Chem. 43, 1151–1160, `DOI ⤶ <https://doi.org/10.1002/jcc.26877>`_
 
.. [#fboreschbruckner] Boresch, S., Bruckner, S. (2011). *Avoiding the van der Waals endpoint problem using serial atomic insertion* J.Comput. Chem. 2011,32, 2449–2458, `DOI ⤶ <https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.21829>`_
 
.. [#fshirts] Shirts, M. R., Mobley, D. L. and Chodera, J. D. (2007). *Alchemical free energy calculations: Ready for prime time?*  Annual Reports in Computational Chemistry, Vol. 3
 
.. [#fjohannes] Karwounopoulos, J., Wieder, M., and Boresch, S. (2022). *Relative binding free energy calculations with Transformato: a molecular dynamics engine-independent tool.* Frontiers, *submitted*.
 
.. [#fMDAnalysis1] Gowers, R. et al., 2016. *MDAnalysis: A Python Package for the Rapid Analysis of Molecular Dynamics Simulations*. In Proceedings of the Python in Science Conference.  SciPy.
 
.. [#fMDAnalysis2] Michaud-Agrawal, N. et al., 2011. *MDAnalysis: A toolkit for the analysis of molecular dynamics simulations*. Journal of Computational Chemistry, 32(10), pp.2319–2327.
 
.. [#fNumpy1] Van Der Walt, S., Colbert, S.C. & Varoquaux, G., 2011. *The NumPy array: a structure for efficient numerical computation*. Computing in Science & Engineering, 13(2), pp.22–30.
 
.. [#fpymbar] Shirts MR and Chodera JD. *Statistically optimal analysis of samples from multiple equilibrium states.* J. Chem. Phys. 129:124105 (2008). `DOI ⤶ <http://dx.doi.org/10.1063/1.2978177>`_
 
.. [#fparmed1]  Michael R. Shirts, Christoph Klein et al., *2016. Lessons learned from comparing molecular dynamics engines on the SAMPL5 dataset*. Journal of Computer-Aided Molecular Design
 
 
 
.. |trafo| replace:: :math:`\texttt{TRANSFORMATO}`
 
 
 

