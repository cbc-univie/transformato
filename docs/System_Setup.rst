System Setup
===============


After successfully installing transformato, you can now setup your first system. For this, you require three items:


+ The output of `CHARMM-GUI <http://www.charmm-gui.org/>`_ 's solution builder for your ligand - solvated once at its position within the protein-ligand complex and once in a pure waterbox.
+ A `config.yaml` which describes your general simulation parameters.
+ A `submit.ipynb` (or just standard .py) that generates your intermediate states.


.. hint:: 
    It is not strictly necessary to use CHARMM - GUI to solvate your system. You may build and solvate your system yourself. In this case however, you will also need to supply proper parameter files along with input scripts.

Finally, to actually analyze your simulations and calculate the resulting free energies, you'll need an `analysis.py` script - this is covered on the *Analysis* page.

File structure and supplying your protein/ligand
#####################################################

Please refer to this convenient and practical flowchart to understand transformato's folder structure:

.. image:: assets/images/transformato-tree.svg
    :alt: transformato-tree

While this may look complicated, all you need to do is to:

#. Create a folder called `your-structure-1` for your first structure
#. Take the CHARMM-GUI output folder (called something like charmm-gui-4842148) for your solvated ligand-protein complex, move it to this folder and rename it into `complex`
#. Take the CHARMM-GUI output folder for your solvated ligand in the waterbox, move it to this folder and rename it inti `waterbox`
#. Do the same with your second, third, fourth... etc. structure
#. **Run the equillibration** part of the provided simulation input. An equillibrated system (and the generated .rst file) are necesary for transformato.

You`ll create the remaining files down below.

.. note:: 
    The submit.ipynb file is covered in :doc:`Running_Simulation`

The config.yaml
#################

The config.yaml is perhaps the most important file you'll create. It contains:

#. Definitions of your structures and under which name to find them
#. Simulation parameters, such as simulation step and length
#. Eventual restraints, both automatic and manual





Restraints
###########

.. warning:: The documentation here is not yet implemented in the main branch yet and more aspirational/for development branches

.. danger:: This section only applies if you are running your simulations with openMM. Should you run your simulations using CHARMM, it will not apply the restraints **and give no warning about it**.



Transformato supports two types of restraints: automatic and manual restraints.

Automatic Restraints
**************************

To activate automatic restraints for your ligand, add 

.. code-block:: yaml

    simulation:
        restraints: "auto"
        
to your `config.yaml`. This wil restrain your ligand in its original position using a harmonic potential of the form :math:`E=0.5k \cdot (r-r_0)²` (where  :math:`r`` is the current distance, :math:`r_0`` the initial distance, and :math:`k`` the force constant), applied within the centers of mass of your ligand and the surrounding protein structure, keeping these two vaguely fixed.

You may specify further keywords:

.. rubric:: Options for automatic restraints
    
restraints:
    :code:`auto`
        Required. Tells transformato to look for and apply automatic restraints
    :code:`k=[int]`
        *optional:* Defines the spring constant used for force calculations. Default is 3
    :code:`extremities=[int]`
        *optional:* If used, transformato will not restraint the entire common core but rather look for [int] extremities. These are then restrained to the surrounding protein carbon-alphas.
    :code:`shape="harmonic" ["harmonic","flatbottom"]`
        *optional:* Defines the shape of the energy potential. Only "harmonic" is currently implemented.
    :code:`scaling`
        *optional:* If present, the k - Value of the force is linearly scaled in the first four intermediate states (0,0.25,0.5,0.75)

A full command might thus look like this:



.. code-block:: yaml

    restraints: "auto k=10 extremities=3 shape=harmonic scaling" 



.. caution:: Be somewhat sure of what your structure looks like, and do a sanity check on the generated restraints before production. As all restraints only act on the common core, setting an arbitrarily high number of extermities can lead to strange results

It should be noted that this means that a small file called `restraints.yaml` is created in your `intst*` - folders.
These have the following structure:


.. code-block:: yaml

    system:
        structure:
            tlc: LIG # same as in the config.yaml, but only one structure (as only one relevant)

    simulation:
        restraints: "auto" # same as in config.yaml
        ccs:  # this represents an array of your common core, upon which restraints can be applied
            - C1
            - C2
            - H2
    intst:
        scaling:0.8 # for non-immediate switches, how far along the scaling is. Only relevant for harmonic potentials.


It is not recommended to manually edit these files, as they are automatically created for each intermediate state.

Manual Restraints
*******************

To activate manual restraints for your ligand, add 

*config.yaml*

.. code-block:: yaml

    simulation:
        restraints: "manual"

to your config.yaml. Below, you may now specify an arbitrary number of restraints using the `MDAnalysis selection syntax <https://docs.mdanalysis.org/stable/documentation_pages/selections.html#simple-selections>`_ :

*config.yaml*

.. code-block:: yaml

    simulation:
        restraints: "manual"
        manualrestraints:
            restraint1:
                shape: "harmonic"
                group1: "resname LIG and type C"
                group2: "protein and type CA"
                k: 30
                r0: 2.41

You may define as many restraints as you like:

Code example with multiple restraints:

*config.yaml*

.. code-block:: yaml

    simulation:
        restraints: "manual"
        manualrestraints:
            restraint1:
                shape: "harmonic"
                group1: "resname LIG and type C"
                group2: "protein and type CA"
            restraint2:
                shape: "flatbottom"
                group1: "resname LIG and type C"
                group2: "protein and type CA"
            restraint3:
                shape: "harmonic"
                group1: "resname LIG and name C14"
                group2: "sphlayer 5 15 name C14 and protein and type CA"

Note that the individual restraints all need to have distinct names (restraint1, restraint2 etc.). It is common that they are numbered, but not required - they simply need to adhere to the yaml syntax.

.. rubric:: Options for manual restraints

manualrestraints
    :code:`shape="harmonic" ["harmonic","flatbottom"]'`
        Shape of the energy potential. Default is "harmonic", "flatbottom" is not yet implemented
    :code:`group1,group2=[MDAnalysis selection string]`
        Defines which Common Core atoms are members of group1 or group2. Please note that group1 **must** be the ligand, and group2 the protein.
    :code:`k=[int]`
        *(optional):* Defines the harmonic force constant. Default is 3.


As with automatic restraints, even manually specified restraints will never act on atoms not in the common core, as this would lead to nonsensical energy calculations.
