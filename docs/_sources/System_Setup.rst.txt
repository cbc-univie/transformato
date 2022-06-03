System Setup
################


After successfully installing transformato, you can now setup your first system. For this, you require three items:

- The output of [CHARMM-GUI](http://www.charmm-gui.org/)'s solution builder for your ligand - solvated once at its position within the protein-ligand complex and once in a pure waterbox.

- A `config.yaml` which describes your general simulation parameters


- A `submit.ipynb` (or just standard .py) that generates your intermediate states

Finally, to actually analyze your simulations and calculate the resulting free energies, you'll need an `analysis.py` script - this is covered on the *Analysis* page.

File structure and supplying your protein/ligand
***************************************************

Please refer to this convenient and practical flowchart to understand transformato's folder structure:
.. image:: https://user-images.githubusercontent.com/79014444/169664208-56cc9d73-bd93-49ed-8f06-14cc15c04f34.svg


While this may look complicated, all you need to do is to:
1. Create a folder called `your-structure-1` for your first structure
2. Take the CHARMM-GUI output folder (called something like charmm-gui-4842148) for your solvated ligand-protein complex, move it to this folder and rename it into `complex`
3. Take the CHARMM-GUI output folder for your solvated ligand in the waterbox, move it to this folder and rename it inti `waterbox`
4. Do the same with your second, third, fourth... etc. structure
5. **Run the equillibration** part of the provided simulation input. An equillibrated system (and the generated .rst file) are necesary for transformato.

You`ll create the remaining files down below.

## The config.yaml

## The submit script



Restraints
************

WARNING: The documentation here is not yet implemented and more aspirational/for development branches

**OpenMM only**

Transformato supports two types of restraints: automatic and manual restraints.

Automatic Restraints
######################

To activate automatic restraints for your ligand, add 

.. code-block:: yaml

    simulation:
        restraints: "auto"
        
to your `config.yaml`. This wil restrain your ligand in its original position using a harmonic potential of the form `E=0.5*k*(r-r0)²` (where  r is the current distance, r0 the initial distance, and k the force constant), applied within the centers of mass of your ligand and the surrounding protein structure, keeping these two vaguely fixed.

You may specify further keywords:

*Keywords for restraints:*
|Keyword|Values|Description|
|-------|-------|----|
|`k=`|`[int]`|*(optional)* Specifies the force constant. Default value is 3. Sensible values are between 0.1 - 500, depending on the system, but tougher restraints have larger impact on system sampling and accuracy of free energy calculations. Default is `k=3`
|`extremities=` |`[int]`|*(optional)* Instead of restraining the entire common core by its center of mass, attaches `n` bonds to the calculated extremities of the ligand and nearby protein structures, thus cutting down on rotational freedom. If not specified, the entire ligand is restrained by its center of mass
|`shape=`|`harmonic` `flatbottom`|*(optional)* Exchanges the harmonic potential underlying the restraints for a flat-bottomed one of the form `step(r-r_0) · (K/2)·(r-r_0)^2` Default is `harmonic`
|`scaling`|Presence toggles|*(optional)* Linearly scales the k value through the first 4 states (0, 0.25, 0.5, 0.75). Default is `False`|

A full command might thus look like this:



.. code-block:: yaml

    restraints: "auto k=10 extremities=3 shape=harmonic scaling" 



CAUTION: Be somewhat sure of what your structure looks like, and do a sanity check on the generated restraints before production. As all restraints only act on the common core, setting an arbitrarily high number of extermities can lead to strange results

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
########################

To activate manual restraints for your ligand, add 

*config.yaml*

.. code-block:: yaml

    simulation:
        restraints: "manual"

to your config.yaml. Below, you may now specify an arbitrary number of restraints using the [MDAnalysis selection syntax](https://docs.mdanalysis.org/stable/documentation_pages/selections.html#simple-selections):

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

Code example with many restraints:

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
                    


|Keyword|Values|Description|
|-------|-------|-------------|
|`shape`|`"harmonic"` `"flatbottom"`|Shape of the energy potential. Default is "harmonic"|
|`group1` `group2`|Any [MDAnalysis selection string](https://docs.mdanalysis.org/stable/documentation_pages/selections.html#simple-selections)|Defines which Common Core atoms are members of group1 or group2. Please note that group1 **must** be the ligand, and group2 the protein.|
|`k`|`[float]`|*(optional)*: Defines the harmonic force constant. Default is 3.|
|`r0`|`[float]`|*(optional)*: **Harmonic Potential:** By default, forces are 0 when the atom groups have their initial distance to each other. By setting this parameter, you may define a custom r0, allowing you to push or pull atom groups.  **Flat-Bottom Potential:** Represents the well size. Defaults to 1.5

As with automatic restraints, even manually specified restraints will never act on atoms not in the common core, as this would lead to nonsensical energy calculations.
