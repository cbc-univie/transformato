

Running Simulations
======================

.. _rst_submitfiledesc:

.. seealso::
    If you are unsure or need inspiration, the `transformato/notebooks <https://github.com/wiederm/transformato/tree/master/notebooks>`_ folder has a number of detailed examples.

After setting up a system, actually running it is comprised of two steps:

+   Creating the intermediate states

    During this stage you may also make final adjustments to the common core and dummy regions of your system.

+   Simulating the system

    Here, you will essentially conducting small MD simulations for complex and waterbox in each of the intermediary states. This will generate the trajcetory files to be analysed later.


Creating intermediary states
#################################


.. tip::
    We recommend using an ipython interpreter like jupyter-notebook or VSCode to run the notebooks, as this will enable you to easier visualize the common-core and make the necessary adjustments.

Loading structures and generating forcefields
**********************************************

.. code-block:: python

    from transformato import load_config_yaml, SystemStructure, IntermediateStateFactory, ProposeMutationRoute
    from transformato.mutate import perform_mutations
    from IPython.display import SVG
    from transformato.utils import run_simulation, postprocessing
    import warnings
    warnings.filterwarnings("ignore", module='parmed')

    configuration = load_config_yaml(config='24to25.yaml',input_dir='../', output_dir='./longrun-norestraints-1/')

During this part, we load the required modules and tell transformato to load our config.yaml. Please note input_dir and output_dir as additional arguments for load_config_yaml: These folders are where to find you structures and where to write your intermediate states, respectively.

Generating the mutation route and common core
**********************************************

.. code-block:: python
    
    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)

These commands tell transformato to create a mutation route. It'll take some time for CGenFF to create parameter files for the proposed route.

.. code-block:: python

    s1_to_s2.propose_common_core()
    SVG(s1_to_s2.show_common_core_on_mol1())
    SVG(s1_to_s2.show_common_core_on_mol2())

This will instruct transformato to create and assign both the common core and dummy regions. Using the SVG() order you can then take a look at your structure and decide if the common core is proper or if adjustments need to be made.

.. code-block:: python

    s1_to_s2.add_idx_to_common_core_of_mol1([34,33,36,38,12,13,37,10,35,8,9,11,7,2])
    s1_to_s2.add_idx_to_common_core_of_mol2([34,33,36,38,12,13,37,10,35,8,9,11,7,2])

This will instruct transformato to add the atoms with the given idxs to the common core. The idxs are identical to the numbers visualized in the SVG. Manually added atoms are not added directly, but staged until finish_common_core() is called.

.. code-block:: python

    s1_to_s2.finish_common_core()

transformato will add staged idxs and start attempting to scale charges. If you want, you may now repeat the SVG() command from above to see the changes you made.

Create mutations and write intermediate states
*************************************************

.. code-block:: python

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    print(mutation_list.keys())
    i = IntermediateStateFactory(
    system=s1,
    configuration=configuration,
    )

This will generate the necessary mutation list from one endstate to the common core.

.. code-block:: python

    perform_mutations(nr_of_mutation_steps_charge=3, configuration=configuration, i=i, mutation_list=mutation_list)

With this command, the actual intermediate state directories are written. After this has finished without errors, you may proceed to actually running the simulation.


Running the simulation
########################

.. note::
    It is *technically* possible to run all simulations locally, as long as you have a CUDA - capable device. However, even with top-end hardware, expect simulation and analysis to take *at least* 15 hours per replicate, depending on timestep, system size, simulation length and especially the number of intermediary states required.

