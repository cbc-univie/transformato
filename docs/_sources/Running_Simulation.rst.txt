

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

If you take  a look at the intst*/ directories now created (located at :code:`project-folder/replicate-folder/combinedstructure/singlestructure/intst*`) you'll find two scripts: :code:`simulation.sh` and :code:`simulation_charmm.sh`

Somewhat unsuprisingly, these are responsible for running the simulation as either openMM or charmm, containing the required information and arguments.

.. important::
    You only need to run *one* of the options below. Please note, however, that CHARMM does not have the same features as openMM. If you need to modify these scripts for all future use in some way, you may find their sources in :code:`transformato/bin`

openMM
********

For openMM, use :code:`simulate.sh`:

.. code-block:: bash

    #!/bin/bash
    #SBATCH -p lgpu
    #SBATCH --gres=gpu


    source ~/anaconda3/etc/profile.d/conda.sh
    conda activate fep

    path=$1

    cd ${path}
    pwd
    hostname

    input=lig_in_complex
    init=lig_in_complex
    pstep=lig_in_complex
    istep=lig_in_complex
    irst=lig_in_complex
    orst=lig_in_complex_rst
    python -u openmm_run.py -i ${input}.inp -t toppar.str -p ${init}.psf -c ${init}.crd -irst ${irst}.rst -orst ${irst} -odcd ${istep}.dcd &> complex_out.log 

    input=lig_in_waterbox
    init=lig_in_waterbox
    pstep=lig_in_waterbox
    istep=lig_in_waterbox
    irst=lig_in_waterbox
    orst=lig_in_waterbox_rst
    python -u openmm_run.py -i ${input}.inp -t toppar.str -p ${init}.psf -c ${init}.crd -irst ${irst}.rst -orst ${irst} -odcd ${istep}.dcd  &> waterbox_out.log 

Essentially, this file just tells openMM what to run and then runs it, providing a bit of debug information along the way. Importantly, it takes **the working directory as argument**.

So, if you want to simulate this intermediate state, you would run:

.. code-block:: bash
    
    ./simulation.sh /absolute/or/relative/path/to/that/intstate

.. note:: Assuming your current working directory is that intstate, you can just supply :code:`./` as argument.


CHARMM
*******

For CHARMM, use :code:`simulate_charmm.sh`

.. code-block:: bash

    #!/bin/bash
    #SBATCH -p lgpu
    #SBATCH --gres=gpu
    path=$1
    SWITCH=$2

    cd ${path}
    pwd
    hostname


    run_complex () {
    input=charmm_run_complex
    OMP_NUM_THREADS=8 ${CHARMM} -i ${input}.inp > log_complex.out
    }

    run_waterbox () {
    input=charmm_run_waterbox
    OMP_NUM_THREADS=8 ${CHARMM} -i ${input}.inp > log_solv.out
    }


    case ${SWITCH} in
    1)
    run_complex
    ;;
    2)
    run_complex
    run_waterbox
    ;;
    esac

Importantly, unlike openMM this takes **two arguments**: The path, same as before, and a CASE statement. You need to run both, so that will always be 2 for you. You also need to have the ${CHARMM} system variable defined.

So, to run the simulation using CHARMM:

.. code-block:: bash

    ./simulate_charmm.sh /absolute/or/relative/path/to/your/intstate 2


Automation and offloading to a cluster
****************************************

There are two methods to automate running simulations:

+   Running it from the submit.ipynb
    
    Has the advantage of no extra step being necessary, but pythons multiprocessing facilities can be capricious. This variant also does not allow inspection of the intermediate states before submission.

    If you'd like this option, add this code to the bottom of your submit.ipynb:

    .. code-block:: python

        import os
        import glob
        import subprocess

        wd = os.getcwd()
        output_files = glob.glob(wd + f"/{folder}/*/*/intst*", recursive = True) 
            #Whether {folder} is necessary depends on your folder setup relative to your working directory -
            # in general, this should point to your intst** folders

        for path in sorted(output_files):
            # because path is object not string
            print(f"Start sampling for: {path}")
            exe = subprocess.Popen(["ssh", "your clustername", "sbatch", f"{str(path)}/simulation.sh", str(path)], text=True, stdout=subprocess.PIPE )
            output = exe.stdout.read()
            print(output)

    You will have to replace 'your-clustername' and modify the paths according to your setup.

+   Running it via script

    Has the advantage of allowing inspection beforehand and being more reliable, as well as allowing the modification of cluster instructions during submit.

    To use, create a script similar to this:

    .. code-block:: bash

        for {i} in ./**/**/intst**/; #This assumes you start from the directory above the initial structure directory - you may change as necessary.
        do cd ${i}; # to prevent issues, it is preferred to switch to the script directory
        sbash simulation.sh; # sbash being the command for the workload manager slurm - you may need to replace as necessary.
        cd ../../..; # switch working directory back to original so next loop starts properly
        done;

    If you run this from the folder containing the replicate directories, any correctly built transformato replicate should be submitted to run. You can modify the loop glob to restrict the script to certain directories.


