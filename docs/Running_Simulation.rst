 
 
Running Simulations
======================
 
.. _rst_submitfiledesc:
 
.. seealso::
   If you are unsure or need inspiration,
   the `transformato/notebooks <https://github.com/wiederm/transformato/tree/master/notebooks>`_
   folder has a number of detailed examples.
 
After setting up a system, actually running it comprises two steps:
 
+   Creating the intermediate states
 
   During this stage, you may also make final adjustments to the common core and dummy regions of your system.
 
+   Simulating the system
 
   Here, you will essentially conduct small MD simulations for complex and
   waterbox in each of the intermediary states. This will generate the trajcetory files to be analysed later.
 
 
Step 3 - Creating intermediary states
######################################
 
We are now following the notebook shown here XXXXXXXXXXXXXXXX
 
.. tip::
   We recommend using an ipython interpreter like jupyter-notebook or VSCode
   to run the notebooks, as this will enable you to easier visualize the common-core
   and make the necessary adjustments.
 
 
#. Loading structures
 
   During this part, we load the required modules and tell transformato to load our config.yaml.
   Please note, input_dir and output_dir as additional arguments for load_config_yaml:
   These folders are where to find your structures and where to write your intermediate states, respectively.
 
   .. code-block:: python
 
       from transformato import load_config_yaml, SystemStructure, IntermediateStateFactory, ProposeMutationRoute
       from transformato.mutate import perform_mutations
       from IPython.display import SVG
       from transformato.utils import run_simulation, postprocessing
       import warnings
       warnings.filterwarnings("ignore", module='parmed')
 
       configuration = load_config_yaml(config='24to25.yaml',input_dir='../', output_dir='./longrun-norestraints-1/')
 
   .. object:: load_config_yaml
 
       These parameters need to be adjusted by the user
      
       .. option:: config
 
           path to your config yaml file
 
       .. option:: input_dir
 
           The path where your input structures are stored
 
       .. option:: output_dir
 
           The path to a folder where you want to write your intermediat states and later run the simulations
 
 
#. Generating the mutation route and common core
 
   Now the two ligands will be processed and the alchemical mutation route connecting the two ligands to their
   common core structure via several intermediate states will be created.
 
 
   .. code-block:: python
      
       s1 = SystemStructure(configuration, "structure1")
       s2 = SystemStructure(configuration, "structure2")
       s1_to_s2 = ProposeMutationRoute(s1, s2)
       s1_to_s2.propose_common_core()
 
 
#. Visualize the Common Core
 
   Now the common core suggested by |trafo| can be shown using the following command
 
   .. code-block:: python
 
           SVG(s1_to_s2.show_common_core_on_mol1())
           SVG(s1_to_s2.show_common_core_on_mol2())
 
   If one does **not** like the suggestion, one can manually interfere and either add atoms to the common core
 
   .. code-block:: python
 
       s1_to_s2.add_idx_to_common_core_of_mol1([idx1,idx2,...])
       s1_to_s2.add_idx_to_common_core_of_mol2([idx1,idx2,...])
 
   or remove atoms:
 
   .. code-block:: python
 
       s1_to_s2.remove_idx_to_common_core_of_mol1([idx1,idx2,...])
       s1_to_s2.remove_idx_to_common_core_of_mol2([idx1,idx2,...])
 
   In both cases, the idx of the respective atoms can be found in the graphic shown beforehand.
 
   .. attention::
 
       If you add or remove atoms for both structures, make sure
       you are using the correct (the same) order for the both structures
 
 
#. Finalize the common core
 
   .. code-block:: python
 
       s1_to_s2.finish_common_core()
 
   |trafo| will add staged idxs and start attempting to scale charges.
   If you want, you may now repeat the SVG() command from above to see
   the changes you made.
 
 
#. Create mutations and write intermediate states
 
 
   .. code-block:: python
 
       mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
       i = IntermediateStateFactory(
       system=s1,
       configuration=configuration,
       )
 
   This will generate the necessary mutation list from one endstate to the common core.
 
   .. code-block:: python
 
       perform_mutations(nr_of_mutation_steps_charge=3, configuration=configuration, i=i, mutation_list=mutation_list)
 
   With this command, the actual intermediate state directories are written. After this has finished without errors, you may proceed to actually running the simulation.
 
   .. object:: perform_mutations
 
       If needed, you can adjust the amount of intermediate steps here
      
       .. option:: nr_of_mutation_steps_charge (default: 5)
 
           how many intermediate states should be created to scale
           the charges of the atoms in the dummy region
 
       .. option:: nr_of_mutation_steps_cc (default: 5)
 
           how many intermediate states should be created to interpolate
           parameters between the two common core regions. This is only necessary for
           ligand 1!
 
#. Running the simulation
 
   Now the simulation can be started locally with:
 
   .. code-block:: python
 
       run_simulation(i.output_files, engine="openMM")
 
   .. object:: engine (["openMM","CHARMM"])
 
       you can decide whether you want to use openMM or CHARMM
 
       .. note::
 
           openMM is already available in the conda ``fep`` environment, CHARMM needs
           to be installed by oneself.
 
 
   .. attention::
 
       It is *technically* possible to run **RBFE** simulations locally but it can
       take a lot of time. For that reason, it is strongly suggested to offloade them to a
       supercomputer. Running **RSFE** on a local machine though, on a GPU can be done easily.
 
   If you take  a look at the intst*/ directories now created
   (located at :code:`project-folder/replicate-folder/combinedstructure/singlestructure/intst*`) you'll
   two scripts: :code:`simulation.sh` and :code:`simulation_charmm.sh`
 
   Somewhat unsurprisingly, these are responsible for running the simulation as
   either openMM or CHARMM, containing the required information and arguments.
 
   .. important::
       You only need to run *one* of the options below. Please note, however,
       that CHARMM does not have the same features as openMM. If you need to
       modify these scripts for all future use in some way, you may find their
       sources in :code:`transformato/bin`
 
   In each intermediate state directory (called intst1, intst2, ...), there is a simulation.sh file for the use with
   openMM and a simulation_charmm.sh file for the use with CHARMM.
 
   **Using openMM:**
 
   For openMM, use :code:`simulation.sh`:
 
   .. code-block:: bash
 
       #!/bin/bash
       #SBATCH -p lgpu
       #SBATCH --gres=gpu
 
 
       source ~/anaconda3/etc/profile.d/conda.sh # you might need to source your conda environment
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
 
   Essentially, this file just tells openMM what to run and then runs it,
   providing a bit of debug information along the way. Importantly, it
   takes **the working directory as argument**.
 
   So, if you want to simulate this intermediate state, you would run:
 
   .. code-block:: bash
      
       ./simulation.sh /absolute/or/relative/path/to/that/intermediatestate
 
   .. note::
      
       Assuming your current working directory is that intermediate state,
       you can just supply :code:`./` as argument.
 
 
   **Using CHARMM:**
 
   For CHARMM, use :code:`simulation_charmm.sh`
 
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
 
   Importantly, unlike openMM this takes **two arguments**: The path, same as before,
   and a CASE statement. You need to run both, so that will always be 2 for you.
   You also need to have the ${CHARMM} system variable defined.
 
   So, to run the simulation using CHARMM:
 
   .. code-block:: bash
 
       ./simulation_charmm.sh /absolute/or/relative/path/to/your/intermediatestate 2
 
 
#. Automation and offloading to a cluster
 
 
   There are two methods to automate running simulations:
 
   +   Running it from the ``submit.ipynb``
      
       Has the advantage of no extra step being necessary, but pythons multiprocessing facilities can be
       capricious. This variant also does not allow inspection of the intermediate states before submission.
 
       If you'd like this option, add this code to the bottom of your ``submit.ipynb``:
 
       .. code-block:: python
 
           import os
           import glob
           import subprocess
 
           wd = os.getcwd()
           output_files = glob.glob(wd + f"/{folder}/*/*/intst*", recursive = True)
               # Whether {folder} is necessary depends on your folder setup
               # relative to your working directory -
               # in general, this should point to your intst** folders
 
           for path in sorted(output_files):
               print(f"Start sampling for: {path}")
               exe = subprocess.Popen(["ssh", "your clustername", "sbatch", f"{str(path)}/simulation.sh", str(path)], text=True, stdout=subprocess.PIPE )
               output = exe.stdout.read()
               print(output)
 
       You will have to replace **your-clustername** and modify the paths according to your setup.
 
   +   Running it via script
 
       Has the advantage of allowing inspection beforehand and being more reliable,
       as well as allowing the modification of cluster instructions during submit.
 
       To use, create a script similar to this:
 
       .. code-block:: bash
 
           for {i} in ./**/**/intst**/; #This assumes you start from the directory above the initial structure directory - you may change as necessary.
           do cd ${i}; # to prevent issues, it is preferred to switch to the script directory
           sbash simulation.sh; # sbash being the command for the workload manager slurm - you may need to replace as necessary.
           cd ../../..; # switch working directory back to original so next loop starts properly
           done;
 
       If you run this from the folder containing the replicate directories,
       any correctly built |trafo| replicate should be submitted to run.
       You can modify the loop glob to restrict the script to certain directories.
 
 
 
.. |trafo| replace:: :math:`\texttt{TRANSFORMATO}`
 

