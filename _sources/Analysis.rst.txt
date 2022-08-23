Trajectory Analysis
=====================

In theory, you could run your analysis directly in your :code:`submit.ipynb`, after the simulations are done.
However, we instead recommend (at least for **RBFE** calculations) putting it into
a separate file (usually called :code:`analysis.py``) that can be aliased and run by shell/script commands.

The *analysis.py* file
***********************

Below you can see the code necessary for analysing your trajetories. You can run this code snippet 
in your jupyter-notebook after all simulation runs finished or run it in a terminal in your ``fep`` 
environment with:

.. code-block:: python

    python analysis.py output_dir path conf

.. note::

    To check whether your run have finished properly, go into the *intstX* directories and check the 
    waterbox_out.log and complex_out.log (for RSFE: vacuum_out.log) files

.. attention::
    This can be done only for **RSFE** in a decent amout of time!
    Analysis, just like simulation, takes a decent amount of computational resources, 
    depending on the number of intermediate states.
    Like simulations, it is therefore recommended to offload it to a cluster (especially for **RBFE**).


.. code-block:: python

    from transformato import load_config_yaml
    from transformato.utils import postprocessing
    import sys

    # comes frome the .sh file only needed when running 
    # the analysis.py as standalone script outside your
    # jupyter-notebook
    output_dir = sys.argv[1]
    path = sys.argv[2]
    conf = sys.argv[3]

    configuration = load_config_yaml(config=conf, input_dir=path, output_dir=output_dir)

    for name in (
        configuration["system"]["structure1"]["name"],
        configuration["system"]["structure2"]["name"],
    ):

        ddG_openMM, dddG, f_openMM = postprocessing(
            configuration,
            name=name,
            engine="openMM",
            max_snapshots=10000,
            num_proc=6,
            analyze_traj_with="mda",
            show_summary=True,
        )
        print(f"Free energy difference: {ddG_openMM} +- {dddG} [kT]")


In the ``postprocessing`` function you can set up several things

.. object:: postprocessing
  
    .. option:: name

        if you want to analyze one ligand only enter the ligand folders name here

    .. option:: engine (["openMM","CHARMM"], default: openMM)

        choose your engine

        .. important::

            Here you need to choose the same engine as you used befor for your 
            simulations


    .. option:: max_snapshots ([int], default: 300)

        The amount of snapshots you want to analyze per intermediate state. The more snapshots
        you choose the longer the analysis will take!

    .. option:: analyze_traj_with (["mda","mdtraj"], default: mdtraj)

        Choose which program you want to use when analyzing your trajectory. The main difference is that
        ``mda`` supports multiprocessing and needs less memory on the GPU. 

.. attention::

    You now get the free energy difference from your ligand to the common core substructure, e.g. for ligand 1
    ::math:`\Delta\Delta G^{bind}_{L1\rightarrow |DL_1| - R_{CC}}`. To calculate the final binding free energy 
    ::math:`\Delta\Delta G^{bind}_{L1\rightarrow L2}` you need to substract the free energy difference of ligand2 from 
    ligand 1: ::math:`\Delta\Delta G^{bind}_{L1\rightarrow L2} = \Delta\Delta G^{bind}_{L1\rightarrow |DL_1| - R_{CC}} - \Delta\Delta G^{bind}_{L2\rightarrow\ |DL_2| - R_{CC}}`
    as explained in :doc:`whatistransformato`. The same applies for RSFE calculations. 


Offloading the analysis to a cluster
****************************************

As for the simulations there are two options to start the analysis:

    + You can run it via your ``submit.ipynb``:

    .. code-block:: python

        with open(f'analysis_{folder}.sh', 'w+') as f:
            f.write(f'#!/bin/bash \n')                       # This are the start
            f.write(f'#SBATCH --gres=gpu \n')                # input lines when using SLURM
            f.write(f'#SBATCH -p lgpu \n')                   # when running directly in the notebook
            f.write(f'#SBATCH -d afterany:"{jobid}"  \n')    # make sure that all simulations have finished
            f.write(f' \n')
            f.write(f' source ~/miniconda3/etc/profile.d/conda.sh \n')  # necessary to use the conda environment
            f.write(f' conda activate fep \n')
            f.write(f' \n')        
            f.write(f'output_dir={folder} \n')
            f.write(f' \n')
            f.write(f'cd {wd} \n')
            f.write(f' \n')
            f.write(f'time python analysis.py {output_dir} {input_dir} {config} > {folder}/analysis.out \n')

        exe = subprocess.Popen(["ssh", "your-clustername", "sbatch", f"{wd}/analysis_{folder}.sh"], text=True, stdout=subprocess.PIPE )
        output = exe.stdout.read()
        print(output)

    The variables ``input_dir`` and ``config`` should be known from the beginning of your ``submit.ipynb`` file
    (see also :doc:`Running_Simulation`).


    + You can offload it manually

    To further simplify matters, create an executable :code:`analysis.sh` file containing the following code:

    .. code-block:: bash

        #!/bin/bash 
        #SBATCH --gres=gpu 
        #SBATCH -p lgpu            
        
        folder=$1 
        input_dir=$2
        config=$3
        
        python analysis.py ./${folder}/ ${input_dir}  ${config} > ./analysis_${folder}.out 

    You can run this script in your terminal

    .. code-block:: bash

        sbatch analysis.sh ./your-output-folder ./your-input-folder ./your-yaml-file


In both cases you will get your results in a file called :code:`analysis.out`, in addition a ``results`` 
folder is created.


.. |trafo| replace:: :math:`\texttt{TRANSFORMATO}`