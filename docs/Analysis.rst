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

        if you only want to analyze one ligand enter the ligand folders name here

    .. option:: engine ({"openMM","CHARMM")

        choose your engine

        




Offloading the analysis to a cluster
****************************************

As for the simulations there are two options to start the analysis:

    + You can run it via your ``submit.ipynb``:


    + You can offload it manually

        To further simplify matters, create an executable :code:`analysis.sh` file containing the following code:

        .. code-block:: bash

            #!/bin/bash 
            #SBATCH --gres=gpu 
            #SBATCH -p lgpu 
            #SBATCH -d afterany:179572  
            #SBATCH --exclude="n00[01-10]" 
            
            
            folder=$1 
            config=$2

            
            time python analysis.py ./${folder}/ ./../  ${config} > ./analysis_${folder}.out 

        The :code:`SBATCH` - lines are cluster instructions and will differ depending on your needs and workload - manager.

        If you have aliased this script as e.g. :code:`transformato-analyse` , you may now simply go to 
        the folder containing your replicates and run:

        .. code-block:: bash

            python transformato-analyse ./replicate-folder/ ./your-config.yaml

        You will get your results in a file called :code:`analysis_replicate-folder.out`



.. |trafo| replace:: :math:`\texttt{TRANSFORMATO}`