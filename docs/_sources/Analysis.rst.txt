Trajectory Analysis
=====================

In theory, you could run your analysis directly in your :code:`submit.ipynb`, after the simulations are done. However, we instead recommend putting it into
a separate file (usually called :code:`analysis.py``) that can be aliased and run by shell/script commands.

.. note::
    Analysis, just like simulation, takes a decent amount of time and computational resources, depending on the number of intermediate states.

    Like simulations, it is therefore recommended to offload it to a cluster.

**analysis.py**

.. code-block:: python

    import subprocess
    from transformato import load_config_yaml
    from transformato.utils import postprocessing
    import sys

    # comes frome the .sh file
    folder = sys.argv[1]
    path = sys.argv[2]
    conf = sys.argv[3]

    # proc = sys.argv[4]
    # frames = sys.argv[5]

    # print(f"Using {proc} CPUS and frames {frames}")

    configuration = load_config_yaml(config=conf, input_dir=path, output_dir=folder)

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

If you have aliased this script as e.g. :code:`transformato-analyse` , you may now simply go to the folder containing your replicates and run:

.. code-block:: bash

    transformato-analyse ./replicate-folder/ ./your-config.yaml

You will get your results in a file called :code:`analysis_replicate-folder.out`