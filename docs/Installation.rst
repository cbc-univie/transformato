Installation
==============
 
Prerequisites
##############
 
#. A working version of `conda/miniconda <https://docs.conda.io/en/latest/>`_
#. A working version of git
#. A CUDA - capable device (locally or node). A CUDA device is not required to generate the intermediate states, but one is required to do simulations or analyses.
 
Getting |trafo|
#########################
 
First, install conda on your machine. Then clone the `repository <https://github.com/wiederm/transformato.git>`_
 
.. code-block:: bash
 
   git clone https://github.com/wiederm/transformato.git
  
Within the newly created ``transformato`` directory, a ``yaml`` file can be found in  ``devtools/conda-envs/``. Use conda to create
an environment called ``fep`` by going into the transformato directory and using conda:
 
.. code-block:: bash
 
   cd transformato
   conda env create --file devtools/conda-envs/fep_env.yaml
 
.. caution::
   If you are doing this the first time, this may take a while (up to 10 minutes!)
 
Activate the environment with ``conda activate fep``.
 
.. note::
   For running simulations on an **NVIDIA** machine using **CUDA**, openMM will install the python
   package cudatoolkit. Make sure the version in your environment is not higher than the
   CUDA version of your NVIDIA driver (you can check by using ``nvidia-smi`` for the driver
   and by ``conda list`` you can see your current version). If necessary, you can install
   the correct cudatoolkit version by using ``conda install cudatoolkit=X.Y -c conda-forge``
 
Now install |trafo| with:
 
.. code-block:: bash
 
   python setup.py install
 
This will install |trafo| in your current (``fep``) conda environment.  
 
Now, you can open ``python`` try the following command:
 
::
  
   import transformato
   print(transformato.__version__)
 
If you output the current transformato version, then congratulations! You are now ready to use transformato.
 
 
 
.. |trafo| replace:: :math:`\texttt{TRANSFORMATO}`
 

