Installation
==============

Prerequisites
##############

#. A working version of `conda/miniconda <https://docs.conda.io/en/latest/>`_ 
#. A working version of git
#. A CUDA - capable device (locally or node). A CUDA device is not required to generate the intermediate states, but one is required to do simulations or analysises.

Getting transformato
#########################

First, install conda on your machine. Then clone the `repository <https://github.com/wiederm/transformato.git>`_

.. code-block:: bash

    git clone `https://github.com/wiederm/transformato.git
    
Within the newly created ``transformato`` directory, a ``yaml`` file can be found in  ``devtools/conda-envs/``. Use conda to create 
an environment called ``fep`` by going into the transformato directory and using conda:

.. code-block:: bash

    cd transformato
    conda env create --file devtools/conda-envs/fep_env.yaml

Activate the environment with ``conda activate fep``. Now install |trafo| with:

.. code-block:: bash

    python setup.py install

This will install |trafo| in your current (``fep``) conda environment.   

Now, you can open ``python`` try the following command:

::
    
    import transformato
    print(transformato.__version__)

If you output the current transformato version, then congratulations! You are now ready to use transformato.



.. |trafo| replace:: :math:`\texttt{TRANSFORMATO}`