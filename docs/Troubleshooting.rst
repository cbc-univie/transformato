Troubleshooting
================



Getting debug information
##############################

For the simulation:

+ in the output file of your worload manager (e.g. `slurm-812311.out`)

+ for the simulation: in `/intstate*/complex_out.log` or `/intstate*/waterbox_out.log` respectively

For the analysis:

+ in the output file you specify when running `analysis.py`, typically some form of `analysis_*.out` 

Error messages and common causes:
####################################


Errors during system setup:
******************************

.. code-block:: python

    load_config_yaml(config, input_dir, output_dir)
         except yaml.YAMLError as exc:
             print(exc)
    if (     settingsMap["simulation"]["parameters"].get("nstep") == None
         or settingsMap["simulation"]["parameters"].get("nstdcd") == None
     ):
         raise KeyError("nsteps or nstdcd is not defined in config file")
    else:

    UnboundLocalError: local variable 'settingsMap' referenced before assignment

This usually means that there is a typo in your yaml file. Remember that yaml requires a space between colon and data!

.. code-block:: python

    RuntimeError: Charge compensation failed. Introducing non integer total charge: -0.584000000000001. 
    Target total charge: 0.

This can have a number of root causes. Ensure that:

#. your dummy regions are sensible
#. your terminal X is where it should be
#. your structures do not require transformato to create new dummy atoms. Put the structure with more atoms first.


Errors during simulation:
**************************

.. code-block:: python

    FileNotFoundError: [Errno 2] No such file or directory: '/site/raid4/student4/production_transformato/tablit_24to25/run1/2oj9_tablit_struc24-2oj9_tablit_struc25-rbfe/2oj9_tablit_struc24//intst1/lig_in_complex.rst'

Occurs if no .rst restart file was outputted during equilibration, either due to the equilibration not running correctly (or at all) or due to no output file specified.


Errors during analysis:
*************************

.. code-block:: python

    Traceback (most recent call last):
    File "/site/raid4/student4/alex/production_transformato/taapdb_24to25/analysis.py", line 23, in <module>
        ddG_openMM, dddG, f_openMM = postprocessing(
    File "/home/student4/anaconda3/envs/fep/lib/python3.9/site-packages/transformato-0.1+129.gee57388-py3.9.egg/transformato/utils.py", line 60, in postprocessing
        f.calculate_dG_to_common_core(
    File "/home/student4/anaconda3/envs/fep/lib/python3.9/site-packages/transformato-0.1+129.gee57388-py3.9.egg/transformato/analysis.py", line 636, in calculate_dG_to_common_core
        self.mbar_results[env] = self._analyse_results_using_mda(
    File "/home/student4/anaconda3/envs/fep/lib/python3.9/site-packages/transformato-0.1+129.gee57388-py3.9.egg/transformato/analysis.py", line 507, in _analyse_results_using_mda
        self.nr_of_states = len(next(os.walk(f"{self.base_path}"))[1])
    StopIteration


This means that it could not find the path to the directory containing the intermediate states. Most commonly, this is caused by a typo when specifying the path.


