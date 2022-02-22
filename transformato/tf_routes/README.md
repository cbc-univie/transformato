tf_routes
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/jalhackl/tf_routes/workflows/CI/badge.svg)](https://github.com/jalhackl/tf_routes/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jalhackl/tf_routes/branch/master/graph/badge.svg)](https://codecov.io/gh/jalhackl/tf_routes/branch/master)


# Mutation Routes for Transformato

In order to use the route algorithms in Transformato,

replace

 **_mol_to_nx(mol: Chem.Mol)** from transformato/system.py
 
 with
 
**tf_routes.preprocessing._mol_to_nx_full_weight(mol)**

and 

**_calculate_order_of_LJ_mutations** from mutate.py

by one of the new mutation algorithms described below.

# tf_routes.preprocessing 
to a large extent slightly modified functions from transformato solely for testing purposes, 

except for _mol_to_nx_full_weight(mol), which is needed for setting up the weighted graph for the new mutation algorithms

# tf_routes.visualization

animated_visualization_3d_v1(mol, mutationl, ccoremol, hits): animated visualization using py3dmol

animated_visualization_3d_v2(mol, mutationl, ccoremol, hits): animated visualization using py3dmol, additionally shows the common core

_show_common_core_gradient( mol, highlight, mutationl, percomponent = False, numbers = False): visualization with colors and numbers

_show_common_core_gradient_write( mol, highlight, mutationl, percomponent = False): visualization with colors and numbers, additionally to the svg returns a png (can be used to subsequently write visualizations in files)

# tf_routes.routes 

## 4 different functions for calculating mutations:
_calculate_order_of_LJ_mutations: naive dfs (as currently in transformato)


_calculate_order_of_LJ_mutations_new: bfs/djikstra-algorithm applied once for route


_calculate_order_of_LJ_mutations_new_iter: bfs/djikstra-algorithm applied iteratively, i.e. after each removal of an atom 


_calculate_order_of_LJ_mutations_new_iter_change: works iteratively, after each removal of an atom the algorithm is chosen depending on current state





### Copyright

Copyright (c) 2022, jalh


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.6.
