#!/usr/bin/env python
# coding: utf-8

#imports
import argparse
from transformato import load_config_yaml, IntermediateStateFactory, SystemStructure
from transformato.mutate import perform_mutations, ProposeMutationRoute
import glob
import subprocess
import os
import warnings

parser=argparse.ArgumentParser()
parser.add_argument("-mol", metavar="Molecule", dest="mol")
args = parser.parse_args()

molecule = f"{args.mol}"


#input directories - beware that transformato has some hardcoded (sub)directory names!
input_dir = "./data/"
config = f"./data/config/{args.mol}.yaml"


#generate configuration dict
configuration = load_config_yaml(config=config, input_dir=input_dir, output_dir=".")

#SAI mutation route
s1 = SystemStructure(configuration, "structure1")
s1_absolute = ProposeMutationRoute(s1)
s1_absolute.propose_common_core()
s1_absolute.finish_common_core()


mutation_list = s1_absolute.generate_mutations_to_common_core_for_mol1()

#generate intermediate states
i = IntermediateStateFactory(
    system=s1,
    configuration=configuration,
    multiple_runs=5,
)

#current implementation scales charge linearly, with step size 1/nr_of_mutation_steps_charge
perform_mutations(
    configuration=configuration,
    nr_of_mutation_steps_lj_of_hydrogens=1,
    nr_of_mutation_steps_charge=4,
    i=i,
    mutation_list=mutation_list,
)