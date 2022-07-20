from venv import create
from transformato import (
    load_config_yaml,
    SystemStructure,
    IntermediateStateFactory,
)

import logging
import pytest
import sys, os

from transformato.annihilation import ProposeMutationRouteASFE
from transformato.mutate import perform_mutations
from transformato.tests.paths import get_test_output_dir
from transformato_testsystems.testsystems import get_testsystems_dir

import warnings

warnings.filterwarnings("ignore", module="parmed")


def create_asfe_system(configuration):

    s1 = SystemStructure(configuration, "structure1")
    s1_absolute = ProposeMutationRouteASFE(s1)
    s1_absolute.propose_common_core()
    s1_absolute.finish_common_core()
    mutation_list = s1_absolute.generate_mutations_to_common_core_for_mol1()

    return s1, mutation_list


def test_create_asfe_system():

    configuration = load_config_yaml(
        config=f"{get_testsystems_dir()}/config/methanol-asfe.yaml",
        input_dir=get_testsystems_dir(),
        output_dir=get_test_output_dir(),
    )

    s1, mutation_list = create_asfe_system(configuration)

    i = IntermediateStateFactory(system=s1, configuration=configuration)

    perform_mutations(
        configuration=configuration,
        nr_of_mutation_steps_charge=2,
        i=i,
        mutation_list=mutation_list,
    )

    assert len(i.output_files) == 7
    assert len((mutation_list)["charge"][0].atoms_to_be_mutated) == 6

    # ddG_openMM, dddG, f_openMM = postprocessing(
    #     configuration,
    #     name="methanol",
    #     engine="openMM",
    #     max_snapshots=10000,
    #     num_proc=5,
    #     analyze_traj_with="mda",
    #     show_summary=True,
    # )
    # print(f"Free energy difference: {ddG_openMM} +- {dddG} [kT]")
