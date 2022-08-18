
from transformato import (
    load_config_yaml,
    SystemStructure,
    IntermediateStateFactory,
)

from transformato.utils import run_simulation, postprocessing
from transformato.annihilation import ProposeMutationRouteASFE
from transformato.mutate import perform_mutations
from transformato.tests.paths import get_test_output_dir
from transformato_testsystems.testsystems import get_testsystems_dir

import numpy as np
import warnings

warnings.filterwarnings("ignore", module="parmed")


def create_asfe_system(configuration):

    s1 = SystemStructure(configuration, "structure1")
    s1_absolute = ProposeMutationRouteASFE(s1)
    s1_absolute.propose_common_core()
    s1_absolute.finish_common_core()
    mutation_list = s1_absolute.generate_mutations_to_common_core_for_mol1()

    return s1, mutation_list


def test_run_asfe_system():

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

    run_simulation(i.output_files, engine="openMM")

def analyse_asfe_with_mda():

    configuration = load_config_yaml(
        config=f"{get_testsystems_dir()}/config/methanol-asfe.yaml",
        input_dir=get_testsystems_dir(),
        output_dir=get_test_output_dir(),
    )

    final_dg = []
    # runs = 1
    # for run in range(1,runs + 1):
    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration,
        name="methanol",
        engine="openMM",
        max_snapshots=50,
        num_proc=2,
        analyze_traj_with="mda",
        show_summary=True,
    )
    print(f"Free energy difference: {ddG_openMM} +- {dddG} [kT]")
    final_dg.append(ddG_openMM)
    # print(f"Final free energy is {round(np.average(final_dg),2)} +- {round(np.std(final_dg))} of the {runs} individual runs {final_dg}")

def analyse_asfe_with_mdtraj():

    configuration = load_config_yaml(
        config=f"{get_testsystems_dir()}/config/methanol-asfe.yaml",
        input_dir=get_testsystems_dir(),
        output_dir=get_test_output_dir(),
    )

    final_dg = []
    runs = 3
    for run in range(1,runs + 1):
        ddG_openMM, dddG, f_openMM = postprocessing(
            configuration,
            name="methanol",
            engine="openMM",
            max_snapshots=100,
            analyze_traj_with="mdtraj",
            show_summary=True,
            consecutive_runs=run,
        )
        print(f"Free energy difference: {ddG_openMM} +- {dddG} [kT]")
        final_dg.append(ddG_openMM)
    print(f"Final free energy is {round(np.average(final_dg),2)} +- {round(np.std(final_dg))} of the {runs} individual runs {final_dg}")