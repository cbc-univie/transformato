from transformato import (
    load_config_yaml,
    SystemStructure,
    IntermediateStateFactory,
)

from transformato.annihilation import ProposeMutationRouteASFE
from transformato.mutate import *
from transformato.utils import *
from transformato.tests.paths import get_test_output_dir

import warnings

warnings.filterwarnings("ignore", module="parmed")


def test_create_asfe_system():

    configuration = load_config_yaml(
        config="/site/raid2/johannes/transformato/transformato/tests/config/test-methanol-asfe.yaml",
        input_dir="data/",
        output_dir="/site/raid4/johannes/test-RSFE/run_4/",
    )

    s1 = SystemStructure(configuration, "structure1")
    s1_absolute = ProposeMutationRouteASFE(s1)
    s1_absolute.propose_common_core()
    s1_absolute.finish_common_core()
    mutation_list = s1_absolute.generate_mutations_to_common_core_for_mol1()

    i = IntermediateStateFactory(system=s1, configuration=configuration)

    perform_mutations(
        configuration=configuration,
        nr_of_mutation_steps_charge=2,
        i=i,
        mutation_list=mutation_list,
    )

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
