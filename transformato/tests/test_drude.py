from transformato.utils import (
    run_simulation,
    postprocessing,
    get_test_output_dir,
    load_config_yaml,
)
from transformato.mutate import ProposeMutationRoute, perform_mutations
from transformato.system import SystemStructure
from transformato.state import IntermediateStateFactory
import pytest
import os


@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_create_asfe_system():
    configuration = load_config_yaml(
        config=f"/site/raid5/marvin/drude_test/data/config/acetamide.yaml",
        input_dir=f"/site/raid5/marvin/drude_test/data/",
        output_dir=f"/site/raid5/marvin/drude_test",
    )

    s1 = SystemStructure(configuration, "structure1")
    s1_absolute = ProposeMutationRoute(s1)
    s1_absolute.propose_common_core()
    s1_absolute.finish_common_core()
    mutation_list = s1_absolute.generate_mutations_to_common_core_for_mol1()

    multiple_runs = 3
    i = IntermediateStateFactory(
        system=s1, multiple_runs=multiple_runs, configuration=configuration
    )

    perform_mutations(
        configuration=configuration,
        nr_of_mutation_steps_charge=5,
        i=i,
        mutation_list=mutation_list,
    )


@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_analyze_drude():

    configuration = load_config_yaml(
        config=f"/site/raid5/marvin/drude_test/data/config/methanol.yaml",
        input_dir=f"/site/raid5/marvin/drude_test/data/",
        output_dir=f"/site/raid5/marvin/drude_test",
    )

    ddG_openMM, dddG, _ = postprocessing(
        configuration,
        name="methanol",
        engine="openMM",
        max_snapshots=50,
        num_proc=6,
        analyze_traj_with="mda",
        show_summary=True,
        multiple_runs=1,
    )

    print(f"Free energy is {ddG_openMM} +- {dddG} kT")
