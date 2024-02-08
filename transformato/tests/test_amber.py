from transformato import (
    load_config_yaml,
    SystemStructure,
    IntermediateStateFactory,
    ProposeMutationRoute,
)
from transformato.mutate import perform_mutations
from transformato.utils import postprocessing, get_test_output_dir
from transformato_testsystems.testsystems import get_testsystems_dir
import pytest
import os


@pytest.mark.rbfe
@pytest.mark.charmm
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_full_amber_test():
    molecule = "ethane-methanol"

    configuration = load_config_yaml(
        config=f"/site/raid3/johannes/amber_tests/data/config/{molecule}_amber.yaml",
        input_dir="/site/raid3/johannes/amber_tests/data/",
        output_dir="/site/raid3/johannes/amber_tests/",
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.propose_common_core()
    s1_to_s2.finish_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )
    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
        nr_of_mutation_steps_charge=2,
        nr_of_mutation_steps_cc=2,
    )

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol2()
    print(mutation_list.keys())
    i = IntermediateStateFactory(
        system=s2,
        configuration=configuration,
    )

    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
        nr_of_mutation_steps_charge=2,
    )


@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_asfe_amber_test():
    molecule = "ethane-methanol_amber"

    configuration = load_config_yaml(
        config=f"/site/raid3/johannes/amber_tests/data/config/{molecule}.yaml",
        input_dir="/site/raid3/johannes/amber_tests/data/",
        output_dir="/site/raid3/johannes/amber_tests/",
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.propose_common_core()
    s1_to_s2.finish_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
        multiple_runs=3,
    )
    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
        nr_of_mutation_steps_charge=3,
    )
    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol2()
    i = IntermediateStateFactory(
        system=s2,
        configuration=configuration,
        multiple_runs=3,
    )
    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
        nr_of_mutation_steps_charge=3,
    )


@pytest.mark.rbfe
@pytest.mark.charmm
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_amber_analysis():
    molecule = "ethanol"
    configuration = load_config_yaml(
        config=f"/site/raid3/johannes/amber_tests/data/config/{molecule}_asfe.yaml",
        input_dir="/site/raid3/johannes/amber_tests/data/",
        output_dir="/site/raid3/johannes/amber_tests/",
    )

    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration,
        name="ethanol_amber",
        engine="openMM",
        max_snapshots=10000,
        num_proc=6,
        analyze_traj_with="mda",
        show_summary=True,
        multiple_runs=2,
    )

    print(f"Free energy is {ddG_openMM} +- {dddG} kT")


@pytest.mark.charmm
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_amber_rbfe():

    molecule = "jmc_28_ejm_31"
    configuration = load_config_yaml(
        config=f"/site/raid3/johannes/amber_tests/data/config/{molecule}.yaml",
        input_dir="/site/raid3/johannes/amber_tests/data/",
        output_dir="/site/raid3/johannes/amber_tests/",
    )

    # molecule = "1_cano-1_U"
    # configuration = load_config_yaml(
    #     config=f"/site/raid3/johannes/amber_tests/data_GU/config/{molecule}.yaml",
    #     input_dir="/site/raid3/johannes/amber_tests/data_GU/",
    #     output_dir="/site/raid3/johannes/amber_tests/",
    # )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)

    s1_to_s2.propose_common_core()

    # # s1_to_s2.add_idx_to_common_core_of_mol1([31])
    # # s1_to_s2.add_idx_to_common_core_of_mol2([30, 31])

    s1_to_s2.finish_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
        multiple_runs=3,
    )
    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
        nr_of_mutation_steps_charge=3,
    )

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol2()
    i = IntermediateStateFactory(
        system=s2,
        configuration=configuration,
        multiple_runs=3,
    )
    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
        nr_of_mutation_steps_charge=3,
    )
