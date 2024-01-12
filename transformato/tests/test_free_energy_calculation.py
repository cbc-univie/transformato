import logging
import os
import warnings
import numpy as np
import pytest


from transformato import (
    IntermediateStateFactory,
    ProposeMutationRoute,
    SystemStructure,
    load_config_yaml,
)
from transformato.utils import postprocessing, run_simulation, get_test_output_dir
from transformato.mutate import perform_mutations
from transformato_testsystems.testsystems import get_testsystems_dir

warnings.filterwarnings("ignore", module="parmed")


@pytest.mark.rsfe
@pytest.mark.full_workflow
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_1a0q_1a07_rsfe_with_openMM(caplog):
    # Test that TF can handel multiple dummy regions
    caplog.set_level(logging.DEBUG)
    conf = f"{get_testsystems_dir()}/config/test-1a0q-1a07-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir=get_testsystems_dir(), output_dir=get_test_output_dir()
    )
    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.propose_common_core()
    s1_to_s2.finish_common_core()
    # generate the mutation list for the original
    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    print(mutation_list.keys())
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )
    perform_mutations(configuration=configuration, i=i, mutation_list=mutation_list)
    run_simulation(i.output_files, engine="openMM")

    # Analysis with mda
    ddG_openMM_mda, dddG_mda, f_openMM_mda = postprocessing(
        configuration,
        name="1a0q",
        engine="openMM",
        analyze_traj_with="mda",
        num_proc=4,
        max_snapshots=50,
    )
    print(f"Free energy difference using mda: {ddG_openMM_mda} +- {dddG_mda} [kT")

    # Analysis with mdtraj
    ddG_openMM_mtraj, dddG_mtraj, f_openMM_mtraj = postprocessing(
        configuration,
        name="1a0q",
        engine="openMM",
        max_snapshots=50,
    )
    print(
        f"Free energy difference using mdtraj: {ddG_openMM_mtraj} +- {dddG_mtraj} [kT"
    )

    assert np.isclose(ddG_openMM_mda, ddG_openMM_mda, rtol=0.2)
    assert np.isclose(dddG_mda, dddG_mda, rtol=0.2)


@pytest.mark.rbfe
@pytest.mark.full_workflow
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_28_1h1q_rbfe_with_openMM():
    # Generating output for a run of the CDK2 Ligand System

    configuration = load_config_yaml(
        config=f"{get_testsystems_dir()}/config/test-28_1h1q_rbfe.yaml",
        input_dir=get_testsystems_dir(),
        output_dir=get_test_output_dir(),
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)

    s1_to_s2.propose_common_core()
    s1_to_s2.finish_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    print(mutation_list.keys())
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

    run_simulation(i.output_files, engine="charmm")
    from transformato.utils import postprocessing

    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration,
        name="cdk2-28",
        engine="openMM",
        max_snapshots=50,
        num_proc=4,
        analyze_traj_with="mda",
    )
    print(f"Free energy difference: {ddG_openMM} +- {dddG} [kT")
