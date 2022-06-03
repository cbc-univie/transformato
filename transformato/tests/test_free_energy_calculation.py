import logging
import os

import numpy as np
import pytest
from transformato import load_config_yaml
from transformato.constants import change_platform_to_test_platform
from transformato.tests.paths import get_test_output_dir
from transformato.utils import postprocessing, run_simulation


@pytest.mark.rsfe
@pytest.mark.full_workflow
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_1a0q_1a07_rsfe_with_openMM(caplog):

    # Test that TF can handel multiple dummy regions
    caplog.set_level(logging.DEBUG)
    import warnings

    from transformato import (IntermediateStateFactory, ProposeMutationRoute,
                              SystemStructure, load_config_yaml)
    from transformato.mutate import perform_mutations

    warnings.filterwarnings("ignore", module="parmed")

    workdir = get_test_output_dir()
    conf = "transformato/tests/config/test-1a0q-1a07-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/test_systems_mutation", output_dir=workdir
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

    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration,
        name="1a0q",
        engine="openMM",
        max_snapshots=10_000,
    )
    print(ddG_openMM, dddG)


@pytest.mark.rsfe
@pytest.mark.full_workflow
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_2oj9_rsfe_with_different_switches(caplog):
    caplog.set_level(logging.WARNING)
    from ..analysis import FreeEnergyCalculator
    from .test_mutation import setup_2OJ9_tautomer_pair_rsfe
    from .test_run_production import run_simulation

    # vfswitch
    conf_path = "transformato/tests/config/test-2oj9-tautomer-pair-rsfe_vfswitch.yaml"
    configuration = load_config_yaml(
        config=conf_path, input_dir="data/", output_dir=get_test_output_dir()
    )  # NOTE: for preprocessing input_dir is the output dir

    # generate samples
    (output_files_t1, output_files_t2), _, _ = setup_2OJ9_tautomer_pair_rsfe(
        configuration=configuration
    )
    run_simulation(output_files_t1)

    f = FreeEnergyCalculator(configuration, "2OJ9-original")

    # 2OJ9-original to tautomer common core
    ddG_charmm, dddG, f_charmm_vfswitch = postprocessing(
        configuration,
        name="2OJ9-original",
        engine="CHARMM",
        max_snapshots=600,
        different_path_for_dcd="data/2OJ9-original-2OJ9-tautomer-rsfe/2OJ9-original",
    )
    ddG_openMM, dddG, f_openMM_vfswitch = postprocessing(
        configuration,
        name="2OJ9-original",
        engine="openMM",
        max_snapshots=600,
        different_path_for_dcd="data/2OJ9-original-2OJ9-tautomer-rsfe/2OJ9-original",
    )

    assert np.isclose(
        f_charmm_vfswitch.vacuum_free_energy_differences[0, -1],
        f_openMM_vfswitch.vacuum_free_energy_differences[0, -1],
        atol=1e-2,
    )

    assert np.isclose(
        f_charmm_vfswitch.waterbox_free_energy_differences[0, -1],
        f_openMM_vfswitch.waterbox_free_energy_differences[0, -1],
        atol=1e-2,
    )

    # switch
    conf_path = "transformato/tests/config/test-2oj9-tautomer-pair-rsfe_vswitch.yaml"
    configuration = load_config_yaml(
        config=conf_path, input_dir="data/", output_dir=get_test_output_dir()
    )  # NOTE: for preprocessing input_dir is the output dir

    run_simulation(output_files_t1)
    configuration = load_config_yaml(
        config=conf_path, input_dir="data/", output_dir=get_test_output_dir()
    )  # NOTE: for preprocessing input_dir is the output dir
    f = FreeEnergyCalculator(configuration, "2OJ9-original")

    # 2OJ9-original to tautomer common core
    ddG_charmm, dddG, f_charmm_switch = postprocessing(
        configuration,
        name="2OJ9-original",
        engine="CHARMM",
        max_snapshots=600,
        different_path_for_dcd="data/2OJ9-original-2OJ9-tautomer-rsfe/2OJ9-original",
    )
    ddG_charmm, dddG, f_openMM_switch = postprocessing(
        configuration,
        name="2OJ9-original",
        engine="openMM",
        max_snapshots=600,
        different_path_for_dcd="data/2OJ9-original-2OJ9-tautomer-rsfe/2OJ9-original",
    )

    assert np.isclose(
        f_charmm_switch.vacuum_free_energy_differences[0, -1],
        f_openMM_switch.vacuum_free_energy_differences[0, -1],
        atol=1e-2,
    )

    assert np.isclose(
        f_charmm_switch.waterbox_free_energy_differences[0, -1],
        f_openMM_switch.waterbox_free_energy_differences[0, -1],
        atol=1e-2,
    )

    assert np.isclose(
        f_charmm_vfswitch.vacuum_free_energy_differences[0, -1],
        f_charmm_switch.vacuum_free_energy_differences[0, -1],
        atol=1e-2,
    )
    assert np.isclose(
        f_charmm_vfswitch.waterbox_free_energy_differences[0, -1],
        f_charmm_switch.waterbox_free_energy_differences[0, -1],
        atol=1e-2,
    )
    assert np.isclose(
        f_openMM_vfswitch.waterbox_free_energy_differences[0, -1],
        f_openMM_switch.waterbox_free_energy_differences[0, -1],
        atol=1e-2,
    )


@pytest.mark.rbfe
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_28_1h1q_rbfe_with_openMM():
    # Generating output for a run of the CDK2 Ligand System
    from transformato import (IntermediateStateFactory, ProposeMutationRoute,
                              SystemStructure, load_config_yaml)
    from transformato.mutate import perform_mutations
    from transformato.utils import postprocessing, run_simulation

    configuration = load_config_yaml(
        config="transformato/tests/config/test-28_1h1q_rbfe.yaml",
        input_dir="data/",
        output_dir='/site/raid4/johannes/test',
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
