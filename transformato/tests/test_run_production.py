"""
Unit and regression test for the transformato package.
"""

import os
import logging

import pytest
import shutil
from transformato import load_config_yaml
from transformato.tests.paths import get_test_output_dir
from transformato.utils import run_simulation
from transformato.constants import change_platform_to_test_platform


@pytest.mark.rsfe
@pytest.mark.skipif(
    os.getenv("CI") == True,
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_1a0q_1a07_rsfe_with_openMM(caplog):
    import logging

    # Test that TF can handel multiple dummy regions
    caplog.set_level(logging.DEBUG)
    import warnings
    from transformato import (
        load_config_yaml,
        SystemStructure,
        IntermediateStateFactory,
        ProposeMutationRoute,
    )
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


@pytest.mark.rsfe
@pytest.mark.skipif(
    os.getenv("CI") == True,
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_toluene_to_methane_cc_rsfe_with_openMM():
    from transformato.testsystems import mutate_toluene_to_methane_cc
    from .test_run_production import run_simulation

    workdir = get_test_output_dir()

    configuration = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-rsfe.yaml",
        input_dir="data/",
        output_dir=workdir,
    )
    change_platform_to_test_platform(configuration=configuration, engine="openMM")
    output_files = mutate_toluene_to_methane_cc(configuration=configuration)

    run_simulation(output_files)


@pytest.mark.skipif(
    os.getenv("CI") == True,
    reason="Skipping tests that cannot pass in github actions",
)
@pytest.mark.rsfe
def test_run_methane_to_methane_cc_with_openMM():
    from transformato.testsystems import mutate_methane_to_methane_cc

    workdir = get_test_output_dir()

    configuration = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-rsfe.yaml",
        input_dir="data/",
        output_dir=workdir,
    )
    change_platform_to_test_platform(configuration=configuration, engine="openMM")
    output_files = mutate_methane_to_methane_cc(configuration=configuration)
    run_simulation(output_files)


@pytest.mark.rsfe
@pytest.mark.skipif(
    os.getenv("CI") == True,
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_methane_to_methane_common_core_with_CHARMM():
    from transformato.testsystems import mutate_methane_to_methane_cc

    workdir = get_test_output_dir()

    configuration = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-rsfe.yaml",
        input_dir="data/",
        output_dir=workdir,
    )
    change_platform_to_test_platform(configuration=configuration, engine="CHARMM")
    output_files = mutate_methane_to_methane_cc(configuration=configuration)
    run_simulation(output_files)
    f = "/".join(output_files[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


@pytest.mark.rsfe
@pytest.mark.skipif(
    os.getenv("CI") == True,
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_acetylacetone_tautomer_pair(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_acetylacetone_tautomer_pair

    workdir = get_test_output_dir()

    configuration = load_config_yaml(
        config="transformato/tests/config/test-acetylacetone-tautomer-rsfe.yaml",
        input_dir="data/",
        output_dir=workdir,
    )
    change_platform_to_test_platform(configuration=configuration, engine="openMM")
    (output_files_t1, output_files_t2), _, _ = setup_acetylacetone_tautomer_pair(
        configuration=configuration
    )
    run_simulation(output_files_t1)
    run_simulation(output_files_t2)


@pytest.mark.skipif(
    os.getenv("CI") == True,
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_acetylacetone_tautomer_pair_only_in_vacuum(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_acetylacetone_tautomer_pair

    workdir = get_test_output_dir()

    conf = "transformato/tests/config/test-acetylacetone-tautomer-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir=workdir
    )  # NOTE: for preprocessing input_dir is the output dir
    change_platform_to_test_platform(configuration=configuration, engine="openMM")

    (output_files_t1, output_files_t2), _, _ = setup_acetylacetone_tautomer_pair(
        configuration=configuration, nr_of_bonded_windows=16
    )

    run_simulation(output_files_t1, only_vacuum=True)
    run_simulation(output_files_t2, only_vacuum=True)


@pytest.mark.rsfe
@pytest.mark.skipif(
    os.getenv("CI") == True,
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_2OJ9_tautomer_pair_with_openMM(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_2OJ9_tautomer_pair_rsfe

    workdir = get_test_output_dir()

    conf = "transformato/tests/config/test-2oj9-tautomer-pair-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir=workdir
    )  # NOTE: for preprocessing input_dir is the output dir
    change_platform_to_test_platform(configuration=configuration, engine="openMM")

    (output_files_t1, output_files_t2), _, _ = setup_2OJ9_tautomer_pair_rsfe(
        configuration=configuration
    )
    run_simulation(output_files_t1)
    run_simulation(output_files_t2)


@pytest.mark.rbfe
@pytest.mark.skipif(
    os.getenv("CI") == True,
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_2OJ9_tautomer_pair_with_openMM(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_2OJ9_tautomer_pair_rbfe

    workdir = get_test_output_dir()

    conf = "transformato/tests/config/test-2oj9-tautomer-pair-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir=workdir
    )  # NOTE: for preprocessing input_dir is the output dir
    change_platform_to_test_platform(configuration=configuration, engine="openMM")

    (output_files_t1, output_files_t2), _, _ = setup_2OJ9_tautomer_pair_rbfe(
        configuration=configuration
    )
    run_simulation(output_files_t1)
    run_simulation(output_files_t2)


@pytest.mark.rsfe
@pytest.mark.skipif(
    os.getenv("CI") == True,
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_2OJ9_tautomer_pair_with_CHARMM(caplog):
    caplog.set_level(logging.WARNING)
    from transformato.constants import change_platform_to
    from .test_mutation import setup_2OJ9_tautomer_pair_rsfe

    workdir = get_test_output_dir()
    conf = "transformato/tests/config/test-2oj9-tautomer-pair-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir=workdir
    )  # NOTE: for preprocessing input_dir is the output dir
    change_platform_to_test_platform(configuration=configuration, engine="CHARMM")

    (output_files_t1, output_files_t2), _, _ = setup_2OJ9_tautomer_pair_rsfe(
        configuration=configuration
    )

    run_simulation(output_files_t1, engine="CHARMM")
    run_simulation(output_files_t2, engine="CHARMM")
