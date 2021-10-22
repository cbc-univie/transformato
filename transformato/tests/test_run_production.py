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


@pytest.mark.slowtest
@pytest.mark.rsfe
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
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

    output_files = mutate_toluene_to_methane_cc(configuration=configuration)

    run_simulation(output_files)
    f = "/".join(output_files[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


@pytest.mark.slowtest
@pytest.mark.rsfe
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_methane_to_methane_cc_rsfe_with_openMM():
    from transformato.testsystems import mutate_methane_to_methane_cc
    workdir = get_test_output_dir()

    configuration = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-rsfe.yaml",
        input_dir="data/",
        output_dir=workdir,
    )

    output_files = mutate_methane_to_methane_cc(configuration=configuration)
    run_simulation(output_files)
    f = "/".join(output_files[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


@pytest.mark.slowtest
@pytest.mark.rsfe
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_methane_to_methane_cc_rsfe_with_CHARMM():
    from transformato.testsystems import mutate_methane_to_methane_cc
    workdir = get_test_output_dir()

    configuration = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-rsfe.yaml",
        input_dir="data/",
        output_dir=workdir,
    )

    output_files = mutate_methane_to_methane_cc(
        configuration=configuration
    )
    run_simulation(output_files)
    f = "/".join(output_files[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


@pytest.mark.slowtest
@pytest.mark.rsfe
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_acetylacetone_tautomer_pair_rsfe(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_acetylacetone_tautomer_pair
    workdir = get_test_output_dir()

    configuration = load_config_yaml(
        config="transformato/tests/config/test-acetylacetone-tautomer-rsfe.yaml",
        input_dir="data/",
        output_dir=workdir,
    )

    (output_files_t1, output_files_t2), _, _ = setup_acetylacetone_tautomer_pair(
        configuration=configuration
    )
    run_simulation(output_files_t1)
    run_simulation(output_files_t2)
    f = "/".join(output_files_t1[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_acetylacetone_tautomer_pair_only_in_vacuum(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_acetylacetone_tautomer_pair
    workdir = get_test_output_dir()

    conf = "transformato/tests/config/test-acetylacetone-tautomer-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir=workdir
    )  # NOTE: for preprocessing input_dir is the output dir

    (output_files_t1, output_files_t2), _, _ = setup_acetylacetone_tautomer_pair(
        configuration=configuration, nr_of_bonded_windows=16
    )

    run_simulation(output_files_t1, only_vacuum=True)
    run_simulation(output_files_t2, only_vacuum=True)
    f = "/".join(output_files_t1[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


@pytest.mark.slowtest
@pytest.mark.rsfe
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_2OJ9_tautomer_pair_rsfe_openMM(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_2OJ9_tautomer_pair_rsfe
    workdir = get_test_output_dir()

    conf = "transformato/tests/config/test-2oj9-tautomer-pair-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir=workdir
    )  # NOTE: for preprocessing input_dir is the output dir

    (output_files_t1, output_files_t2), _, _ = setup_2OJ9_tautomer_pair_rsfe(
        configuration=configuration
    )
    run_simulation(output_files_t1)
    run_simulation(output_files_t2)
    f = "/".join(output_files_t1[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


@pytest.mark.slowtest
@pytest.mark.rbfe
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_2OJ9_tautomer_pair_rbfe_openMM(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_2OJ9_tautomer_pair_rbfe
    workdir = get_test_output_dir()

    conf = "transformato/tests/config/test-2oj9-tautomer-pair-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir=workdir
    )  # NOTE: for preprocessing input_dir is the output dir

    (output_files_t1, output_files_t2), _, _ = setup_2OJ9_tautomer_pair_rbfe(
        configuration=configuration
    )
    run_simulation(output_files_t1)
    run_simulation(output_files_t2)
    f = "/".join(output_files_t1[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


@pytest.mark.slowtest
@pytest.mark.rsfe
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_2OJ9_tautomer_pair_charmm(caplog):
    caplog.set_level(logging.WARNING)
    from transformato.constants import change_platform_to 
    from .test_mutation import setup_2OJ9_tautomer_pair_rsfe
    workdir = get_test_output_dir()
    conf = "transformato/tests/config/test-2oj9-tautomer-pair-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir=workdir
    )  # NOTE: for preprocessing input_dir is the output dir
    change_platform_to(configuration, change_to='cpu')

    (output_files_t1, output_files_t2), _, _ = setup_2OJ9_tautomer_pair_rsfe(
        configuration=configuration
    )

    run_simulation(output_files_t1, engine="CHARMM")
    run_simulation(output_files_t2, engine="CHARMM")
    f = "/".join(output_files_t1[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_get_free_energy_acetylaceton_tautomer_pair(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_acetylacetone_tautomer_pair
    workdir = get_test_output_dir()

    conf = "transformato/tests/config/test-2oj9-tautomer-pair-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir=workdir
    )  # NOTE: for preprocessing input_dir is the output dir

    (output_files_t1, output_files_t2), conf, _ = setup_acetylacetone_tautomer_pair(
        configuration=configuration
    )
    run_simulation(output_files_t1, only_vacuum=True)
    run_simulation(output_files_t2, only_vacuum=True)
    f = "/".join(output_files_t1[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)
