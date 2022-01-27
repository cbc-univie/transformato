"""
Unit and regression test for the transformato package.
"""

# Import package, test suite, and other packages as needed
import logging
import os
import warnings
from io import StringIO

import parmed as pm
import pytest
# read in specific topology with parameters
# read in specific topology with parameters
from transformato import SystemStructure, load_config_yaml, psf_correction
from transformato.constants import loeffler_testsystems_dir
from transformato.tests.paths import get_test_output_dir

warnings.filterwarnings("ignore", module="parmed")

def test_read_yaml():
    """Sample test, will check ability to read yaml files"""
    settingsMap = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-rsfe.yaml",
        input_dir=".",
        output_dir="data/",
    )

    assert settingsMap["system"]["name"] == "toluene-methane-rsfe"
    assert settingsMap["system"]["structure1"]["tlc"] == "UNL"


def test_io_psf_files():
    from simtk.openmm.app import CharmmPsfFile
    from transformato.testsystems import mutate_toluene_to_methane_cc

    from .test_run_production import run_simulation

    configuration = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-rsfe.yaml",
        input_dir="data/",
        output_dir=get_test_output_dir(),
    )

    output_files = mutate_toluene_to_methane_cc(configuration=configuration)
    output_path = output_files[0]
    print(output_path)
    CharmmPsfFile(f"{output_path}/lig_in_waterbox.psf")
    CharmmPsfFile(f"{output_path}/lig_in_waterbox_corr.psf")


def test_psf_files():
    test_psf = pm.charmm.psf.CharmmPsfFile("transformato/tests/config/test_input.psf")
    output = StringIO()
    test_psf.write_psf(output)
    corrected_psf = psf_correction(output)
    correction_on = False
    for line in corrected_psf.split("\n"):  # split on newline charactar
        if "!NATOM" in line:  # if !NATOM is found start correction mode
            correction_on = True
            continue

        if "!NBOND" in line:  # if !NBOND is found exit correction mode
            correction_on = False

        if (
            correction_on == True
        ):  # if in correction mode take the string, split on whitespace and put the values in a newly formated string
            if len(line) == 0:
                pass
            else:
                assert len(line) == 118
                values = line.split()
                assert len(values) == 11


def test_initialize_systems(caplog):
    caplog.set_level(logging.DEBUG)

    configuration = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-rsfe.yaml",
        input_dir="data/",
        output_dir=get_test_output_dir(),
    )

    s1 = SystemStructure(configuration, "structure1")
    assert int(s1.offset["waterbox"]) == 0
    assert int(s1.offset["vacuum"]) == 0

    s2 = SystemStructure(configuration, "structure2")
    assert int(s2.offset["waterbox"]) == 0
    assert int(s2.offset["vacuum"]) == 0

    assert "vacuum" in s1.envs and "vacuum" in s2.envs
    assert "waterbox" in s1.envs and "waterbox" in s2.envs

    configuration = load_config_yaml(
        config="transformato/tests/config/test-2oj9-tautomer-pair-rsfe.yaml",
        input_dir="data/",
        output_dir=get_test_output_dir(),
    )

    s1 = SystemStructure(configuration, "structure1")
    assert int(s1.offset["waterbox"]) == 0
    assert int(s1.offset["vacuum"]) == 0

    s2 = SystemStructure(configuration, "structure2")
    assert int(s2.offset["waterbox"]) == 0
    assert int(s2.offset["vacuum"]) == 0


@pytest.mark.rsfe
@pytest.mark.requires_loeffler_systems
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_setup_system_for_methane_common_core():
    from transformato.testsystems import mutate_methane_to_methane_cc

    configuration = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-rsfe.yaml",
        input_dir=loeffler_testsystems_dir,
        output_dir=get_test_output_dir(),
    )
    output_files = mutate_methane_to_methane_cc(configuration=configuration)

    assert len(output_files) == 3


@pytest.mark.rsfe
@pytest.mark.requires_loeffler_systems
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_setup_system_for_methane_common_core_with_HMR():
    from transformato.testsystems import mutate_methane_to_methane_cc

    configuration = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-rsfe-HMR.yaml",
        input_dir=loeffler_testsystems_dir,
        output_dir=get_test_output_dir(),
    )
    output_files = mutate_methane_to_methane_cc(configuration=configuration)

    assert len(output_files) == 3
