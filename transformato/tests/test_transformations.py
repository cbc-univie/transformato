"""
Unit and regression test for the transformato package.
"""

import sys
import shutil

# Import package, test suite, and other packages as needed
import logging

# read in specific topology with parameters
from transformato import (
    load_config_yaml,
)
from transformato.tests.paths import get_test_output_dir
from transformato.constants import loeffler_testsystems_dir


def test_transformato_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "transformato" in sys.modules


def test_generate_output_for_acetylacetone_tautomer_pair(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_acetylacetone_tautomer_pair

    conf = "transformato/tests/config/test-acetylacetone-tautomer-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir=get_test_output_dir()
    )

    (output_files_t1, output_files_t2), _, _ = setup_acetylacetone_tautomer_pair(
        configuration=configuration, nr_of_bonded_windows=8
    )
    f = "/".join(output_files_t1[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


def test_generate_output_for_methane_cc_rsfe():
    from transformato.testsystems import mutate_methane_to_methane_cc

    configuration = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-rsfe.yaml",
        input_dir=loeffler_testsystems_dir,
        output_dir=get_test_output_dir(),
    )
    output_files = mutate_methane_to_methane_cc(configuration=configuration)

    assert len(output_files) == 3
    f = "/".join(output_files[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


def test_generate_output_for_toluene_cc_rsfe():
    from transformato.testsystems import mutate_toluene_to_methane_cc

    configuration = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-rsfe.yaml",
        input_dir=loeffler_testsystems_dir,
        output_dir=get_test_output_dir(),
    )

    output_files = mutate_toluene_to_methane_cc(configuration=configuration)

    f = "/".join(output_files[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


def test_generate_output_for_2MIN_cc_rsfe():
    from transformato.testsystems import mutate_2_methylindole_to_methane_cc

    configuration = load_config_yaml(
        config="transformato/tests/config/test-2MIN-methane-rsfe.yaml",
        input_dir=loeffler_testsystems_dir,
        output_dir=get_test_output_dir(),
    )

    output_files = mutate_2_methylindole_to_methane_cc(configuration=configuration)
    print(output_files)
    f = "/".join(output_files[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


def test_generate_output_for_2MFN_cc_rsfe():
    from transformato.testsystems import mutate_2_methylfuran_to_methane_cc

    configuration = load_config_yaml(
        config="transformato/tests/config/test-2MFN-methane-rsfe.yaml",
        input_dir=loeffler_testsystems_dir,
        output_dir=get_test_output_dir(),
    )

    output_files = mutate_2_methylfuran_to_methane_cc(configuration=configuration)
    print(output_files)
    f = "/".join(output_files[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


def test_generate_output_for_neopentane_cc_rsfe():
    from transformato.testsystems import mutate_neopentane_to_methane_cc

    configuration = load_config_yaml(
        config="transformato/tests/config/test-neopentane-methane-rsfe.yaml",
        input_dir=loeffler_testsystems_dir,
        output_dir=get_test_output_dir(),
    )

    output_files = mutate_neopentane_to_methane_cc(configuration=configuration)
    print(output_files)
    f = "/".join(output_files[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


def test_generate_output_for_methanol_cc_rsfe():
    from transformato.testsystems import mutate_methanol_to_methane_cc

    configuration = load_config_yaml(
        config="transformato/tests/config/test-methanol-methane-rsfe.yaml",
        input_dir=loeffler_testsystems_dir,
        output_dir=get_test_output_dir(),
    )

    output_files = mutate_methanol_to_methane_cc(configuration=configuration)
    print(output_files)
    f = "/".join(output_files[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


def test_generate_output_for_2_CPI_rsfe():
    from transformato.testsystems import mutate_2_CPI_to_7_CPI_cc

    configuration = load_config_yaml(
        config="transformato/tests/config/test-7-CPI-2-CPI-rsfe.yaml",
        input_dir=loeffler_testsystems_dir,
        output_dir=get_test_output_dir(),
    )

    output_files = mutate_2_CPI_to_7_CPI_cc(configuration=configuration)
    print(output_files)
    f = "/".join(output_files[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


def test_generate_output_for_7_CPI_rsfe():
    from transformato.testsystems import mutate_7_CPI_to_2_CPI_cc

    configuration = load_config_yaml(
        config="transformato/tests/config/test-7-CPI-2-CPI-rsfe.yaml",
        input_dir=loeffler_testsystems_dir,
        output_dir=get_test_output_dir(),
    )

    output_files = mutate_7_CPI_to_2_CPI_cc(configuration=configuration)
    print(output_files)
    f = "/".join(output_files[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)
