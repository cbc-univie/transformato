"""
Unit and regression test for the transformato package.
"""

import sys
import shutil

# Import package, test suite, and other packages as needed
import logging
import pytest

# read in specific topology with parameters
from transformato import (
    load_config_yaml,
)
from transformato.tests.paths import get_test_output_dir
from transformato.constants import loeffler_testsystems_dir
from transformato.utils import run_simulation


def test_transformato_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "transformato" in sys.modules


@pytest.mark.rsfe
def test_rsfe_mutate_acetylaceton_methyl_common_core():
    from ..testsystems import mutate_acetylaceton_methyl_common_core

    conf = "transformato/tests/config/test-acetylacetone-tautomer-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir=get_test_output_dir()
    )  # NOTE: for preprocessing input_dir is the output dir

    mutate_acetylaceton_methyl_common_core(configuration=configuration)


@pytest.mark.rbfe
def test_rbfe_mutate_2oj9():
    from ..mutate import mutate_pure_tautomers
    from transformato import (
        ProposeMutationRoute,
        SystemStructure,
        load_config_yaml,
    )

    conf_path = "transformato/tests/config/test-2oj9-tautomer-pair-rbfe.yaml"

    configuration = load_config_yaml(
        config=conf_path, input_dir="data/", output_dir=get_test_output_dir()
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()
    return (
        mutate_pure_tautomers(
            s1_to_s2,
            s1,
            s2,
            configuration,
            nr_of_bonded_windows=4,
        ),
        configuration,
        s1_to_s2,
    )


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


def test_generate_output_for_1a0q_1a07_rsfe(caplog):
    import logging

    # Test that TF can handel multiple dummy regions
    caplog.set_level(logging.INFO)
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
    assert i.current_step == 24