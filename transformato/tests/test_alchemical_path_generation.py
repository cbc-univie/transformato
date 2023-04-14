"""
Unit and regression test for the transformato package.
"""

# Import package, test suite, and other packages as needed
import logging
import os
import sys
import warnings

import pytest

# read in specific topology with parameters
from transformato import (
    IntermediateStateFactory,
    ProposeMutationRoute,
    SystemStructure,
    load_config_yaml,
)
from transformato.mutate import perform_mutations
from transformato.tests.paths import get_test_output_dir
from transformato_testsystems.testsystems import get_testsystems_dir

warnings.filterwarnings("ignore", module="parmed")


def test_transformato_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "transformato" in sys.modules


@pytest.mark.rsfe
def test_generate_alchemical_path_acetylaceton_methyl_common_core():
    from transformato_testsystems.testsystems import (
        mutate_acetylaceton_methyl_common_core,
    )

    conf = f"{get_testsystems_dir()}/config/test-acetylacetone-tautomer-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf,
        input_dir=get_testsystems_dir(),
        output_dir=get_test_output_dir(),
    )  # NOTE: for preprocessing input_dir is the output dir

    mutate_acetylaceton_methyl_common_core(configuration=configuration)


@pytest.mark.rbfe
def test_rbfe_mutate_2oj9():
    from transformato import ProposeMutationRoute, SystemStructure, load_config_yaml

    from ..mutate import mutate_pure_tautomers

    conf = f"{get_testsystems_dir()}/config/test-2oj9-tautomer-pair-rbfe.yaml"

    configuration = load_config_yaml(
        config=conf,
        input_dir=get_testsystems_dir(),
        output_dir=get_test_output_dir(),
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


@pytest.mark.rsfe
def test_generate_alchemical_path_for_acetylacetone_tautomer_pair(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_acetylacetone_tautomer_pair

    conf = f"{get_testsystems_dir()}/config/test-acetylacetone-tautomer-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir=get_testsystems_dir(), output_dir=get_test_output_dir()
    )

    (output_files_t1, output_files_t2), _, _ = setup_acetylacetone_tautomer_pair(
        configuration=configuration, nr_of_bonded_windows=8
    )


@pytest.mark.rsfe
def test_generate_alchemical_path_for_toluene_commmon_core():
    from transformato_testsystems.testsystems import mutate_toluene_to_methane_cc

    configuration = load_config_yaml(
        config=f"{get_testsystems_dir()}/config/test-toluene-methane-rsfe.yaml",
        input_dir=get_testsystems_dir(),
        output_dir=get_test_output_dir(),
    )

    output_files = mutate_toluene_to_methane_cc(configuration=configuration)
    assert len(output_files) == 16

    with open(
        f"{get_test_output_dir()}/toluene-methane-rsfe/toluene/intst1/openmm_run.py",
        "r",
    ) as f:
        i = 0
        for line in f.readlines():
            if "# Set platform" in line and i == 0:
                i += 1
            elif i == 1:
                assert "CPU" in line
                print(line)
                i += 1


@pytest.mark.rsfe
def test_generate_alchemical_path_for_toluene_commmon_core_with_CUDA():
    from transformato_testsystems.testsystems import mutate_toluene_to_methane_cc

    configuration = load_config_yaml(
        config=f"{get_testsystems_dir()}/config/test-toluene-methane-rsfe-CUDA.yaml",
        input_dir=get_testsystems_dir(),
        output_dir=get_test_output_dir(),
    )

    output_files = mutate_toluene_to_methane_cc(configuration=configuration)
    assert len(output_files) == 16

    with open(
        f"{get_test_output_dir()}/toluene-methane-rsfe/toluene/intst1/openmm_run.py",
        "r",
    ) as f:
        i = 0
        for line in f.readlines():
            if "# Set platform" in line and i == 0:
                i += 1
            elif i == 1:
                assert "CUDA" in line
                print(line)
                i += 1
            elif i == 2:
                print(line)
                assert "dict(CudaPrecision='mixed')" in line
                i += 1


@pytest.mark.rsfe
def test_generate_alchemical_path_for_2MIN_common_core():
    from transformato_testsystems.testsystems import mutate_2_methylindole_to_methane_cc

    configuration = load_config_yaml(
        config=f"{get_testsystems_dir()}/config/test-2MIN-methane-rsfe.yaml",
        input_dir=get_testsystems_dir(),
        output_dir=get_test_output_dir(),
    )

    output_files = mutate_2_methylindole_to_methane_cc(configuration=configuration)
    assert len(output_files) == 21


@pytest.mark.rsfe
def test_generate_alchemical_path_for_2MFN_common_core():
    from transformato_testsystems.testsystems import mutate_2_methylfuran_to_methane_cc

    configuration = load_config_yaml(
        config=f"{get_testsystems_dir()}/config/test-2MFN-methane-rsfe.yaml",
        input_dir=get_testsystems_dir(),
        output_dir=get_test_output_dir(),
    )

    output_files = mutate_2_methylfuran_to_methane_cc(configuration=configuration)
    assert len(output_files) == 17


@pytest.mark.rsfe
def test_generate_alchemical_path_for_neopentane_common_core():
    from transformato_testsystems.testsystems import mutate_neopentane_to_methane_cc

    configuration = load_config_yaml(
        config=f"{get_testsystems_dir()}/config/test-neopentane-methane-rsfe.yaml",
        input_dir=get_testsystems_dir(),
        output_dir=get_test_output_dir(),
    )

    output_files = mutate_neopentane_to_methane_cc(configuration=configuration)
    assert len(output_files) == 16


@pytest.mark.rsfe
def test_generate_alchemical_path_for_methanol_common_core():
    from transformato_testsystems.testsystems import perform_generic_mutation

    configuration = load_config_yaml(
        config=f"{get_testsystems_dir()}/config/test-methanol-methane-rsfe.yaml",
        input_dir=get_testsystems_dir(),
        output_dir=get_test_output_dir(),
    )
    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )
    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
        nr_of_mutation_steps_charge=3,
    )
    assert len(i.output_files) == 11


@pytest.mark.rsfe
def test_generate_alchemical_path_for_2_CPI_to_common_core():
    from transformato_testsystems.testsystems import mutate_2_CPI_to_7_CPI_cc

    configuration = load_config_yaml(
        config=f"{get_testsystems_dir()}/config/test-7-CPI-2-CPI-rsfe.yaml",
        input_dir=get_testsystems_dir(),
        output_dir=get_test_output_dir(),
    )

    output_files = mutate_2_CPI_to_7_CPI_cc(configuration=configuration)
    assert len(output_files) == 17


@pytest.mark.rsfe
def test_generate_alchemical_path_for_7_CPI_to_common_core():
    from transformato_testsystems.testsystems import mutate_7_CPI_to_2_CPI_cc

    configuration = load_config_yaml(
        config=f"{get_testsystems_dir()}/config/test-7-CPI-2-CPI-rsfe.yaml",
        input_dir=get_testsystems_dir(),
        output_dir=get_test_output_dir(),
    )

    output_files = mutate_7_CPI_to_2_CPI_cc(configuration=configuration)
    assert len(output_files) == 12


@pytest.mark.rsfe
def test_generate_alchemical_path_for_1a0q_1a07(caplog):
    # Test that TF can handel multiple dummy regions
    caplog.set_level(logging.INFO)
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
    assert i.current_step == 24


@pytest.mark.rbfe
def test_generate_path_for_ppar_system():
    conf = f"{get_testsystems_dir()}/config/ppar-cpd31_cpd25.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir=get_testsystems_dir(), output_dir=get_test_output_dir()
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.propose_common_core()
    s1_to_s2.finish_common_core()

    assert (
        len(s1_to_s2.get_common_core_idx_mol1())
        == len(s1_to_s2.get_common_core_idx_mol2())
        == 46
    )

    assert s1_to_s2.get_idx_not_in_common_core_for_mol1() == [31, 32, 33, 49, 50, 51]
    assert s1_to_s2.get_idx_not_in_common_core_for_mol2() == [0, 1, 40, 49, 50, 51]

    assert s1_to_s2.terminal_real_atom_cc1 == [19, 12]
    assert s1_to_s2.terminal_real_atom_cc2 == [21, 14]

    assert s1_to_s2.dummy_region_cc1.lj_default == [31, 33]
    assert s1_to_s2.dummy_region_cc2.lj_default == [40, 1]
    assert s1_to_s2.dummy_region_cc2.connected_dummy_regions == [
        [49, 50, 51, 0, 1],
        [40],
    ]
    assert s1_to_s2.dummy_region_cc1.connected_dummy_regions == [
        [33],
        [51, 50, 49, 32, 31],
    ]

    assert s1_to_s2.matching_terminal_atoms_between_cc == [(12, 14), (19, 21)]
