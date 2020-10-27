"""
Unit and regression test for the transformato package.
"""

import copy
import logging
import os
import pathlib
import shutil
import subprocess
import sys

import numpy as np
import parmed as pm
import pytest

from io import StringIO
import filecmp

# Import package, test suite, and other packages as needed
import transformato

# read in specific topology with parameters
from parmed.charmm.parameters import CharmmParameterSet
from transformato import (
    IntermediateStateFactory,
    ProposeMutationRoute,
    SystemStructure,
    load_config_yaml,
    psf_correction,
)


def read_params(output_file_base):
    extlist = ["rtf", "prm", "str"]
    print(output_file_base)
    parFiles = ()
    toppar_file = f"{output_file_base}/toppar.str"
    for line in open(toppar_file, "r"):
        if "!" in line:
            line = line.split("!")[0]
        parfile = line.strip()
        if len(parfile) != 0:
            ext = parfile.lower().split(".")[-1]
            if not ext in extlist:
                continue
            parFiles += (f"{output_file_base}/{parfile}",)

    params = CharmmParameterSet(*parFiles)
    return params


def generate_psf(output_file_base, env):
    parms = read_params(output_file_base)
    target_psf = pm.charmm.CharmmPsfFile(f"{output_file_base}/lig_in_{env}.psf")
    target_psf.load_parameters(parms)
    return target_psf


def test_transformato_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "transformato" in sys.modules


def test_read_yaml():
    """Sample test, will check ability to read yaml files"""
    settingsMap = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
        input_dir=".",
        output_dir="data/",
    )

    assert settingsMap["system"]["name"] == "toluene-methane-solvation-free-energy"


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


def test_initialize_systems():
    configuration = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
        input_dir="data/",
        output_dir=".",
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
        config="transformato/tests/config/test-2oj9-binding-free-energy.yaml",
        input_dir="data/",
        output_dir=".",
    )

    s1 = SystemStructure(configuration, "structure1")
    assert int(s1.offset["waterbox"]) == 0
    assert int(s1.offset["complex"]) == 4811

    s2 = SystemStructure(configuration, "structure2")
    assert int(s2.offset["waterbox"]) == 0
    assert int(s2.offset["complex"]) == 4692

    assert "complex" in s1.envs and "complex" in s2.envs
    assert "waterbox" in s1.envs and "waterbox" in s2.envs


def test_proposed_mutation_mcs():

    from rdkit.Chem import rdFMCS

    for conf in [
        "transformato/tests/config/test-2oj9-solvation-free-energy.yaml",
        "transformato/tests/config/test-2oj9-binding-free-energy.yaml",
    ]:
        configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        a = ProposeMutationRoute(s1, s2)
        # find mcs
        a._find_mcs("m1", "m2")

        assert str(a.s1_tlc) == "BMI"
        assert str(a.s2_tlc) == "UNK"

        cc1 = set(
            [
                0,
                3,
                6,
                5,
                4,
                14,
                24,
                23,
                26,
                22,
                25,
                17,
                16,
                28,
                27,
                29,
                46,
                47,
                48,
                45,
                41,
                44,
                2,
                7,
                11,
                9,
                8,
                10,
                13,
                12,
                39,
                38,
                36,
                37,
                34,
                35,
                30,
                32,
                33,
                31,
            ]
        )
        cc2 = set(
            [
                0,
                3,
                6,
                5,
                4,
                14,
                19,
                18,
                21,
                17,
                20,
                16,
                15,
                23,
                22,
                24,
                43,
                44,
                45,
                42,
                40,
                41,
                2,
                7,
                11,
                9,
                8,
                10,
                13,
                12,
                38,
                37,
                35,
                36,
                33,
                34,
                29,
                31,
                32,
                30,
            ]
        )
        assert set(a.get_common_core_idx_mol1()) == cc1
        assert set(a.get_common_core_idx_mol2()) == cc2

        a.bondCompare = rdFMCS.BondCompare.CompareOrder
        # find mcs
        a._find_mcs("m1", "m2")

        cc1 = set(
            [
                0,
                3,
                6,
                33,
                31,
                14,
                24,
                23,
                26,
                22,
                25,
                17,
                16,
                28,
                27,
                29,
                46,
                47,
                48,
                45,
                41,
                44,
            ]
        )
        assert set(a.get_common_core_idx_mol1()) == cc1
        cc2 = set(
            [
                0,
                3,
                6,
                32,
                30,
                14,
                19,
                18,
                21,
                17,
                20,
                16,
                15,
                23,
                22,
                24,
                43,
                44,
                45,
                42,
                40,
                41,
            ]
        )
        assert set(a.get_common_core_idx_mol2()) == cc2

    for conf in [
        "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml"
    ]:
        configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        a = ProposeMutationRoute(s1, s2)
        # find mcs
        a._find_mcs("m1", "m2")

        assert str(a.s1_tlc) == "UNL"
        assert str(a.s2_tlc) == "LIG"

        print(set(a.get_common_core_idx_mol1()))
        print(set(a.get_common_core_idx_mol2()))
        cc1 = set([8, 5, 6, 7])
        cc2 = set([0, 1, 2, 3])
        assert set(a.get_common_core_idx_mol1()) == cc1
        assert set(a.get_common_core_idx_mol2()) == cc2

        a.bondCompare = rdFMCS.BondCompare.CompareOrder
        # find mcs
        a._find_mcs("m1", "m2")

        cc1 = set([8, 5, 6, 7])
        cc2 = set([0, 1, 2, 3])
        assert set(a.get_common_core_idx_mol1()) == cc1
        assert set(a.get_common_core_idx_mol2()) == cc2


def test_proposed_mutation_terminal_dummy_real_atom_match():
    from rdkit.Chem import rdFMCS

    for conf in [
        "transformato/tests/config/test-7-CPI-2-CPI-solvation-free-energy.yaml",
    ]:
        configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        a = ProposeMutationRoute(s1, s2)
        # find mcs
        a.bondCompare = rdFMCS.BondCompare.CompareOrderExact
        a.completeRingsOnly = True
        a._find_mcs("m1", "m2")
        a.remove_idx_from_common_core_of_mol1([14])
        a.remove_idx_from_common_core_of_mol2([6])

        # find terminal dummy/real atoms
        a._set_common_core_parameters()
        # match terminal real/dummy atoms
        print(a.terminal_real_atom_cc1)
        print(a.terminal_real_atom_cc2)
        assert set(a.terminal_real_atom_cc1) == set([15])
        assert set(a.terminal_dummy_atom_cc1) == set([14])
        assert set(a.terminal_real_atom_cc2) == set([15])
        assert set(a.terminal_dummy_atom_cc2) == set([6])

        match_terminal_atoms_cc1 = a._match_terminal_real_and_dummy_atoms_for_mol1()
        match_terminal_atoms_cc2 = a._match_terminal_real_and_dummy_atoms_for_mol2()

        # are the correct terminal common core atoms identified?
        assert match_terminal_atoms_cc1[15] == set([14])
        assert match_terminal_atoms_cc2[15] == set([6])

        # terminal atoms match between the two common cores
        assert a.matching_terminal_atoms_between_cc[0] == (15, 15)
        # INFO     transformato.mutate:mutate.py:139 Matching terminal atoms from cc1 to cc2. cc1: 0 : cc2: 0
        # INFO     transformato.mutate:mutate.py:139 Matching terminal atoms from cc1 to cc2. cc1: 16 : cc2: 15


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_find_connected_dummy_regions1():

    from rdkit.Chem import rdFMCS

    conf = "transformato/tests/config/test-7-CPI-2-CPI-solvation-free-energy.yaml"
    configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    a = ProposeMutationRoute(s1, s2)
    # find mcs
    a.bondCompare = rdFMCS.BondCompare.CompareOrderExact
    a.completeRingsOnly = True
    a._find_mcs("m1", "m2")
    a.remove_idx_from_common_core_of_mol1([14])
    a.remove_idx_from_common_core_of_mol2([6])

    # find terminal dummy/real atoms
    a._set_common_core_parameters()
    # match terminal real/dummy atoms
    match_terminal_atoms_cc1 = a._match_terminal_real_and_dummy_atoms_for_mol1()
    match_terminal_atoms_cc2 = a._match_terminal_real_and_dummy_atoms_for_mol2()
    print(match_terminal_atoms_cc1)
    print(match_terminal_atoms_cc2)

    # find connected dummy regions
    connected_dummy_regions_cc1 = a._find_connected_dummy_regions(
        "m1", match_terminal_atoms_cc1
    )
    connected_dummy_regions_cc2 = a._find_connected_dummy_regions(
        "m2", match_terminal_atoms_cc2
    )

    lj_default_cc1, lj_default_cc2 = a._match_terminal_dummy_atoms_between_common_cores(
        match_terminal_atoms_cc1,
        match_terminal_atoms_cc2,
        connected_dummy_regions_cc1,
        connected_dummy_regions_cc2,
    )

    print(connected_dummy_regions_cc1)
    print(connected_dummy_regions_cc2)


def test_find_connected_dummy_regions2():

    from rdkit.Chem import rdFMCS

    ##################################################
    conf = "transformato/tests/config/test-2oj9-solvation-free-energy.yaml"
    configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    a = ProposeMutationRoute(s1, s2)
    # find mcs
    a._find_mcs("m1", "m2")
    a._set_common_core_parameters()
    # match the real/dummy atoms
    match_terminal_atoms_cc1 = a._match_terminal_real_and_dummy_atoms_for_mol1()
    match_terminal_atoms_cc2 = a._match_terminal_real_and_dummy_atoms_for_mol2()

    # find connected dummy regions
    connected_dummy_regions_cc1 = a._find_connected_dummy_regions(
        "m1", match_terminal_atoms_cc1
    )
    connected_dummy_regions_cc2 = a._find_connected_dummy_regions(
        "m2", match_terminal_atoms_cc2
    )

    lj_default_cc1, lj_default_cc2 = a._match_terminal_dummy_atoms_between_common_cores(
        match_terminal_atoms_cc1,
        match_terminal_atoms_cc2,
        connected_dummy_regions_cc1,
        connected_dummy_regions_cc2,
    )

    dummy_region_m1 = transformato.mutate.DummyRegion(
        "m1",
        match_terminal_atoms_cc1,
        connected_dummy_regions_cc1,
        s1.tlc,
        lj_default=lj_default_cc1,
    )

    print(connected_dummy_regions_cc1)

    dummy_region_m2 = transformato.mutate.DummyRegion(
        "m2",
        match_terminal_atoms_cc2,
        connected_dummy_regions_cc2,
        s2.tlc,
        lj_default=lj_default_cc2,
    )

    print(connected_dummy_regions_cc2)

    # match
    assert dummy_region_m1.connected_dummy_regions[0] == {
        40,
        42,
        43,
        15,
        18,
        19,
        20,
        21,
    }

    assert dummy_region_m1.connected_dummy_regions[1] == {1}
    assert dummy_region_m2.connected_dummy_regions[0] == {1, 26, 27, 28}
    assert dummy_region_m2.connected_dummy_regions[1] == {25}
    assert dummy_region_m2.connected_dummy_regions[2] == {39}

    # return real atom that connects dummy region to mol
    assert (
        dummy_region_m1.return_connecting_real_atom(
            dummy_region_m1.connected_dummy_regions[0]
        )
        == 16
    )

    assert (
        dummy_region_m1.return_connecting_real_atom(
            dummy_region_m1.connected_dummy_regions[1]
        )
        == 0
    )

    print(f"Matched dummy region: {dummy_region_m1.lj_default}")
    print(f"Matched dummy region: {dummy_region_m2.lj_default}")
    assert dummy_region_m1.lj_default == [1, 18]
    assert dummy_region_m2.lj_default == [1, 39]


def test_common_core_system1():
    from rdkit.Chem import rdFMCS

    for conf in [
        "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
        "transformato/tests/config/test-neopentane-methane-solvation-free-energy.yaml",
        "transformato/tests/config/test-ethane-methanol-solvation-free-energy.yaml",
    ]:
        configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        a = ProposeMutationRoute(s1, s2)
        if (
            conf
            == "transformato/tests/config/test-neopentane-methane-solvation-free-energy.yaml"
        ):
            a.propose_common_core()
            a.finish_common_core(
                connected_dummy_regions_cc1=[
                    {0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
                ]
            )
        else:
            # find mcs, find terminal dummy/real atoms, generate charge compensated psfs
            a.calculate_common_core()


def test_mutation_list():
    from ..utils import print_mutations
    from rdkit.Chem import rdFMCS

    for conf, system_name in zip(
        [
            "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
            "transformato/tests/config/test-neopentane-methane-solvation-free-energy.yaml",
            "transformato/tests/config/test-ethane-methanol-solvation-free-energy.yaml",
        ],
        ["toluene-methane", "neopentane-methane", "ethane-ethanol"],
    ):

        configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        s1_to_s2 = ProposeMutationRoute(s1, s2)

        if system_name == "ethane-methanol":
            s1_to_s2.calculate_common_core()
            mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
            assert set(mutation_list.keys()) == set(
                ["charge", "terminal-lj", "transform"]
            )
        if system_name == "toluene-methane":
            s1_to_s2.calculate_common_core()
            mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
            assert set(mutation_list.keys()) == set(
                ["charge", "hydrogen-lj", "lj", "transform", "terminal-lj"]
            )

        if system_name == "neopentane-methane":
            s1_to_s2.propose_common_core()
            s1_to_s2.finish_common_core(
                connected_dummy_regions_cc1=[
                    {0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
                ]
            )
            mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
            assert set(mutation_list.keys()) == set(
                ["charge", "hydrogen-lj", "lj", "transform", "terminal-lj"]
            )


def test_endpoint_mutation():
    # test physical endpoint systems
    from ..utils import print_mutations
    from ..mutate import MutationDefinition

    for conf in [
        "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
        "transformato/tests/config/test-ethane-methanol-solvation-free-energy.yaml",
    ]:
        configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        s1_to_s2 = ProposeMutationRoute(s1, s2)
        s2_to_s1 = ProposeMutationRoute(s2, s1)
        try:
            for a, system in zip([s1_to_s2, s2_to_s1], [s1, s2]):
                a.calculate_common_core()
                i = IntermediateStateFactory(
                    system=system,
                    configuration=configuration,
                )
                output_file_base = i.write_state(mutation_conf=[], intst_nr=0)

        finally:
            pass
            # shutil.rmtree(output_file_base)


def test_charge_mutation_test_system1():
    from ..utils import print_mutations
    from ..mutate import MutationDefinition

    for conf, system_name in zip(
        [
            "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
            "transformato/tests/config/test-ethane-methanol-solvation-free-energy.yaml",
        ],
        ["toluene-methane" "ethane-ethanol"],
    ):

        # try writing endstate in all directions
        configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        s1_to_s2 = ProposeMutationRoute(s1, s2)
        s2_to_s1 = ProposeMutationRoute(s2, s1)
        for a, system in zip([s1_to_s2, s2_to_s1], [s1, s2]):
            i = IntermediateStateFactory(
                system=system,
                configuration=configuration,
            )
            output_file_base = i.write_state(mutation_conf=[], intst_nr=0)

            # original psfs without charge change
            original_psf = {}
            for env in system.envs:
                original_psf[env] = copy.deepcopy(system.psfs[env])

            try:
                for env in system.envs:
                    offset = system.offset[env]

                    mutated_psf = generate_psf(output_file_base, env)
                    for atom in system.mol.GetAtoms():
                        idx = atom.GetIdx()
                        assert np.isclose(
                            original_psf[env].atoms[idx + offset].charge,
                            mutated_psf.atoms[idx + offset].charge,
                            rtol=1e-3,
                        )
            finally:
                shutil.rmtree(output_file_base)


def test_charge_mutation_test_system2():
    from ..utils import print_mutations
    from ..mutate import MutationDefinition

    for conf, system_name in zip(
        [
            "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
            "transformato/tests/config/test-neopentane-methane-solvation-free-energy.yaml",
            "transformato/tests/config/test-ethane-methanol-solvation-free-energy.yaml",
        ],
        ["toluene-methane", "neopentane-methane", "ethane-methanol"],
    ):

        configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
        # scale charges with 0.5
        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        s1_to_s2 = ProposeMutationRoute(s1, s2)
        # original psfs without charge change
        original_psf = {}
        for env in s1.envs:
            original_psf[env] = copy.deepcopy(s1.psfs[env])

        for lambda_charge in [1.0, 0.5, 0.25, 0.0]:
            for a, system in zip([s1_to_s2], [s1]):
                if system_name == "neopentane-methane":
                    a.propose_common_core()
                    a.finish_common_core(
                        connected_dummy_regions_cc1=[
                            {0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
                        ]
                    )
                else:
                    a.calculate_common_core()

                i = IntermediateStateFactory(
                    system=system,
                    configuration=configuration,
                )

                mutation_list = a.generate_mutations_to_common_core_for_mol1()

                try:
                    charges = mutation_list["charge"]
                    output_file_base = i.write_state(
                        mutation_conf=charges,
                        lambda_value_electrostatic=lambda_charge,
                        intst_nr=0,
                    )
                    for env in system.envs:
                        offset = system.offset[env]
                        # read in newly generated psf
                        mutated_psf = generate_psf(output_file_base, env)
                        for idx in charges[0].atoms_to_be_mutated:
                            assert np.isclose(
                                original_psf[env].atoms[idx + offset].charge
                                * lambda_charge,
                                mutated_psf.atoms[idx + offset].charge,
                                rtol=1e-03,
                            )
                finally:
                    shutil.rmtree(output_file_base)


def test_vdw_mutation_for_hydrogens_system1():
    # testing terminal lj
    from rdkit.Chem import rdFMCS

    for conf, system_name in zip(
        [
            "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
            "transformato/tests/config/test-neopentane-methane-solvation-free-energy.yaml",
            "transformato/tests/config/test-ethane-methanol-solvation-free-energy.yaml",
        ],
        ["toluene-methane", "neopentane-methane", "ethane-methanol"],
    ):
        configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        s1_to_s2 = ProposeMutationRoute(s1, s2)
        if system_name == "neopentane-methane":
            s1_to_s2.propose_common_core()
            s1_to_s2.finish_common_core(
                connected_dummy_regions_cc1=[
                    {0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
                ]
            )
        else:
            s1_to_s2.calculate_common_core()

        mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
        i = IntermediateStateFactory(
            system=s1,
            configuration=configuration,
        )

        original_psf = {}
        for env in s1.envs:
            original_psf[env] = copy.deepcopy(s1.psfs[env])

        try:
            terminal_lj_mutations = mutation_list["terminal-lj"]
            output_file_base = i.write_state(
                mutation_conf=terminal_lj_mutations,
                lambda_value_vdw=0.0,
                intst_nr=0,
            )
            print("Set epsilon/rmin to zero for selected atoms")

            for env in s1.envs:
                new_psf = generate_psf(output_file_base, env)
                # are the terminal lj parameter correctly set?
                offset = s1.offset[env]
                for terminal_lj in terminal_lj_mutations:
                    for idx in terminal_lj.vdw_atom_idx:
                        idxo = idx + offset
                        assert np.isclose(
                            -0.15,
                            new_psf.atoms[idxo].epsilon,
                            rtol=1e-3,
                        )
                        assert np.isclose(
                            1.5,
                            new_psf.atoms[idxo].rmin,
                            rtol=1e-3,
                        )

                # make sure that all other idx are not touched
                for idx in range(len(original_psf[env].atoms)):
                    idxo = idx - offset  # NOTE: the '-'
                    if idxo not in terminal_lj.vdw_atom_idx:
                        assert np.isclose(
                            original_psf[env].atoms[idx].epsilon,
                            new_psf.atoms[idx].epsilon,
                            rtol=1e-3,
                        )
                        assert np.isclose(
                            original_psf[env].atoms[idx].rmin,
                            new_psf.atoms[idx].rmin,
                            rtol=1e-3,
                        )

        finally:
            pass
            # shutil.rmtree(output_file_base)


def test_vdw_mutation_for_hydrogens_system2():
    from rdkit.Chem import rdFMCS

    for conf, system_name in zip(
        [
            "transformato/tests/config/test-7-CPI-2-CPI-solvation-free-energy.yaml",
        ],
        ["7-CPI-2-CPI"],
    ):
        configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        s1_to_s2 = ProposeMutationRoute(s1, s2)
        s1_to_s2.bondCompare = rdFMCS.BondCompare.CompareOrderExact
        s1_to_s2.completeRingsOnly = True
        s1_to_s2.calculate_common_core()

        mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
        i = IntermediateStateFactory(
            system=s1,
            configuration=configuration,
        )

        original_psf = {}
        for env in s1.envs:
            original_psf[env] = copy.deepcopy(s1.psfs[env])

        for lambda_value in [1.0, 0.5, 0.25]:
            try:
                # lj hydrogens scaling
                hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
                output_file_base = i.write_state(
                    mutation_conf=hydrogen_lj_mutations,
                    lambda_value_vdw=lambda_value,
                    intst_nr=0,
                )
                print("Set epsilon/rmin for selected atoms")

                for env in s1.envs:
                    new_psf = generate_psf(output_file_base, env)
                    # are the terminal lj parameter correctly set?
                    offset = s1.offset[env]
                    for hydrogen_lj in hydrogen_lj_mutations:
                        for idx in hydrogen_lj.vdw_atom_idx:
                            idxo = idx + offset
                            assert np.isclose(
                                original_psf[env].atoms[idxo].epsilon * lambda_value,
                                new_psf.atoms[idxo].epsilon,
                                rtol=1e-3,
                            )
                            assert np.isclose(
                                original_psf[env].atoms[idxo].rmin * lambda_value,
                                new_psf.atoms[idxo].rmin,
                                rtol=1e-3,
                            )

                    # make sure that all other idx are not touched
                    for idx in range(len(original_psf[env].atoms)):
                        idxo = idx - offset  # NOTE: the '-'
                        if idxo not in hydrogen_lj.vdw_atom_idx:
                            assert np.isclose(
                                original_psf[env].atoms[idx].epsilon,
                                new_psf.atoms[idx].epsilon,
                                rtol=1e-3,
                            )
                            assert np.isclose(
                                original_psf[env].atoms[idx].rmin,
                                new_psf.atoms[idx].rmin,
                                rtol=1e-3,
                            )

            finally:
                shutil.rmtree(output_file_base)


def test_vdw_mutation_for_hydrogens_and_heavy_atoms():
    from rdkit.Chem import rdFMCS

    for conf, system_name in zip(
        [
            "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
            "transformato/tests/config/test-neopentane-methane-solvation-free-energy.yaml",
            "transformato/tests/config/test-ethane-methanol-solvation-free-energy.yaml",
        ],
        ["toluene-methane", "neopentane-methane", "ethane-methanol"],
    ):
        configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        s1_to_s2 = ProposeMutationRoute(s1, s2)
        if system_name == "neopentane-methane":
            s1_to_s2.propose_common_core()
            s1_to_s2.finish_common_core(
                connected_dummy_regions_cc1=[
                    {0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
                ]
            )
        elif system_name == "7-CPI-2-CPI":
            # find mcs
            s1_to_s2.bondCompare = rdFMCS.BondCompare.CompareOrderExact
            s1_to_s2.completeRingsOnly = True
            s1_to_s2.calculate_common_core()
        else:
            s1_to_s2.calculate_common_core()

        mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
        i = IntermediateStateFactory(
            system=s1,
            configuration=configuration,
        )

        original_psf = {}
        for env in s1.envs:
            original_psf[env] = copy.deepcopy(s1.psfs[env])

        #
        for lambda_vdw in [1.0, 0.0]:
            print(f"Lambda: {lambda_vdw}")
            output_files = []

            intst_nr = 0
            try:

                all_atoms_for_which_lj_turned_off = []
                print(mutation_list.keys())
                # turn off hydrogen lj
                hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
                print(
                    f"Turn off lj for hydrogen atoms : {[e.atoms_to_be_mutated for e in hydrogen_lj_mutations]}"
                )
                output_file_base = i.write_state(
                    mutation_conf=hydrogen_lj_mutations,
                    lambda_value_vdw=lambda_vdw,
                    intst_nr=intst_nr,
                )
                output_files.append(output_file_base)
                for mutation in hydrogen_lj_mutations:
                    all_atoms_for_which_lj_turned_off.extend(mutation.vdw_atom_idx)

                # turn off heavy atom lj
                for mutation in mutation_list["lj"]:
                    intst_nr += 1
                    print(
                        f"Turn off lj for heavy atom : {mutation.atoms_to_be_mutated}"
                    )

                    output_file_base = i.write_state(
                        mutation_conf=[mutation],
                        lambda_value_vdw=lambda_vdw,
                        intst_nr=intst_nr,
                    )
                    all_atoms_for_which_lj_turned_off.extend(mutation.vdw_atom_idx)
                    output_files.append(output_file_base)

                # change to default lj
                intst_nr += 1
                terminal_lj_mutations = mutation_list["terminal-lj"]
                terminal_idx = []
                for mutation in terminal_lj_mutations:
                    all_atoms_for_which_lj_turned_off.extend(mutation.vdw_atom_idx)
                    terminal_idx.extend(mutation.vdw_atom_idx)
                print(
                    f"Turn off lj for terminal atom : {all_atoms_for_which_lj_turned_off}"
                )

                output_file_base = i.write_state(
                    mutation_conf=terminal_lj_mutations,
                    lambda_value_vdw=0.0,
                    intst_nr=intst_nr,
                )
                output_files.append(output_file_base)

                print(f"Set epsilon/rmin to base * {lambda_vdw} for selected atoms")

                print(all_atoms_for_which_lj_turned_off)
                for env in s1.envs:
                    print(env)
                    # read in generated psf at last mutation step
                    new_psf = generate_psf(output_file_base, env)
                    offset = s1.offset[env]
                    for idx in all_atoms_for_which_lj_turned_off:
                        idxo = idx + offset
                        # match terminal jl
                        if idx in terminal_idx:
                            assert np.isclose(
                                -0.15,
                                new_psf.atoms[idxo].epsilon,
                                rtol=1e-3,
                            )
                            assert np.isclose(
                                1.5,
                                new_psf.atoms[idxo].rmin,
                                rtol=1e-3,
                            )
                        else:
                            # match all other lj
                            assert np.isclose(
                                original_psf[env].atoms[idxo].epsilon * lambda_vdw,
                                new_psf.atoms[idxo].epsilon,
                                rtol=1e-3,
                            )
                            assert np.isclose(
                                original_psf[env].atoms[idxo].rmin * lambda_vdw,
                                new_psf.atoms[idxo].rmin,
                                rtol=1e-3,
                            )

                    # make sure that all other idx are not touched
                    for idx in range(s1.mol.GetNumAtoms()):
                        if idx not in all_atoms_for_which_lj_turned_off:
                            print(idx)
                            assert np.isclose(
                                original_psf[env].atoms[idx].epsilon,
                                new_psf.atoms[idx].epsilon,
                                rtol=1e-3,
                            )
                            assert np.isclose(
                                original_psf[env].atoms[idx].rmin,
                                new_psf.atoms[idx].rmin,
                                rtol=1e-3,
                            )

            finally:
                for output_file_base in output_files:
                    shutil.rmtree(output_file_base)


def test_full_mutation_system1():
    from rdkit.Chem import rdFMCS

    for conf, system_name in zip(
        [
            "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
            "transformato/tests/config/test-neopentane-methane-solvation-free-energy.yaml",
            "transformato/tests/config/test-ethane-methanol-solvation-free-energy.yaml",
        ],
        ["toluene-methane", "neopentane-methane", "ethane-methanol"],
    ):
        configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        s1_to_s2 = ProposeMutationRoute(s1, s2)
        if system_name == "neopentane-methane":
            s1_to_s2.propose_common_core()
            s1_to_s2.finish_common_core(
                connected_dummy_regions_cc1=[
                    {0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
                ]
            )
        else:
            s1_to_s2.calculate_common_core()

        mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
        i = IntermediateStateFactory(
            system=s1,
            configuration=configuration,
        )

        intst = 0
        charges = mutation_list["charge"]
        try:
            for lambda_value in np.linspace(0, 1, 5):
                print(lambda_value)
                intst += 1
                # turn off charges
                output_file_base = i.write_state(
                    mutation_conf=charges,
                    lambda_value_electrostatic=1 - lambda_value,
                    intst_nr=intst,
                )
        finally:
            shutil.rmtree(output_file_base)

        original_psf = {}
        for env in s1.envs:
            original_psf[env] = copy.deepcopy(s1.psfs[env])

        #
        lambda_vdw = 0.0
        print(f"Lambda: {lambda_vdw}")
        output_files = []
        try:
            intst += 1

            all_atoms_for_which_lj_turned_off = []
            print(mutation_list.keys())
            # Turn of hydrogens
            terminal_lj_mutations = mutation_list["terminal-lj"]
            output_file_base = i.write_state(
                mutation_conf=terminal_lj_mutations,
                lambda_value_vdw=lambda_vdw,
                intst_nr=intst,
            )
            output_files.append(output_file_base)
            for mutation in terminal_lj_mutations:
                all_atoms_for_which_lj_turned_off.extend(mutation.vdw_atom_idx)

        finally:
            for output_file_base in output_files:
                shutil.rmtree(output_file_base)


def test_full_mutation_system2():
    from rdkit.Chem import rdFMCS

    for conf, system_name in zip(
        [
            "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
            "transformato/tests/config/test-neopentane-methane-solvation-free-energy.yaml",
        ],
        ["toluene-methane", "neopentane-methane"],
    ):
        configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        s1_to_s2 = ProposeMutationRoute(s1, s2)
        if system_name == "neopentane-methane":
            s1_to_s2.propose_common_core()
            s1_to_s2.finish_common_core(
                connected_dummy_regions_cc1=[
                    {0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
                ]
            )
        else:
            s1_to_s2.calculate_common_core()

        mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
        i = IntermediateStateFactory(
            system=s1,
            configuration=configuration,
        )

        intst = 0
        charges = mutation_list["charge"]
        try:
            for lambda_value in np.linspace(0, 1, 5):
                intst += 1
                # turn off charges
                output_file_base = i.write_state(
                    mutation_conf=charges,
                    lambda_value_electrostatic=1 - lambda_value,
                    intst_nr=intst,
                )
        finally:
            shutil.rmtree(output_file_base)

        original_psf = {}
        for env in s1.envs:
            original_psf[env] = copy.deepcopy(s1.psfs[env])

        #
        lambda_vdw = 0.0
        print(f"Lambda: {lambda_vdw}")
        output_files = []
        try:
            intst += 1

            all_atoms_for_which_lj_turned_off = []
            print(mutation_list.keys())
            # Turn of hydrogens
            hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
            output_file_base = i.write_state(
                mutation_conf=hydrogen_lj_mutations,
                lambda_value_vdw=lambda_vdw,
                intst_nr=intst,
            )
            output_files.append(output_file_base)
            for mutation in hydrogen_lj_mutations:
                all_atoms_for_which_lj_turned_off.extend(mutation.vdw_atom_idx)

            # turn off heavy atoms
            for mutation in mutation_list["lj"]:
                intst += 1

                output_file_base = i.write_state(
                    mutation_conf=[mutation],
                    lambda_value_vdw=lambda_vdw,
                    intst_nr=intst,
                )
                all_atoms_for_which_lj_turned_off.extend(mutation.vdw_atom_idx)
                output_files.append(output_file_base)

            new_psf = generate_psf(output_file_base, env)
            print(f"Set epsilon/rmin to base * {lambda_vdw} for selected atoms")

            print(all_atoms_for_which_lj_turned_off)
            offset = s1.offset[env]
            for idx in all_atoms_for_which_lj_turned_off:
                idxo = idx + offset
                assert np.isclose(
                    original_psf[env].atoms[idxo].epsilon * lambda_vdw,
                    new_psf.atoms[idxo].epsilon,
                    rtol=1e-3,
                )
                assert np.isclose(
                    original_psf[env].atoms[idxo].rmin * lambda_vdw,
                    new_psf.atoms[idxo].rmin,
                    rtol=1e-3,
                )

            # make sure that all other idx are not touched
            for idx in range(len(original_psf[env].atoms)):
                idxo = idx - offset  # NOTE: the '-'
                if idxo not in all_atoms_for_which_lj_turned_off:
                    assert np.isclose(
                        original_psf[env].atoms[idx].epsilon,
                        new_psf.atoms[idx].epsilon,
                        rtol=1e-3,
                    )
                    assert np.isclose(
                        original_psf[env].atoms[idx].rmin,
                        new_psf.atoms[idx].rmin,
                        rtol=1e-3,
                    )

        finally:
            for output_file_base in output_files:
                shutil.rmtree(output_file_base)


@pytest.mark.slowtest
def test_bonded_mutation():

    for conf in [
        "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
    ]:
        configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        s1_to_s2 = ProposeMutationRoute(s1, s2)

        original_psf = {}
        original_psf = {}
        output_files = []
        for env in s1.envs:
            original_psf[env] = copy.deepcopy(s1.psfs[env])

        s1_to_s2.calculate_common_core()

        mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
        i = IntermediateStateFactory(
            system=s1,
            configuration=configuration,
        )

        # mutate everything else before touching bonded terms
        intst = 0
        charges = mutation_list["charge"]
        intst += 1
        # turn off charges
        output_file_base = i.write_state(
            mutation_conf=charges,
            lambda_value_electrostatic=0.0,
            intst_nr=intst,
        )
        output_files.append(output_file_base)

        # Turn of hydrogens
        intst += 1
        hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
        output_file_base = i.write_state(
            mutation_conf=hydrogen_lj_mutations,
            lambda_value_vdw=0.0,
            intst_nr=intst,
        )
        output_files.append(output_file_base)

        # turn off heavy atoms
        intst += 1

        output_file_base = i.write_state(
            mutation_conf=mutation_list["lj"],
            lambda_value_vdw=0.0,
            intst_nr=intst,
        )
        output_files.append(output_file_base)

        # generate terminal lj
        intst += 1

        output_file_base = i.write_state(
            mutation_conf=mutation_list["terminal-lj"],
            lambda_value_vdw=0.0,
            intst_nr=intst,
        )
        output_files.append(output_file_base)

        m = mutation_list["transform"]
        for lambda_value in np.linspace(0.25, 1, 3):
            intst += 1
            print(lambda_value)
            # turn off charges
            output_file_base = i.write_state(
                mutation_conf=m,
                common_core_transformation=1 - lambda_value,
                intst_nr=intst,
            )
            output_files.append(output_file_base)


def _mutate_toluene_to_methane_cc():
    conf = "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml"
    configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )

    output_files = []
    # mutate everything else before touching bonded terms
    intst = 0
    charges = mutation_list["charge"]
    intst += 1
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # Turn of hydrogens
    intst += 1
    hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
    output_file_base = i.write_state(
        mutation_conf=hydrogen_lj_mutations,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=mutation_list["lj"],
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # generate terminal lj
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=mutation_list["terminal-lj"],
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    m = mutation_list["transform"]
    for lambda_value in np.linspace(0.25, 1, 3):
        intst += 1
        print(lambda_value)
        # turn off charges
        output_file_base = i.write_state(
            mutation_conf=m,
            common_core_transformation=1 - lambda_value,
            intst_nr=intst,
        )
        output_files.append(output_file_base)
    return output_files, configuration


def _mutate_methane_to_methane_cc():
    conf = "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml"
    configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol2()
    i = IntermediateStateFactory(
        system=s2,
        configuration=configuration,
    )

    output_files = []
    # mutate everything else before touching bonded terms
    intst = 1
    charges = mutation_list["charge"]
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # generate terminal lj
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=mutation_list["terminal-lj"],
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    return output_files, configuration


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_toluene_to_methane_cc_solvation_free_energy_with_openMM():
    from transformato import FreeEnergyCalculator

    output_files, configuration = _mutate_toluene_to_methane_cc()

    for path in sorted(output_files):
        # because path is object not string
        print(f"Start sampling for: {path}")
        try:
            exe = subprocess.run(
                ["bash", f"{str(path)}/simulation.sh", str(path)],
                check=True,
                capture_output=True,
                text=True,
            )
        except TypeError:
            exe = subprocess.run(
                ["bash", f"{str(path)}/simulation.sh", str(path)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
        print(exe.stdout)
        print("Capture stderr")
        print(exe.stderr)

    f = FreeEnergyCalculator(configuration, "toluene")
    f.load_trajs(thinning=1)
    f.calculate_dG_to_common_core()
    ddG, dddG = f.end_state_free_energy_difference
    print(f"Free energy difference: {ddG}")
    print(f"Uncertanty: {dddG}")
    np.isclose(ddG, 8.9984, rtol=1e-2)
    f.show_summary()


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_methane_to_methane_cc_solvation_free_energy_with_openMM():
    from transformato import FreeEnergyCalculator

    output_files, configuration = _mutate_methane_to_methane_cc()

    for path in sorted(output_files):
        # because path is object not string
        print(f"Start sampling for: {path}")
        try:
            exe = subprocess.run(
                ["bash", f"{str(path)}/simulation.sh", str(path)],
                check=True,
                capture_output=True,
                text=True,
            )
        except TypeError:
            exe = subprocess.run(
                ["bash", f"{str(path)}/simulation.sh", str(path)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
        print(exe.stdout)
        print("Capture stderr")
        print(exe.stderr)

    f = FreeEnergyCalculator(configuration, "methane")
    f.load_trajs(thinning=1)
    f.calculate_dG_to_common_core()
    ddG, dddG = f.end_state_free_energy_difference
    print(f"Free energy difference: {ddG}")
    print(f"Uncertanty: {dddG}")
    np.isclose(ddG, 8.9984, rtol=1e-2)
    f.show_summary()


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_methane_to_methane_cc_solvation_free_energy_with_CHARMM():
    from transformato import FreeEnergyCalculator

    output_files, configuration = _mutate_methane_to_methane_cc()

    for path in sorted(output_files):
        # because path is object not string
        print(f"Start sampling for: {path}")
        try:
            exe = subprocess.run(
                ["bash", f"{str(path)}/simulation_charmm.sh", str(path)],
                check=True,
                capture_output=True,
                text=True,
            )
        except TypeError:
            exe = subprocess.run(
                ["bash", f"{str(path)}/simulation_charmm.sh", str(path)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
        print(exe.stdout)
        print("Capture stderr")
        print(exe.stderr)

    # f = FreeEnergyCalculator(configuration, "methane")
    # f.load_trajs(thinning=1)
    # f.calculate_dG_to_common_core()
    # ddG, dddG = f.end_state_free_energy_difference
    # print(f"Free energy difference: {ddG}")
    # print(f"Uncertanty: {dddG}")
    # np.isclose(ddG, 8.9984, rtol=1e-2)
    # f.show_summary()


# @pytest.mark.slowtest
# @pytest.mark.skipif(
#     os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
# )
# def test_run_example2_systems_solvation_free_energy():
#     from transformato import FreeEnergyCalculator

#     for conf in [
#         "transformato/tests/config/test-ethane-ethanol-solvation-free-energy.yaml"
#     ]:
#         configuration = load_config_yaml(
#             config=conf, input_dir="data/", output_dir="data/"
#         )

#         # load systems
#         s1 = SystemStructure(configuration, "structure1")
#         s2 = SystemStructure(configuration, "structure2")
#         a = ProposeMutationRoute(s1, s2)

#         # manually matching Oxygen (0) with Hydrogen (4)
#         a.add_idx_to_common_core_of_mol1(4)
#         a.add_idx_to_common_core_of_mol2(0)

#         # generate mutation route
#         mutation_list = a.generate_mutations_to_common_core_for_mol1(
#             nr_of_steps_for_electrostatic=5, nr_of_steps_for_cc_transformation=5
#         )

#         try:
#             # write intermediate states for systems
#             i = IntermediateStateFactory(
#                 system=s1, mutation_list=mutation_list, configuration=configuration
#             )
#             i.generate_intermediate_states()
#             paths = pathlib.Path(i.path).glob("**/*.sh")
#             for path in sorted(paths):
#                 run_dir = path.parent
#                 # because path is object not string
#                 print(f"Start sampling for: {path}")
#                 print(f"In directory: {run_dir}")
#                 try:
#                     exe = subprocess.run(
#                         ["bash", str(path), str(run_dir)],
#                         check=True,
#                         capture_output=True,
#                         text=True,
#                     )
#                 except TypeError:
#                     exe = subprocess.run(
#                         ["bash", str(path), str(run_dir)],
#                         check=True,
#                         stdout=subprocess.PIPE,
#                         stderr=subprocess.PIPE,
#                     )
#                 print(exe.stdout)
#                 print("Capture stderr")
#                 print(exe.stderr)

#             # generate mutation route
#             mutation_list = a.generate_mutations_to_common_core_for_mol2(
#                 nr_of_steps_for_electrostatic=5
#             )
#             # write intermediate states
#             i = IntermediateStateFactory(
#                 system=s2, mutation_list=mutation_list, configuration=configuration
#             )
#             i.generate_intermediate_states()

#             paths = pathlib.Path(i.path).glob("**/*.sh")
#             for path in sorted(paths):
#                 run_dir = path.parent
#                 # because path is object not string
#                 print(f"Start sampling for: {path}")
#                 print(f"In directory: {run_dir}")
#                 try:
#                     exe = subprocess.run(
#                         ["bash", str(path), str(run_dir)],
#                         check=True,
#                         capture_output=True,
#                         text=True,
#                     )
#                 except TypeError:
#                     exe = subprocess.run(
#                         ["bash", str(path), str(run_dir)],
#                         check=True,
#                         stdout=subprocess.PIPE,
#                         stderr=subprocess.PIPE,
#                     )
#                 print(exe.stdout)
#                 print("Capture stderr")
#                 print(exe.stderr)

#             f = FreeEnergyCalculator(configuration, "ethane")
#             f.load_trajs(thinning=2)
#             f.calculate_dG_to_common_core()
#             ddG, dddG = f.end_state_free_energy_difference
#             print(f"Free energy difference: {ddG}")
#             print(f"Uncertanty: {dddG}")
#             # assert(ddG == 10.0)
#             f.show_summary()

#             f = FreeEnergyCalculator(configuration, "ethanol")
#             f.load_trajs(thinning=2)
#             f.calculate_dG_to_common_core()
#             ddG, dddG = f.end_state_free_energy_difference
#             print(f"Free energy difference: {ddG}")
#             print(f"Uncertanty: {dddG}")
#             f.show_summary()
#         finally:
#             pass
#             # shutil.rmtree(pathlib.Path(i.path).parent)
