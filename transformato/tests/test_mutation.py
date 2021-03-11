"""
Unit and regression test for the transformato package.
"""

import copy
import os
import shutil
import logging

import numpy as np
import parmed as pm
import pytest

# Import package, test suite, and other packages as needed
import transformato

# read in specific topology with parameters
from parmed.charmm.parameters import CharmmParameterSet
from transformato import (
    IntermediateStateFactory,
    ProposeMutationRoute,
    SystemStructure,
    load_config_yaml,
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


def test_proposed_mutation_mcs():

    from rdkit.Chem import rdFMCS

    for conf in [
        "transformato/tests/config/test-2oj9-solvation-free-energy.yaml",
    ]:
        configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        a = ProposeMutationRoute(s1, s2)
        # find mcs
        a._find_mcs("m1", "m2")

        assert str(a.s1_tlc) == "UNK"
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


def setup_systems(conf):
    configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()
    i_s1 = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )
    i_s2 = IntermediateStateFactory(
        system=s2,
        configuration=configuration,
    )

    mutation_list_mol1 = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    mutation_list_mol2 = s1_to_s2.generate_mutations_to_common_core_for_mol2()

    return (configuration, mutation_list_mol1, mutation_list_mol2, i_s1, i_s2)


def test_mutation_list():

    for conf, system_name in zip(
        [
            "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
            "transformato/tests/config/test-neopentane-methane-solvation-free-energy.yaml",
            "transformato/tests/config/test-ethane-methanol-solvation-free-energy.yaml",
        ],
        ["toluene-methane", "neopentane-methane", "ethane-ethanol"],
    ):

        if system_name == "ethane-methanol":
            (
                configuration,
                mutation_list_mol1,
                mutation_list_mol2,
                i_s1,
                i_s2,
            ) = setup_systems(conf)
            assert set(mutation_list_mol1.keys()) == set(
                ["charge", "default-lj", "transform"]
            )
        if system_name == "toluene-methane":
            (
                configuration,
                mutation_list_mol1,
                mutation_list_mol2,
                i_s1,
                i_s2,
            ) = setup_systems(conf)
            assert set(mutation_list_mol1.keys()) == set(
                ["charge", "hydrogen-lj", "lj", "transform", "default-lj"]
            )

        if system_name == "neopentane-methane":
            configuration = load_config_yaml(
                config=conf, input_dir="data/", output_dir="."
            )
            s1 = SystemStructure(configuration, "structure1")
            s2 = SystemStructure(configuration, "structure2")

            s1_to_s2 = ProposeMutationRoute(s1, s2)
            s1_to_s2.propose_common_core()
            s1_to_s2.finish_common_core(
                connected_dummy_regions_cc1=[
                    {0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
                ]
            )
            mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
            assert set(mutation_list.keys()) == set(
                ["charge", "hydrogen-lj", "lj", "transform", "default-lj"]
            )


def test_endpoint_mutation():
    # test physical endpoint systems
    from ..utils import print_mutations
    from ..mutate import MutationDefinition

    for conf in [
        "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
        "transformato/tests/config/test-ethane-methanol-solvation-free-energy.yaml",
    ]:
        (
            configuration,
            mutation_list_mol1,
            mutation_list_mol2,
            i_s1,
            i_s2,
        ) = setup_systems(conf)

        output_file_base, _ = i_s1.write_state(mutation_conf=[], intst_nr=0)
        shutil.rmtree(output_file_base)

        output_file_base, _ = i_s2.write_state(mutation_conf=[], intst_nr=0)
        shutil.rmtree(output_file_base)


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
        (
            configuration,
            mutation_list_mol1,
            mutation_list_mol2,
            i_s1,
            i_s2,
        ) = setup_systems(conf)
        for i in [i_s1, i_s2]:
            output_file_base, _ = i.write_state(mutation_conf=[], intst_nr=0)

            # original psfs without charge change
            original_psf = {}
            for env in i.system.envs:
                original_psf[env] = copy.deepcopy(i.system.psfs[env])

            for env in i.system.envs:
                offset = i.system.offset[env]

                mutated_psf = generate_psf(output_file_base, env)
                for atom in i.system.mol.GetAtoms():
                    idx = atom.GetIdx()
                    assert np.isclose(
                        original_psf[env].atoms[idx + offset].charge,
                        mutated_psf.atoms[idx + offset].charge,
                        rtol=1e-3,
                    )
            shutil.rmtree(output_file_base)


def test_setup_dual_junction_system():

    conf = "transformato/tests/config/test-2oj9-solvation-free-energy.yaml"
    configuration, mutation_list_mol1, mutation_list_mol2, i_s1, i_s2 = setup_systems(
        conf
    )
    # write out endpoint
    output_files = []
    output_file_base, intst = i_s1.write_state(mutation_conf=[], intst_nr=0)
    output_files.append(output_file_base)

    charges = mutation_list_mol1["charge"]
    # start with charges
    # turn off charges
    output_file_base, intst = i_s1.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # Turn off hydrogens
    hydrogen_lj_mutations = mutation_list_mol1["hydrogen-lj"]
    output_file_base, intst = i_s1.write_state(
        mutation_conf=hydrogen_lj_mutations,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)
    shutil.rmtree(output_files)


def test_charge_mutation_test_system2():

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
                charges = mutation_list["charge"]
                output_file_base, intst_nr = i.write_state(
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
        terminal_lj_mutations = mutation_list["default-lj"]
        output_file_base, intst_nr = i.write_state(
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

        shutil.rmtree(output_file_base)


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
            # lj hydrogens scaling
            hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
            output_file_base, _ = i.write_state(
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

            shutil.rmtree(output_file_base)


@pytest.mark.slowtest
def test_bonded_mutation():

    for conf in [
        "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
    ]:

        original_psf = {}
        output_files = []
        for env in s1.envs:
            original_psf[env] = copy.deepcopy(s1.psfs[env])

        # mutate everything else before touching bonded terms
        charges = mutation_list["charge"]
        # turn off charges
        output_file_base, intst_nr = i.write_state(
            mutation_conf=charges,
            lambda_value_electrostatic=0.0,
            intst_nr=0,
        )
        output_files.append(output_file_base)

        # Turn of hydrogens
        hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
        output_file_base, intst_nr = i.write_state(
            mutation_conf=hydrogen_lj_mutations,
            lambda_value_vdw=0.0,
            intst_nr=intst_nr,
        )
        output_files.append(output_file_base)

        # turn off heavy atoms
        output_file_base, intst_nr = i.write_state(
            mutation_conf=mutation_list["lj"],
            lambda_value_vdw=0.0,
            intst_nr=intst_nr,
        )
        output_files.append(output_file_base)

        # generate terminal lj
        output_file_base, intst_nr = i.write_state(
            mutation_conf=mutation_list["default-lj"],
            lambda_value_vdw=0.0,
            intst_nr=intst_nr,
        )
        output_files.append(output_file_base)

        m = mutation_list["transform"]
        for lambda_value in np.linspace(0.75, 0, 3):
            print(lambda_value)
            # turn off charges
            output_file_base, intst_nr = i.write_state(
                mutation_conf=m,
                common_core_transformation=lambda_value,
                intst_nr=intst_nr,
            )
            output_files.append(output_file_base)
        shutil.rmtree(output_file_base)


@pytest.mark.slowtest
def test_bonded_mutation_atoms(caplog):
    caplog.set_level(logging.CRITICAL)

    from .test_mutation import setup_2OJ9_tautomer_pair

    (output_files_t1, output_files_t2), _, p = setup_2OJ9_tautomer_pair()
    psf_at_endstate_t1 = generate_psf(output_files_t1[0], "vacuum")
    psf_at_t1_cc = generate_psf(output_files_t1[-1], "vacuum")
    psf_at_t2_cc = generate_psf(output_files_t2[-1], "vacuum")
    psf_at_endstate_t2 = generate_psf(output_files_t2[0], "vacuum")
    print("@@@@@@@@@@@@@@@@@@")

    cc1_idx, cc2_idx = p.get_common_core_idx_mol1(), p.get_common_core_idx_mol2()
    for atom_t1_idx, atom_t2_idx in zip(cc1_idx, cc2_idx):
        atom_t1 = psf_at_t1_cc.atoms[atom_t1_idx]
        atom_t2 = psf_at_t2_cc.atoms[atom_t2_idx]
        assert atom_t1.charge == atom_t2.charge
        assert atom_t1.sigma == atom_t2.sigma

    for atom_t1_idx, atom_t2_idx in zip(cc1_idx, cc2_idx):
        atom_t1 = psf_at_endstate_t1.atoms[atom_t1_idx]
        atom_t2 = psf_at_endstate_t2.atoms[atom_t2_idx]
        try:
            assert atom_t1.charge == atom_t2.charge
            assert atom_t1.sigma == atom_t2.sigma
        except AssertionError:
            pass


@pytest.mark.slowtest
def test_bonded_mutation_bonds(caplog):
    caplog.set_level(logging.CRITICAL)

    from .test_mutation import setup_2OJ9_tautomer_pair

    (output_files_t1, output_files_t2), _, p = setup_2OJ9_tautomer_pair()
    psf_at_endstate_t1 = generate_psf(output_files_t1[0], "vacuum")
    psf_at_t1_cc = generate_psf(output_files_t1[-1], "vacuum")
    psf_at_t2_cc = generate_psf(output_files_t2[-1], "vacuum")
    psf_at_endstate_t2 = generate_psf(output_files_t2[0], "vacuum")
    print("@@@@@@@@@@@@@@@@@@")

    cc1_idx, cc2_idx = p.get_common_core_idx_mol1(), p.get_common_core_idx_mol2()
    # compare cc enstates
    for bond_t1 in psf_at_t1_cc.bonds:
        atom1_t1_idx = bond_t1.atom1.idx
        atom2_t1_idx = bond_t1.atom2.idx
        if atom1_t1_idx not in cc1_idx or atom2_t1_idx not in cc1_idx:
            continue

        # get index in common core
        idx1 = cc1_idx.index(atom1_t1_idx)
        idx2 = cc1_idx.index(atom2_t1_idx)

        atom1_t2 = psf_at_t2_cc[cc2_idx[idx1]]
        atom2_t2 = psf_at_t2_cc[cc2_idx[idx2]]
        bond_t2 = None
        for bond in psf_at_t2_cc.bonds:
            if atom1_t2 in bond and atom2_t2 in bond:
                bond_t2 = bond

        assert bond_t1.type.k == bond_t2.type.k
        assert bond_t1.type.req == bond_t2.type.req

    # compare physical endstates
    for bond_t1 in psf_at_endstate_t1.bonds:
        atom1_t1_idx = bond_t1.atom1.idx
        atom2_t1_idx = bond_t1.atom2.idx
        if atom1_t1_idx not in cc1_idx or atom2_t1_idx not in cc1_idx:
            continue

        # get index in common core
        idx1 = cc1_idx.index(atom1_t1_idx)
        idx2 = cc1_idx.index(atom2_t1_idx)

        atom1_t2 = psf_at_endstate_t2[cc2_idx[idx1]]
        atom2_t2 = psf_at_endstate_t2[cc2_idx[idx2]]
        bond_t2 = None
        for bond in psf_at_endstate_t2.bonds:
            if atom1_t2 in bond and atom2_t2 in bond:
                bond_t2 = bond

        try:
            assert bond_t1.type.k == bond_t2.type.k
            assert bond_t1.type.req == bond_t2.type.req
        except AssertionError:
            pass


@pytest.mark.slowtest
def test_bonded_mutation_angles(caplog):
    caplog.set_level(logging.CRITICAL)

    from .test_mutation import setup_2OJ9_tautomer_pair

    (output_files_t1, output_files_t2), _, p = setup_2OJ9_tautomer_pair()
    psf_at_endstate_t1 = generate_psf(output_files_t1[0], "vacuum")
    psf_at_t1_cc = generate_psf(output_files_t1[-1], "vacuum")
    psf_at_t2_cc = generate_psf(output_files_t2[-1], "vacuum")
    psf_at_endstate_t2 = generate_psf(output_files_t2[0], "vacuum")
    print("@@@@@@@@@@@@@@@@@@")

    cc1_idx, cc2_idx = p.get_common_core_idx_mol1(), p.get_common_core_idx_mol2()
    for angle_t1 in psf_at_t1_cc.angles:
        atom1_t1_idx = angle_t1.atom1.idx
        atom2_t1_idx = angle_t1.atom2.idx
        atom3_t1_idx = angle_t1.atom3.idx
        if (
            atom1_t1_idx not in cc1_idx
            or atom2_t1_idx not in cc1_idx
            or atom3_t1_idx not in cc1_idx
        ):
            continue

        # get index in common core
        idx1 = cc1_idx.index(atom1_t1_idx)
        idx2 = cc1_idx.index(atom2_t1_idx)
        idx3 = cc1_idx.index(atom3_t1_idx)

        atom1_t2 = psf_at_t2_cc[cc2_idx[idx1]]
        atom2_t2 = psf_at_t2_cc[cc2_idx[idx2]]
        atom3_t2 = psf_at_t2_cc[cc2_idx[idx3]]
        angle_t2 = None
        for angle in psf_at_t2_cc.angles:
            if atom1_t2 in angle and atom2_t2 in angle and atom3_t2 in angle:
                angle_t2 = angle

        assert angle_t1.type.k == angle_t2.type.k
        assert angle_t1.type.theteq == angle_t2.type.theteq


@pytest.mark.slowtest
def test_bonded_mutation_dihedrals(caplog):
    caplog.set_level(logging.CRITICAL)

    from .test_mutation import setup_2OJ9_tautomer_pair

    (output_files_t1, output_files_t2), _, p = setup_2OJ9_tautomer_pair()
    psf_at_endstate_t1 = generate_psf(output_files_t1[0], "vacuum")
    psf_at_t1_cc = generate_psf(output_files_t1[-1], "vacuum")
    print(output_files_t1[-1])
    psf_at_t2_cc = generate_psf(output_files_t2[-1], "vacuum")
    psf_at_endstate_t2 = generate_psf(output_files_t2[0], "vacuum")
    print("@@@@@@@@@@@@@@@@@@")

    assert len(psf_at_t1_cc.dihedrals) == len(psf_at_t2_cc.dihedrals)
    # compare the common core parameters at the cc state
    cc1_idx, cc2_idx = p.get_common_core_idx_mol1(), p.get_common_core_idx_mol2()
    for dihedral_t1 in psf_at_t1_cc.dihedrals:
        atom1_t1_idx = dihedral_t1.atom1.idx
        atom2_t1_idx = dihedral_t1.atom2.idx
        atom3_t1_idx = dihedral_t1.atom3.idx
        atom4_t1_idx = dihedral_t1.atom4.idx
        if (
            atom1_t1_idx not in cc1_idx
            or atom2_t1_idx not in cc1_idx
            or atom3_t1_idx not in cc1_idx
            or atom4_t1_idx not in cc1_idx
        ):
            print("Not present:")
            print(dihedral_t1)

            continue

        # get index in common core
        idx1 = cc1_idx.index(atom1_t1_idx)
        idx2 = cc1_idx.index(atom2_t1_idx)
        idx3 = cc1_idx.index(atom3_t1_idx)
        idx4 = cc1_idx.index(atom4_t1_idx)

        atom1_t2 = psf_at_t2_cc[cc2_idx[idx1]]
        atom2_t2 = psf_at_t2_cc[cc2_idx[idx2]]
        atom3_t2 = psf_at_t2_cc[cc2_idx[idx3]]
        atom4_t2 = psf_at_t2_cc[cc2_idx[idx4]]

        dihedral_t2 = None
        for dihedral in psf_at_t2_cc.dihedrals:
            if (
                atom1_t2 in dihedral
                and atom2_t2 in dihedral
                and atom3_t2 in dihedral
                and atom4_t2 in dihedral
            ):
                dihedral_t2 = dihedral
                break

        assert dihedral_t2 != None
        print(dihedral_t1)
        print(dihedral_t2)
        for t1, t2 in zip(dihedral_t1.type, dihedral_t2.type):
            assert t1.phi_k == t2.phi_k
            assert t1.phase == t2.phase
    assert False
    # # compare the endstates
    # cc1_idx, cc2_idx = p.get_common_core_idx_mol1(), p.get_common_core_idx_mol2()
    # for dihedral_t1 in psf_at_endstate_t1.dihedrals:
    #     atom1_t1_idx = dihedral_t1.atom1.idx
    #     atom2_t1_idx = dihedral_t1.atom2.idx
    #     atom3_t1_idx = dihedral_t1.atom3.idx
    #     atom4_t1_idx = dihedral_t1.atom4.idx
    #     if (
    #         atom1_t1_idx not in cc1_idx
    #         or atom2_t1_idx not in cc1_idx
    #         or atom3_t1_idx not in cc1_idx
    #         or atom4_t1_idx not in cc1_idx
    #     ):
    #         continue

    #     # get index in common core
    #     idx1 = cc1_idx.index(atom1_t1_idx)
    #     idx2 = cc1_idx.index(atom2_t1_idx)
    #     idx3 = cc1_idx.index(atom3_t1_idx)
    #     idx4 = cc1_idx.index(atom4_t1_idx)

    #     atom1_t2 = psf_at_t2_cc[cc2_idx[idx1]]
    #     atom2_t2 = psf_at_t2_cc[cc2_idx[idx2]]
    #     atom3_t2 = psf_at_t2_cc[cc2_idx[idx3]]
    #     atom4_t2 = psf_at_t2_cc[cc2_idx[idx4]]

    #     dihedral_t2 = None
    #     for dihedral in psf_at_endstate_t2.dihedrals:
    #         if (
    #             atom1_t2 in dihedral
    #             and atom2_t2 in dihedral
    #             and atom3_t2 in dihedral
    #             and atom4_t2 in dihedral
    #         ):
    #             dihedral_t2 = dihedral
    #             break

    #     assert dihedral_t2 != None
    #     print(dihedral_t1)
    #     print(dihedral_t2)
    #     for t1, t2 in zip(dihedral_t1.type, dihedral_t2.type):
    #         assert t1.phi_k == t2.phi_k
    #         assert t1.phase == t2.phase


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
        intst_nr = 1
        for lambda_vdw in [1.0, 0.0]:
            print(f"Lambda: {lambda_vdw}")
            output_files = []
            all_atoms_for_which_lj_turned_off = []
            print(mutation_list.keys())
            # turn off hydrogen lj
            hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
            print(
                f"Turn off lj for hydrogen atoms : {[e.atoms_to_be_mutated for e in hydrogen_lj_mutations]}"
            )
            output_file_base, intst_nr = i.write_state(
                mutation_conf=hydrogen_lj_mutations,
                lambda_value_vdw=lambda_vdw,
                intst_nr=intst_nr,
            )
            output_files.append(output_file_base)
            for mutation in hydrogen_lj_mutations:
                all_atoms_for_which_lj_turned_off.extend(mutation.vdw_atom_idx)

            # turn off heavy atom lj
            for mutation in mutation_list["lj"]:
                print(f"Turn off lj for heavy atom : {mutation.atoms_to_be_mutated}")

                output_file_base, intst_nr = i.write_state(
                    mutation_conf=[mutation],
                    lambda_value_vdw=lambda_vdw,
                    intst_nr=intst_nr,
                )
                all_atoms_for_which_lj_turned_off.extend(mutation.vdw_atom_idx)
                output_files.append(output_file_base)

            # change to default lj
            terminal_lj_mutations = mutation_list["default-lj"]
            terminal_idx = []
            for mutation in terminal_lj_mutations:
                all_atoms_for_which_lj_turned_off.extend(mutation.vdw_atom_idx)
                terminal_idx.extend(mutation.vdw_atom_idx)
            print(
                f"Turn off lj for terminal atom : {all_atoms_for_which_lj_turned_off}"
            )

            output_file_base, intst_nr = i.write_state(
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

        shutil.rmtree(f"{system_name}-solvation-free-energy")


def setup_2OJ9_tautomer_pair():
    from ..mutate import mutate_pure_tautomers
    from ..constants import check_platform

    conf_path = (
        "transformato/tests/config/test-2oj9-tautomer-pair-solvation-free-energy.yaml"
    )
    configuration = load_config_yaml(
        config=conf_path, input_dir="data/", output_dir="."
    )
    check_platform(configuration)

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()
    return (
        mutate_pure_tautomers(s1_to_s2, s1, s2, configuration),
        configuration,
        s1_to_s2,
    )


def setup_acetylacetone_tautomer_pair():
    from ..mutate import mutate_pure_tautomers
    from ..constants import check_platform

    conf_path = "transformato/tests/config/test-acetylaceton-tautomer-solvation-free-energy.yaml"
    configuration = load_config_yaml(
        config=conf_path, input_dir="data/", output_dir="."
    )
    check_platform(configuration)

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()
    return mutate_pure_tautomers(s1_to_s2, s1, s2, configuration), configuration


def test_acetylacetone_tautomer_pair(caplog):
    caplog.set_level(logging.DEBUG)
    setup_acetylacetone_tautomer_pair()
    shutil.rmtree("acetylacetone-keto-acetylacetone-enol-solvation-free-energy")


def test_2OJ9_tautomer_pair(caplog):
    caplog.set_level(logging.DEBUG)
    setup_2OJ9_tautomer_pair()
    shutil.rmtree("2OJ9-original-2OJ9-tautomer-solvation-free-energy")


def test_full_mutation_system1(caplog):
    caplog.set_level(logging.WARNING)

    for conf, system_name in zip(
        [
            "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
            "transformato/tests/config/test-neopentane-methane-solvation-free-energy.yaml",
            "transformato/tests/config/test-ethane-methanol-solvation-free-energy.yaml",
        ],
        ["toluene-methane", "neopentane-methane", "ethane-methanol"],
    ):
        print(system_name)
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
        intst = 1
        charges = mutation_list["charge"]
        for lambda_value in np.linspace(0, 1, 5):
            # turn off charges
            output_file_base, intst = i.write_state(
                mutation_conf=charges,
                lambda_value_electrostatic=1 - lambda_value,
                intst_nr=intst,
            )

        original_psf = {}
        for env in s1.envs:
            original_psf[env] = copy.deepcopy(s1.psfs[env])

        #
        lambda_vdw = 0.0

        all_atoms_for_which_lj_turned_off = []
        print(mutation_list.keys())
        # Turn of hydrogens
        terminal_lj_mutations = mutation_list["default-lj"]
        _, intst = i.write_state(
            mutation_conf=terminal_lj_mutations,
            lambda_value_vdw=lambda_vdw,
            intst_nr=intst,
        )
        for mutation in terminal_lj_mutations:
            all_atoms_for_which_lj_turned_off.extend(mutation.vdw_atom_idx)

        shutil.rmtree(f"{system_name}-solvation-free-energy")


def test_full_mutation_system2():

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

        intst = 1
        charges = mutation_list["charge"]
        for lambda_value in np.linspace(0, 1, 5):
            # turn off charges
            output_file_base, intst = i.write_state(
                mutation_conf=charges,
                lambda_value_electrostatic=1 - lambda_value,
                intst_nr=intst,
            )

        original_psf = {}
        for env in s1.envs:
            original_psf[env] = copy.deepcopy(s1.psfs[env])

        #
        lambda_vdw = 0.0
        print(f"Lambda: {lambda_vdw}")
        all_atoms_for_which_lj_turned_off = []
        print(mutation_list.keys())
        # Turn of hydrogens
        hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
        output_file_base, intst = i.write_state(
            mutation_conf=hydrogen_lj_mutations,
            lambda_value_vdw=lambda_vdw,
            intst_nr=intst,
        )

        for mutation in hydrogen_lj_mutations:
            all_atoms_for_which_lj_turned_off.extend(mutation.vdw_atom_idx)

        # turn off heavy atoms
        for mutation in mutation_list["lj"]:

            output_file_base, intst = i.write_state(
                mutation_conf=[mutation],
                lambda_value_vdw=lambda_vdw,
                intst_nr=intst,
            )
            all_atoms_for_which_lj_turned_off.extend(mutation.vdw_atom_idx)

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

        shutil.rmtree(f"{system_name}-solvation-free-energy")
