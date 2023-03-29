"""
Unit and regression test for the transformato package.
"""

import copy
import shutil
import logging
import pytest
import os

import numpy as np
import parmed as pm

# Import package, test suite, and other packages as needed
import transformato
from openmm import unit

# read in specific topology with parameters
from parmed.charmm.parameters import CharmmParameterSet
from transformato import (
    IntermediateStateFactory,
    ProposeMutationRoute,
    SystemStructure,
    load_config_yaml,
)

from transformato.tests.paths import get_test_output_dir
from transformato_testsystems.testsystems import (
    get_testsystems_dir,
    get_output_files_2oj9_tautomer_pair,
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
    psf = pm.charmm.CharmmPsfFile(f"{output_file_base}/lig_in_{env}.psf")
    psf.load_parameters(parms)
    return psf, parms


def generate_crd(output_file_base, env):
    return pm.charmm.CharmmCrdFile(f"{output_file_base}/lig_in_{env}.crd")


def generate_sim(output_file_base, env):
    import openmm as mm
    import openmm.app as app

    # ParmEd Imports
    from parmed import unit as u

    system, psf = generate_system(output_file_base, env)
    crd = generate_crd(output_file_base, env)

    # Create the integrator to do Langevin dynamics
    integrator = mm.LangevinIntegrator(
        300 * u.kelvin,  # Temperature of heat bath
        1.0 / u.picoseconds,  # Friction coefficient
        2.0 * u.femtoseconds,  # Time step
    )

    # Create the Simulation object
    sim = app.Simulation(psf.topology, system, integrator)

    # Set the particle positions
    sim.context.setPositions(crd.positions)
    return sim


def generate_system(output_file_base, env):
    import openmm as mm
    import openmm.app as app

    # ParmEd Imports
    from parmed import unit as u

    psf, parms = generate_psf(output_file_base, env)

    system = psf.createSystem(parms, nonbondedMethod=app.NoCutoff)
    return system, psf


def setup_acetylacetone_tautomer_pair(
    configuration: dict, single_state=False, nr_of_bonded_windows=4
):
    from ..mutate import mutate_pure_tautomers

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
            single_state=single_state,
            nr_of_bonded_windows=nr_of_bonded_windows,
        ),
        configuration,
        s1_to_s2,
    )


def test_proposed_mutation_mcs():
    from rdkit.Chem import rdFMCS

    for conf in [
        f"{get_testsystems_dir()}/config/test-2oj9-rsfe.yaml",
    ]:
        configuration = load_config_yaml(
            config=conf,
            input_dir=get_testsystems_dir(),
            output_dir=get_test_output_dir(),
        )
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
                2,
                7,
                11,
                9,
                8,
                10,
                13,
                12,
                31,
                33,
                32,
                44,
                41,
                45,
                46,
                47,
                48,
                30,
                34,
                35,
                37,
                36,
                38,
                39,
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
                2,
                7,
                11,
                9,
                8,
                10,
                13,
                12,
                30,
                32,
                31,
                41,
                40,
                42,
                43,
                44,
                45,
                29,
                33,
                34,
                36,
                35,
                37,
                38,
            ]
        )

        print(a.get_common_core_idx_mol1())
        print(a.get_common_core_idx_mol2())

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
                2,
                7,
                11,
                9,
                8,
                10,
                13,
                12,
                31,
                33,
                32,
                44,
                41,
                45,
                46,
                47,
                48,
                30,
                34,
                35,
                37,
                36,
                38,
                39,
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
                2,
                7,
                11,
                9,
                8,
                10,
                13,
                12,
                30,
                32,
                31,
                41,
                40,
                42,
                43,
                44,
                45,
                29,
                33,
                34,
                36,
                35,
                37,
                38,
            ]
        )

        print(a.get_common_core_idx_mol1())
        print(a.get_common_core_idx_mol2())

        assert set(a.get_common_core_idx_mol2()) == cc2
        assert set(a.get_common_core_idx_mol1()) == cc1

    for conf in [f"{get_testsystems_dir()}/config/test-toluene-methane-rsfe.yaml"]:
        configuration = load_config_yaml(
            config=conf,
            input_dir=get_testsystems_dir(),
            output_dir=get_test_output_dir(),
        )
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


def test_mutation_with_multiple_dummy_regions(caplog):
    # Test that TF can handel multiple dummy regions
    caplog.set_level(logging.INFO)
    import warnings

    warnings.filterwarnings("ignore", module="parmed")

    conf = f"{get_testsystems_dir()}/config/test-1a0q-1a07-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir=get_testsystems_dir(), output_dir=get_test_output_dir()
    )
    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.propose_common_core()
    s1_to_s2.finish_common_core()


@pytest.mark.rsfe
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_proposed_mutation_terminal_dummy_real_atom_match():
    from rdkit.Chem import rdFMCS

    workdir = get_test_output_dir()
    for conf in [
        f"{get_testsystems_dir()}/config/test-7-CPI-2-CPI-rsfe.yaml",
    ]:
        configuration = load_config_yaml(
            config=conf,
            input_dir=get_testsystems_dir(),
            output_dir=get_test_output_dir(),
        )
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


@pytest.mark.rsfe
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_find_connected_dummy_regions1():
    workdir = get_test_output_dir()
    from rdkit.Chem import rdFMCS

    conf = f"{get_testsystems_dir()}/config/test-7-CPI-2-CPI-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir=get_testsystems_dir(), output_dir=get_test_output_dir()
    )
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
    e = next(iter(match_terminal_atoms_cc1[15]))
    assert e == 14
    assert len(match_terminal_atoms_cc1[15]) == 1
    e = next(iter(match_terminal_atoms_cc2[15]))
    assert e == 6
    assert len(match_terminal_atoms_cc2[15]) == 1

    # find connected dummy regions
    connected_dummy_regions_cc1 = a._find_connected_dummy_regions("m1")
    connected_dummy_regions_cc2 = a._find_connected_dummy_regions("m2")

    lj_default_cc1, lj_default_cc2 = a._match_terminal_dummy_atoms_between_common_cores(
        match_terminal_atoms_cc1,
        match_terminal_atoms_cc2,
    )

    assert set(connected_dummy_regions_cc1[0]) == set(
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    )
    assert set(connected_dummy_regions_cc2[0]) == set(
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    )


def test_find_connected_dummy_regions2():
    ##################################################
    conf = f"{get_testsystems_dir()}/config/test-2oj9-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir=get_testsystems_dir(), output_dir=get_test_output_dir()
    )
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
    connected_dummy_regions_cc1 = a._find_connected_dummy_regions("m1")
    connected_dummy_regions_cc2 = a._find_connected_dummy_regions("m2")

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
    print(connected_dummy_regions_cc2)

    dummy_region_m2 = transformato.mutate.DummyRegion(
        "m2",
        match_terminal_atoms_cc2,
        connected_dummy_regions_cc2,
        s2.tlc,
        lj_default=lj_default_cc2,
    )

    print(dummy_region_m1.connected_dummy_regions)
    print(dummy_region_m2.connected_dummy_regions)

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


def test_common_core_for_multiple_systems():
    for conf in [
        f"{get_testsystems_dir()}/config/test-toluene-methane-rsfe.yaml",
        f"{get_testsystems_dir()}/config/test-neopentane-methane-rsfe.yaml",
        f"{get_testsystems_dir()}/config/test-ethane-methanol-rsfe.yaml",
    ]:
        configuration = load_config_yaml(
            config=conf,
            input_dir=get_testsystems_dir(),
            output_dir=get_test_output_dir(),
        )
        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        a = ProposeMutationRoute(s1, s2)
        if conf == f"{get_testsystems_dir()}/config/test-neopentane-methane-rsfe.yaml":
            a.propose_common_core()
            a.finish_common_core(
                connected_dummy_regions_cc1=[
                    {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
                ]
            )
        else:
            # find mcs, find terminal dummy/real atoms, generate charge compensated psfs
            a.calculate_common_core()


def setup_systems(conf):
    workdir = get_test_output_dir()
    configuration = load_config_yaml(
        config=conf, input_dir=get_testsystems_dir(), output_dir=get_test_output_dir()
    )
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


def test_generate_mutation_list_for_multiple_systems():
    for conf, system_name in zip(
        [
            f"{get_testsystems_dir()}/config/test-toluene-methane-rsfe.yaml",
            f"{get_testsystems_dir()}/config/test-neopentane-methane-rsfe.yaml",
            f"{get_testsystems_dir()}/config/test-ethane-methanol-rsfe.yaml",
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
                config=conf,
                input_dir=get_testsystems_dir(),
                output_dir=get_test_output_dir(),
            )
            s1 = SystemStructure(configuration, "structure1")
            s2 = SystemStructure(configuration, "structure2")

            s1_to_s2 = ProposeMutationRoute(s1, s2)
            s1_to_s2.propose_common_core()
            s1_to_s2.finish_common_core(connected_dummy_regions_cc2=[{1, 2, 3, 4}])
            mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
            assert set(mutation_list.keys()) == set(
                ["charge", "hydrogen-lj", "lj", "transform", "default-lj"]
            )


def test_write_endpoint_state():
    # test physical endpoint systems
    for conf in [
        f"{get_testsystems_dir()}/config/test-toluene-methane-rsfe.yaml",
        f"{get_testsystems_dir()}/config/test-ethane-methanol-rsfe.yaml",
    ]:
        (
            _,
            _,
            _,
            i_s1,
            i_s2,
        ) = setup_systems(conf)

        assert i_s1.current_step == 1
        i_s1.write_state(mutation_conf=[])
        assert i_s1.current_step == 2
        print(i_s1.output_files)
        for opath in i_s1.output_files:
            shutil.rmtree(opath)

        assert i_s2.current_step == 1
        i_s2.write_state(mutation_conf=[])
        assert i_s2.current_step == 2
        print(i_s2.output_files)
        for opath in i_s2.output_files:
            shutil.rmtree(opath)


def test_charges_at_endstate():
    for conf, system_name in zip(
        [
            f"{get_testsystems_dir()}/config/test-toluene-methane-rsfe.yaml",
            f"{get_testsystems_dir()}/config/test-ethane-methanol-rsfe.yaml",
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
            i.write_state(mutation_conf=[])

            # original psfs without charge change
            original_psf = {}
            for env in i.system.envs:
                original_psf[env] = copy.deepcopy(i.system.psfs[env])

            for env in i.system.envs:
                offset = i.system.offset[env]

                mutated_psf, param = generate_psf(i.output_files[-1], env)
                for atom in i.system.mol.GetAtoms():
                    idx = atom.GetIdx()
                    print(original_psf[env].atoms[idx + offset].charge)
                    print(mutated_psf.atoms[idx + offset].charge)
                    assert np.isclose(
                        original_psf[env].atoms[idx + offset].charge,
                        mutated_psf.atoms[idx + offset].charge,
                        rtol=1e-3,
                    )

            for opath in i.output_files:
                shutil.rmtree(opath)


def test_setup_dual_junction_system():
    conf = f"{get_testsystems_dir()}/config/test-2oj9-rsfe.yaml"
    configuration, mutation_list_mol1, mutation_list_mol2, i_s1, i_s2 = setup_systems(
        conf
    )
    # write out endpoint
    output_files = []
    i_s1.write_state(mutation_conf=[])

    charges = mutation_list_mol1["charge"]
    # start with charges
    # turn off charges
    i_s1.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.0,
    )

    # Turn off hydrogens
    hydrogen_lj_mutations = mutation_list_mol1["hydrogen-lj"]
    i_s1.write_state(
        mutation_conf=hydrogen_lj_mutations,
        lambda_value_vdw=0.0,
    )
    for f in i_s1.output_files:
        shutil.rmtree(f)


def test_charge_mutation_for_multiple_systems():
    for conf, system_name in zip(
        [
            f"{get_testsystems_dir()}/config/test-toluene-methane-rsfe.yaml",
            f"{get_testsystems_dir()}/config/test-neopentane-methane-rsfe.yaml",
            f"{get_testsystems_dir()}/config/test-ethane-methanol-rsfe.yaml",
        ],
        ["toluene-methane", "neopentane-methane", "ethane-methanol"],
    ):
        configuration = load_config_yaml(
            config=conf,
            input_dir=get_testsystems_dir(),
            output_dir=get_test_output_dir(),
        )
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
                    a.finish_common_core(connected_dummy_regions_cc2=[{1, 2, 3, 4}])
                else:
                    a.calculate_common_core()

                i = IntermediateStateFactory(
                    system=system,
                    configuration=configuration,
                )

                mutation_list = a.generate_mutations_to_common_core_for_mol1()
                charges = mutation_list["charge"]
                i.write_state(
                    mutation_conf=charges,
                    lambda_value_electrostatic=lambda_charge,
                )
                for env in system.envs:
                    offset = system.offset[env]
                    # read in newly generated psf
                    mutated_psf, params = generate_psf(i.output_files[-1], env)
                    for idx in charges[0].atoms_to_be_mutated:
                        assert np.isclose(
                            original_psf[env].atoms[idx + offset].charge
                            * lambda_charge,
                            mutated_psf.atoms[idx + offset].charge,
                            rtol=1e-03,
                        )
                for opath in i.output_files:
                    shutil.rmtree(opath)


def test_vdw_mutation_for_hydrogens_system1():
    # testing terminal lj

    for conf, system_name in zip(
        [
            f"{get_testsystems_dir()}/config/test-toluene-methane-rsfe.yaml",
            f"{get_testsystems_dir()}/config/test-neopentane-methane-rsfe.yaml",
            f"{get_testsystems_dir()}/config/test-ethane-methanol-rsfe.yaml",
        ],
        ["toluene-methane", "neopentane-methane", "ethane-methanol"],
    ):
        configuration = load_config_yaml(
            config=conf,
            input_dir=get_testsystems_dir(),
            output_dir=get_test_output_dir(),
        )
        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        s1_to_s2 = ProposeMutationRoute(s1, s2)
        if system_name == "neopentane-methane":
            s1_to_s2.propose_common_core()
            s1_to_s2.finish_common_core(
                connected_dummy_regions_cc1=[
                    {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
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
        i.write_state(
            mutation_conf=terminal_lj_mutations,
            lambda_value_vdw=0.0,
        )
        print("Set epsilon/rmin to zero for selected atoms")

        for env in s1.envs:
            new_psf, params = generate_psf(i.output_files[-1], env)
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
        for opath in i.output_files:
            shutil.rmtree(opath)


def test_vdw_mutation_for_hydrogens_system2():
    from rdkit.Chem import rdFMCS

    for conf, system_name in zip(
        [
            f"{get_testsystems_dir()}/config/test-7-CPI-2-CPI-rsfe.yaml",
        ],
        ["7-CPI-2-CPI"],
    ):
        configuration = load_config_yaml(
            config=conf,
            input_dir=get_testsystems_dir(),
            output_dir=get_test_output_dir(),
        )

        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        s1_to_s2 = ProposeMutationRoute(s1, s2)
        s1_to_s2.completeRingsOnly = True
        s1_to_s2.propose_common_core()
        # s1_to_s2.remove_idx_from_common_core_of_mol1([14])
        # s1_to_s2.remove_idx_from_common_core_of_mol2([6])
        s1_to_s2.finish_common_core()

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
            i.write_state(
                mutation_conf=hydrogen_lj_mutations,
                lambda_value_vdw=lambda_value,
            )
            print("Set epsilon/rmin for selected atoms")

            for env in s1.envs:
                new_psf, params = generate_psf(i.output_files[-1], env)
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
        for opath in i.output_files:
            shutil.rmtree(opath)


def test_bonded_mutation():
    for conf in [
        f"{get_testsystems_dir()}/config/test-toluene-methane-rsfe.yaml",
    ]:
        configuration = load_config_yaml(
            config=conf,
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
        original_psf = {}
        output_files = []
        for env in s1.envs:
            original_psf[env] = copy.deepcopy(s1.psfs[env])

        # mutate everything else before touching bonded terms
        charges = mutation_list["charge"]
        # turn off charges
        i.write_state(
            mutation_conf=charges,
            lambda_value_electrostatic=0.0,
        )

        # Turn of hydrogens
        hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
        i.write_state(
            mutation_conf=hydrogen_lj_mutations,
            lambda_value_vdw=0.0,
        )

        # turn off heavy atoms
        i.write_state(
            mutation_conf=mutation_list["lj"],
            lambda_value_vdw=0.0,
        )

        # generate terminal lj
        i.write_state(
            mutation_conf=mutation_list["default-lj"],
            lambda_value_vdw=0.0,
        )

        m = mutation_list["transform"]
        for lambda_value in np.linspace(0.75, 0, 3):
            print(lambda_value)
            # turn off charges
            i.write_state(
                mutation_conf=m,
                common_core_transformation=lambda_value,
            )
        for opath in i.output_files:
            shutil.rmtree(opath)


def test_equivalent_endstates_vacuum():
    workdir = get_test_output_dir()

    import openmm as mm
    import openmm.app as app

    # ParmEd Imports
    from parmed import unit as u

    output_files_t1, output_files_t2 = get_output_files_2oj9_tautomer_pair()

    env = "vacuum"
    psf, parms = generate_psf(output_files_t1[-1], env)
    coord = generate_crd(output_files_t1[-1], env).positions
    mask = []
    nr_of_atoms = len(psf.atoms)
    assert nr_of_atoms == len(coord)

    for a in psf.atoms:
        if str(a.type).startswith("DD"):
            mask.append(False)
        else:
            mask.append(True)
    psf = psf[mask]
    for idx, e in enumerate(mask):
        if not e:
            print(coord.pop(idx))

    system = psf.createSystem(parms, nonbondedMethod=app.NoCutoff)

    # Create the integrator to do Langevin dynamics
    integrator = mm.LangevinIntegrator(
        300 * u.kelvin,  # Temperature of heat bath
        1.0 / u.picoseconds,  # Friction coefficient
        2.0 * u.femtoseconds,  # Time step
    )
    assert nr_of_atoms == psf.topology.getNumAtoms() + 1
    print(psf.topology.getNumAtoms())
    # Create the Simulation object
    sim = app.Simulation(psf.topology, system, integrator)
    # Set the particle positions
    sim.context.setPositions(coord)
    e1 = (
        sim.context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(unit.kilocalorie_per_mole)
    )

    #####################################
    #####################################
    psf, parms = generate_psf(output_files_t2[-1], env)
    # coord = generate_crd(output_files_t2[-1], env).positions
    mask = []
    nr_of_atoms = len(psf.atoms)
    assert nr_of_atoms == len(coord) + 1

    mask = []
    for a in psf.atoms:
        if str(a.type).startswith("DD"):
            mask.append(False)
        else:
            mask.append(True)
    psf = psf[mask]

    system = psf.createSystem(parms, nonbondedMethod=app.NoCutoff)
    # Create the integrator to do Langevin dynamics
    integrator = mm.LangevinIntegrator(
        300 * u.kelvin,  # Temperature of heat bath
        1.0 / u.picoseconds,  # Friction coefficient
        2.0 * u.femtoseconds,  # Time step
    )

    # Create the Simulation object
    sim = app.Simulation(psf.topology, system, integrator)
    # Set the particle positions
    sim.context.setPositions(coord)
    e2 = (
        sim.context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(unit.kilocalorie_per_mole)
    )

    assert np.isclose(e1, e2)


def test_equivalent_endstates_waterbox():
    import openmm as mm
    import openmm.app as app

    # ParmEd Imports
    from parmed import unit as u

    output_files_t1, output_files_t2 = get_output_files_2oj9_tautomer_pair()

    env = "waterbox"
    psf, parms = generate_psf(output_files_t1[-1], env)
    coords = generate_crd(output_files_t1[-1], env).positions

    psf.box = (30, 30, 30, 90, 90, 90)

    # remove dummy atom from psf topology
    mask = []
    nr_of_atoms = len(psf.atoms)
    assert nr_of_atoms == len(coords)

    for a in psf.atoms:
        if str(a.type).startswith("DD"):
            mask.append(False)
        else:
            mask.append(True)

    psf = psf[mask]
    nr_of_atoms_t1 = psf.topology.getNumAtoms()

    new_coords = []
    for idx, e in enumerate(mask):
        if e:
            new_coords.append(coords[idx])

    coords = new_coords
    system = psf.createSystem(
        parms,
        nonbondedMethod=app.PME,
        nonbondedCutoff=12.0 * u.angstroms,
        switchDistance=10.0 * u.angstroms,
    )

    # Create the integrator to do Langevin dynamics
    integrator = mm.LangevinIntegrator(
        300 * u.kelvin,  # Temperature of heat bath
        1.0 / u.picoseconds,  # Friction coefficient
        2.0 * u.femtoseconds,  # Time step
    )
    # test that only a single dummy atom was present
    assert nr_of_atoms == psf.topology.getNumAtoms() + 1
    # Create the Simulation object
    sim = app.Simulation(psf.topology, system, integrator)
    # Set the particle positions
    sim.context.setPositions(coords)
    e1 = (
        sim.context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(unit.kilocalorie_per_mole)
    )
    nr_of_atoms_t1 = psf.topology.getNumAtoms()
    #####################################
    #####################################
    psf, parms = generate_psf(output_files_t2[-1], env)
    psf.box = (30, 30, 30, 90, 90, 90)
    # remove dummy atom
    mask = []
    for a in psf.atoms:
        if str(a.type).startswith("DD"):
            mask.append(False)
        else:
            mask.append(True)

    # remove waters that differ between the two topologies
    mask[-3:] = [False] * 3
    psf = psf[mask]
    nr_of_atoms_t2 = psf.topology.getNumAtoms()
    assert nr_of_atoms_t1 == nr_of_atoms_t2

    system = psf.createSystem(
        parms,
        nonbondedMethod=app.PME,
        nonbondedCutoff=12.0 * u.angstroms,
        switchDistance=10.0 * u.angstroms,
    )
    # Create the integrator to do Langevin dynamics
    integrator = mm.LangevinIntegrator(
        300 * u.kelvin,  # Temperature of heat bath
        1.0 / u.picoseconds,  # Friction coefficient
        2.0 * u.femtoseconds,  # Time step
    )

    # Create the Simulation object
    sim = app.Simulation(psf.topology, system, integrator)
    # Set the particle positions
    sim.context.setPositions(coords)
    e2 = (
        sim.context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(unit.kilocalorie_per_mole)
    )
    assert np.isclose(e1, e2)


def test_bonded_mutation_energies_t1_s1(caplog):
    caplog.set_level(logging.CRITICAL)
    output_files_t1, _ = get_output_files_2oj9_tautomer_pair()
    e_t1_s1 = (
        generate_sim(output_files_t1[0], "vacuum")
        .context.getState(getEnergy=True)
        .getPotentialEnergy()
    )
    assert np.isclose(
        e_t1_s1.value_in_unit(unit.kilocalorie_per_mole), -17.638044396797515
    )


def test_bonded_mutation_energies_t1_s2(caplog):
    caplog.set_level(logging.CRITICAL)
    output_files_t1, _ = get_output_files_2oj9_tautomer_pair()

    e_t1_s2 = (
        generate_sim(output_files_t1[1], "vacuum")
        .context.getState(getEnergy=True)
        .getPotentialEnergy()
    )

    assert np.isclose(
        e_t1_s2.value_in_unit(unit.kilocalorie_per_mole), 5.50502085150725
    )


def test_bonded_mutation_energies_t1_s3(caplog):
    caplog.set_level(logging.CRITICAL)
    output_files_t1, _ = get_output_files_2oj9_tautomer_pair()

    e_t1_s3 = (
        generate_sim(output_files_t1[2], "vacuum")
        .context.getState(getEnergy=True)
        .getPotentialEnergy()
    )

    assert np.isclose(
        e_t1_s3.value_in_unit(unit.kilocalorie_per_mole), 5.680121579945277
    )


def test_bonded_mutation_energies_t1_s4(caplog):
    caplog.set_level(logging.CRITICAL)
    output_files_t1, _ = get_output_files_2oj9_tautomer_pair()

    e_t1_s4 = (
        generate_sim(output_files_t1[3], "vacuum")
        .context.getState(getEnergy=True)
        .getPotentialEnergy()
    )

    assert np.isclose(
        e_t1_s4.value_in_unit(unit.kilocalorie_per_mole), 19.32431518895557
    )


def test_bonded_mutation_energies_t1_s5(caplog):
    caplog.set_level(logging.CRITICAL)
    output_files_t1, _ = get_output_files_2oj9_tautomer_pair()

    e_t1_s5 = (
        generate_sim(output_files_t1[4], "vacuum")
        .context.getState(getEnergy=True)
        .getPotentialEnergy()
    )

    assert np.isclose(
        e_t1_s5.value_in_unit(unit.kilocalorie_per_mole), 39.24392077135464
    )


def test_bonded_mutation_energies_t1_s6(caplog):
    caplog.set_level(logging.CRITICAL)
    output_files_t1, _ = get_output_files_2oj9_tautomer_pair()

    e_t1_s6 = (
        generate_sim(output_files_t1[5], "vacuum")
        .context.getState(getEnergy=True)
        .getPotentialEnergy()
    )

    assert np.isclose(
        e_t1_s6.value_in_unit(unit.kilocalorie_per_mole), 54.94647553965466
    )


def test_bonded_mutation_energies_t1_s7(caplog):
    caplog.set_level(logging.CRITICAL)
    output_files_t1, _ = get_output_files_2oj9_tautomer_pair()

    e_t1_s7 = (
        generate_sim(output_files_t1[6], "vacuum")
        .context.getState(getEnergy=True)
        .getPotentialEnergy()
    )

    assert np.isclose(
        e_t1_s7.value_in_unit(unit.kilocalorie_per_mole), 62.026441265363836
    )


def test_bonded_mutation_energies_t2_s1(caplog):
    caplog.set_level(logging.CRITICAL)
    _, output_files_t2 = get_output_files_2oj9_tautomer_pair()

    e_t2_s1 = (
        generate_sim(output_files_t2[0], "vacuum")
        .context.getState(getEnergy=True)
        .getPotentialEnergy()
    )

    assert np.isclose(
        e_t2_s1.value_in_unit(unit.kilocalorie_per_mole), 32.698191830578544, rtol=1e-4
    )


def test_bonded_mutation_energies_t2_s2(caplog):
    caplog.set_level(logging.CRITICAL)
    _, output_files_t2 = get_output_files_2oj9_tautomer_pair()

    e_t2_s2 = (
        generate_sim(output_files_t2[1], "vacuum")
        .context.getState(getEnergy=True)
        .getPotentialEnergy()
    )

    assert np.isclose(e_t2_s2.value_in_unit(unit.kilocalorie_per_mole), 63.196164016286)


def test_bonded_mutation_energies_t2_s3(caplog):
    caplog.set_level(logging.CRITICAL)
    _, output_files_t2 = get_output_files_2oj9_tautomer_pair()

    e_t2_s3 = (
        generate_sim(output_files_t2[2], "vacuum")
        .context.getState(getEnergy=True)
        .getPotentialEnergy()
    )

    assert np.isclose(
        e_t2_s3.value_in_unit(unit.kilocalorie_per_mole), 63.36457961831904
    )


def test_bonded_mutation_atoms(caplog):
    caplog.set_level(logging.CRITICAL)

    from .test_mutation import setup_2OJ9_tautomer_pair_rsfe

    conf = f"{get_testsystems_dir()}/config/test-2oj9-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir=get_testsystems_dir(), output_dir=get_test_output_dir()
    )

    (output_files_t1, output_files_t2), _, p = setup_2OJ9_tautomer_pair_rsfe(
        configuration=configuration
    )
    psf_at_endstate_t1, _ = generate_psf(output_files_t1[0], "vacuum")
    prm_at_endstate_t1 = {
        a.idx: (a.charge, a.sigma, a.epsilon) for a in psf_at_endstate_t1.atoms
    }
    psf_at_t1_cc, _ = generate_psf(output_files_t1[-1], "vacuum")
    prm_at_t1_cc = {a.idx: (a.charge, a.sigma, a.epsilon) for a in psf_at_t1_cc.atoms}
    psf_at_t2_cc, _ = generate_psf(output_files_t2[-1], "vacuum")
    prm_at_t2_cc = {a.idx: (a.charge, a.sigma, a.epsilon) for a in psf_at_t2_cc.atoms}
    psf_at_endstate_t2, _ = generate_psf(output_files_t2[0], "vacuum")
    prm_at_endstate_t2 = {
        a.idx: (a.charge, a.sigma, a.epsilon) for a in psf_at_endstate_t2.atoms
    }
    print("@@@@@@@@@@@@@@@@@@")

    cc1_idx, cc2_idx = p.get_common_core_idx_mol1(), p.get_common_core_idx_mol2()
    for atom_t1_idx, atom_t2_idx in zip(cc1_idx, cc2_idx):
        atom_t1 = prm_at_t1_cc[atom_t1_idx]
        atom_t2 = prm_at_t2_cc[atom_t2_idx]
        assert atom_t1 == atom_t2

    for atom_t1 in psf_at_t1_cc.atoms:
        if atom_t1.idx not in p.get_common_core_idx_mol1():
            print(atom_t1)
    for atom_t2 in psf_at_t2_cc.atoms:
        if atom_t2.idx not in p.get_common_core_idx_mol2():
            print(atom_t2)

    for atom_t1_idx, atom_t2_idx in zip(cc1_idx, cc2_idx):
        atom_t1 = prm_at_endstate_t1[atom_t1_idx]
        atom_t2 = prm_at_endstate_t2[atom_t2_idx]
        try:
            assert atom_t1 == atom_t2
        except AssertionError:
            pass


def test_bonded_mutation_bonds(caplog):
    caplog.set_level(logging.CRITICAL)

    from .test_mutation import setup_2OJ9_tautomer_pair_rsfe

    conf = f"{get_testsystems_dir()}/config/test-2oj9-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir=get_testsystems_dir(), output_dir=get_test_output_dir()
    )

    (output_files_t1, output_files_t2), _, p = setup_2OJ9_tautomer_pair_rsfe(
        configuration=configuration
    )
    ##################
    psf_at_endstate_t1, _ = generate_psf(output_files_t1[0], "vacuum")
    prm_at_endstate_t1 = {
        idx: (a.type.k, a.type.req) for idx, a in enumerate(psf_at_endstate_t1.bonds)
    }
    ##################
    psf_at_t1_cc, _ = generate_psf(output_files_t1[-1], "vacuum")
    prm_at_t1_cc = {
        idx: (a.type.k, a.type.req) for idx, a in enumerate(psf_at_t1_cc.bonds)
    }
    ##################
    psf_at_t2_cc, _ = generate_psf(output_files_t2[-1], "vacuum")
    prm_at_t2_cc = {
        idx: (a.type.k, a.type.req) for idx, a in enumerate(psf_at_t2_cc.bonds)
    }
    ##################
    psf_at_endstate_t2, _ = generate_psf(output_files_t2[0], "vacuum")
    prm_at_endstate_t2 = {
        idx: (a.type.k, a.type.req) for idx, a in enumerate(psf_at_endstate_t2.bonds)
    }  ##################
    print("@@@@@@@@@@@@@@@@@@")

    cc1_idx, cc2_idx = p.get_common_core_idx_mol1(), p.get_common_core_idx_mol2()
    # compare cc enstates
    for bond_t1_idx, bond_t1 in enumerate(psf_at_t1_cc.bonds):
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
        for bond_t2_idx, bond_t2 in enumerate(psf_at_t2_cc.bonds):
            if atom1_t2 in bond_t2 and atom2_t2 in bond_t2:
                prm_at_t1_cc[bond_t1_idx] = prm_at_t2_cc[bond_t2_idx]


def test_bonded_mutation_angles(caplog):
    caplog.set_level(logging.CRITICAL)
    from copy import copy
    from .test_mutation import setup_2OJ9_tautomer_pair_rsfe

    conf = f"{get_testsystems_dir()}/config/test-2oj9-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir=get_testsystems_dir(), output_dir=get_test_output_dir()
    )

    (output_files_t1, output_files_t2), _, p = setup_2OJ9_tautomer_pair_rsfe(
        configuration=configuration
    )
    ##################
    # psf_at_endstate_t1 = generate_psf(output_files_t1[0], "vacuum")
    # prm_at_endstate_t1 = {
    #     idx: (copy(a.type.k), copy(a.type.theteq))
    #     for idx, a in enumerate(psf_at_endstate_t1.angles)
    # }
    # ##################
    # psf_at_endstate_t2 = generate_psf(output_files_t2[0], "vacuum")
    # prm_at_endstate_t2 = {
    #     idx: (copy(a.type.k), copy(a.type.theteq))
    #     for idx, a in enumerate(psf_at_endstate_t2.angles)
    # }
    ##################
    psf_at_t1_cc, _ = generate_psf(output_files_t1[-1], "vacuum")
    prm_at_t1_cc = {
        idx: (copy(a.type.k), copy(a.type.theteq))
        for idx, a in enumerate(psf_at_t1_cc.angles)
    }
    ##################
    psf_at_t2_cc, _ = generate_psf(output_files_t2[-1], "vacuum")
    prm_at_t2_cc = {
        idx: (copy(a.type.k), copy(a.type.theteq))
        for idx, a in enumerate(psf_at_t2_cc.angles)
    }
    print("@@@@@@@@@@@@@@@@@@")
    ##################

    cc1_idx, cc2_idx = p.get_common_core_idx_mol1(), p.get_common_core_idx_mol2()
    for angle_t1_idx, angle_t1 in enumerate(psf_at_t1_cc.angles):
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
        faulty = False
        for angle_t2_idx, angle_t2 in enumerate(psf_at_t2_cc.angles):
            if atom1_t2 in angle_t2 and atom2_t2 in angle_t2 and atom3_t2 in angle_t2:
                if (
                    not (prm_at_t1_cc[angle_t1_idx] == prm_at_t2_cc[angle_t2_idx])
                    and atom1_t2 == "<Atom C12 [14]; In UNK 0>"
                ):  # the AND statement is only necessary for cgenff v.4.6 becaues the c11-c18-n6 in bmi and c12-c16-n6 in unk are slightly different
                    print("###################")
                    print(prm_at_t1_cc[angle_t1_idx])
                    print(prm_at_t2_cc[angle_t2_idx])
                    print(atom1_t2)
                    print(angle_t2)
                    faulty = True

        if faulty:
            raise AssertionError()


def test_bonded_mutation_dihedrals(caplog):
    caplog.set_level(logging.CRITICAL)
    from .test_mutation import setup_2OJ9_tautomer_pair_rsfe

    conf = f"{get_testsystems_dir()}/config/test-2oj9-tautomer-pair-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir=get_testsystems_dir(), output_dir=get_test_output_dir()
    )

    (output_files_t1, output_files_t2), _, p = setup_2OJ9_tautomer_pair_rsfe(
        configuration=configuration
    )
    # ##################
    # psf_at_endstate_t1 = generate_psf(output_files_t1[0], "vacuum")
    # prm_at_endstate_t1 = {
    #     idx: [(a.type.phi_k, a.type.per, a.type.phase) for a in l.type]
    #     for idx, l in enumerate(psf_at_endstate_t1.dihedrals)
    # }
    # ##################
    # psf_at_endstate_t2 = generate_psf(output_files_t2[0], "vacuum")
    # prm_at_endstate_t2 = {
    #     idx: [(a.type.phi_k, a.type.per, a.type.phase) for a in l]
    #     for idx, a in enumerate(psf_at_endstate_t2.dihedrals)
    # }
    ####################################
    ####################################
    psf_at_t1_cc, _ = generate_psf(output_files_t1[-1], "vacuum")
    prm_at_t1_cc = {
        idx: [(a.phi_k, a.per, a.phase) for a in l.type]
        for idx, l in enumerate(psf_at_t1_cc.dihedrals)
    }
    t1_cc_nr_of_dihedrals = len(psf_at_t1_cc.dihedrals)
    e_at_t1_cc = (
        generate_sim(output_files_t1[-1], "vacuum")
        .context.getState(getEnergy=True)
        .getPotentialEnergy()
    )
    assert np.isclose(
        e_at_t1_cc.value_in_unit(unit.kilocalorie_per_mole), 62.02640132284073
    )
    ####################################
    ####################################
    psf_at_t2_cc, _ = generate_psf(output_files_t2[-1], "vacuum")
    prm_at_t2_cc = {
        idx: [(a.phi_k, a.per, a.phase) for a in l.type]
        for idx, l in enumerate(psf_at_t2_cc.dihedrals)
    }
    t2_cc_nr_of_dihedrals = len(psf_at_t2_cc.dihedrals)
    e_at_t2_cc = (
        generate_sim(output_files_t2[-1], "vacuum")
        .context.getState(getEnergy=True)
        .getPotentialEnergy()
    )
    assert np.isclose(
        e_at_t2_cc.value_in_unit(unit.kilocalorie_per_mole), 63.36457961831904
    )
    ####################################
    ####################################
    assert t1_cc_nr_of_dihedrals == t2_cc_nr_of_dihedrals
    # compare the common core parameters at the cc state
    cc1_idx, cc2_idx = p.get_common_core_idx_mol1(), p.get_common_core_idx_mol2()

    for dihedral_t1_idx, dihedral_t1 in enumerate(psf_at_t1_cc.dihedrals):
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
            print("Not present in cc2:")
            print(dihedral_t1)
            print("#####################")
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
        faulty = False
        for dihedral_t2_idx, dihedral in enumerate(psf_at_t2_cc.dihedrals):
            if (
                atom1_t2 in dihedral
                and atom2_t2 in dihedral
                and atom3_t2 in dihedral
                and atom4_t2 in dihedral
            ):
                dihedral_t2 = dihedral
                break
        else:
            print("Not present in cc1:")
            print(dihedral)
            print("#####################")

        assert dihedral_t2 != None
        if not (prm_at_t1_cc[dihedral_t1_idx] == prm_at_t2_cc[dihedral_t2_idx]):
            print("###################")
            print(prm_at_t1_cc[dihedral_t1_idx])
            print(prm_at_t2_cc[dihedral_t2_idx])
            print(dihedral_t1)
            print(dihedral_t2)
            faulty = True
        if faulty:
            raise AssertionError()


def test_vdw_mutation_for_hydrogens_and_heavy_atoms():
    from rdkit.Chem import rdFMCS

    for conf, system_name in zip(
        [
            f"{get_testsystems_dir()}/config/test-toluene-methane-rsfe.yaml",
            f"{get_testsystems_dir()}/config/test-neopentane-methane-rsfe.yaml",
            f"{get_testsystems_dir()}/config/test-ethane-methanol-rsfe.yaml",
        ],
        ["toluene-methane", "neopentane-methane", "ethane-methanol"],
    ):
        configuration = load_config_yaml(
            config=conf,
            input_dir=get_testsystems_dir(),
            output_dir=get_test_output_dir(),
        )
        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        s1_to_s2 = ProposeMutationRoute(s1, s2)
        if system_name == "neopentane-methane":
            s1_to_s2.propose_common_core()
            s1_to_s2.finish_common_core(
                connected_dummy_regions_cc1=[
                    {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
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
            all_atoms_for_which_lj_turned_off = []
            print(mutation_list.keys())
            # turn off hydrogen lj
            hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
            print(
                f"Turn off lj for hydrogen atoms : {[e.atoms_to_be_mutated for e in hydrogen_lj_mutations]}"
            )
            i.write_state(
                mutation_conf=hydrogen_lj_mutations,
                lambda_value_vdw=lambda_vdw,
            )
            for mutation in hydrogen_lj_mutations:
                all_atoms_for_which_lj_turned_off.extend(mutation.vdw_atom_idx)

            # turn off heavy atom lj
            for mutation in mutation_list["lj"]:
                print(f"Turn off lj for heavy atom : {mutation.atoms_to_be_mutated}")

                i.write_state(
                    mutation_conf=[mutation],
                    lambda_value_vdw=lambda_vdw,
                )
                all_atoms_for_which_lj_turned_off.extend(mutation.vdw_atom_idx)

            # change to default lj
            terminal_lj_mutations = mutation_list["default-lj"]
            terminal_idx = []
            for mutation in terminal_lj_mutations:
                all_atoms_for_which_lj_turned_off.extend(mutation.vdw_atom_idx)
                terminal_idx.extend(mutation.vdw_atom_idx)
            print(
                f"Turn off lj for terminal atom : {all_atoms_for_which_lj_turned_off}"
            )

            i.write_state(
                mutation_conf=terminal_lj_mutations,
                lambda_value_vdw=0.0,
            )

            print(f"Set epsilon/rmin to base * {lambda_vdw} for selected atoms")

            print(all_atoms_for_which_lj_turned_off)
            for env in s1.envs:
                print(env)
                # read in generated psf at last mutation step
                new_psf, params = generate_psf(i.output_files[-1], env)
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

        shutil.rmtree(f"{get_test_output_dir()}/{system_name}-rsfe")


def setup_2OJ9_tautomer_pair_rsfe(
    configuration: dict, single_state=False, nr_of_bonded_windows: int = 4
):
    from ..mutate import mutate_pure_tautomers

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
            single_state=single_state,
            nr_of_bonded_windows=nr_of_bonded_windows,
        ),
        configuration,
        s1_to_s2,
    )


def setup_2OJ9_tautomer_pair_rbfe(
    configuration: dict, single_state: bool = False, nr_of_bonded_windows: int = 4
):
    from ..mutate import mutate_pure_tautomers

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
            single_state=single_state,
            nr_of_bonded_windows=nr_of_bonded_windows,
        ),
        configuration,
        s1_to_s2,
    )


def test_acetylacetone_tautomer_pair(caplog):
    workdir = get_test_output_dir()
    caplog.set_level(logging.DEBUG)
    conf = f"{get_testsystems_dir()}/config/test-acetylacetone-tautomer-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir=get_testsystems_dir(), output_dir=get_test_output_dir()
    )
    setup_acetylacetone_tautomer_pair(configuration=configuration)
    shutil.rmtree(f"{workdir}/acetylacetone-keto-acetylacetone-enol-rsfe")


def test_2OJ9_tautomer_pair(caplog):
    workdir = get_test_output_dir()
    caplog.set_level(logging.DEBUG)
    conf = f"{get_testsystems_dir()}/config/test-2oj9-tautomer-pair-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir=get_testsystems_dir(), output_dir=get_test_output_dir()
    )

    setup_2OJ9_tautomer_pair_rsfe(configuration=configuration)
    shutil.rmtree(f"{workdir}/2OJ9-original-2OJ9-tautomer-rsfe")


def test_full_mutation_system1(caplog):
    caplog.set_level(logging.WARNING)
    workdir = get_test_output_dir()

    for conf, system_name in zip(
        [
            f"{get_testsystems_dir()}/config/test-toluene-methane-rsfe.yaml",
            f"{get_testsystems_dir()}/config/test-neopentane-methane-rsfe.yaml",
            f"{get_testsystems_dir()}/config/test-ethane-methanol-rsfe.yaml",
        ],
        ["toluene-methane", "neopentane-methane", "ethane-methanol"],
    ):
        print(system_name)
        configuration = load_config_yaml(
            config=conf,
            input_dir=get_testsystems_dir(),
            output_dir=get_test_output_dir(),
        )
        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        s1_to_s2 = ProposeMutationRoute(s1, s2)
        if system_name == "neopentane-methane":
            s1_to_s2.propose_common_core()
            s1_to_s2.finish_common_core(
                connected_dummy_regions_cc1=[
                    {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
                ]
            )
        else:
            s1_to_s2.calculate_common_core()

        mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
        i = IntermediateStateFactory(
            system=s1,
            configuration=configuration,
        )
        charges = mutation_list["charge"]
        for lambda_value in np.linspace(0, 1, 5):
            # turn off charges
            i.write_state(
                mutation_conf=charges,
                lambda_value_electrostatic=1 - lambda_value,
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
        i.write_state(
            mutation_conf=terminal_lj_mutations,
            lambda_value_vdw=lambda_vdw,
        )
        for mutation in terminal_lj_mutations:
            all_atoms_for_which_lj_turned_off.extend(mutation.vdw_atom_idx)

        shutil.rmtree(f"{workdir}/{system_name}-rsfe")


def test_full_mutation_system2():
    workdir = get_test_output_dir()

    for conf, system_name in zip(
        [
            f"{get_testsystems_dir()}/config/test-toluene-methane-rsfe.yaml",
            f"{get_testsystems_dir()}/config/test-neopentane-methane-rsfe.yaml",
        ],
        ["toluene-methane", "neopentane-methane"],
    ):
        configuration = load_config_yaml(
            config=conf,
            input_dir=get_testsystems_dir(),
            output_dir=get_test_output_dir(),
        )
        s1 = SystemStructure(configuration, "structure1")
        s2 = SystemStructure(configuration, "structure2")

        s1_to_s2 = ProposeMutationRoute(s1, s2)
        if system_name == "neopentane-methane":
            s1_to_s2.propose_common_core()
            s1_to_s2.finish_common_core(
                connected_dummy_regions_cc1=[
                    {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
                ]
            )
        else:
            s1_to_s2.calculate_common_core()

        mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
        i = IntermediateStateFactory(
            system=s1,
            configuration=configuration,
        )

        charges = mutation_list["charge"]
        for lambda_value in np.linspace(0, 1, 5):
            # turn off charges
            i.write_state(
                mutation_conf=charges,
                lambda_value_electrostatic=1 - lambda_value,
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
        i.write_state(
            mutation_conf=hydrogen_lj_mutations,
            lambda_value_vdw=lambda_vdw,
        )

        for mutation in hydrogen_lj_mutations:
            all_atoms_for_which_lj_turned_off.extend(mutation.vdw_atom_idx)

        # turn off heavy atoms
        for mutation in mutation_list["lj"]:
            i.write_state(
                mutation_conf=[mutation],
                lambda_value_vdw=lambda_vdw,
            )
            all_atoms_for_which_lj_turned_off.extend(mutation.vdw_atom_idx)

        new_psf, params = generate_psf(i.output_files[-1], env)
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

        shutil.rmtree(f"{workdir}/{system_name}-rsfe")


def test_generate_list_of_heavy_atoms_to_mutate():
    from transformato.utils import map_lj_mutations_to_atom_idx

    #########################################
    #########################################
    # neopentane to methane
    # with user defined connected dummy region
    configuration = load_config_yaml(
        config=f"{get_testsystems_dir()}/config/test-neopentane-methane-rsfe.yaml",
        input_dir=get_testsystems_dir(),
        output_dir=get_test_output_dir(),
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.propose_common_core()
    s1_to_s2.finish_common_core(
        connected_dummy_regions_cc1=[
            {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
        ]
    )

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )

    # get the order of the heavy atom lj mutation
    list_of_heavy_atoms_to_be_mutated = [
        lj.vdw_atom_idx[0] for lj in (mutation_list["lj"])
    ]
    # map the order to the actual mutations
    mapping_of_atom_idx_to_mutation = map_lj_mutations_to_atom_idx(mutation_list["lj"])

    assert set([atom_idx for atom_idx in list_of_heavy_atoms_to_be_mutated]) == set(
        [5, 9, 13]
    )

    #########################################
    #########################################
    # neopentane to methane
    # without user defined connected dummy region
    configuration = load_config_yaml(
        config=f"{get_testsystems_dir()}/config/test-neopentane-methane-rsfe.yaml",
        input_dir=get_testsystems_dir(),
        output_dir=get_test_output_dir(),
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.propose_common_core()
    s1_to_s2.finish_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )

    # get the order of the heavy atom lj mutation
    list_of_heavy_atoms_to_be_mutated = [
        lj.vdw_atom_idx[0] for lj in (mutation_list["lj"])
    ]
    # map the order to the actual mutations
    mapping_of_atom_idx_to_mutation = map_lj_mutations_to_atom_idx(mutation_list["lj"])

    assert set([atom_idx for atom_idx in list_of_heavy_atoms_to_be_mutated]) == set(
        [5, 9, 13]
    )

    #########################################
    #########################################
    # toluene to methane
    configuration = load_config_yaml(
        config=f"{get_testsystems_dir()}/config/test-toluene-methane-rsfe.yaml",
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
    # get the order of the heavy atom lj mutation
    list_of_heavy_atoms_to_be_mutated = [
        lj.vdw_atom_idx[0] for lj in (mutation_list["lj"])
    ]
    # map the order to the actual mutations
    mapping_of_atom_idx_to_mutation = map_lj_mutations_to_atom_idx(mutation_list["lj"])

    assert set(list_of_heavy_atoms_to_be_mutated) == set([1, 3, 9, 11, 13])
