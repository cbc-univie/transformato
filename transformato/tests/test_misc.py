"""
Unit and regression test for the transformato package.
"""
from transformato.constants import temperature as T
from openmm import unit
import numpy as np
import mdtraj as md
from transformato.utils import load_config_yaml
from transformato.analysis import return_reduced_potential
from transformato.tests.paths import get_test_output_dir
from transformato_testsystems.testsystems import (
    get_testsystems_dir,
    get_output_files_2oj9_tautomer_pair,
)


def test_reduced_energy():

    # with openMM generated traj evaluated with openMM
    e = -41264.39524669979 * unit.kilojoule_per_mole
    rE = return_reduced_potential(e, volume=0, temperature=T)
    print(rE)
    assert np.isclose(rE, -16371.30301422)

    # with openMM generated traj evaluated with CHARMM
    e = -1099.41855 * unit.kilocalorie_per_mole
    rE = return_reduced_potential(e, volume=0, temperature=T)
    print(rE)
    assert np.isclose(rE, -1824.9982986086145)

    # energy term in CHARMM traj
    # DYNA>        0      0.00000  -7775.74490   1914.51007  -9690.25498    377.14828
    e = -9690.25498 * unit.kilocalorie_per_mole  # ENERGgy
    rE = return_reduced_potential(e, volume=0, temperature=T)
    print(rE)
    assert np.isclose(rE, -16085.501605902184)


def test_convert_to_kT():

    kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
    beta = 1.0 / (kB * T)

    e = -41264.39524669979 * unit.kilojoule_per_mole
    rE = e * beta
    assert np.isclose(rE, -16371.30301422)


def test_change_platform():
    from transformato.constants import (
        change_platform_to_test_platform,
        test_platform_openMM,
        test_platform_CHARMM,
    )

    configuration = load_config_yaml(
        config=f"{get_testsystems_dir()}/config/test-toluene-methane-rsfe.yaml",
        input_dir=".",
        output_dir=get_test_output_dir(),
    )

    change_platform_to_test_platform(configuration, engine="openMM")
    print(configuration["simulation"]["GPU"])
    print(test_platform_openMM)
    if test_platform_openMM.upper() == "CPU":
        assert configuration["simulation"]["GPU"] == False
    else:
        assert configuration["simulation"]["GPU"] == True

    change_platform_to_test_platform(configuration, engine="CHARMM")
    print(configuration["simulation"]["GPU"])
    print(test_platform_openMM)
    if test_platform_CHARMM.upper() == "CPU":
        assert configuration["simulation"]["GPU"] == False
    else:
        assert configuration["simulation"]["GPU"] == True


def test_scaling():

    for i in np.linspace(1, 0, 11):
        f = max((1 - ((1 - i) * 2)), 0.0)
        print(f"{i}:{f}")

    print("##########################")
    for i in np.linspace(1, 0, 11):
        f = 1 - min((i) * 2, 1.0)
        print(f"{i}:{f}")

    print("##########################")


def test_old_scaling():

    for i in np.linspace(1, 0, 11):
        f = 1 - (1 - i) * 2
        print(f"{i}:{f}")
    print("##########################")

    for i in np.linspace(1, 0, 11):
        f = 1 - (1 - (1 - i)) * 2
        print(f"{i}:{f}")


def test_reading_of_coords():

    env = "vacuum"
    output_files_t1, _ = get_output_files_2oj9_tautomer_pair()

    conf = f"{get_testsystems_dir()}/config/test-2oj9-tautomer-pair-rsfe.yaml"

    configuration = load_config_yaml(
        config=conf, input_dir=get_testsystems_dir(), output_dir=get_test_output_dir()
    )  # NOTE: for preprocessing input_dir is the output dir

    b = output_files_t1[0]
    print(b)
    traj_load = md.load_dcd(
        f"{b}/lig_in_{env}.dcd",
        f"{b}/lig_in_{env}.psf",
    )
    print(traj_load.xyz[0])

    traj_open = md.open(f"{b}/lig_in_{env}.dcd")
    xyz, unitcell_lengths, _ = traj_open.read()
    xyz = xyz / 10
    print(xyz[0])
    assert np.allclose(xyz[0], traj_load.xyz[0])
