"""
Unit and regression test for the transformato package.
"""

from transformato.utils import load_config_yaml
from transformato.analysis import return_reduced_potential
from simtk import unit
import numpy as np
# read in specific topology with parameters
from transformato import (
    load_config_yaml,
)
from .test_mutation import (
    _set_output_files_2oj9_tautomer_pair,
)


def test_reduced_energy():
    from transformato.constants import temperature as T

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
    from transformato.constants import temperature as T

    kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
    beta = 1.0 / (kB * T)

    e = -41264.39524669979 * unit.kilojoule_per_mole
    rE = e * beta
    assert np.isclose(rE, -16371.30301422)


def test_change_platform():
    from ..constants import change_platform, platform
    from ..utils import load_config_yaml

    configuration = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-rsfe.yaml",
        input_dir=".",
        output_dir=".",
    )

    change_platform(configuration)
    print(configuration["simulation"]["GPU"])
    print(platform)
    if platform.upper() == "CPU":
        assert configuration["simulation"]["GPU"] == False
    else:
        assert configuration["simulation"]["GPU"] == True


def test_scaling():
    import numpy as np

    for i in np.linspace(1, 0, 11):
        f = max((1 - ((1 - i) * 2)), 0.0)
        print(f"{i}:{f}")

    print("##########################")
    for i in np.linspace(1, 0, 11):
        f = 1 - min((i) * 2, 1.0)
        print(f"{i}:{f}")

    print("##########################")


def test_old_scaling():
    import numpy as np

    for i in np.linspace(1, 0, 11):
        f = 1 - (1 - i) * 2
        print(f"{i}:{f}")
    print("##########################")

    for i in np.linspace(1, 0, 11):
        f = 1 - (1 - (1 - i)) * 2
        print(f"{i}:{f}")

def test_reading_of_coords():

    import mdtraj as md

    env = "vacuum"
    base = "data/2OJ9-original-2OJ9-tautomer-rsfe/2OJ9-original/"
    output_files_t1, _ = _set_output_files_2oj9_tautomer_pair()

    conf = "transformato/tests/config/test-2oj9-tautomer-pair-rsfe.yaml"

    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    b = output_files_t1[0]
    traj_load = md.load_dcd(
        f"{b}/lig_in_{env}.dcd",
        f"{b}/lig_in_{env}.psf",
    )
    print(traj_load.xyz[0])
        
    traj_open = md.open(f"{b}/lig_in_{env}.dcd")
    xyz, unitcell_lengths, _ = traj_open.read()
    xyz = xyz/10
    print(xyz[0])
    assert np.allclose(xyz[0], traj_load.xyz[0])