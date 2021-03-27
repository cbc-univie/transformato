"""
Unit and regression test for the transformato package.
"""

import logging
from transformato.utils import load_config_yaml
from transformato.analysis import return_reduced_potential
import pytest
from simtk import unit
import numpy as np


def test_reduced_energy():
    from transformato.constants import temperature as T

    e = -41264.39524669979 * unit.kilojoule_per_mole
    rE = return_reduced_potential(e, volume=0, temperature=T)
    print(rE)
    assert np.isclose(rE, -16371.30301422)

    e = -1099.41855 * unit.kilocalorie_per_mole
    rE = return_reduced_potential(e, volume=0, temperature=T)
    print(rE)
    assert np.isclose(rE, -1824.9982986086145)


def test_convert_to_kT():
    from transformato.constants import temperature as T

    kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
    beta = 1.0 / (kB * T)

    e = -41264.39524669979 * unit.kilojoule_per_mole
    rE = e * beta
    assert np.isclose(rE, -16371.30301422)


def test_change_platform():
    from ..constants import check_platform, platform
    from ..utils import load_config_yaml

    configuration = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
        input_dir=".",
        output_dir=".",
    )

    check_platform(configuration)
    print(configuration)
    if platform == "CPU":
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
