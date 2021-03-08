"""
Unit and regression test for the transformato package.
"""

import logging
from transformato.utils import load_config_yaml

import pytest


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
