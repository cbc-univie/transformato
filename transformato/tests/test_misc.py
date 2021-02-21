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
