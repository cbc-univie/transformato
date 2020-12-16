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


def test_read_yaml():
    """Sample test, will check ability to read yaml files"""
    settingsMap = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
        input_dir=".",
        output_dir="data/",
    )

    assert settingsMap["system"]["name"] == "toluene-methane-solvation-free-energy"
    assert settingsMap["system"]["structure1"]["tlc"] == "UNL"

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
