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


def test_transformato_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "transformato" in sys.modules


def test_run_toluene_to_methane_cc_solvation_free_energy_with_openMM():
    from transformato import FreeEnergyCalculator

    conf = (
        "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
    )

    configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")

    f = FreeEnergyCalculator(configuration, "toluene")


def test_generate_output_for_methane_cc_solvation_free_energy():
    from transformato.loeffler_systems import mutate_methane_to_methane_cc

    output_files, configuration = mutate_methane_to_methane_cc()


def test_generate_output_for_toluene_cc_solvation_free_energy():
    from transformato.loeffler_systems import mutate_toluene_to_methane_cc

    output_files, configuration = mutate_toluene_to_methane_cc()


def test_generate_output_for_toluene_cc_solvation_free_energy_with_test_conf():
    from transformato.loeffler_systems import mutate_toluene_to_methane_cc

    output_files, configuration = mutate_toluene_to_methane_cc(
        conf="transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml"
    )
    print(output_files)
