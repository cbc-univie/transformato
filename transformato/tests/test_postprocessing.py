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


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_methane_to_methane_cc_solvation_free_energy_with_CHARMM_postprocessing():
    from transformato import FreeEnergyCalculator

    conf = "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    f = FreeEnergyCalculator(configuration, "methane")
    f.load_trajs(nr_of_max_snapshots=300)
    f.calculate_dG_to_common_core(engine="CHARMM")
    ddG, dddG = f.end_state_free_energy_difference
    print(f"Free energy difference: {ddG}")
    print(f"Uncertanty: {dddG}")
    np.isclose(ddG, -1.2102764838282152, rtol=1e-8)
    f.show_summary()


def test_run_methane_to_methane_cc_solvation_free_energy_with_openMM_postprocessing():
    from transformato import FreeEnergyCalculator

    conf = "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    f = FreeEnergyCalculator(configuration, "methane")
    f.load_trajs(nr_of_max_snapshots=300)
    f.calculate_dG_to_common_core()
    ddG, dddG = f.end_state_free_energy_difference
    print(f"Free energy difference: {ddG}")
    print(f"Uncertanty: {dddG}")
    print(ddG)
    np.isclose(ddG, -1.2102764838282152, rtol=1e-8)
    f.show_summary()


def test_postprocessing_thinning():
    from transformato import FreeEnergyCalculator

    conf = "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    f = FreeEnergyCalculator(configuration, "methane")
    f.load_trajs(nr_of_max_snapshots=1000)
    assert len(f.snapshots.keys()) == 2  # entry for vacuum and waterbox
    assert f.nr_of_states == 3  # nr of lambda states

    print(f.snapshots["vacuum"])
    assert (
        len(f.snapshots["vacuum"]) == 2250
    )  # total:3000 frames, 75%: 2250 ---> for the individual traj: 1000 frames, 75% are 750, and we take all of these
    assert len(f.snapshots["waterbox"]) == 2250  # total:3000 frames, 75%: 2250

    f = FreeEnergyCalculator(configuration, "methane")
    f.load_trajs(nr_of_max_snapshots=500)

    print(f.snapshots["vacuum"])
    assert (
        len(f.snapshots["vacuum"]) == 1500
    )  # input had length of 1000, 25% removed gives 750 frames, taking only 500 of these we end up with 1500 frames
    assert len(f.snapshots["waterbox"]) == 1500  #
