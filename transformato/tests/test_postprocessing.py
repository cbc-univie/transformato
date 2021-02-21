"""
Unit and regression test for the transformato package.
"""

import os

import numpy as np
import pytest


# read in specific topology with parameters
from transformato import (
    load_config_yaml,
)


def postprocessing(configuration, name="methane", engine="openMM", max_snapshots=300):
    from transformato import FreeEnergyCalculator

    f = FreeEnergyCalculator(configuration, name)
    f.load_trajs(nr_of_max_snapshots=max_snapshots)
    f.calculate_dG_to_common_core(engine=engine)
    ddG, dddG = f.end_state_free_energy_difference
    print(f"Free energy difference: {ddG}")
    print(f"Uncertanty: {dddG}")
    f.show_summary()
    return ddG, dddG


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

    ddG, dddG = postprocessing(
        configuration, name="methane", engine="CHARMM", max_snapshots=300
    )
    np.isclose(ddG, -1.2102764838282152, rtol=1e-8)


def test_run_methane_to_methane_cc_solvation_free_energy_with_openMM_postprocessing():

    conf = "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir
    ddG, dddG = postprocessing(
        configuration, name="methane", engine="openMM", max_snapshots=300
    )
    np.isclose(ddG, -1.2102764838282152, rtol=1e-8)


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
