"""
Unit and regression test for the transformato package.
"""

import os

import numpy as np
import pytest
import logging

# read in specific topology with parameters
from transformato import (
    load_config_yaml,
)


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_compare_energies_2OJ9_tautomer_pair_vacuum(caplog):
    caplog.set_level(logging.WARNING)
    from transformato import FreeEnergyCalculator
    import mdtraj as md

    base = "/home/mwieder/Work/Projects/transformato/data/2OJ9-original-2OJ9-tautomer-solvation-free-energy/2OJ9-original/"
    # run_simulation(output_files_t1[:2])
    # print(output_files_t1)
    output_files_t1 = [
        f"{base}/intst1/",
        f"{base}/intst2/",
        f"{base}/intst3/",
        f"{base}/intst4/",
        f"{base}/intst5/",
        f"{base}/intst6/",
        f"{base}/intst7/",
    ]

    conf = (
        "transformato/tests/config/test-2oj9-tautomer-pair-solvation-free-energy.yaml"
    )

    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    f = FreeEnergyCalculator(configuration, "2OJ9-original")
    for idx, b in enumerate(output_files_t1):
        traj = md.load_dcd(
            f"{b}/lig_in_vacuum.dcd",
            f"{b}/lig_in_vacuum.psf",
        )
        traj.save_dcd(f"{base}/traj.dcd")
        l_charmm = f._evaluated_e_on_all_snapshots_CHARMM(traj, idx + 1, "vacuum")
        l_openMM = f._evaluated_e_on_all_snapshots_openMM(traj, idx + 1, "vacuum")

        assert len(l_charmm) == len(l_openMM)
        s = abs(np.array(l_charmm) - np.array(l_openMM))
        mae = np.sum(s) / len(s)
        assert mae < 0.005
        for e_charmm, e_openMM in zip(l_charmm, l_openMM):
            assert np.isclose(e_charmm, e_openMM, rtol=1e-2)


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_compare_energies_2OJ9_tautomer_pair_waterbox(caplog):
    caplog.set_level(logging.WARNING)
    from transformato import FreeEnergyCalculator
    import mdtraj as md

    base = "/home/mwieder/Work/Projects/transformato/data/2OJ9-original-2OJ9-tautomer-solvation-free-energy/2OJ9-original/"
    output_files = [
        f"{base}/intst1/",
        f"{base}/intst2/",
        f"{base}/intst3/",
        f"{base}/intst4/",
        f"{base}/intst5/",
        f"{base}/intst6/",
        f"{base}/intst7/",
    ]

    conf = (
        "transformato/tests/config/test-2oj9-tautomer-pair-solvation-free-energy.yaml"
    )

    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir
    f = FreeEnergyCalculator(configuration, "2OJ9-original")
    for idx, b in enumerate(output_files[:1]):
        traj = md.load_dcd(
            f"{b}/lig_in_waterbox.dcd",
            f"{b}/lig_in_waterbox.psf",
        )
        traj.save_dcd(f"{base}/traj.dcd")
        l_charmm = f._evaluated_e_on_all_snapshots_CHARMM(traj, idx + 1, "waterbox")
        l_openMM = f._evaluated_e_on_all_snapshots_openMM(traj, idx + 1, "waterbox")
        assert len(l_charmm) == len(l_openMM)
        s = abs(np.array(l_charmm) - np.array(l_openMM))
        mae = np.sum(s) / len(s)
        print(mae)
        print(l_charmm)
        print(l_openMM)

        # assert mae < 0.005
        # for e_charmm, e_openMM in zip(l_charmm, l_openMM):
        #    assert np.isclose(e_charmm, e_openMM, rtol=1e-1)
    assert False


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
