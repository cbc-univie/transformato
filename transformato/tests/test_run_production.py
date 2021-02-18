"""
Unit and regression test for the transformato package.
"""

import os
import subprocess
import logging

import numpy as np
import pytest


def run_simulation(output_files):
    for path in sorted(output_files):
        # because path is object not string
        print(f"Start sampling for: {path}")
        exe = subprocess.run(
            ["bash", f"{str(path)}/simulation.sh", str(path)],
            check=True,
            capture_output=True,
            text=True,
        )
        print(exe.stdout)
        print("Capture stderr")
        print(exe.stderr)


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_toluene_to_methane_cc_solvation_free_energy_with_openMM():
    from transformato import FreeEnergyCalculator
    from transformato.loeffler_systems import mutate_toluene_to_methane_cc
    from .test_run_production import run_simulation
    from .test_postprocessing import postprocessing

    output_files, configuration = mutate_toluene_to_methane_cc(
        conf="transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml"
    )

    run_simulation(output_files)
    ddG, dddG = postprocessing(
        configuration, name="toluene", engine="openMM", max_snapshots=100
    )
    np.isclose(ddG, 8.9984, rtol=1e-2)


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_methane_to_methane_cc_solvation_free_energy_with_openMM():
    from transformato import FreeEnergyCalculator
    from transformato.loeffler_systems import mutate_methane_to_methane_cc
    from .test_run_production import run_simulation
    from .test_postprocessing import postprocessing

    output_files, configuration = mutate_methane_to_methane_cc(
        conf="transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml"
    )
    run_simulation(output_files)
    ddG, dddG = postprocessing(
        configuration, name="methane", engine="openMM", max_snapshots=300
    )
    np.isclose(ddG, 8.9984, rtol=1e-2)


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_methane_to_methane_cc_solvation_free_energy_with_openMM_generate_trajs():
    from transformato.loeffler_systems import mutate_methane_to_methane_cc
    from .test_run_production import run_simulation

    output_files, configuration = mutate_methane_to_methane_cc(
        conf="transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml"
    )
    run_simulation(output_files)


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_methane_to_methane_cc_solvation_free_energy_with_CHARMM_generate_trajs():
    from transformato.loeffler_systems import mutate_methane_to_methane_cc
    from .test_run_production import run_simulation
    from .test_postprocessing import postprocessing

    output_files, configuration = mutate_methane_to_methane_cc(
        conf="transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
        output_dir=".",
    )
    run_simulation(output_files)
    postprocessing(configuration, name="methane", engine="CHARMM")


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_2OJ_tautomer_pair(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_2OJ9_tautomer_pair
    from .test_run_production import run_simulation

    output_files_t1, output_files_t2 = setup_2OJ9_tautomer_pair()
    run_simulation(output_files_t1)
    run_simulation(output_files_t2)
