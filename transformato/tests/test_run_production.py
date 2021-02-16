"""
Unit and regression test for the transformato package.
"""

import os
import subprocess
import logging

import numpy as np
import pytest
from transformato import (
    ProposeMutationRoute,
    SystemStructure,
    load_config_yaml,
)


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_toluene_to_methane_cc_solvation_free_energy_with_openMM():
    from transformato import FreeEnergyCalculator
    from transformato.loeffler_systems import mutate_toluene_to_methane_cc

    output_files, configuration = mutate_toluene_to_methane_cc(
        conf="transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml"
    )

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

    f = FreeEnergyCalculator(configuration, "toluene")
    f.load_trajs(nr_of_max_snapshots=100)
    f.calculate_dG_to_common_core()
    ddG, dddG = f.end_state_free_energy_difference
    print(f"Free energy difference: {ddG}")
    print(f"Uncertanty: {dddG}")
    np.isclose(ddG, 8.9984, rtol=1e-2)
    f.show_summary()


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_methane_to_methane_cc_solvation_free_energy_with_openMM():
    from transformato import FreeEnergyCalculator
    from transformato.loeffler_systems import mutate_methane_to_methane_cc

    output_files, configuration = mutate_methane_to_methane_cc(
        conf="transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml"
    )

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

    f = FreeEnergyCalculator(configuration, "methane")
    f.load_trajs(nr_of_max_snapshots=300)
    f.calculate_dG_to_common_core()
    ddG, dddG = f.end_state_free_energy_difference
    print(f"Free energy difference: {ddG}")
    print(f"Uncertanty: {dddG}")
    np.isclose(ddG, 8.9984, rtol=1e-2)
    f.show_summary()


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_methane_to_methane_cc_solvation_free_energy_with_openMM_generate_trajs():
    from transformato.loeffler_systems import mutate_methane_to_methane_cc

    output_files, configuration = mutate_methane_to_methane_cc(
        conf="transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml"
    )
    print(output_files)
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
def test_run_methane_to_methane_cc_solvation_free_energy_with_CHARMM_generate_trajs():
    from transformato.loeffler_systems import mutate_methane_to_methane_cc
    from transformato import FreeEnergyCalculator

    output_files, configuration = mutate_methane_to_methane_cc(
        conf="transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
        output_dir=".",
    )

    for path in sorted(output_files):
        # because path is object not string
        print(f"Start sampling for: {path}")
        exe = subprocess.run(
            ["bash", f"{str(path)}/simulation.sh", str(path)],
            check=True,
            capture_output=True,
            text=True,
        )

    f = FreeEnergyCalculator(configuration, "methane")
    f.load_trajs(nr_of_max_snapshots=300)
    f.calculate_dG_to_common_core(engine="CHARMM")
    ddG, dddG = f.end_state_free_energy_difference
    print(f"Free energy difference: {ddG}")
    print(f"Uncertanty: {dddG}")
    f.show_summary()


@pytest.mark.slowtest
def test_run_2OJ_tautomer_pair(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_2OJ9_tautomer_pair

    output_files = setup_2OJ9_tautomer_pair()
    output_files_t1 = output_files[0]
    for path in sorted(output_files_t1):
        # because path is object not string
        print(f"Start sampling for: {path}")
        exe = subprocess.run(
            ["bash", f"{str(path)}/simulation.sh", str(path)],
            check=True,
            capture_output=True,
            text=True,
        )
        print(exe)