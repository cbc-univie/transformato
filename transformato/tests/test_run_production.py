"""
Unit and regression test for the transformato package.
"""

import os
import subprocess
import logging

import numpy as np
import pytest
import shutil


def run_simulation(output_files, engine="openMM"):
    for path in sorted(output_files):
        # because path is object not string
        print(f"Start sampling for: {path}")
        runfile = "simulation.sh"
        if engine.upper() == "CHARMM":
            runfile = "simulation_charmm.sh"

        exe = subprocess.run(
            ["bash", f"{str(path)}/{runfile}", str(path)],
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
    from transformato.loeffler_systems import mutate_toluene_to_methane_cc
    from .test_run_production import run_simulation
    from .test_postprocessing import postprocessing

    output_files, configuration = mutate_toluene_to_methane_cc(
        conf="transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml"
    )

    run_simulation(output_files)
    f = "/".join(output_files[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_methane_to_methane_cc_solvation_free_energy_with_openMM():
    from transformato.loeffler_systems import mutate_methane_to_methane_cc
    from .test_run_production import run_simulation
    from .test_postprocessing import postprocessing

    output_files, configuration = mutate_methane_to_methane_cc(
        conf="transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml"
    )
    run_simulation(output_files)
    f = "/".join(output_files[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_methane_to_methane_cc_solvation_free_energy_with_CHARMM():
    from transformato.loeffler_systems import mutate_methane_to_methane_cc
    from .test_run_production import run_simulation
    from .test_postprocessing import postprocessing

    output_files, configuration = mutate_methane_to_methane_cc(
        conf="transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml",
        output_dir=".",
    )
    run_simulation(output_files)
    f = "/".join(output_files[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_acetylacetone_tautomer_pair(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_acetylacetone_tautomer_pair
    from .test_run_production import run_simulation

    (output_files_t1, output_files_t2), _ = setup_acetylacetone_tautomer_pair()
    run_simulation(output_files_t1)
    run_simulation(output_files_t2)
    f = "/".join(output_files_t1[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_2OJ9_tautomer_pair_openMM(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_2OJ9_tautomer_pair
    from .test_run_production import run_simulation

    conf_path = (
        "transformato/tests/config/test-2oj9-tautomer-pair-solvation-free-energy.yaml"
    )

    (output_files_t1, output_files_t2), _, _ = setup_2OJ9_tautomer_pair(
        conf_path=conf_path
    )
    run_simulation(output_files_t1)
    run_simulation(output_files_t2)
    f = "/".join(output_files_t1[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_2OJ9_tautomer_pair_charmm(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_2OJ9_tautomer_pair
    from .test_run_production import run_simulation

    conf_path = (
        "transformato/tests/config/test-2oj9-tautomer-pair-solvation-free-energy.yaml"
    )

    (output_files_t1, output_files_t2), _, _ = setup_2OJ9_tautomer_pair(
        conf_path=conf_path
    )
    run_simulation(output_files_t1, engine="CHARMM")
    run_simulation(output_files_t2, engine="CHARMM")
    f = "/".join(output_files_t1[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_get_free_energy_acetylaceton_tautomer_pair(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_acetylacetone_tautomer_pair
    from .test_run_production import run_simulation
    from .test_postprocessing import postprocessing

    (output_files_t1, output_files_t2), conf = setup_acetylacetone_tautomer_pair()
    run_simulation(output_files_t1)
    run_simulation(output_files_t2)
    f = "/".join(output_files_t1[0].split("/")[:-3])
    print(f)
    shutil.rmtree(f)
