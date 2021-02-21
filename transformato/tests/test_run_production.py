"""
Unit and regression test for the transformato package.
"""

import os
import subprocess
import logging

import numpy as np
import pytest
import shutil


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
    np.isclose(ddG, 10.074696575528037, rtol=1e-5)
    print(ddG)
    shutil.rmtree("toluene-methane-solvation-free-energy")


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
def test_run_acetylacetone_tautomer_pair(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_acetylacetone_tautomer_pair
    from .test_run_production import run_simulation

    (output_files_t1, output_files_t2), _ = setup_acetylacetone_tautomer_pair()
    run_simulation(output_files_t1)
    run_simulation(output_files_t2)
    shutil.rmtree("acetylacetone-keto-acetylacetone-enol-solvation-free-energy")


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_2OJ0_tautomer_pair(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_2OJ9_tautomer_pair
    from .test_run_production import run_simulation

    (output_files_t1, output_files_t2), _ = setup_2OJ9_tautomer_pair()
    run_simulation(output_files_t1)
    run_simulation(output_files_t2)
    shutil.rmtree("2OJ9-original-2OJ9-tautomer-solvation-free-energy")


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_get_free_energy_2OJ9_tautomer_pair(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_2OJ9_tautomer_pair
    from .test_run_production import run_simulation
    from .test_postprocessing import postprocessing

    (output_files_t1, output_files_t2), conf = setup_2OJ9_tautomer_pair()
    run_simulation(output_files_t1)
    run_simulation(output_files_t2)
    ddG1, dddG1 = postprocessing(
        conf, name="2OJ9-original", engine="openMM"
    )  # TODO: correct namings
    ddG2, dddG2 = postprocessing(
        conf, name="2OJ9-tautomer", engine="openMM"
    )  # TODO: correct namings

    ddG = ddG2 - ddG1
    assert np.isclose(ddG, 1.12)


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
    ddG1, dddG1 = postprocessing(
        conf, name="acetylacetone-keto", engine="openMM"
    )  # TODO: correct namings
    ddG2, dddG2 = postprocessing(
        conf, name="acetylacetone-enol", engine="openMM"
    )  # TODO: correct namings

    ddG = ddG2 - ddG1
    assert np.isclose(ddG, 1.12)
    shutil.rmtree("acetylacetone-keto-acetylacetone-enol-solvation-free-energy")
