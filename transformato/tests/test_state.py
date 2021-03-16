"""
Unit and regression test for the transformato package.
"""

import os
import logging

import pytest


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_2OJ9_tautomer_pair_vfswitch(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_2OJ9_tautomer_pair

    conf_path = "transformato/tests/config/test-2oj9-tautomer-pair-solvation-free-energy_vfswitch.yaml"

    (output_files_t1, output_files_t2), _, _ = setup_2OJ9_tautomer_pair(
        single_state=True, conf_path=conf_path
    )


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_2OJ9_tautomer_pair_vswitch(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_2OJ9_tautomer_pair

    conf_path = "transformato/tests/config/test-2oj9-tautomer-pair-solvation-free-energy_vswitch.yaml"

    (output_files_t1, output_files_t2), _, _ = setup_2OJ9_tautomer_pair(
        single_state=True, conf_path=conf_path
    )


@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_run_2OJ9_tautomer_pair_no_switch(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_2OJ9_tautomer_pair

    conf_path = "transformato/tests/config/test-2oj9-tautomer-pair-solvation-free-energy_no-switch.yaml"

    (output_files_t1, output_files_t2), _, _ = setup_2OJ9_tautomer_pair(
        single_state=True, conf_path=conf_path
    )