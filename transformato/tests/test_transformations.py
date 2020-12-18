"""
Unit and regression test for the transformato package.
"""

import sys

import numpy as np
import parmed as pm
import pytest


# Import package, test suite, and other packages as needed
import transformato

# read in specific topology with parameters
from transformato import (
    load_config_yaml,
)


def test_transformato_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "transformato" in sys.modules


def test_generate_output_for_methane_cc_solvation_free_energy():
    from transformato.loeffler_systems import mutate_methane_to_methane_cc

    output_files, configuration = mutate_methane_to_methane_cc(
        conf="transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml"
    )
    
    assert len(output_files) == 3


def test_generate_output_for_toluene_cc_solvation_free_energy():
    from transformato.loeffler_systems import mutate_toluene_to_methane_cc

    output_files, configuration = mutate_toluene_to_methane_cc(
        conf="transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml"
    )
    print(output_files)


def test_generate_output_for_neopentane_cc_solvation_free_energy():
    from transformato.loeffler_systems import mutate_neopentane_to_methane_cc

    output_files, configuration = mutate_neopentane_to_methane_cc(
        conf="transformato/tests/config/test-neopentane-methane-solvation-free-energy.yaml"
    )
    print(output_files)


def test_generate_output_for_methanol_cc_solvation_free_energy():
    from transformato.loeffler_systems import mutate_methanol_to_methane_cc

    output_files, _ = mutate_methanol_to_methane_cc(
        conf="transformato/tests/config/test-methanol-methane-solvation-free-energy.yaml"
    )
    print(output_files)
