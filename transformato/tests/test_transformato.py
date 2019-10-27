"""
Unit and regression test for the transformato package.
"""

# Import package, test suite, and other packages as needed
import transformato
import pytest
import sys
import logging

def test_transformato_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "transformato" in sys.modules


def test_read_yaml():

    from transformato.utils import load_config_yaml
    settingsMap = load_config_yaml(config='config/2oj9-test.yaml',
                       data_dir='.', output_dir='data/')

    assert(settingsMap['system']['name'] == '2OJ9-2OJ9-mod')


def test_initialize_systems():
    from transformato.utils import load_config_yaml
    conf = load_config_yaml(config='config/2oj9-test.yaml',
                       data_dir='data/', output_dir='.')
    from transformato.system import SystemStructure

    s1 = SystemStructure(conf, 'structure1')
    assert(int(s1.waterbox_offset) == 0)
    assert(int(s1.complex_offset) == 4811)

    s2 = SystemStructure(conf, 'structure2')
    assert(int(s2.waterbox_offset) == 0)
    assert(int(s2.complex_offset) == 4692)

    assert(s1.envs[0] == 'complex')
    assert(s1.envs[1] == 'waterbox')

def test_proposed_mutation():
    from transformato import load_config_yaml, SystemStructure, ProposeMutationRoute
    conf = load_config_yaml(config='config/2oj9-test.yaml',
                       data_dir='data/', output_dir='.')
    s1 = SystemStructure(conf, 'structure1')
    s2 = SystemStructure(conf, 'structure2')

    a = ProposeMutationRoute(s1, s2)
    cc1 = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 41, 44, 45, 46, 47, 48]
    assert(a.get_common_core_idx_mol1() == cc1)
    assert(str(a.s1_tlc) == 'BMI')
    

