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
    settingsMap = load_config_yaml('2oj9-test')

    assert(settingsMap['system']['name'] == '2OJ9-2OJ9-mod')


def test_initialize_systems():
    from transformato.utils import load_config_yaml
    conf = load_config_yaml('2oj9-test')
    from transformato.system import SystemStructure

    s1 = SystemStructure(conf, 'structure1')
    assert(int(s1.waterbox_offset) == 0)
    assert(int(s1.complex_offset) == 4811)

    s2 = SystemStructure(conf, 'structure2')
    assert(int(s2.waterbox_offset) == 0)
    assert(int(s2.complex_offset) == 4692)

    assert(s1.envs[0] == 'complex')
    assert(s1.envs[1] == 'waterbox')

