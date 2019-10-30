"""
Unit and regression test for the transformato package.
"""

# Import package, test suite, and other packages as needed
import transformato
import pytest
import sys
import logging

def read_params(filename):
    from parmed.charmm.parameters import CharmmParameterSet
    extlist = ['rtf', 'prm', 'str']
    print(filename)
    parFiles = ()
    toppar_file = f"{filename}/toppar.str"
    for line in open(toppar_file, 'r'):
        if '!' in line: line = line.split('!')[0]
        parfile = line.strip()
        if len(parfile) != 0:
            ext = parfile.lower().split('.')[-1]
            if not ext in extlist: continue
            parFiles += ( f"{filename}/{parfile}", )

    params = CharmmParameterSet( *parFiles )
    return params


def test_transformato_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "transformato" in sys.modules


def test_read_yaml():

    from transformato.utils import load_config_yaml
    settingsMap = load_config_yaml(config='config/2oj9-test.yaml',
                       input_dir='.', output_dir='data/')

    assert(settingsMap['system']['name'] == '2OJ9-2OJ9-mod')


def test_initialize_systems():
    from transformato.utils import load_config_yaml
    configuration = load_config_yaml(config='config/2oj9-test.yaml',
                       input_dir='data/', output_dir='.')
    from transformato.system import SystemStructure

    s1 = SystemStructure(configuration, 'structure1')
    assert(int(s1.waterbox_offset) == 0)
    assert(int(s1.complex_offset) == 4811)

    s2 = SystemStructure(configuration, 'structure2')
    assert(int(s2.waterbox_offset) == 0)
    assert(int(s2.complex_offset) == 4692)

    assert(s1.envs[0] == 'complex')
    assert(s1.envs[1] == 'waterbox')

def test_proposed_mutation():
    from transformato import load_config_yaml, SystemStructure, ProposeMutationRoute
    configuration = load_config_yaml(config='config/2oj9-test.yaml',
                       input_dir='data/', output_dir='.')
    s1 = SystemStructure(configuration, 'structure1')
    s2 = SystemStructure(configuration, 'structure2')

    a = ProposeMutationRoute(s1, s2)
    cc1 = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 41, 44, 45, 46, 47, 48]
    assert(a.get_common_core_idx_mol1() == cc1)
    assert(str(a.s1_tlc) == 'BMI')
    

def test_specific_mutations():
    from transformato import load_config_yaml, SystemStructure, ProposeMutationRoute, IntermediateStateFactory
    configuration = load_config_yaml(config='config/2oj9-test.yaml',
                       input_dir='data/', output_dir='.')
    s1 = SystemStructure(configuration, 'structure1')
    s2 = SystemStructure(configuration, 'structure2')

    a = ProposeMutationRoute(s1, s2)
    mutation_list = a.generate_mutations_to_common_core_for_mol1(nr_of_steps_for_el=2, nr_of_steps_for_bonded_parameters=2)
    i = IntermediateStateFactory(system=s1, mutation_list=mutation_list, configuration=configuration)
    m = mutation_list[0]
    output_file_base = i.generate_specific_intermediate_state(m, 0)
    parms = read_params(output_file_base)
