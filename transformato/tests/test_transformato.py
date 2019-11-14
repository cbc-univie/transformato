"""
Unit and regression test for the transformato package.
"""

# Import package, test suite, and other packages as needed
import transformato
import pytest
import sys
import logging
import shutil
from transformato import load_config_yaml, SystemStructure, ProposeMutationRoute, IntermediateStateFactory
import parmed as pm
import copy
import numpy as np
# read in specific topology with parameters
from parmed.charmm.parameters import CharmmParameterSet


def read_params(output_file_base):
    extlist = ['rtf', 'prm', 'str']
    print(output_file_base)
    parFiles = ()
    toppar_file = f"{output_file_base}/toppar.str"
    for line in open(toppar_file, 'r'):
        if '!' in line: line = line.split('!')[0]
        parfile = line.strip()
        if len(parfile) != 0:
            ext = parfile.lower().split('.')[-1]
            if not ext in extlist: continue
            parFiles += ( f"{output_file_base}/{parfile}", )

    params = CharmmParameterSet( *parFiles )
    return params

def generate_psf(output_file_base, env):
    parms = read_params(output_file_base)  
    target_psf = pm.charmm.CharmmPsfFile(f"{output_file_base}/lig_in_{env}.psf")
    target_psf.load_parameters(parms)
    return target_psf


def test_transformato_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "transformato" in sys.modules


def test_read_yaml():
    settingsMap = load_config_yaml(config='config/2oj9-test.yaml',
                       input_dir='.', output_dir='data/')

    assert(settingsMap['system']['name'] == '2OJ9-test1-2OJ9-test2')


def test_initialize_systems():
    configuration = load_config_yaml(config='config/2oj9-test.yaml',
                       input_dir='data/', output_dir='.')

    s1 = SystemStructure(configuration, 'structure1')
    assert(int(s1.waterbox_offset) == 0)
    assert(int(s1.complex_offset) == 4811)

    s2 = SystemStructure(configuration, 'structure2')
    assert(int(s2.waterbox_offset) == 0)
    assert(int(s2.complex_offset) == 4692)

    assert(s1.envs[0] == 'complex')
    assert(s1.envs[1] == 'waterbox')

def test_proposed_mutation():
    configuration = load_config_yaml(config='config/2oj9-test.yaml',
                       input_dir='data/', output_dir='.')
    s1 = SystemStructure(configuration, 'structure1')
    s2 = SystemStructure(configuration, 'structure2')

    a = ProposeMutationRoute(s1, s2)
    cc1 = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 41, 44, 45, 46, 47, 48]
    assert(a.get_common_core_idx_mol1() == cc1)
    assert(str(a.s1_tlc) == 'BMI')
    

def test_endpoint():

    configuration = load_config_yaml(config='config/2oj9-test.yaml',
                    input_dir='data/', output_dir='.')
    s1 = SystemStructure(configuration, 'structure1')
    s2 = SystemStructure(configuration, 'structure2')

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s2_to_s1 = ProposeMutationRoute(s2, s1)
    for a, system in zip([s1_to_s2, s2_to_s1], [s1, s2]):
        mutation_list = a.generate_mutations_to_common_core_for_mol1(nr_of_steps_for_el=3, nr_of_steps_for_bonded_parameters=3)
        i = IntermediateStateFactory(system=system, mutation_list=mutation_list, configuration=configuration)

        # extract first mutation
        m = mutation_list[0]
        current_step = 0
        # write out structure
        output_file_base = i.generate_specific_intermediate_state(m, current_step)

        # get a copy of the original structure
        original_psf_waterbox = copy.deepcopy(system.waterbox_psf)
        original_psf_complex = copy.deepcopy(system.complex_psf)

        # test for waterbox and complex
        for env in ['waterbox', 'complex']:
            # read in psf and parameters
            new_psf = generate_psf(output_file_base, env)
            if env == 'waterbox':
                original_psf = original_psf_waterbox
            else:
                original_psf = original_psf_complex
            offset = min([a.idx for a in original_psf.view[f":{system.tlc.upper()}"].atoms])


            # test if atom parameters are the same
            for idx in m.atom_idx:
                assert(original_psf.atoms[idx+offset].epsilon  == new_psf.atoms[idx+offset].epsilon)
                assert(original_psf.atoms[idx+offset].charge  == new_psf.atoms[idx+offset].charge)
                assert(original_psf.atoms[idx+offset].rmin  == new_psf.atoms[idx+offset].rmin)
                
                # test all bond parameters
                for o_bond, n_bond in zip(original_psf.atoms[idx+offset].bonds, new_psf.atoms[idx+offset].bonds):
                    assert(o_bond.type.k  == n_bond.type.k)
                
                # test all angle parameters
                for o_angle, n_angle in zip(original_psf.atoms[idx+offset].angles, new_psf.atoms[idx+offset].angles):
                    assert(o_angle.type.k  == n_angle.type.k)
        
        shutil.rmtree(output_file_base) 
                
def test_charge_mutation():
    
    configuration = load_config_yaml(config='config/2oj9-test.yaml',
                    input_dir='data/', output_dir='.')
    s1 = SystemStructure(configuration, 'structure1')
    s2 = SystemStructure(configuration, 'structure2')

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s2_to_s1 = ProposeMutationRoute(s2, s1)
    for a, system in zip([s1_to_s2, s2_to_s1], [s1, s2]):
        mutation_list = a.generate_mutations_to_common_core_for_mol1(nr_of_steps_for_el=4, nr_of_steps_for_bonded_parameters=4)
        i = IntermediateStateFactory(system=system, mutation_list=mutation_list, configuration=configuration)

        # test endpoint - charge not scaled
        m = mutation_list[0]
        current_step = 0
        output_file_base = i.generate_specific_intermediate_state(m, current_step)
        # get a copy of the original structure
        original_psf_waterbox = copy.deepcopy(system.waterbox_psf)
        original_psf_complex = copy.deepcopy(system.complex_psf)
        
        for env in ['waterbox', 'complex']:
            if env == 'waterbox':
                original_psf = original_psf_waterbox
            else:
                original_psf = original_psf_complex
            offset = min([a.idx for a in original_psf.view[f":{system.tlc.upper()}"].atoms])

            new_psf = generate_psf(output_file_base, env)
            scale = current_step/(m.nr_of_steps)
            for idx in m.atom_idx:
                assert(np.isclose(original_psf.atoms[idx+offset].charge * (1 - scale), new_psf.atoms[idx+offset].charge))
                assert(np.isclose(original_psf.atoms[idx+offset].charge, new_psf.atoms[idx+offset].charge))

        shutil.rmtree(output_file_base) 

        # scale charges with 0.25
        current_step = 1
        output_file_base = i.generate_specific_intermediate_state(m, current_step)
        for env in ['waterbox', 'complex']:
            if env == 'waterbox':
                original_psf = original_psf_waterbox
            else:
                original_psf = original_psf_complex

            offset = min([a.idx for a in original_psf.view[f":{system.tlc.upper()}"].atoms])
            new_psf = generate_psf(output_file_base, env)
            scale = current_step/(m.nr_of_steps)
            print('Scaling charges with: {}'.format((1 - scale)))
            for idx in m.atom_idx:
                assert(np.isclose(original_psf.atoms[idx+offset].charge * (1 - scale), new_psf.atoms[idx+offset].charge, rtol=1e-03))
        shutil.rmtree(output_file_base) 

        # scale charges with 1.0
        current_step = 4
        output_file_base = i.generate_specific_intermediate_state(m, current_step)
        for env in ['waterbox', 'complex']:
            if env == 'waterbox':
                original_psf = original_psf_waterbox
            else:
                original_psf = original_psf_complex

            offset = min([a.idx for a in original_psf.view[f":{system.tlc.upper()}"].atoms])
            new_psf = generate_psf(output_file_base, env)
            scale = current_step/(m.nr_of_steps)
            print('Scaling charges with: {}'.format((1 - scale)))
            for idx in m.atom_idx:
                idxo = idx+offset
                assert(np.isclose(original_psf.atoms[idxo].charge * (1 - scale), new_psf.atoms[idxo].charge, rtol=1e-03))
                assert(np.isclose(0.0, new_psf.atoms[idxo].charge))
        shutil.rmtree(output_file_base) 

def test_vdw_mutation():
    
    configuration = load_config_yaml(config='config/2oj9-test.yaml',
                    input_dir='data/', output_dir='.')
    s1 = SystemStructure(configuration, 'structure1')
    s2 = SystemStructure(configuration, 'structure2')

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s2_to_s1 = ProposeMutationRoute(s2, s1)
    for a, system in zip([s1_to_s2, s2_to_s1], [s1, s2]):
        mutation_list = a.generate_mutations_to_common_core_for_mol1(nr_of_steps_for_el=4, nr_of_steps_for_bonded_parameters=4)
        i = IntermediateStateFactory(system=system, mutation_list=mutation_list, configuration=configuration)

        m1 = mutation_list[1]
        current_step = 0
        output_file_base = i.generate_specific_intermediate_state(m1, current_step)
        # get a copy of the original structure
        original_psf_waterbox = copy.deepcopy(system.waterbox_psf)
        original_psf_complex = copy.deepcopy(system.complex_psf)

        for env in ['waterbox', 'complex']:
            
            if env == 'waterbox':
                original_psf = original_psf_waterbox
            else:
                original_psf = original_psf_complex
            offset = min([a.idx for a in original_psf.view[f":{system.tlc.upper()}"].atoms])
            
            new_psf = generate_psf(output_file_base, env)
            print('Set epsilon/rmin to zero for selected atoms')
            for idx in m1.atom_idx:
                idxo = idx + offset
                assert(np.isclose(original_psf.atoms[idxo].epsilon * 0, new_psf.atoms[idxo].epsilon))
                assert(np.isclose(original_psf.atoms[idxo].rmin * 0, new_psf.atoms[idxo].rmin))

            # make sure that all other idx are not touched
            for idx in range(len(original_psf.atoms)):
                idxo = idx - offset # NOTE: the '-'
                if idxo not in m1.atom_idx:
                    assert(np.isclose(original_psf.atoms[idx].epsilon, new_psf.atoms[idx].epsilon))
                    assert(np.isclose(original_psf.atoms[idx].rmin, new_psf.atoms[idx].rmin))

        shutil.rmtree(output_file_base) 

        m2 = mutation_list[2]
        current_step = 1
        output_file_base = i.generate_specific_intermediate_state(m2, current_step)
        for env in ['waterbox', 'complex']:
            if env == 'waterbox':
                original_psf = original_psf_waterbox
            else:
                original_psf = original_psf_complex
            
            offset = min([a.idx for a in original_psf.view[f":{system.tlc.upper()}"].atoms])
            new_psf = generate_psf(output_file_base, env)
            print('Set epsilon/rmin to zero for selected atoms')
            for idx in m2.atom_idx:
                idxo = idx + offset
                assert(np.isclose(original_psf.atoms[idxo].epsilon * 0, new_psf.atoms[idxo].epsilon))
                assert(np.isclose(original_psf.atoms[idxo].rmin * 0, new_psf.atoms[idxo].rmin))

            # make sure that all other idx are not touched
            print(m2.atom_idx + m1.atom_idx)
            for idx in range(len(original_psf.atoms)):
                idxo = idx - offset # NOTE: the '-'
                if idxo not in m2.atom_idx + m1.atom_idx:
                    assert(np.isclose(original_psf.atoms[idx].epsilon, new_psf.atoms[idx].epsilon))
                    assert(np.isclose(original_psf.atoms[idx].rmin, new_psf.atoms[idx].rmin))

        shutil.rmtree(output_file_base) 


def test_bonded_mutation():

    configuration = load_config_yaml(config='config/2oj9-test.yaml',
                    input_dir='data/', output_dir='.')
    s1 = SystemStructure(configuration, 'structure1')
    s2 = SystemStructure(configuration, 'structure2')

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s2_to_s1 = ProposeMutationRoute(s2, s1)


    for a, system, template in zip([s1_to_s2, s2_to_s1], [SystemStructure(configuration, 'structure1'), SystemStructure(configuration, 'structure2')], 
                                [SystemStructure(configuration, 'structure2'), SystemStructure(configuration, 'structure1')]):
        
        original_psf_waterbox = copy.deepcopy(system.waterbox_psf)
        original_psf_complex = copy.deepcopy(system.complex_psf)
        template_psf_waterbox = copy.deepcopy(template.waterbox_psf)
        template_psf_complex = copy.deepcopy(template.complex_psf)
        
        mutation_list = a.generate_mutations_to_common_core_for_mol1(nr_of_steps_for_el=4, nr_of_steps_for_bonded_parameters=4)
        i = IntermediateStateFactory(system=system, mutation_list=mutation_list, configuration=configuration)

        print(mutation_list)
        # we select the bonded mutation which should be at postition [-2]
        m1 = mutation_list[-2]
        assert(type(m1) == transformato.mutate.BondedParameterMutation)
        current_step = 1
        output_file_base = i.generate_specific_intermediate_state(m1, current_step)

        
        for env in ['waterbox', 'complex']:
            print('Env : {}'.format(env))
            if env == 'waterbox':
                original_psf = original_psf_waterbox
                template_psf = template_psf_waterbox
            else:
                original_psf = original_psf_complex
                template_psf = template_psf_complex
            
            new_psf = generate_psf(output_file_base, env)
            
            
            # atom names are the same in new and orginal psf
            for natom, oatom in zip(new_psf.view[f":{system.tlc.upper()}"].atoms, 
                                    original_psf.view[f":{system.tlc.upper()}"].atoms):
                
                assert(natom.name== oatom.name)

            # changes are only in the ligand
            for natom, oatom in zip(new_psf.view[f"!(:{system.tlc.upper()})"].atoms, 
                                    original_psf.view[f"!(:{system.tlc.upper()})"].atoms):
                
                assert(np.isclose(natom.charge, oatom.charge))
                assert(np.isclose(natom.epsilon, oatom.epsilon))
            
            
            # get mapping between original/new and template psf
            match_atom_names_cc1_to_cc2 = dict()
            cc1_offset = min([a.idx for a in original_psf.view[f":{m1.tlc_cc1.upper()}"].atoms])
            cc2_offset = min([a.idx for a in template_psf.view[f":{m1.tlc_cc2.upper()}"].atoms])
            scale = current_step/(m1.nr_of_steps)

            # test atom parameters
            for cc1, cc2 in zip(m1.cc1_idx, m1.cc2_idx):
                # did atom type change? if not continue

                cc1_oidx = cc1 + cc1_offset
                cc2_oidx = cc2 + cc2_offset
                cc1_a = original_psf[cc1_oidx]
                cc2_a = template_psf[cc2_oidx]
                if cc1_a.type == cc2_a.type:
                    continue

                match_atom_names_cc1_to_cc2[cc1_a.name] = cc2_a.name
                                    
                
                assert(np.isclose((1.0 - scale) * cc1_a.epsilon + scale * cc2_a.epsilon, 
                                new_psf[cc1_oidx].epsilon))
                assert(np.isclose((1.0 - scale) * cc1_a.sigma + scale * cc2_a.sigma, 
                                new_psf[cc1_oidx].sigma, rtol=1e-03))
                
            
            # get mapping between original/new and template psf
            for cc1_bond, new_bond in zip(original_psf.view[f":{m1.tlc_cc1.upper()}"].bonds, 
                                        new_psf.view[f":{m1.tlc_cc1.upper()}"].bonds):
                cc1_a1 = cc1_bond.atom1.name
                cc1_a2 = cc1_bond.atom2.name
                # all atoms of the bond must be in cc
                # everything outside the cc are bonded terms between dummies or 
                # between real atoms and dummies and we can ignore them for now
                if not all(elem in match_atom_names_cc1_to_cc2 for elem in [cc1_a1, cc1_a2]):
                    assert(np.isclose(cc1_bond.type.k, new_bond.type.k))
                    continue

                for cc2_bond in template_psf.view[f":{m1.tlc_cc2.upper()}"].bonds:
                    cc2_a1 = cc2_bond.atom1.name
                    cc2_a2 = cc2_bond.atom2.name
                    # all atoms of the bond must be in cc
                    if not all(elem in match_atom_names_cc1_to_cc2.values() for elem in [cc2_a1, cc2_a2]):
                        continue

                    # match the two bonds
                    if sorted([match_atom_names_cc1_to_cc2[e] for e in [cc1_a1, cc1_a2]]) == sorted([cc2_a1, cc2_a2]):
                        scaled = (1.0 - scale) * cc1_bond.type.k + scale * cc2_bond.type.k
                        assert(np.isclose(scaled, new_bond.type.k))
            
            # make sure everything else in not changed
            for cc1_bond, new_bond in zip(original_psf.view[f"!:{m1.tlc_cc1.upper()}"].bonds, 
                                        new_psf.view[f"!:{m1.tlc_cc1.upper()}"].bonds):
                assert(np.isclose(cc1_bond.type.k, new_bond.type.k))

            
            
            # get mapping between original/new and template psf
            for cc1_angle, new_angle in zip(original_psf.view[f":{m1.tlc_cc1.upper()}"].angles, 
                                        new_psf.view[f":{m1.tlc_cc1.upper()}"].angles):
                cc1_a1 = cc1_angle.atom1.name
                cc1_a2 = cc1_angle.atom2.name
                cc1_a3 = cc1_angle.atom3.name
                # all atoms of the bond must be in cc
                # everything outside the cc are bonded terms between dummies or 
                # between real atoms and dummies and we can ignore them for now
                if not all(elem in match_atom_names_cc1_to_cc2 for elem in [cc1_a1, cc1_a2, cc1_a3]):
                    assert(np.isclose(cc1_angle.type.k, new_angle.type.k))
                    continue
                
                for cc2_angle in template_psf.view[f":{m1.tlc_cc2.upper()}"].angles:
                    cc2_a1 = cc2_angle.atom1.name
                    cc2_a2 = cc2_angle.atom2.name
                    cc2_a3 = cc2_angle.atom3.name
                    # all atoms of the bond must be in cc
                    if not all(elem in match_atom_names_cc1_to_cc2.values() for elem in [cc2_a1, cc2_a2, cc2_a3]):
                        continue

                    # match the two bonds
                    if sorted([match_atom_names_cc1_to_cc2[e] for e in [cc1_a1, cc1_a2, cc1_a3]]) == sorted([cc2_a1, cc2_a2, cc2_a3]):
                        scaled = (1.0 - scale) * cc1_angle.type.k + scale * cc2_angle.type.k
                        assert(np.isclose(scaled, new_angle.type.k))

    shutil.rmtree(output_file_base) 
