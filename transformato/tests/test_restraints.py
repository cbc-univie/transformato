"""
Unit, integration and end-to-end testing for restraints.py in transformato

Uses the following markers:

@pytest.mark.restraints_unittest
@pytest.mark.restraints_integrationtest
@pytest.mark.restraints_endtoendtest


Tests are mainly based on 2OJ9
"""
import pytest
import transformato.restraints as tfrs
import transformato.utils as tfut
import yaml
import simtk.openmm
import numpy as np
import os
import sys

TRAFO_DIR="./transformato/"
PATH_2OJ9=f"{TRAFO_DIR}/../data/2OJ9-original/complex/openmm/step3_input.pdb"

@pytest.mark.restraints
@pytest.mark.restraints_unittest
def test_createRestraintsFromConfig():
    
    with open(f"{TRAFO_DIR}/tests/config/test-2oj9-restraints.yaml","r") as stream:
        config=yaml.safe_load(stream)
    
    assert type(config)==dict # checks if config yaml is properly loaded

    restraints=tfrs.CreateRestraintsFromConfig(config,PATH_2OJ9)

    assert type(restraints)==list

    for restraint in restraints:
        assert isinstance(restraint,tfrs.Restraint)


@pytest.mark.restraints
@pytest.mark.restraints_unittest
def test_Restraints():

    testrestraint=tfrs.Restraint("resname BMI and type C","protein and name CA",PATH_2OJ9,14)

    assert isinstance(testrestraint, tfrs.Restraint)

    testrestraint.createForce(["C14","C12","C11","C9"]) # intentionally wrong CC
    print(testrestraint.g1_openmm)
    assert [comp  in testrestraint.g1_openmm for comp in [4822,4825,4826,2828]] #Test proper selection and translation
    

    assert isinstance(testrestraint.force, simtk.openmm.CustomCentroidBondForce)

    assert isinstance(testrestraint.get_force(), simtk.openmm.CustomCentroidBondForce)

@pytest.mark.restraints
@pytest.mark.restraints_unittest
def test_3DDistance():
    assert(tfrs.get3DDistance(np.asarray([1,0,0]),np.asarray([0,0,0])))==1

@pytest.mark.restraints
@pytest.mark.restraints_unittest
def test_write_yaml(tmp_path):

    class MockSystem():
        def __init__(self):
            self.tlc="LIG"
            self.structure="structure2"
    sys.modules["transformato.mutate"].cc_names_struc1=["C1","C2"]
    sys.modules["transformato.mutate"].cc_names_struc2=["C1","C2"]
    path=tmp_path/"test-config.yaml"
    config=tfut.load_config_yaml(f"{TRAFO_DIR}/tests/config/test-2oj9-rsfe-restraints.yaml",".","./tmp/")
    system=MockSystem()
    current_step=4
    tfrs.write_restraints_yaml(path,system,config,current_step)

    assert os.path.exists(path)
