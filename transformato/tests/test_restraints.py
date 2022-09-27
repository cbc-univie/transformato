"""
Unit, integration and end-to-end testing for restraints.py in transformato

Uses the following markers:

@pytest.mark.restraints_unittest
@pytest.mark.restraints_integrationtest
@pytest.mark.restraints_endtoendtest


Tests are mainly based on 2OJ9
"""
import glob
import logging
import os
import sys
import warnings

import numpy as np
import openmm
import pytest
import transformato.restraints as tfrs
import transformato.utils as tfut
import yaml
from transformato_testsystems.testsystems import get_testsystems_dir

PATH_2OJ9 = f"{get_testsystems_dir()}/2OJ9-original/complex/openmm/step3_input.pdb"
PATH_2OJ9_DIR = f"{get_testsystems_dir()}/2OJ9-original/complex/openmm/"

# Disable useless parmed warnings for structure generation

warnings.filterwarnings("ignore", module="parmed")
logger = logging.getLogger(__name__)


@pytest.mark.restraints
@pytest.mark.restraints_unittest
def test_create_restraints_from_config():

    with open(
        f"{get_testsystems_dir()}/config/test-2oj9-restraints.yaml", "r"
    ) as stream:
        config = yaml.safe_load(stream)

    assert type(config) == dict  # checks if config yaml is properly loaded
    restraints = tfrs.create_restraints_from_config(config, PATH_2OJ9)
    assert type(restraints) == list
    for restraint in restraints:
        assert isinstance(restraint, tfrs.Restraint)


@pytest.mark.restraints
@pytest.mark.restraints_unittest
def test_restraints():

    testrestraint = tfrs.Restraint(
        "resname BMI and type C", "protein and name CA", PATH_2OJ9, 14
    )

    testrestraint_fb = tfrs.Restraint(
        "resname BMI and type C",
        "protein and name CA",
        PATH_2OJ9,
        14,
        shape="flatbottom",
        wellsize=0.12,
    )
    assert isinstance(testrestraint, tfrs.Restraint)
    assert isinstance(testrestraint, tfrs.Restraint)

    testrestraint.createForce(["C14", "C12", "C11", "C9"])  # intentionally wrong CC
    testrestraint_fb.createForce(["C14", "C12", "C11", "C9"])
    print(testrestraint.g1_openmm)
    assert [
        comp in testrestraint.g1_openmm for comp in [4822, 4825, 4826, 2828]
    ]  # Test proper selection and translation

    assert isinstance(testrestraint.force, openmm.CustomCentroidBondForce)
    assert isinstance(testrestraint.get_force(), openmm.CustomCentroidBondForce)
    assert isinstance(testrestraint_fb.force, openmm.CustomCentroidBondForce)
    assert isinstance(testrestraint_fb.get_force(), openmm.CustomCentroidBondForce)


@pytest.mark.restraints
@pytest.mark.restraints_unittest
def test_3D_distance():
    assert (tfrs.get3DDistance(np.asarray([1, 0, 0]), np.asarray([0, 0, 0]))) == 1


@pytest.mark.restraints
@pytest.mark.restraints_unittest
def test_write_yaml(tmp_path):
    class MockSystem:
        def __init__(self):
            self.tlc = "LIG"
            self.structure = "structure2"

    sys.modules["transformato.mutate"].cc_names_struc1 = ["C1", "C2"]
    sys.modules["transformato.mutate"].cc_names_struc2 = ["C1", "C2"]
    path = tmp_path / "test-config.yaml"
    config = tfut.load_config_yaml(
        f"{get_testsystems_dir()}/config/test-2oj9-rsfe-restraints.yaml", ".", "./tmp/"
    )
    system = MockSystem()
    current_step = 4
    tfrs.write_restraints_yaml(path, system, config, current_step)

    assert os.path.exists(path)


@pytest.mark.restraints
@pytest.mark.restraints_integrationtests
def test_integration():
    """
    Full scale integration test of automatic and manual restraints, including an openMM test system.
    Essentially a modified openmm_run.py
    """
    from openmm import LangevinIntegrator, Platform, CustomCentroidBondForce
    from openmm.app import CharmmCrdFile, CharmmParameterSet, CharmmPsfFile, Simulation
    from openmm.unit import kelvin, mole, nanometers, picoseconds, kilojoule
    from .restraint_helper_functions import (
        read_inputs,
        gen_box,
        vfswitch,
        barostat,
        restraints,
    )
    import warnings

    warnings.filterwarnings("ignore", module="parmed")

    inputs = read_inputs(f"{PATH_2OJ9_DIR}step5_production.inp")
    params = CharmmParameterSet(
        *(
            glob.glob(f"{PATH_2OJ9_DIR}../*/*.prm")
            + glob.glob(f"{PATH_2OJ9_DIR}../*/*.str")
            + glob.glob(f"{PATH_2OJ9_DIR}../*/*.rtf")
            + glob.glob(f"{PATH_2OJ9_DIR}../toppar/*.prm")
        )
    )
    print(params)
    psf = CharmmPsfFile(f"{PATH_2OJ9_DIR}step3_input.psf")
    crd = CharmmCrdFile(f"{PATH_2OJ9_DIR}step3_input.crd")

    top = gen_box(psf, crd)

    # Build system
    nboptions = dict(
        nonbondedMethod=inputs.coulomb,
        nonbondedCutoff=inputs.r_off * nanometers,
        constraints=inputs.cons,
        ewaldErrorTolerance=inputs.ewald_Tol,
    )

    if inputs.vdw == "Switch":
        nboptions["switchDistance"] = inputs.r_on * nanometers
    system = psf.createSystem(params, **nboptions)
    if inputs.vdw == "Force-switch":
        system = vfswitch(system, psf, inputs)

    system = barostat(system, inputs)
    if inputs.rest == "yes":
        system = restraints(system, crd, inputs)
    integrator = LangevinIntegrator(
        inputs.temp * kelvin, inputs.fric_coeff / picoseconds, inputs.dt * picoseconds
    )

    # Set platform

    platform = Platform.getPlatformByName("CPU")

    prop = dict()
    pdbpath = PATH_2OJ9
    with open(
        f"{get_testsystems_dir()}/config/test-2oj9-restraints.yaml", "r"
    ) as stream:
        try:
            configuration = yaml.safe_load(stream)
            logger.debug(configuration)
        except yaml.YAMLError as exc:
            print(exc)

    cc_names = configuration["system"]["structure"]["ccs"]

    # Add forces via transformato.restraints

    if not os.path.exists(pdbpath):
        raise FileNotFoundError(
            f"Couldnt find {pdbpath} necessary for Restraint Analysis"
        )

    restraintList = tfrs.create_restraints_from_config(configuration, pdbpath)

    # Test an additional, simple restraint
    logger.debug("generating simple selection")
    selstr = tfrs.generate_simple_selection(configuration, pdbpath)
    tlc = configuration["system"]["structure"]["tlc"]
    restraintList.append(tfrs.Restraint(f"resname {tlc} and type C", selstr, pdbpath))

    logger.debug(
        "****************** ALL RESTRAINTS CREATED SUCCESSFULLY ***************************"
    )
    num_standard_forces = len(system.getForces())
    for restraint in restraintList:
        restraint.createForce(cc_names)
    logger.debug("FORCES CREATED")
    for restraint in restraintList:
        restraint.applyForce(system)
    logger.debug("FORCES APPLIED")
    for i, f in enumerate(system.getForces()):
        f.setForceGroup(i)

    forcesinsystem = system.getForces()

    assert len(forcesinsystem) == num_standard_forces + len(restraintList)

    logger.debug(f"Number of forces in system: {len(forcesinsystem)}")
    simulation = Simulation(top.topology, system, integrator, platform, prop)
    simulation.context.setPositions(crd.positions)

    logger.info(
        f"Potential energy {simulation.context.getState(getEnergy=True).getPotentialEnergy()}"
    )

    simulation.minimizeEnergy(
        tolerance=inputs.mini_Tol * kilojoule / mole, maxIterations=10
    )
    logger.info(
        f"Potential energy after 1-iter minimization {simulation.context.getState(getEnergy=True).getPotentialEnergy()}"
    )
    for i, f in enumerate(system.getForces()):
        if isinstance(f, CustomCentroidBondForce):
            state = simulation.context.getState(getEnergy=True, groups={i})
            logger.info("Force contributions before steps:")
            logger.info(f"{f}::{state.getPotentialEnergy()}")
            if (
                f.getEnergyFunction()
                == "step(abs(distance(g1,g2)-r0)-w)*k*(distance(g1,g2)-r0)^2"
            ):  # flatbottom-twosided, has set r0
                assert (
                    0.9 in f.getBondParameters(0)[1]
                )  # assert the manually set r0 value is found in the bond parameters

    logger.info("Simulation stepping")
    simulation.step(9)
    temp_results = []
    for i, f in enumerate(system.getForces()):
        if isinstance(f, CustomCentroidBondForce):
            state = simulation.context.getState(getEnergy=True, groups={i})
            logger.info("Force contributions after 10 steps:")
            logger.info(f"{f}::{state.getPotentialEnergy()}")
            temp_results.append(state.getPotentialEnergy())
    simulation.step(10)
    for i, f in enumerate(system.getForces()):
        if isinstance(f, CustomCentroidBondForce):
            state = simulation.context.getState(getEnergy=True, groups={i})
            logger.info("Force contributions after 20 steps:")
            logger.info(f"{f}::{state.getPotentialEnergy()}")

    assert all(
        [results < 0.1 * kilojoule / mole for results in temp_results]
    )  # assert all force contributions at this point are very small
    assert (
        0 * kilojoule / mole in temp_results
    )  # assert there is at least one flatbottom with 0 kJ/mol

    # Change position of r0 and check force contributions now
    # Both harmonic and flatbottom potentials should be off the charts

    logger.info("Changing Bond Parameters to ensure restraint activation")

    for i, f in enumerate(system.getForces()):
        if isinstance(f, CustomCentroidBondForce):
            logger.debug(f.getBondParameters(0))
            if (
                f.getEnergyFunction() == "0.5*k*(distance(g1,g2)-r0)^2"
            ):  # harmonic shape.
                f.setBondParameters(0, [0, 1], [25, 5])

            elif (
                f.getEnergyFunction()
                == "step(distance(g1,g2)-r0) * (k/2)*(distance(g1,g2)-r0)^2"
            ):  # flatbottom-oneside-soft shape:
                f.setBondParameters(0, [0, 1], [250, 0.01])

            elif (
                f.getEnergyFunction()
                == "step(distance(g1,g2)-r0) * (k/2)*(distance(g1,g2))^2"
            ):  # flatbottom-oneside-sharp:
                f.setBondParameters(0, [0, 1], [250, 0.01])

            elif (
                f.getEnergyFunction()
                == "step(abs(distance(g1,g2)-r0)-w)*k*(distance(g1,g2)-r0)^2"
            ):  # flatbottom-twosided, has an additional wellsize parameter:
                f.setBondParameters(0, [0, 1], [250, 12, 1])

            f.updateParametersInContext(simulation.context)

    for i, f in enumerate(system.getForces()):
        if isinstance(f, CustomCentroidBondForce):
            logger.debug(f.getBondParameters(0))

    temp_results = []
    simulation.step(10)
    for i, f in enumerate(system.getForces()):
        if isinstance(f, CustomCentroidBondForce):
            state = simulation.context.getState(getEnergy=True, groups={i})
            logger.info("Force contributions after 20 steps:")
            logger.info(f"{f}::{state.getPotentialEnergy()}")
            temp_results.append(state.getPotentialEnergy())

    assert all(
        [results > 50 * kilojoule / mole for results in temp_results]
    )  # assert all force contributions at this point are rather large
    assert all(
        [results != 0 * kilojoule / mole for results in temp_results]
    )  # assert there is no restraint without force contribution

    logger.info("Test complete")
