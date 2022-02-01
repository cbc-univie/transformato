"""
Unit and regression test for the transformato package.
"""

# Import package, test suite, and other packages as needed
import logging
import os
import warnings
from io import StringIO

import parmed as pm
import pytest

# read in specific topology with parameters
# read in specific topology with parameters
from transformato import SystemStructure, load_config_yaml, psf_correction
from transformato.constants import loeffler_testsystems_dir
from transformato.tests.paths import get_test_output_dir

warnings.filterwarnings("ignore", module="parmed")


def test_read_yaml():
    """Sample test, will check ability to read yaml files"""
    settingsMap = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-rsfe.yaml",
        input_dir=".",
        output_dir="data/",
    )

    assert settingsMap["system"]["name"] == "toluene-methane-rsfe"
    assert settingsMap["system"]["structure1"]["tlc"] == "UNL"


def test_io_psf_files():
    from simtk.openmm.app import CharmmPsfFile
    from transformato.testsystems import mutate_toluene_to_methane_cc

    from .test_run_production import run_simulation

    configuration = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-rsfe.yaml",
        input_dir="data/",
        output_dir=get_test_output_dir(),
    )

    output_files = mutate_toluene_to_methane_cc(configuration=configuration)
    output_path = output_files[0]
    print(output_path)
    CharmmPsfFile(f"{output_path}/lig_in_waterbox.psf")
    CharmmPsfFile(f"{output_path}/lig_in_waterbox_corr.psf")


def test_psf_files():
    test_psf = pm.charmm.psf.CharmmPsfFile("transformato/tests/config/test_input.psf")
    output = StringIO()
    test_psf.write_psf(output)
    corrected_psf = psf_correction(output)
    correction_on = False
    for line in corrected_psf.split("\n"):  # split on newline charactar
        if "!NATOM" in line:  # if !NATOM is found start correction mode
            correction_on = True
            continue

        if "!NBOND" in line:  # if !NBOND is found exit correction mode
            correction_on = False

        if (
            correction_on == True
        ):  # if in correction mode take the string, split on whitespace and put the values in a newly formated string
            if len(line) == 0:
                pass
            else:
                assert len(line) == 118
                values = line.split()
                assert len(values) == 11


def test_initialize_systems(caplog):
    caplog.set_level(logging.DEBUG)

    configuration = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-rsfe.yaml",
        input_dir="data/",
        output_dir=get_test_output_dir(),
    )

    s1 = SystemStructure(configuration, "structure1")
    assert int(s1.offset["waterbox"]) == 0
    assert int(s1.offset["vacuum"]) == 0

    s2 = SystemStructure(configuration, "structure2")
    assert int(s2.offset["waterbox"]) == 0
    assert int(s2.offset["vacuum"]) == 0

    assert "vacuum" in s1.envs and "vacuum" in s2.envs
    assert "waterbox" in s1.envs and "waterbox" in s2.envs

    configuration = load_config_yaml(
        config="transformato/tests/config/test-2oj9-tautomer-pair-rsfe.yaml",
        input_dir="data/",
        output_dir=get_test_output_dir(),
    )

    s1 = SystemStructure(configuration, "structure1")
    assert int(s1.offset["waterbox"]) == 0
    assert int(s1.offset["vacuum"]) == 0

    s2 = SystemStructure(configuration, "structure2")
    assert int(s2.offset["waterbox"]) == 0
    assert int(s2.offset["vacuum"]) == 0


@pytest.mark.rsfe
@pytest.mark.requires_loeffler_systems
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_setup_system_for_methane_common_core():
    from transformato.testsystems import mutate_methane_to_methane_cc

    configuration = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-rsfe.yaml",
        input_dir=loeffler_testsystems_dir,
        output_dir=get_test_output_dir(),
    )
    output_files = mutate_methane_to_methane_cc(configuration=configuration)

    assert len(output_files) == 3


@pytest.mark.rsfe
@pytest.mark.requires_loeffler_systems
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_setup_system_for_toluene_common_core_with_HMR():
    from transformato.testsystems import mutate_toluene_to_methane_cc

    configuration = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-rsfe-HMR.yaml",
        input_dir=loeffler_testsystems_dir,
        output_dir=get_test_output_dir(),
    )
    output_files = mutate_toluene_to_methane_cc(configuration=configuration)
    assert len(output_files) == 16
    # set up openMM system
    base = f"{get_test_output_dir()}/toluene-methane-rsfe/toluene/intst10/"
    import sys

    # sys.path is a list of absolute path strings
    sys.path.append(base)
    system, psf = generate_openMM_system_using_cgui_scripts(base)

    # enumearte and iterate over all constraints in the system
    for cons_idx in range(system.getNumConstraints()):
        idx1, idx2, value = system.getConstraintParameters(cons_idx)
        print(
            f"{psf.atom_list[idx1]}: {system.getParticleMass(idx1)}; {psf.atom_list[idx2]}: {system.getParticleMass(idx2)}"
        )
        print(value)
    print("Finsihed.")
    assert system.getNumConstraints() == 2552


@pytest.mark.rsfe
@pytest.mark.requires_loeffler_systems
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_setup_system_for_methane_common_core_with_HMR():
    from transformato.testsystems import mutate_methane_to_methane_cc

    configuration = load_config_yaml(
        config="transformato/tests/config/test-toluene-methane-rsfe-HMR.yaml",
        input_dir=loeffler_testsystems_dir,
        output_dir=get_test_output_dir(),
    )
    output_files = mutate_methane_to_methane_cc(configuration=configuration)

    assert len(output_files) == 3

    # set up openMM system
    base = f"{get_test_output_dir()}/toluene-methane-rsfe/methane/intst3/"
    import sys

    # sys.path is a list of absolute path strings
    sys.path.append(base)
    system, psf = generate_openMM_system_using_cgui_scripts(base)

    # enumearte and iterate over all constraints in the system
    for cons_idx in range(system.getNumConstraints()):
        idx1, idx2, value = system.getConstraintParameters(cons_idx)
        print(
            f"{psf.atom_list[idx1]}: {system.getParticleMass(idx1)}; {psf.atom_list[idx2]}: {system.getParticleMass(idx2)}"
        )
        print(value)
    print("Finsihed.")
    assert system.getNumConstraints() == 2557


def generate_openMM_system_using_cgui_scripts(base: str):

    # change working directory
    current_dir = os.curdir
    os.chdir(base)
    # imports
    from omm_readinputs import read_inputs
    from omm_readparams import read_params, read_psf, read_crd, gen_box
    from simtk.openmm import unit

    # Load parameters
    print("Loading parameters")
    inputs = read_inputs(f"{base}/lig_in_waterbox.inp")
    params = read_params(f"{base}/toppar.str")
    psf = read_psf(f"{base}/lig_in_waterbox.psf")
    crd = read_crd(f"{base}/lig_in_waterbox.crd")
    psf = gen_box(psf, crd)

    # Build system
    system = psf.createSystem(
        params,
        nonbondedMethod=inputs.coulomb,
        nonbondedCutoff=inputs.r_off * unit.nanometers,
        constraints=inputs.cons,
        ewaldErrorTolerance=inputs.ewald_Tol,
        hydrogenMass=3.0 * unit.atom_mass_units,
    )
    print(inputs.cons)
    os.chdir(current_dir)
    return system, psf


@pytest.mark.rbfe
def test_lonepairs():

    from transformato import (
        load_config_yaml,
        SystemStructure,
        IntermediateStateFactory,
        ProposeMutationRoute,
    )
    from transformato.mutate import perform_mutations
    from IPython.display import SVG
    from transformato.utils import run_simulation, postprocessing
    import glob
    import subprocess
    import warnings
    import os
    import sys

    warnings.filterwarnings("ignore", module="parmed")

    molecule = "17124_18631"

    folder = "/site/raid3/johannes/run_1"

    ### should be the same within each system ###
    input_dir = "/site/raid3/student4/alex/rbfe-data/jnk1/"
    config = "/site/raid3/johannes/rbfe-data/config/jnk1/{}.yaml".format(molecule)

    configuration = load_config_yaml(
        config=config, input_dir=input_dir, output_dir=folder
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)

    s1_to_s2.propose_common_core()

    odered_connected_dummy_regions_cc1 = [[46, 22], [45, 44, 43, 26, 25]]
    s1_to_s2.finish_common_core(
        odered_connected_dummy_regions_cc1=odered_connected_dummy_regions_cc1
    )

    # generate the mutation list for the original
    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )
    # nr_of_mutation_steps_lj_of_hydrogens =3, nr_of_mutation_steps_lj_of_heavy_atoms = 2,
    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
        nr_of_mutation_steps_charge=2,
        nr_of_mutation_steps_cc=2,
    )

    # # Not necessary in this case!
    # mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol2()
    # print(mutation_list.keys())
    # i = IntermediateStateFactory(
    #     system=s2,
    #     configuration=configuration,
    # )
    # # nr_of_mutation_steps_lj_of_hydrogens =3,nr_of_mutation_steps_lj_of_heavy_atoms = 2,
    # perform_mutations(
    #     configuration=configuration,
    #     i=i,
    #     mutation_list=mutation_list,
    #     nr_of_mutation_steps_charge=2,
    # )
