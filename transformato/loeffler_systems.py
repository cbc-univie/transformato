import numpy as np

import transformato
from transformato.mutate import ProposeMutationRoute
from transformato.state import IntermediateStateFactory
from transformato.system import SystemStructure
from transformato.utils import load_config_yaml

transformato_systems_dir = "/scratch/braunsfeld/transformato-systems"


def mutate_methane_to_methane_cc(conf: str = "", output_dir: str = "."):

    configuration = load_config_yaml(
        config=conf, input_dir=transformato_systems_dir, output_dir=output_dir
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol2()
    i = IntermediateStateFactory(
        system=s2,
        configuration=configuration,
    )

    # write out endpoint
    output_files = []
    intst = 1
    output_file_base = i.write_state(mutation_conf=[], intst_nr=intst)
    output_files.append(output_file_base)

    charges = mutation_list["charge"]
    intst += 1
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # generate terminal lj
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=mutation_list["terminal-lj"],
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    return output_files, configuration


def mutate_toluene_to_methane_cc(conf: str = "", output_dir: str = "."):

    configuration = load_config_yaml(
        config=conf, input_dir=transformato_systems_dir, output_dir=output_dir
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )

    # write out endpoint
    output_files = []
    intst = 1
    output_file_base = i.write_state(mutation_conf=[], intst_nr=intst)
    output_files.append(output_file_base)

    charges = mutation_list["charge"]
    intst += 1
    # start with charges
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # Turn off hydrogens
    intst += 1
    hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
    output_file_base = i.write_state(
        mutation_conf=hydrogen_lj_mutations,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    d = transformato.utils.map_lj_mutations_to_atom_idx(mutation_list["lj"])
    m = [d[(13,)]]
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    m = [d[(11,)]]
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    m = [d[(9,)]]
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    m = [d[(3,)]]
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    m = [d[(1,)]]
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # generate terminal lj
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=mutation_list["terminal-lj"],
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    m = mutation_list["transform"]
    for lambda_value in np.linspace(0.75, 0, 4):
        intst += 1
        print(lambda_value)
        # turn off charges
        output_file_base = i.write_state(
            mutation_conf=m,
            common_core_transformation=lambda_value,
            intst_nr=intst,
        )
        output_files.append(output_file_base)
    return output_files, configuration


def mutate_ethane_to_methane_cc(conf: str = "", output_dir: str = "."):

    configuration = load_config_yaml(
        config=conf, input_dir=transformato_systems_dir, output_dir=output_dir
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )

    # write out endpoint
    output_files = []
    intst = 1
    output_file_base = i.write_state(mutation_conf=[], intst_nr=intst)
    output_files.append(output_file_base)

    charges = mutation_list["charge"]
    intst += 1
    # start with charges
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # Turn off hydrogens
    intst += 1
    hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
    output_file_base = i.write_state(
        mutation_conf=hydrogen_lj_mutations,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # generate terminal lj
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=mutation_list["terminal-lj"],
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    m = mutation_list["transform"]
    for lambda_value in np.linspace(0.75, 0, 4):
        intst += 1
        print(lambda_value)
        # turn off charges
        output_file_base = i.write_state(
            mutation_conf=m,
            common_core_transformation=lambda_value,
            intst_nr=intst,
        )
        output_files.append(output_file_base)
    return output_files, configuration


def mutate_methanol_to_methane_cc(conf: str = "", output_dir: str = "."):

    configuration = load_config_yaml(
        config=conf, input_dir=transformato_systems_dir, output_dir=output_dir
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )

    # write out endpoint
    output_files = []
    intst = 1
    output_file_base = i.write_state(mutation_conf=[], intst_nr=intst)
    output_files.append(output_file_base)

    charges = mutation_list["charge"]
    intst += 1
    # start with charges
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.8,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    intst += 1
    # start with charges
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.6,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    intst += 1
    # start with charges
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.4,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    intst += 1
    # start with charges
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.2,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    intst += 1
    # start with charges
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # Turn off hydrogens
    intst += 1
    hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
    output_file_base = i.write_state(
        mutation_conf=hydrogen_lj_mutations,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # generate terminal lj
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=mutation_list["terminal-lj"],
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    m = mutation_list["transform"]
    for lambda_value in np.linspace(0.75, 0, 4):
        intst += 1
        print(lambda_value)
        # turn off charges
        output_file_base = i.write_state(
            mutation_conf=m,
            common_core_transformation=lambda_value,
            intst_nr=intst,
        )
        output_files.append(output_file_base)
    return output_files, configuration


def mutate_ethane_to_methanol_cc(
    conf: str = "", output_dir: str = "."
):  # not a loeffler system

    configuration = load_config_yaml(
        config=conf, input_dir=transformato_systems_dir, output_dir=output_dir
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )

    # write out endpoint
    output_files = []
    intst = 1
    output_file_base = i.write_state(mutation_conf=[], intst_nr=intst)
    output_files.append(output_file_base)

    charges = mutation_list["charge"]
    intst += 1
    # start with charges
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # Turn off hydrogens
    intst += 1
    hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
    output_file_base = i.write_state(
        mutation_conf=hydrogen_lj_mutations,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # generate terminal lj
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=mutation_list["terminal-lj"],
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    m = mutation_list["transform"]
    for lambda_value in np.linspace(0.75, 0, 4):
        intst += 1
        print(lambda_value)
        # turn off charges
        output_file_base = i.write_state(
            mutation_conf=m,
            common_core_transformation=lambda_value,
            intst_nr=intst,
        )
        output_files.append(output_file_base)
    return output_files, configuration


def mutate_2_CPI_7_CPI_cc(
    conf: str = "", output_dir: str = "."
):  # will be tested later on

    configuration = load_config_yaml(
        config=conf, input_dir=transformato_systems_dir, output_dir=output_dir
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )
    # write out endpoint
    output_files = []
    intst = 1
    output_file_base = i.write_state(mutation_conf=[], intst_nr=intst)
    output_files.append(output_file_base)

    charges = mutation_list["charge"]
    intst += 1
    # start with charges
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # Turn off hydrogens
    intst += 1
    hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
    output_file_base = i.write_state(
        mutation_conf=hydrogen_lj_mutations,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # generate terminal lj
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=mutation_list["terminal-lj"],
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    m = mutation_list["transform"]
    for lambda_value in np.linspace(0.75, 0, 4):
        intst += 1
        print(lambda_value)
        # turn off charges
        output_file_base = i.write_state(
            mutation_conf=m,
            common_core_transformation=lambda_value,
            intst_nr=intst,
        )
        output_files.append(output_file_base)
    return output_files, configuration


def mutate_2_methylfuran_to_methane_cc(conf: str = "", output_dir: str = "."):

    configuration = load_config_yaml(
        config=conf, input_dir=transformato_systems_dir, output_dir=output_dir
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )

    # write out endpoint
    output_files = []
    intst = 1
    output_file_base = i.write_state(mutation_conf=[], intst_nr=intst)
    output_files.append(output_file_base)

    charges = mutation_list["charge"]
    intst += 1
    # start with charges
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # Turn off hydrogens
    intst += 1
    hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
    output_file_base = i.write_state(
        mutation_conf=hydrogen_lj_mutations,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    d = transformato.utils.map_lj_mutations_to_atom_idx(mutation_list["lj"])

    # turn off heavy atoms
    m = [d[(10,)]]
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    m = [d[(8,)]]
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    m = [d[(2,)]]
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    m = [d[(0,)]]
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # generate terminal lj
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=mutation_list["terminal-lj"],
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    m = mutation_list["transform"]
    for lambda_value in np.linspace(0.75, 0, 4):
        intst += 1
        print(lambda_value)
        # turn off charges
        output_file_base = i.write_state(
            mutation_conf=m,
            common_core_transformation=lambda_value,
            intst_nr=intst,
        )
        output_files.append(output_file_base)
    return output_files, configuration


def mutate_neopentane_to_methane_cc(conf: str = "", output_dir: str = "."):

    configuration = load_config_yaml(
        config=conf, input_dir=transformato_systems_dir, output_dir=output_dir
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.propose_common_core()
    s1_to_s2.finish_common_core(
        connected_dummy_regions_cc1=[{0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}]
    )

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )

    # write out endpoint
    output_files = []
    intst = 1
    output_file_base = i.write_state(mutation_conf=[], intst_nr=intst)
    output_files.append(output_file_base)

    # start with charges
    # turn off charges
    charges = mutation_list["charge"]
    nr_of_mutation_steps_charge = (
        5  # defines the number of mutation steps for charge mutation
    )
    for lambda_value in np.linspace(1, 0, nr_of_mutation_steps_charge + 1)[1:]:
        intst += 1
        print(lambda_value)
        output_file_base = i.write_state(
            mutation_conf=charges,
            lambda_value_electrostatic=lambda_value,
            intst_nr=intst,
        )
        output_files.append(output_file_base)

    # Turn off hydrogens
    intst += 1
    hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
    output_file_base = i.write_state(
        mutation_conf=hydrogen_lj_mutations,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    d = transformato.utils.map_lj_mutations_to_atom_idx(mutation_list["lj"])

    # turn off heavy atoms
    m = [d[(13,)]]
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    m = [d[(9,)]]
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    m = [d[(5,)]]
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # generate terminal lj
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=mutation_list["terminal-lj"],
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    m = mutation_list["transform"]
    nr_of_mutation_steps_cc_bonded_terms = 5

    for lambda_value in np.linspace(1, 0, nr_of_mutation_steps_cc_bonded_terms + 1)[1:]:
        intst += 1
        print(lambda_value)
        # turn off charges
        output_file_base = i.write_state(
            mutation_conf=m,
            common_core_transformation=lambda_value,
            intst_nr=intst,
        )
        output_files.append(output_file_base)
    return output_files, configuration


def mutate_2_methylindole_to_methane_cc(conf: str = "", output_dir: str = "."):

    configuration = load_config_yaml(
        config=conf, input_dir=transformato_systems_dir, output_dir=output_dir
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )
    # write out endpoint
    output_files = []
    intst = 1
    output_file_base = i.write_state(mutation_conf=[], intst_nr=intst)
    output_files.append(output_file_base)

    charges = mutation_list["charge"]
    intst += 1
    # start with charges
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.5,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    intst += 1
    # start with charges
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # Turn off hydrogens
    intst += 1
    hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
    output_file_base = i.write_state(
        mutation_conf=hydrogen_lj_mutations,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    d = transformato.utils.map_lj_mutations_to_atom_idx(mutation_list["lj"])

    # turn off heavy atoms
    m = [d[(17,)]]
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    m = [d[(15,)]]
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    m = [d[(9,)]]
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    m = [d[(7,)]]
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    m = [d[(3,)]]
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    m = [d[(2,)]]
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    m = [d[(5,)]]
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off heavy atoms
    m = [d[(0,)]]
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # generate terminal lj
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=mutation_list["terminal-lj"],
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    m = mutation_list["transform"]
    for lambda_value in np.linspace(0.75, 0, 4):
        intst += 1
        print(lambda_value)
        # turn off charges
        output_file_base = i.write_state(
            mutation_conf=m,
            common_core_transformation=lambda_value,
            intst_nr=intst,
        )
        output_files.append(output_file_base)
    return output_files, configuration