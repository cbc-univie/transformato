import numpy as np

import transformato
from transformato.mutate import ProposeMutationRoute
from transformato.state import IntermediateStateFactory
from transformato.system import SystemStructure
from transformato.utils import load_config_yaml

transformato_systems_dir = "/home/master/transformato-systems"


def mutate_methane_to_methane_cc(conf: str = "", modifier: str = ""):

    if conf:
        configuration = load_config_yaml(
            config=conf, input_dir=transformato_systems_dir, output_dir="."
        )
    else:
        conf = f"{transformato_systems_dir}/config/toluene-methane-solvation-free-energy.yaml"
        configuration = load_config_yaml(
            config=conf, input_dir=transformato_systems_dir, output_dir="."
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

    # using different switching functinos for CHARMM
    if modifier:
        if "switch" in modifier:
            configuration["simulation"]["parameters"]["switch"] = modifier
        i.path += f"-{modifier}"

    output_files = []
    # mutate everything else before touching bonded terms
    intst = 1
    charges = mutation_list["charge"]
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


def mutate_toluene_to_methane_cc(conf: str = "", modifier: str = ""):

    if conf:
        configuration = load_config_yaml(
            config=conf, input_dir=transformato_systems_dir, output_dir="."
        )
    else:
        conf = f"{transformato_systems_dir}/config/toluene-methane-solvation-free-energy.yaml"
        configuration = load_config_yaml(
            config=conf, input_dir=transformato_systems_dir, output_dir="."
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
    # using different switching functinos for CHARMM
    if modifier:
        if "switch" in modifier:
            configuration["simulation"]["parameters"]["switch"] = modifier
        i.path += f"-{modifier}"

    # start with charges
    output_files = []
    # mutate everything else before touching bonded terms
    intst = 0
    charges = mutation_list["charge"]
    intst += 1
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # Turn of hydrogens
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

    # turn off heavy atoms
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    m = [d[(11,)], d[(9,)]]

    # turn off heavy atoms
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    m = [d[(3,)], d[(1,)]]

    # turn off heavy atoms
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
    for lambda_value in np.linspace(0.25, 1, 3):
        intst += 1
        print(lambda_value)
        # turn off charges
        output_file_base = i.write_state(
            mutation_conf=m,
            common_core_transformation=1 - lambda_value,
            intst_nr=intst,
        )
        output_files.append(output_file_base)
    return output_files, configuration

def mutate_ethane_to_methane_cc(conf: str = "", modifier: str = ""):

    if conf:
        configuration = load_config_yaml(
            config=conf, input_dir=transformato_systems_dir, output_dir="."
        )
    else:
        conf = f"{transformato_systems_dir}/config/ethane-methane-solvation-free-energy.yaml"
        configuration = load_config_yaml(
            config=conf, input_dir=transformato_systems_dir, output_dir="."
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
    # using different switching functinos for CHARMM
    if modifier:
        if "switch" in modifier:
            configuration["simulation"]["parameters"]["switch"] = modifier
        i.path += f"-{modifier}"

    # start with charges
    output_files = []
    # mutate everything else before touching bonded terms
    intst = 0
    charges = mutation_list["charge"]
    intst += 1
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # Turn of hydrogens
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
    for lambda_value in np.linspace(0.25, 1, 3):
        intst += 1
        print(lambda_value)
        # turn off charges
        output_file_base = i.write_state(
            mutation_conf=m,
            common_core_transformation=1 - lambda_value,
            intst_nr=intst,
        )
        output_files.append(output_file_base)
    return output_files, configuration

def mutate_methanol_to_methane_cc(conf: str = "", modifier: str = ""):

    if conf:
        configuration = load_config_yaml(
            config=conf, input_dir=transformato_systems_dir, output_dir="."
        )
    else:
        conf = f"{transformato_systems_dir}/config/methanol-methane-solvation-free-energy.yaml"
        configuration = load_config_yaml(
            config=conf, input_dir=transformato_systems_dir, output_dir="."
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
    # using different switching functinos for CHARMM
    if modifier:
        if "switch" in modifier:
            configuration["simulation"]["parameters"]["switch"] = modifier
        i.path += f"-{modifier}"

    # start with charges
    output_files = []
    # mutate everything else before touching bonded terms
    intst = 0
    charges = mutation_list["charge"]
    intst += 1
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # Turn of hydrogens
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
    for lambda_value in np.linspace(0.25, 1, 3):
        intst += 1
        print(lambda_value)
        # turn off charges
        output_file_base = i.write_state(
            mutation_conf=m,
            common_core_transformation=1 - lambda_value,
            intst_nr=intst,
        )
        output_files.append(output_file_base)
    return output_files, configuration

def mutate_ethane_to_methanol_cc(conf: str = "", modifier: str = ""):

    if conf:
        configuration = load_config_yaml(
            config=conf, input_dir=transformato_systems_dir, output_dir="."
        )
    else:
        conf = f"{transformato_systems_dir}/config/ethane-methanol-solvation-free-energy.yaml"
        configuration = load_config_yaml(
            config=conf, input_dir=transformato_systems_dir, output_dir="."
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
    # using different switching functinos for CHARMM
    if modifier:
        if "switch" in modifier:
            configuration["simulation"]["parameters"]["switch"] = modifier
        i.path += f"-{modifier}"

    # start with charges
    output_files = []
    # mutate everything else before touching bonded terms
    intst = 0
    charges = mutation_list["charge"]
    intst += 1
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # Turn of hydrogens
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
    for lambda_value in np.linspace(0.25, 1, 3):
        intst += 1
        print(lambda_value)
        # turn off charges
        output_file_base = i.write_state(
            mutation_conf=m,
            common_core_transformation=1 - lambda_value,
            intst_nr=intst,
        )
        output_files.append(output_file_base)
    return output_files, configuration

def mutate_2_CPI_7_CPI_cc(conf: str = "", modifier: str = ""):

    if conf:
        configuration = load_config_yaml(
            config=conf, input_dir=transformato_systems_dir, output_dir="."
        )
    else:
        conf = f"{transformato_systems_dir}/config/toluene-methane-solvation-free-energy.yaml"
        configuration = load_config_yaml(
            config=conf, input_dir=transformato_systems_dir, output_dir="."
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
    # using different switching functinos for CHARMM
    if modifier:
        if "switch" in modifier:
            configuration["simulation"]["parameters"]["switch"] = modifier
        i.path += f"-{modifier}"

    # start with charges
    output_files = []
    # mutate everything else before touching bonded terms
    intst = 0
    charges = mutation_list["charge"]
    intst += 1
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # Turn of hydrogens
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
    for lambda_value in np.linspace(0.25, 1, 3):
        intst += 1
        print(lambda_value)
        # turn off charges
        output_file_base = i.write_state(
            mutation_conf=m,
            common_core_transformation=1 - lambda_value,
            intst_nr=intst,
        )
        output_files.append(output_file_base)
    return output_files, configuration

def mutate_2_methylfuran_to_methane_cc(conf: str = "", modifier: str = ""):

    if conf:
        configuration = load_config_yaml(
            config=conf, input_dir=transformato_systems_dir, output_dir="."
        )
    else:
        conf = f"{transformato_systems_dir}/config/2-methylfuran-methane-solvation-free-energy.yaml"
        configuration = load_config_yaml(
            config=conf, input_dir=transformato_systems_dir, output_dir="."
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
    # using different switching functinos for CHARMM
    if modifier:
        if "switch" in modifier:
            configuration["simulation"]["parameters"]["switch"] = modifier
        i.path += f"-{modifier}"

    # start with charges
    output_files = []
    # mutate everything else before touching bonded terms
    intst = 0
    charges = mutation_list["charge"]
    intst += 1
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # Turn of hydrogens
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
    m = [d[(10,)], d[(8,)]]

    # turn off heavy atoms
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    m = [d[(2,)], d[(0,)]]

    # turn off heavy atoms
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
    for lambda_value in np.linspace(0.25, 1, 3):
        intst += 1
        print(lambda_value)
        # turn off charges
        output_file_base = i.write_state(
            mutation_conf=m,
            common_core_transformation=1 - lambda_value,
            intst_nr=intst,
        )
        output_files.append(output_file_base)
    return output_files, configuration

def mutate_neopentane_to_methane_cc(conf: str = "", modifier: str = ""):

    if conf:
        configuration = load_config_yaml(
            config=conf, input_dir=transformato_systems_dir, output_dir="."
        )
    else:
        conf = f"{transformato_systems_dir}/config/neopentane-methane-solvation-free-energy.yaml"
        configuration = load_config_yaml(
            config=conf, input_dir=transformato_systems_dir, output_dir="."
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
    # using different switching functinos for CHARMM
    if modifier:
        if "switch" in modifier:
            configuration["simulation"]["parameters"]["switch"] = modifier
        i.path += f"-{modifier}"

    # start with charges
    output_files = []
    # mutate everything else before touching bonded terms
    intst = 0
    charges = mutation_list["charge"]
    intst += 1
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # Turn of hydrogens
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

    # turn off heavy atoms
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    m = [d[(9,)], d[(5,)]]

    # turn off heavy atoms
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
    for lambda_value in np.linspace(0.25, 1, 3):
        intst += 1
        print(lambda_value)
        # turn off charges
        output_file_base = i.write_state(
            mutation_conf=m,
            common_core_transformation=1 - lambda_value,
            intst_nr=intst,
        )
        output_files.append(output_file_base)
    return output_files, configuration

def mutate_2_methylindole_to_methane_cc(conf: str = "", modifier: str = ""):

    if conf:
        configuration = load_config_yaml(
            config=conf, input_dir=transformato_systems_dir, output_dir="."
        )
    else:
        conf = f"{transformato_systems_dir}/config/2-methylindole-methane-solvation-free-energy.yaml"
        configuration = load_config_yaml(
            config=conf, input_dir=transformato_systems_dir, output_dir="."
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
    # using different switching functinos for CHARMM
    if modifier:
        if "switch" in modifier:
            configuration["simulation"]["parameters"]["switch"] = modifier
        i.path += f"-{modifier}"

    # start with charges
    output_files = []
    # mutate everything else before touching bonded terms
    intst = 0
    charges = mutation_list["charge"]
    intst += 1
    # turn off charges
    output_file_base = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # Turn of hydrogens
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
    m = [d[(17,)], d[(15,)]]

    # turn off heavy atoms
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    m = [d[(9,)], d[(7,)]]

    # turn off heavy atoms
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)
    
    m = [d[(3,)], d[(2,)]]

    # turn off heavy atoms
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=m,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)
    
    m = [d[(5,)], d[(0,)]]

    # turn off heavy atoms
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
    for lambda_value in np.linspace(0.25, 1, 3):
        intst += 1
        print(lambda_value)
        # turn off charges
        output_file_base = i.write_state(
            mutation_conf=m,
            common_core_transformation=1 - lambda_value,
            intst_nr=intst,
        )
        output_files.append(output_file_base)
    return output_files, configuration