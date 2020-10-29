def _mutate_methane_to_methane_cc(modifier: str = ""):
    conf = "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml"
    configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol2()
    i = IntermediateStateFactory(
        system=s2,
        configuration=configuration,
    )

    if modifier:
        # building whole file
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




def _mutate_toluene_to_methane_cc():
    conf = "transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml"
    configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )

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
    intst += 1

    output_file_base = i.write_state(
        mutation_conf=mutation_list["lj"],
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
