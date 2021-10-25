import numpy as np

import transformato
from transformato.mutate import ProposeMutationRoute, perform_mutations
from transformato.state import IntermediateStateFactory
from transformato.system import SystemStructure
from transformato.constants import change_platform


loeffler_testsystems_dir = "/home/mwieder/Work/Projects/transformato-systems/"


def mutate_methane_to_methane_cc(configuration: dict):

    change_platform(configuration)

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol2()
    i = IntermediateStateFactory(
        system=s2,
        configuration=configuration,
    )

    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
        nr_of_mutation_steps_charge=1,
    )

    return i.output_files


def mutate_toluene_to_methane_cc(configuration: dict):

    change_platform(configuration)
    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )

    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
        list_of_heavy_atoms_to_be_mutated=[13, (11, 9), (3, 1)],
    )

    return i.output_files


def mutate_ethane_to_methane_cc(configuration: dict):

    change_platform(configuration)

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )

    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
    )

    return i.output_files


def mutate_methanol_to_methane_cc(configuration: dict):

    change_platform(configuration)

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )

    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
    )

    return i.output_files


def mutate_ethane_to_methanol_cc(configuration: dict):  # not a loeffler system

    change_platform(configuration)

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )
    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
    )

    return i.output_files


def mutate_7_CPI_to_2_CPI_cc(configuration: dict):  # will be tested later on

    change_platform(configuration)

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.completeRingsOnly = True
    s1_to_s2.propose_common_core()
    s1_to_s2.remove_idx_from_common_core_of_mol1([14])
    s1_to_s2.remove_idx_from_common_core_of_mol2([6])
    s1_to_s2.finish_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol2()

    i = IntermediateStateFactory(
        system=s2,
        configuration=configuration,
    )
    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
        list_of_heavy_atoms_to_be_mutated=[13, (11, 8), 0, (2, 10), (4, 7)],
    )
    return i.output_files


def mutate_2_CPI_to_7_CPI_cc(configuration: dict):  # will be tested later on

    change_platform(configuration)

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.completeRingsOnly = True
    s1_to_s2.propose_common_core()
    s1_to_s2.remove_idx_from_common_core_of_mol1([14])
    s1_to_s2.remove_idx_from_common_core_of_mol2([6])
    s1_to_s2.finish_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()

    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )
    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
        list_of_heavy_atoms_to_be_mutated=[(2, 4), (0, 6), 11, 8, (12, 9)],
    )

    return i.output_files


def mutate_2_methylfuran_to_methane_cc(configuration: dict):

    change_platform(configuration)

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
    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
        list_of_heavy_atoms_to_be_mutated=[10, 8, 2, 0],
    )

    return i.output_files


def mutate_neopentane_to_methane_cc(configuration: dict):

    change_platform(configuration)

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
    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
        list_of_heavy_atoms_to_be_mutated=[13, 9, 5],
    )

    return i.output_files


def mutate_2_methylindole_to_methane_cc(configuration: dict):

    change_platform(configuration)

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")

    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )
    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
        list_of_heavy_atoms_to_be_mutated=[17, 15, 9, 7, 3, 2, 5, 0],
    )

    return i.output_files


def mutate_acetylaceton_methyl_common_core(configuration: dict):
    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)

    # modify MCS matching
    from rdkit.Chem import rdFMCS

    s1_to_s2.matchValences = True
    s1_to_s2.bondCompare = rdFMCS.BondCompare.CompareOrderExact
    s1_to_s2.propose_common_core()
    s1_to_s2.remove_idx_from_common_core_of_mol1([2, 10])
    s1_to_s2.remove_idx_from_common_core_of_mol2([2, 10])

    # manually set the dummy region
    s1_to_s2.finish_common_core(
        connected_dummy_regions_cc1=[{11, 2, 10, 3, 6, 4, 13, 14, 12}]
    )

    ###############################
    ########### KETO ##############
    ###############################
    # generate the mutation list for keto acetylaceton
    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    print(mutation_list.keys())
    output_files = []

    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )

    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
        list_of_heavy_atoms_to_be_mutated=[4, 6, 3],
    )

    ###############################
    ########### ENOL ##############
    ###############################

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol2()

    i = IntermediateStateFactory(
        system=s2,
        configuration=configuration,
    )

    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
        list_of_heavy_atoms_to_be_mutated=[0, 5, 1],
    )


def mutate_bmi_small_common_core(configuration: dict):

    print("Setup system ...")
    print(configuration)
    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    print("Propose mutation ...")
    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.propose_common_core()
    print("Remove atom idx ...")

    s1_to_s2.remove_idx_from_common_core_of_mol1(
        [2, 3, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
    )
    s1_to_s2.remove_idx_from_common_core_of_mol2(
        [2, 3, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
    )

    # manually set the dummy region
    print("Manually specify dummy regions ...")

    s1_to_s2.finish_common_core(
        connected_dummy_regions_cc1=[
            {2, 3, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23},
            {38},
        ],
        connected_dummy_regions_cc2=[
            {2, 3, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23},
            {40},
        ],
    )

    ###############################
    ####### 2OJ9 - original #######
    ###############################
    # generate the mutation list for the original
    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    print(mutation_list.keys())

    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )

    # perform mutations
    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
        list_of_heavy_atoms_to_be_mutated=[(18, 22), 14, (21, 16), 20, 11],
    )

    ###############################
    ###### 2OJ9 - tautomer ########
    ###############################

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol2()

    i = IntermediateStateFactory(
        system=s2,
        configuration=configuration,
    )

    # perform mutations
    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
        list_of_heavy_atoms_to_be_mutated=[(18, 22, 14), (21, 16), 20, 11],
    )


def mutate_2ra0_l51a_l51b(configuration):
    print("Setup system ...")
    print(configuration)
    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    print("Propose mutation ...")
    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.propose_common_core()
    print("Remove atom idx ...")

    s1_to_s2.remove_idx_from_common_core_of_mol1(
        [2, 3, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
    )
    s1_to_s2.remove_idx_from_common_core_of_mol2(
        [2, 3, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
    )

    # manually set the dummy region
    print("Manually specify dummy regions ...")

    s1_to_s2.finish_common_core(
        connected_dummy_regions_cc1=[
            {2, 3, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23},
            {38},
        ],
        connected_dummy_regions_cc2=[
            {2, 3, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23},
            {40},
        ],
    )
