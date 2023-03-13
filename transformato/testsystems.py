from transformato.mutate import ProposeMutationRoute, perform_mutations
from transformato.state import IntermediateStateFactory
from transformato.system import SystemStructure


def perform_generic_mutation(configuration: dict):
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
