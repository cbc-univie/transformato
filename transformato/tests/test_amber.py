from transformato import (
    load_config_yaml,
    SystemStructure,
    IntermediateStateFactory,
    ProposeMutationRoute,
)
from transformato.mutate import perform_mutations
from transformato.tests.paths import get_test_output_dir
from transformato_testsystems.testsystems import get_testsystems_dir


def test_full_amber_test():
    molecule = "ethane-methanol"

    configuration = load_config_yaml(
        config=f"/site/raid3/johannes/amber_tests/data/config/{molecule}_amber.yaml",
        input_dir="/site/raid3/johannes/amber_tests/data/",
        output_dir="/site/raid3/johannes/amber_tests/",
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.propose_common_core()
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
        nr_of_mutation_steps_charge=2,
        nr_of_mutation_steps_cc=2,
    )

    # mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol2()
    # print(mutation_list.keys())
    # i = IntermediateStateFactory(
    #     system=s2,
    #     configuration=configuration,
    #     multiple_runs=3,
    # )

    # perform_mutations(
    #     configuration=configuration,
    #     i=i,
    #     mutation_list=mutation_list,
    #     nr_of_mutation_steps_charge=2,
    # )