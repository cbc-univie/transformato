import pytest
from transformato import load_config_yaml


@pytest.mark.rsfe
def test_rsfe_mutate_acetylaceton_methyl_common_core():
    from ..testsystems import mutate_acetylaceton_methyl_common_core

    conf = "transformato/tests/config/test-acetylacetone-tautomer-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="."
    )  # NOTE: for preprocessing input_dir is the output dir

    mutate_acetylaceton_methyl_common_core(configuration=configuration)


@pytest.mark.rsfe
def test_rsfe_mutate_bmi_small_common_core():
    from ..testsystems import mutate_bmi_small_common_core

    conf = "transformato/tests/config/test-2oj9-tautomer-pair-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="."
    )  # NOTE: for preprocessing input_dir is the output dir

    mutate_bmi_small_common_core(configuration=configuration)


@pytest.mark.rbfe
def test_rbfe_mutate_2oj9():
    from ..mutate import mutate_pure_tautomers
    from transformato import (
        ProposeMutationRoute,
        SystemStructure,
        load_config_yaml,
    )

    conf_path = "transformato/tests/config/test-2oj9-tautomer-pair-rbfe.yaml"

    configuration = load_config_yaml(
        config=conf_path, input_dir="data/", output_dir="."
    )
    # change_platform(configuration)

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.calculate_common_core()
    return (
        mutate_pure_tautomers(
            s1_to_s2,
            s1,
            s2,
            configuration,
            nr_of_bonded_windows=4,
        ),
        configuration,
        s1_to_s2,
    )
