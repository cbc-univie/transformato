def test_mutate_acetylaceton_methyl_common_core():
    from ..testsystems import mutate_acetylaceton_methyl_common_core

    mutate_acetylaceton_methyl_common_core(
        conf_path="transformato/tests/config/test-acetylacetone-tautomer-solvation-free-energy.yaml",
        input_dir="data/",
        output_dir=".",
    )
