from transformato import (
    ProposeMutationRoute,
    SystemStructure,
    load_config_yaml,
)
import numpy as np


conf = (
    "data/test-2ra0-l51a-l51b-rbfe.yaml"
)
configuration = load_config_yaml(config=conf, input_dir="data/", output_dir=".")
s1 = SystemStructure(configuration, "structure1")
s2 = SystemStructure(configuration, "structure2")
s1_to_s2 = ProposeMutationRoute(s1, s2)

s1_to_s2.calculate_common_core()
