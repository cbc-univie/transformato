"""
Unit and regression test for the transformato package.
"""

import os

import numpy as np
import pytest
import logging
from transformato.constants import (
    initialize_NUM_PROC,
)
from transformato.tests.paths import get_test_output_dir
from transformato.utils import postprocessing

# read in specific topology with parameters
from transformato import (
    load_config_yaml,
)
from .test_mutation import (
    _set_output_files_2oj9_tautomer_pair,
    _set_output_files_acetylaceton_tautomer_pair,
    _set_output_files_toluene_methane_pair,
)


rbfe_test_systemes_generated = os.path.isdir("data/2OJ9-original-2OJ9-tautomer-rbfe")


###########################################
# 2OJ9-tautomer system
###########################################
@pytest.mark.system_2oj9
def test_compare_energies_2OJ9_original_vacuum(caplog):
    caplog.set_level(logging.WARNING)
    from transformato import FreeEnergyCalculator
    import mdtraj as md

    env = "vacuum"

    base = "data/2OJ9-original-2OJ9-tautomer-rsfe/2OJ9-original/"
    output_files_t1, _ = _set_output_files_2oj9_tautomer_pair()

    conf = "transformato/tests/config/test-2oj9-tautomer-pair-rsfe.yaml"

    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    f = FreeEnergyCalculator(configuration, "2OJ9-original")
    for idx, b in enumerate(output_files_t1):
        # used load_dcd for CHARMM
        traj = md.load_dcd(
            f"{b}/lig_in_{env}.dcd",
            f"{b}/lig_in_{env}.psf",
        )
        traj.save_dcd(f"{base}/traj.dcd")
        l_charmm = f._evaluate_e_on_all_snapshots_CHARMM(traj, idx + 1, env)
        # load dcd with openMM
        traj = md.open(f"{b}/lig_in_{env}.dcd")
        xyz, unitcell_lengths, _ = traj.read()
        xyz = xyz / 10  # correct the conversion
        l_openMM = f._evaluate_e_on_all_snapshots_openMM(
            xyz, unitcell_lengths, idx + 1, env
        )

        assert len(l_charmm) == len(l_openMM)
        s = abs(np.array(l_charmm) - np.array(l_openMM))
        print(s)
        for e_charmm, e_openMM in zip(l_charmm, l_openMM):
            assert np.isclose(e_charmm, e_openMM, rtol=0.2)
        mae = np.sum(s) / len(s)
        assert mae < 0.005


def test_lazy_eval():
    import mdtraj as md

    base_path = f"data/2OJ9-original-2OJ9-tautomer-rsfe/2OJ9-original/intst1/"
    dcd_path = f"{base_path}/lig_in_waterbox.dcd"
    psf_file = f"{base_path}/lig_in_waterbox.psf"
    traj = md.load(
        f"{dcd_path}",
        top=f"{psf_file}",
    )

    with md.open(dcd_path) as f:
        f.seek(10)
        xyz, unitcell_lengths, unitcell_angles = f.read()

    assert (xyz.shape) == (190, 2530, 3)
    assert len(xyz) == 190

    print(unitcell_lengths[0])
    print(len(unitcell_lengths))
    assert len(unitcell_lengths) == 190
    print(unitcell_angles[0])


@pytest.mark.system_2oj9
def test_compare_energies_2OJ9_original_waterbox(caplog):
    caplog.set_level(logging.WARNING)
    from transformato import FreeEnergyCalculator
    import mdtraj as md

    env = "waterbox"

    base = "data/2OJ9-original-2OJ9-tautomer-rsfe/2OJ9-original/"
    output_files_t1, _ = _set_output_files_2oj9_tautomer_pair()

    conf = "transformato/tests/config/test-2oj9-tautomer-pair-rsfe.yaml"

    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    f = FreeEnergyCalculator(configuration, "2OJ9-original")
    for idx, b in enumerate(output_files_t1):
        traj = md.load_dcd(
            f"{b}/lig_in_{env}.dcd",
            f"{b}/lig_in_{env}.psf",
        )
        traj.save_dcd(f"{base}/traj.dcd", force_overwrite=True)
        l_charmm = f._evaluate_e_on_all_snapshots_CHARMM(traj, idx + 1, env)
        traj = md.open(f"{b}/lig_in_{env}.dcd")
        xyz, unitcell_lengths, _ = traj.read()
        xyz = xyz / 10  # correct the conversion
        unitcell_lengths = unitcell_lengths / 10
        l_openMM = f._evaluate_e_on_all_snapshots_openMM(
            xyz, unitcell_lengths, idx + 1, env
        )

        assert len(l_charmm) == len(l_openMM)
        s = abs(np.array(l_charmm) - np.array(l_openMM))
        print(s)
        mae = np.sum(s) / len(s)
        assert mae < 1.0
        for e_charmm, e_openMM in zip(l_charmm, l_openMM):
            assert np.isclose(e_charmm, e_openMM, rtol=0.1)


@pytest.mark.system_2oj9
def test_compare_energies_2OJ9_tautomer_vacuum(caplog):
    caplog.set_level(logging.WARNING)
    from transformato import FreeEnergyCalculator
    import mdtraj as md

    env = "vacuum"
    base = "data/2OJ9-original-2OJ9-tautomer-rsfe/2OJ9-tautomer/"
    _, output_files_t2 = _set_output_files_2oj9_tautomer_pair()

    conf = "transformato/tests/config/test-2oj9-tautomer-pair-rsfe.yaml"

    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    f = FreeEnergyCalculator(configuration, "2OJ9-tautomer")
    for idx, b in enumerate(output_files_t2):
        # used load_dcd for CHARMM
        traj = md.load_dcd(
            f"{b}/lig_in_{env}.dcd",
            f"{b}/lig_in_{env}.psf",
        )
        traj.save_dcd(f"{base}/traj.dcd", force_overwrite=True)
        l_charmm = f._evaluate_e_on_all_snapshots_CHARMM(traj, idx + 1, env)
        # load dcd with openMM
        traj = md.open(f"{b}/lig_in_{env}.dcd")
        xyz, unitcell_lengths, _ = traj.read()
        xyz = xyz / 10  # correct the conversion
        l_openMM = f._evaluate_e_on_all_snapshots_openMM(
            xyz, unitcell_lengths, idx + 1, env
        )

        assert len(l_charmm) == len(l_openMM)
        s = abs(np.array(l_charmm) - np.array(l_openMM))
        print(s)
        mae = np.sum(s) / len(s)
        assert mae < 0.005
        for e_charmm, e_openMM in zip(l_charmm, l_openMM):
            assert np.isclose(e_charmm, e_openMM, rtol=1e-1)


@pytest.mark.system_2oj9
def test_compare_energies_2OJ9_tautomer_waterbox(caplog):
    caplog.set_level(logging.WARNING)
    from transformato import FreeEnergyCalculator
    import mdtraj as md

    env = "waterbox"
    base = "data/2OJ9-original-2OJ9-tautomer-rsfe/2OJ9-tautomer/"
    _, output_files_t2 = _set_output_files_2oj9_tautomer_pair()

    conf = "transformato/tests/config/test-2oj9-tautomer-pair-rsfe.yaml"

    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir
    f = FreeEnergyCalculator(configuration, "2OJ9-tautomer")
    for idx, b in enumerate(output_files_t2):
        print(b)
        traj = md.load_dcd(
            f"{b}/lig_in_{env}.dcd",
            f"{b}/lig_in_{env}.psf",
        )
        # used load_dcd for CHARMM
        traj.save_dcd(f"{base}/traj.dcd", force_overwrite=True)
        l_charmm = f._evaluate_e_on_all_snapshots_CHARMM(traj, idx + 1, env)
        # load dcd with openMM
        traj = md.open(f"{b}/lig_in_{env}.dcd")
        xyz, unitcell_lengths, _ = traj.read()
        xyz = xyz / 10  # correct the conversion
        unitcell_lengths = unitcell_lengths / 10
        l_openMM = f._evaluate_e_on_all_snapshots_openMM(
            xyz, unitcell_lengths, idx + 1, env
        )
        assert len(l_charmm) == len(l_openMM)
        s = abs(np.array(l_charmm) - np.array(l_openMM))
        mae = np.sum(s) / len(s)
        print(mae)
        assert mae < 1.0
        for e_charmm, e_openMM in zip(l_charmm, l_openMM):
            print(f"{e_charmm}, {e_openMM}: {e_charmm - e_openMM}")
            assert np.isclose(e_charmm, e_openMM, rtol=0.1)


# @pytest.mark.system_2oj9
# @pytest.mark.slowtest
# @pytest.mark.skipif(
#     os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
# )
# def test_2oj9_calculate_rbfe_with_openMM_only_complex():

#     initialize_NUM_PROC(1)

#     conf = "transformato/tests/config/test-2oj9-tautomer-pair-rbfe.yaml"
#     configuration = load_config_yaml(
#         config=conf, input_dir="data/", output_dir="data"
#     )  # NOTE: for preprocessing input_dir is the output dir

#     # 2OJ9-original to tautomer common core
#     ddG_openMM, dddG, f_openMM = postprocessing(
#         configuration,
#         name="2OJ9-original",
#         engine="openMM",
#         max_snapshots=10_000,
#         only_single_state='complex'
#     )


@pytest.mark.system_2oj9
@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_2oj9_calculate_rsfe_with_openMM(caplog):
    caplog.set_level(logging.DEBUG)

    conf = "transformato/tests/config/test-2oj9-tautomer-pair-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    # 2OJ9-original to tautomer common core
    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration,
        name="2OJ9-original",
        engine="openMM",
        max_snapshots=600,
    )

    assert np.isclose(ddG_openMM, 1.661582340891158)


@pytest.mark.system_2oj9
@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_2oj9_calculate_rsfe_with_different_engines():

    conf = "transformato/tests/config/test-2oj9-tautomer-pair-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    # 2OJ9-original to tautomer common core
    ddG_charmm, dddG, f_charmm = postprocessing(
        configuration, name="2OJ9-original", engine="CHARMM", max_snapshots=600
    )

    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration, name="2OJ9-original", engine="openMM", max_snapshots=600
    )

    assert np.isclose(
        f_charmm.vacuum_free_energy_differences[0, -1],
        f_openMM.vacuum_free_energy_differences[0, -1],
        rtol=1e-4,
    )

    assert np.isclose(
        f_charmm.waterbox_free_energy_differences[0, -1],
        f_openMM.waterbox_free_energy_differences[0, -1],
        rtol=1e-2,
    )

    assert np.isclose(ddG_openMM, ddG_charmm, rtol=1e-2)
    assert np.isclose(ddG_openMM, 1.6615831380721744)
    assert np.isclose(ddG_charmm, 1.6579906682671464)
    print(ddG_openMM, ddG_charmm)
    # 2OJ9-tautomer to tautomer common core

    ddG_charmm, dddG, f_charmm = postprocessing(
        configuration, name="2OJ9-tautomer", engine="CHARMM", max_snapshots=600
    )
    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration, name="2OJ9-tautomer", engine="openMM", max_snapshots=600
    )

    assert np.isclose(
        f_charmm.vacuum_free_energy_differences[0, -1],
        f_openMM.vacuum_free_energy_differences[0, -1],
        atol=1e-2,
    )

    assert np.isclose(
        f_charmm.waterbox_free_energy_differences[0, -1],
        f_openMM.waterbox_free_energy_differences[0, -1],
        rtol=1e-2,
    )

    assert np.isclose(ddG_openMM, ddG_charmm, rtol=1e-2)
    print(ddG_openMM, ddG_charmm)

    assert np.isclose(ddG_openMM, -0.4446282802464623)
    assert np.isclose(ddG_charmm, -0.4459798541540181)


@pytest.mark.system_2oj9
@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_2oj9_calculate_rsfe_with_different_switches(caplog):
    caplog.set_level(logging.WARNING)
    from .test_mutation import setup_2OJ9_tautomer_pair_rsfe
    from .test_run_production import run_simulation
    from ..analysis import FreeEnergyCalculator

    # vfswitch
    conf_path = "transformato/tests/config/test-2oj9-tautomer-pair-rsfe_vfswitch.yaml"
    configuration = load_config_yaml(
        config=conf_path, input_dir="data/", output_dir=get_test_output_dir()
    )  # NOTE: for preprocessing input_dir is the output dir

    # generate samples
    (output_files_t1, output_files_t2), _, _ = setup_2OJ9_tautomer_pair_rsfe(
        configuration=configuration
    )
    run_simulation(output_files_t1)

    f = FreeEnergyCalculator(configuration, "2OJ9-original")

    # 2OJ9-original to tautomer common core
    ddG_charmm, dddG, f_charmm_vfswitch = postprocessing(
        configuration,
        name="2OJ9-original",
        engine="CHARMM",
        max_snapshots=600,
        different_path_for_dcd="data/2OJ9-original-2OJ9-tautomer-rsfe/2OJ9-original",
    )
    ddG_openMM, dddG, f_openMM_vfswitch = postprocessing(
        configuration,
        name="2OJ9-original",
        engine="openMM",
        max_snapshots=600,
        different_path_for_dcd="data/2OJ9-original-2OJ9-tautomer-rsfe/2OJ9-original",
    )

    assert np.isclose(
        f_charmm_vfswitch.vacuum_free_energy_differences[0, -1],
        f_openMM_vfswitch.vacuum_free_energy_differences[0, -1],
        atol=1e-2,
    )

    assert np.isclose(
        f_charmm_vfswitch.waterbox_free_energy_differences[0, -1],
        f_openMM_vfswitch.waterbox_free_energy_differences[0, -1],
        atol=1e-2,
    )

    # switch
    conf_path = "transformato/tests/config/test-2oj9-tautomer-pair-rsfe_vswitch.yaml"
    configuration = load_config_yaml(
        config=conf_path, input_dir="data/", output_dir=get_test_output_dir()
    )  # NOTE: for preprocessing input_dir is the output dir

    run_simulation(output_files_t1)
    configuration = load_config_yaml(
        config=conf_path, input_dir="data/", output_dir=get_test_output_dir()
    )  # NOTE: for preprocessing input_dir is the output dir
    f = FreeEnergyCalculator(configuration, "2OJ9-original")

    # 2OJ9-original to tautomer common core
    ddG_charmm, dddG, f_charmm_switch = postprocessing(
        configuration,
        name="2OJ9-original",
        engine="CHARMM",
        max_snapshots=600,
        different_path_for_dcd="data/2OJ9-original-2OJ9-tautomer-rsfe/2OJ9-original",
    )
    ddG_charmm, dddG, f_openMM_switch = postprocessing(
        configuration,
        name="2OJ9-original",
        engine="openMM",
        max_snapshots=600,
        different_path_for_dcd="data/2OJ9-original-2OJ9-tautomer-rsfe/2OJ9-original",
    )

    assert np.isclose(
        f_charmm_switch.vacuum_free_energy_differences[0, -1],
        f_openMM_switch.vacuum_free_energy_differences[0, -1],
        atol=1e-2,
    )

    assert np.isclose(
        f_charmm_switch.waterbox_free_energy_differences[0, -1],
        f_openMM_switch.waterbox_free_energy_differences[0, -1],
        atol=1e-2,
    )

    assert np.isclose(
        f_charmm_vfswitch.vacuum_free_energy_differences[0, -1],
        f_charmm_switch.vacuum_free_energy_differences[0, -1],
        atol=1e-2,
    )
    assert np.isclose(
        f_charmm_vfswitch.waterbox_free_energy_differences[0, -1],
        f_charmm_switch.waterbox_free_energy_differences[0, -1],
        atol=1e-2,
    )
    assert np.isclose(
        f_openMM_vfswitch.waterbox_free_energy_differences[0, -1],
        f_openMM_switch.waterbox_free_energy_differences[0, -1],
        atol=1e-2,
    )


@pytest.mark.system_2oj9
@pytest.mark.slowtest
@pytest.mark.skipif(
    rbfe_test_systemes_generated == False,
    reason="These calculations need outside testsystems.",
)
def test_2oj9_calculate_rbfe_with_openMM():

    initialize_NUM_PROC(1)

    conf = "transformato/tests/config/test-2oj9-tautomer-pair-rbfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    # 2OJ9-original to tautomer common core
    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration, name="2OJ9-original", engine="openMM", max_snapshots=200
    )


###########################################
# acetylacetone-tautomer system
###########################################
@pytest.mark.system_acetylacetone
def test_compare_energies_acetylacetone_enol_vacuum(caplog):
    caplog.set_level(logging.WARNING)
    from transformato import FreeEnergyCalculator
    import mdtraj as md

    env = "vacuum"

    base = "data/acetylacetone-keto-acetylacetone-enol-rsfe/acetylacetone-enol/"
    (
        output_files_enol,
        output_files_keto,
    ) = _set_output_files_acetylaceton_tautomer_pair()

    conf = "transformato/tests/config/test-acetylacetone-tautomer-rsfe.yaml"

    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    f = FreeEnergyCalculator(configuration, "acetylacetone-enol")
    for idx, b in enumerate(output_files_enol):
        traj = md.load_dcd(
            f"{b}/lig_in_{env}.dcd",
            f"{b}/lig_in_{env}.psf",
        )
        traj.save_dcd(f"{base}/traj.dcd")
        l_charmm = f._evaluate_e_on_all_snapshots_CHARMM(traj, idx + 1, env)
        # load dcd with openMM
        traj = md.open(f"{b}/lig_in_{env}.dcd")
        xyz, unitcell_lengths, _ = traj.read()
        xyz = xyz / 10  # correct the conversion
        l_openMM = f._evaluate_e_on_all_snapshots_openMM(
            xyz, unitcell_lengths, idx + 1, env
        )
        assert len(l_charmm) == len(l_openMM)
        s = abs(np.array(l_charmm) - np.array(l_openMM))
        mae = np.sum(s) / len(s)
        print(mae)
        assert mae < 1.0

        for e_charmm, e_openMM in zip(l_charmm, l_openMM):
            assert np.isclose(e_charmm, e_openMM, rtol=1e-2)


@pytest.mark.system_acetylacetone
def test_compare_energies_acetylacetone_enol_waterbox(caplog):
    caplog.set_level(logging.WARNING)
    from transformato import FreeEnergyCalculator
    import mdtraj as md

    env = "waterbox"

    base = "data/acetylacetone-keto-acetylacetone-enol-rsfe/acetylacetone-enol/"
    (
        output_files_enol,
        output_files_keto,
    ) = _set_output_files_acetylaceton_tautomer_pair()

    conf = "transformato/tests/config/test-acetylacetone-tautomer-rsfe.yaml"

    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    f = FreeEnergyCalculator(configuration, "acetylacetone-enol")
    for idx, b in enumerate(output_files_enol):
        traj = md.load_dcd(
            f"{b}/lig_in_{env}.dcd",
            f"{b}/lig_in_{env}.psf",
        )
        # used load_dcd for CHARMM
        traj.save_dcd(f"{base}/traj.dcd", force_overwrite=True)
        l_charmm = f._evaluate_e_on_all_snapshots_CHARMM(traj, idx + 1, env)
        # load dcd with openMM
        traj = md.open(f"{b}/lig_in_{env}.dcd")
        xyz, unitcell_lengths, _ = traj.read()
        xyz = xyz / 10  # correct the conversion
        unitcell_lengths = unitcell_lengths / 10
        l_openMM = f._evaluate_e_on_all_snapshots_openMM(
            xyz, unitcell_lengths, idx + 1, env
        )
        assert len(l_charmm) == len(l_openMM)
        s = abs(np.array(l_charmm) - np.array(l_openMM))
        mae = np.sum(s) / len(s)
        print(mae)
        assert mae < 1.0

        for e_charmm, e_openMM in zip(l_charmm, l_openMM):
            assert np.isclose(e_charmm, e_openMM, rtol=1e-2)


@pytest.mark.system_acetylacetone
def test_compare_energies_acetylacetone_keto_vacuum(caplog):
    caplog.set_level(logging.WARNING)
    from transformato import FreeEnergyCalculator
    import mdtraj as md

    env = "vacuum"
    base = "data/acetylacetone-keto-acetylacetone-enol-rsfe/acetylacetone-keto/"
    (
        output_files_enol,
        output_files_keto,
    ) = _set_output_files_acetylaceton_tautomer_pair()

    conf = "transformato/tests/config/test-acetylacetone-tautomer-rsfe.yaml"

    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    f = FreeEnergyCalculator(configuration, "acetylacetone-keto")
    for idx, b in enumerate(output_files_keto):
        traj = md.load_dcd(
            f"{b}/lig_in_{env}.dcd",
            f"{b}/lig_in_{env}.psf",
        )
        # used load_dcd for CHARMM
        traj.save_dcd(f"{base}/traj.dcd", force_overwrite=True)
        l_charmm = f._evaluate_e_on_all_snapshots_CHARMM(traj, idx + 1, env)
        # load dcd with openMM
        traj = md.open(f"{b}/lig_in_{env}.dcd")
        xyz, unitcell_lengths, _ = traj.read()
        xyz = xyz / 10  # correct the conversion
        l_openMM = f._evaluate_e_on_all_snapshots_openMM(
            xyz, unitcell_lengths, idx + 1, env
        )
        assert len(l_charmm) == len(l_openMM)
        s = abs(np.array(l_charmm) - np.array(l_openMM))
        mae = np.sum(s) / len(s)
        print(mae)
        assert mae < 1.0

        for e_charmm, e_openMM in zip(l_charmm, l_openMM):
            assert np.isclose(e_charmm, e_openMM, rtol=1e-2)


@pytest.mark.system_acetylacetone
def test_compare_energies_acetylacetone_keto_waterbox(caplog):
    caplog.set_level(logging.WARNING)
    from transformato import FreeEnergyCalculator
    import mdtraj as md

    env = "waterbox"
    base = "data/acetylacetone-keto-acetylacetone-enol-rsfe/acetylacetone-keto/"
    (
        output_files_enol,
        output_files_keto,
    ) = _set_output_files_acetylaceton_tautomer_pair()

    conf = "transformato/tests/config/test-acetylacetone-tautomer-rsfe.yaml"

    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    f = FreeEnergyCalculator(configuration, "acetylacetone-keto")
    for idx, b in enumerate(output_files_keto):
        traj = md.load_dcd(
            f"{b}/lig_in_{env}.dcd",
            f"{b}/lig_in_{env}.psf",
        )
        # used load_dcd for CHARMM
        traj.save_dcd(f"{base}/traj.dcd", force_overwrite=True)
        l_charmm = f._evaluate_e_on_all_snapshots_CHARMM(traj, idx + 1, env)
        # load dcd with openMM
        traj = md.open(f"{b}/lig_in_{env}.dcd")
        xyz, unitcell_lengths, _ = traj.read()
        xyz = xyz / 10  # correct the conversion
        unitcell_lengths = unitcell_lengths / 10
        l_openMM = f._evaluate_e_on_all_snapshots_openMM(
            xyz, unitcell_lengths, idx + 1, env
        )
        assert len(l_charmm) == len(l_openMM)
        s = abs(np.array(l_charmm) - np.array(l_openMM))
        mae = np.sum(s) / len(s)
        print(mae)
        assert mae < 1.0

        for e_charmm, e_openMM in zip(l_charmm, l_openMM):
            assert np.isclose(e_charmm, e_openMM, rtol=1e-2)


@pytest.mark.system_acetylacetone
@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_acetylacetone_calculate_rsfe_with_different_engines():

    conf = "transformato/tests/config/test-acetylacetone-tautomer-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    # enol
    ddG_charmm, dddG, f_charmm = postprocessing(
        configuration, name="acetylacetone-enol", engine="CHARMM", max_snapshots=500
    )
    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration, name="acetylacetone-enol", engine="openMM", max_snapshots=500
    )
    assert np.isclose(
        f_charmm.vacuum_free_energy_differences[0, -1],
        f_openMM.vacuum_free_energy_differences[0, -1],
        rtol=1e-4,
    )
    assert np.isclose(
        f_charmm.waterbox_free_energy_differences[0, -1],
        f_openMM.waterbox_free_energy_differences[0, -1],
        rtol=1e-4,
    )

    assert np.isclose(ddG_openMM, ddG_charmm, rtol=1e-1)
    print(ddG_openMM)
    print(ddG_charmm)
    assert np.isclose(ddG_openMM, -0.6532604462065663, rtol=1e-3)
    assert np.isclose(ddG_charmm, -0.6591611245563769, rtol=1e-3)

    # keto
    ddG_charmm, dddG, f_charmm = postprocessing(
        configuration, name="acetylacetone-keto", engine="CHARMM", max_snapshots=600
    )
    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration, name="acetylacetone-keto", engine="openMM", max_snapshots=600
    )
    assert np.isclose(
        f_charmm.vacuum_free_energy_differences[0, -1],
        f_openMM.vacuum_free_energy_differences[0, -1],
        rtol=1e-4,
    )
    assert np.isclose(
        f_charmm.waterbox_free_energy_differences[0, -1],
        f_openMM.waterbox_free_energy_differences[0, -1],
        rtol=0.01,
    )
    assert np.isclose(ddG_openMM, ddG_charmm, rtol=1e-1)
    print(ddG_openMM)
    print(ddG_charmm)
    assert np.isclose(ddG_openMM, 2.691482858438775, rtol=0.01)
    assert np.isclose(ddG_charmm, 2.699116266252545, rtol=0.01)


@pytest.mark.system_acetylacetone
@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_acetylacetone_calculate_rsfe_with_different_engines_only_vacuum():
    from ..constants import kT
    from simtk import unit

    conf = "transformato/tests/config/test-acetylacetone-tautomer-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    # enol
    ddG_charmm, dddG, f_charmm = postprocessing(
        configuration,
        name="acetylacetone-enol",
        engine="CHARMM",
        max_snapshots=10_000,
        only_single_state="vacuum",
    )
    print(ddG_charmm, dddG)

    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration,
        name="acetylacetone-enol",
        engine="openMM",
        max_snapshots=10_000,
        only_single_state="vacuum",
    )
    print(ddG_charmm, dddG)

    assert np.isclose(
        f_charmm.vacuum_free_energy_differences[0, -1],
        f_openMM.vacuum_free_energy_differences[0, -1],
        rtol=1e-4,
    )

    # keto
    ddG_charmm, dddG, f_charmm = postprocessing(
        configuration,
        name="acetylacetone-keto",
        engine="CHARMM",
        max_snapshots=10_000,
        only_single_state="vacuum",
    )
    print(ddG_charmm, dddG)
    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration,
        name="acetylacetone-keto",
        engine="openMM",
        max_snapshots=10_000,
        only_single_state="vacuum",
    )
    print(ddG_charmm, dddG)
    assert np.isclose(
        f_charmm.vacuum_free_energy_differences[0, -1],
        f_openMM.vacuum_free_energy_differences[0, -1],
        rtol=1e-4,
    )

    print(
        f"results-openMM: {(f_openMM.vacuum_free_energy_differences[0, -1] * kT).value_in_unit(unit.kilocalorie_per_mole)}"
    )
    print(
        f"results-CHARMM: {(f_charmm.vacuum_free_energy_differences[0, -1] * kT).value_in_unit(unit.kilocalorie_per_mole)}"
    )


###########################################
# toluene-methane system
###########################################
@pytest.mark.system_toluene_methane
def test_compare_energies_methane_vacuum(caplog):
    caplog.set_level(logging.WARNING)
    from transformato import FreeEnergyCalculator
    import mdtraj as md

    env = "vacuum"
    base = "data/toluene-methane-rsfe/methane/"
    (
        output_files_methane,
        output_files_toluene,
    ) = _set_output_files_toluene_methane_pair()

    conf = "transformato/tests/config/test-toluene-methane-rsfe.yaml"

    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    f = FreeEnergyCalculator(configuration, "methane")
    for idx, b in enumerate(output_files_methane):
        traj = md.load_dcd(
            f"{b}/lig_in_{env}.dcd",
            f"{b}/lig_in_{env}.psf",
        )
        # used load_dcd for CHARMM
        traj.save_dcd(f"{base}/traj.dcd", force_overwrite=True)
        l_charmm = f._evaluate_e_on_all_snapshots_CHARMM(traj, idx + 1, env)
        # load dcd with openMM
        traj = md.open(f"{b}/lig_in_{env}.dcd")
        xyz, unitcell_lengths, _ = traj.read()
        xyz = xyz / 10  # correct the conversion
        unitcell_lengths = unitcell_lengths
        l_openMM = f._evaluate_e_on_all_snapshots_openMM(
            xyz, unitcell_lengths, idx + 1, env
        )
        assert len(l_charmm) == len(l_openMM)
        s = abs(np.array(l_charmm) - np.array(l_openMM))
        mae = np.sum(s) / len(s)
        print(mae)
        assert mae < 0.005
        for e_charmm, e_openMM in zip(l_charmm, l_openMM):
            assert np.isclose(e_charmm, e_openMM, rtol=1e-2)


@pytest.mark.system_toluene_methane
def test_compare_energies_methane_waterbox(caplog):
    caplog.set_level(logging.WARNING)
    from transformato import FreeEnergyCalculator
    import mdtraj as md

    env = "waterbox"
    base = "data/toluene-methane-rsfe/methane/"
    (
        output_files_methane,
        output_files_toluene,
    ) = _set_output_files_toluene_methane_pair()

    conf = "transformato/tests/config/test-toluene-methane-rsfe.yaml"

    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    f = FreeEnergyCalculator(configuration, "methane")
    for idx, b in enumerate(output_files_methane):
        traj = md.load_dcd(
            f"{b}/lig_in_{env}.dcd",
            f"{b}/lig_in_{env}.psf",
        )
        traj.save_dcd(f"{base}/traj.dcd", force_overwrite=True)
        l_charmm = f._evaluate_e_on_all_snapshots_CHARMM(traj, idx + 1, env)
        # load dcd with openMM
        traj = md.open(f"{b}/lig_in_{env}.dcd")
        xyz, unitcell_lengths, _ = traj.read()
        xyz = xyz / 10  # correct the conversion
        unitcell_lengths = unitcell_lengths / 10
        l_openMM = f._evaluate_e_on_all_snapshots_openMM(
            xyz, unitcell_lengths, idx + 1, env
        )
        assert len(l_charmm) == len(l_openMM)
        s = abs(np.array(l_charmm) - np.array(l_openMM))
        mae = np.sum(s) / len(s)
        assert mae < 0.7
        for e_charmm, e_openMM in zip(l_charmm, l_openMM):
            assert np.isclose(e_charmm, e_openMM, rtol=1e-2)


@pytest.mark.system_toluene_methane
def test_compare_energies_toluene_vacuum(caplog):
    caplog.set_level(logging.WARNING)
    from transformato import FreeEnergyCalculator
    import mdtraj as md

    env = "vacuum"
    base = "data/toluene-methane-rsfe/toluene/"
    (
        output_files_methane,
        output_files_toluene,
    ) = _set_output_files_toluene_methane_pair()

    conf = "transformato/tests/config/test-toluene-methane-rsfe.yaml"

    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    f = FreeEnergyCalculator(configuration, "toluene")
    for idx, b in enumerate(output_files_toluene):
        traj = md.load_dcd(
            f"{b}/lig_in_{env}.dcd",
            f"{b}/lig_in_{env}.psf",
        )
        # used load_dcd for CHARMM
        traj.save_dcd(f"{base}/traj.dcd", force_overwrite=True)
        l_charmm = f._evaluate_e_on_all_snapshots_CHARMM(traj, idx + 1, env)
        # load dcd with openMM
        traj = md.open(f"{b}/lig_in_{env}.dcd")
        xyz, unitcell_lengths, _ = traj.read()
        xyz = xyz / 10  # correct the conversion
        unitcell_lengths = unitcell_lengths
        l_openMM = f._evaluate_e_on_all_snapshots_openMM(
            xyz, unitcell_lengths, idx + 1, env
        )
        assert len(l_charmm) == len(l_openMM)
        s = abs(np.array(l_charmm) - np.array(l_openMM))
        mae = np.sum(s) / len(s)
        assert mae < 0.005
        for e_charmm, e_openMM in zip(l_charmm, l_openMM):
            assert np.isclose(e_charmm, e_openMM, rtol=1e-2)


@pytest.mark.system_toluene_methane
def test_compare_energies_toluene_waterbox(caplog):
    caplog.set_level(logging.WARNING)
    from transformato import FreeEnergyCalculator
    import mdtraj as md

    env = "waterbox"
    base = "data/toluene-methane-rsfe/toluene/"
    (
        output_files_methane,
        output_files_toluene,
    ) = _set_output_files_toluene_methane_pair()

    conf = "transformato/tests/config/test-toluene-methane-rsfe.yaml"

    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    f = FreeEnergyCalculator(configuration, "toluene")
    for idx, b in enumerate(output_files_toluene):
        traj = md.load_dcd(
            f"{b}/lig_in_{env}.dcd",
            f"{b}/lig_in_{env}.psf",
        )
        # used load_dcd for CHARMM
        traj.save_dcd(f"{base}/traj.dcd", force_overwrite=True)
        l_charmm = f._evaluate_e_on_all_snapshots_CHARMM(traj, idx + 1, env)
        # load dcd with openMM
        traj = md.open(f"{b}/lig_in_{env}.dcd")
        xyz, unitcell_lengths, _ = traj.read()
        xyz = xyz / 10  # correct the conversion
        unitcell_lengths = unitcell_lengths / 10
        l_openMM = f._evaluate_e_on_all_snapshots_openMM(
            xyz, unitcell_lengths, idx + 1, env
        )
        assert len(l_charmm) == len(l_openMM)
        s = abs(np.array(l_charmm) - np.array(l_openMM))
        mae = np.sum(s) / len(s)
        print(mae)
        assert mae < 0.8
        for e_charmm, e_openMM in zip(l_charmm, l_openMM):
            assert np.isclose(e_charmm, e_openMM, rtol=1e-2)


@pytest.mark.system_toluene_methane
@pytest.mark.slowtest
@pytest.mark.skipif(
    os.environ.get("TRAVIS", None) == "true", reason="Skip slow test on travis."
)
def test_toluene_to_methane_calculate_rsfe_with_different_engines():

    conf = "transformato/tests/config/test-toluene-methane-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir
    # methane
    ddG_charmm, dddG, f_charmm = postprocessing(
        configuration, name="methane", engine="CHARMM", max_snapshots=600
    )
    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration, name="methane", engine="openMM", max_snapshots=600
    )
    assert np.isclose(
        f_charmm.vacuum_free_energy_differences[0, -1],
        f_openMM.vacuum_free_energy_differences[0, -1],
        rtol=1e-4,
    )
    assert np.isclose(
        f_charmm.waterbox_free_energy_differences[0, -1],
        f_openMM.waterbox_free_energy_differences[0, -1],
        rtol=1e-3,
    )

    assert np.isclose(ddG_openMM, ddG_charmm, rtol=1e-1)
    print(f_charmm.vacuum_free_energy_differences[0, -1])
    print(f_openMM.vacuum_free_energy_differences[0, -1])
    assert np.isclose(ddG_openMM, -1.3681988336807516)
    assert np.isclose(ddG_charmm, -1.3674494407692004)

    # toluene
    ddG_charmm, dddG, f_charmm = postprocessing(
        configuration, name="toluene", engine="CHARMM", max_snapshots=600
    )
    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration, name="toluene", engine="openMM", max_snapshots=600
    )
    assert np.isclose(
        f_charmm.vacuum_free_energy_differences[0, -1],
        f_openMM.vacuum_free_energy_differences[0, -1],
        rtol=1e-4,
    )
    assert np.isclose(
        f_charmm.waterbox_free_energy_differences[0, -1],
        f_openMM.waterbox_free_energy_differences[0, -1],
        rtol=1e-4,
    )

    assert np.isclose(ddG_openMM, ddG_charmm, rtol=1e-1)
    print(f_charmm.vacuum_free_energy_differences[0, -1])
    print(f_openMM.vacuum_free_energy_differences[0, -1])
    assert np.isclose(ddG_openMM, 5.651522313536532, rtol=1e-3)
    assert np.isclose(ddG_charmm, 5.651678173410401, rtol=1e-3)


def test_postprocessing_thinning():
    from transformato import FreeEnergyCalculator

    conf = "transformato/tests/config/test-toluene-methane-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/", output_dir="data"
    )  # NOTE: for preprocessing input_dir is the output dir

    f = FreeEnergyCalculator(configuration, "methane")
    f.load_trajs(nr_of_max_snapshots=300)
    assert len(f.snapshots.keys()) == 2  # entry for vacuum and waterbox
    assert f.nr_of_states == 3  # nr of lambda states

    print(f.snapshots["vacuum"])
    assert len(f.snapshots["vacuum"]) == 900
    assert len(f.snapshots["waterbox"]) == 900

    f = FreeEnergyCalculator(configuration, "methane")
    f.load_trajs(nr_of_max_snapshots=200)

    print(f.snapshots["vacuum"])
    assert len(f.snapshots["vacuum"]) == 600
    assert len(f.snapshots["waterbox"]) == 600

    f = FreeEnergyCalculator(configuration, "toluene")
    f.load_trajs(nr_of_max_snapshots=200)
    assert len(f.snapshots.keys()) == 2  # entry for vacuum and waterbox
    assert f.nr_of_states == 13  # nr of lambda states

    print(f.snapshots["vacuum"])
    assert len(f.snapshots["vacuum"]) == 1950
    assert len(f.snapshots["waterbox"]) == 1950
