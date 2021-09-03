import numpy as np
from transformato.constants import initialize_NUM_PROC

# read in specific topology with parameters
from transformato import load_config_yaml, FreeEnergyCalculator
from time import sleep
import multiprocessing as mp
import numpy as np
from itertools import repeat
from multiprocessing import shared_memory


def postprocessing(
    configuration: dict,
    name: str,
    engine: str = "openMM",
    max_snapshots: int = 300,
    show_summary: bool = False,
    different_path_for_dcd: str = "",
    only_vacuum: bool = False,
):
    f = FreeEnergyCalculator(configuration, name)
    if only_vacuum:
        f.envs = ["vacuum"]

    if different_path_for_dcd:
        # this is needed if the trajectories are stored at a different location than the
        # potential definitions
        path = f.base_path
        f.base_path = different_path_for_dcd
        f.load_trajs(nr_of_max_snapshots=max_snapshots)
        f.base_path = path
    else:
        f.load_trajs(nr_of_max_snapshots=max_snapshots)

    f.calculate_dG_to_common_core(engine=engine)
    if only_vacuum:
        return -1, -1, f
    else:
        ddG, dddG = f.end_state_free_energy_difference
        print(f"Free energy difference: {ddG}")
        print(f"Uncertanty: {dddG}")
        if show_summary:
            f.show_summary()
        return ddG, dddG, f


def calculate_rsfe_mp():

    conf = "config/test-2oj9-tautomer-pair-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="../../data/", output_dir="../../data"
    )  # NOTE: for preprocessing input_dir is the output dir

    # 2OJ9-original to tautomer common core
    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration, name="2OJ9-original", engine="openMM", max_snapshots=800
    )


def calculate_rbfe_mp():

    conf = "config/test-2oj9-tautomer-pair-rbfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="../../data/", output_dir="../../data"
    )  # NOTE: for preprocessing input_dir is the output dir

    # 2OJ9-original to tautomer common core
    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration, name="2OJ9-original", engine="openMM", max_snapshots=1000
    )


def procft(snapshots):
    print(snapshots.n_frames)
    for i in range(snapshots.n_frames):
        j = i * i
        print(
            j,
        )
        sleep(1)
    return j


def procsh(shr_name, dims, i):
    print(shr_name)
    print(i)
    print(dims)
    existing_shm = shared_memory.SharedMemory(name=shr_name)
    np_array = np.ndarray(dims, dtype=np.float32, buffer=existing_shm.buf)

    for i in range(len(np_array))[:50]:
        j = i * i
        print(
            j,
        )
        sleep(1)
    existing_shm.close()


def create_data():

    name = "2OJ9-original"
    conf = "config/test-2oj9-tautomer-pair-rbfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="../../data/", output_dir="../../data"
    )  # NOTE: for preprocessing input_dir is the output dir

    f = FreeEnergyCalculator(configuration, name)
    f.load_trajs(nr_of_max_snapshots=1_000)
    snapshots = f.snapshots["complex"]
    print(snapshots.n_frames)
    xyz = np.asarray(snapshots.xyz)
    dims = xyz.shape
    print(dims)
    print(xyz.dtype)
    shm = shared_memory.SharedMemory(create=True, size=xyz.nbytes)
    shm_xyz = np.ndarray(xyz.shape, dtype=xyz.dtype, buffer=shm.buf)
    shm_xyz[:] = xyz[:]
    return shm, shm_xyz, dims


def load_traj_span_mp(shm, shm_xyz, dims):

    ctx = mp.get_context("spawn")
    pool = ctx.Pool(processes=4)
    # pool.map(proct, [i for i in range(snapshots.n_frames)])

    # pool.map(procf, repeat(snapshots))
    # map(procf, zip(repeat(snapshots), [i for i in range(7)]))
    # r = map(procf, (snapshots, 1))
    # print(list(r))
    # r = map(procft, [snapshots])
    # print(list(r))
    # pool.map(procft, [snapshots,snapshots,snapshots,snapshots])
    pool.starmap(procsh, zip(repeat(shm.name), repeat(dims), [i for i in range(4)]))


if __name__ == "__main__":
    initialize_NUM_PROC(1)
    # calculate_rsfe_mp()
    calculate_rbfe_mp()
    # shm, shm_xyz, dims = create_data()
    # load_traj_span_mp(shm, shm_xyz, dims)
    # shm.close()
    # shm.unlink()
