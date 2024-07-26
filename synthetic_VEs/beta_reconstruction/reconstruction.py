from typing import List, Tuple
import pathlib

import numpy as np
import networkx as nx
from tqdm.auto import tqdm
from defdap import ebsd
from defdap.quat import Quat

from beta_reconstruction.crystal_relations import (
    unq_hex_syms, hex_syms, unq_cub_syms, burg_trans
)


def calc_beta_oris(alpha_ori: Quat) -> List[Quat]:
    """Calculate the possible beta orientations for a given alpha orientation.

    Uses the Burgers relation and crystal symmetries to calculate beta orientations.

    Parameters
    ----------
    alpha_ori
        Orientation of an alpha grain

    Returns
    -------
    list of Quat
        List of possible beta orientations
    """
    beta_oris = []

    for sym in unq_hex_syms:
        beta_oris.append(burg_trans * sym.conjugate * alpha_ori)

    return beta_oris


def beta_oris_from_cub_sym(
    alpha_ori: Quat,
    unq_cub_sym_idx: int,
    hex_sym_idx: int
) -> List[Quat]:
    """

    Parameters
    ----------
    alpha_ori
        The orientation of the grain in the alpha phase.
    unq_cub_sym_idx

    hex_sym_idx


    Returns
    -------
    list of Quat
        Possible beta orientations from given symmetries

    """
    if not (0 <= unq_cub_sym_idx <= 11):
        raise ValueError("unq_cub_sym_idx must be between 0 and 11 inclusive")
    if not (0 <= hex_sym_idx <= 11):
        raise ValueError("hex_sym_idx must be between 0 and 11 inclusive")

    beta_oris = []

    beta_ori_base = hex_syms[hex_sym_idx].conjugate * alpha_ori

    # all cases have one possible beta orientation
    beta_oris.append(burg_trans * beta_ori_base)

    if unq_cub_sym_idx == 9:
        # two extra possible beta orientations
        # B - unq_hex_syms[1] is C^+_3z:
        beta_oris.append(
            burg_trans * unq_hex_syms[1].conjugate * beta_ori_base
        )
        # C - unq_hex_syms[2] is C^+_6z:
        beta_oris.append(
            burg_trans * unq_hex_syms[2].conjugate * beta_ori_base
        )

    if unq_cub_sym_idx > 9:
        # one extra possible beta orientations
        # D - unq_hex_syms[4] is C'_22:
        beta_oris.append(
            burg_trans * unq_hex_syms[4].conjugate * beta_ori_base
        )

    return beta_oris


def calc_misori_of_variants(
    alpha_ori_inv: Quat,
    neighbour_ori: Quat,
    unq_cub_sym_comps: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Calculate possible symmetry variants between two orientations.

    Calculate all possible sym variants for misorientation between two
    orientations undergoing a Burgers type transformation. Then calculate
    the misorientation to the nearest cubic symmetry, this is the deviation
    to a perfect Burgers transformation.

    Parameters
    ----------
    alpha_ori_inv
        Inverse of first orientation
    neighbour_ori
        Second orientation
    unq_cub_sym_comps
        Components of the unique cubic symmetries

    Returns
    -------
    min_misoris : np.ndarray 
       The minimum misorientation for each of the possible beta variants - shape (12, 12)

    min_cub_sym_idx : np.ndarray
       The minimum cubic symmetry index for each of the possible variants - shape (12, 12)

    """
    # calculate all possible S^B_m (eqn 11. from [1]) from the
    # measured misorientation from 2 neighbour alpha grains
    # for each S^B_m calculate the 'closest' cubic symmetry
    # (from reduced subset) and the deviation from this symmetry

    # Vectorised calculation of:
    # hex_sym[j].inv * ((neighbour_ori * alpha_ori_inv) * hex_sym[i])
    # labelled: d = h2.inv * (c * h1)
    hex_sym_comps = Quat.extract_quat_comps(hex_syms)
    c = (neighbour_ori * alpha_ori_inv).quatCoef
    h1 = np.repeat(hex_sym_comps, 12, axis=1)  # outer loop
    h2 = np.tile(hex_sym_comps, (1, 12))  # inner loop
    d = np.zeros_like(h1)

    c_dot_h1 = c[1]*h1[1] + c[2]*h1[2] + c[3]*h1[3]
    c_dot_h2 = c[1]*h2[1] + c[2]*h2[2] + c[3]*h2[3]
    h1_dot_h2 = h1[1]*h2[1] + h1[2]*h2[2] + h1[3]*h2[3]

    d[0] = (c[0]*h1[0]*h2[0] - h2[0]*c_dot_h1 +
            c[0]*h1_dot_h2 + h1[0]*c_dot_h2 +
            h2[1] * (c[2]*h1[3] - c[3]*h1[2]) +
            h2[2] * (c[3]*h1[1] - c[1]*h1[3]) +
            h2[3] * (c[1]*h1[2] - c[2]*h1[1]))
    d[1] = (c[0]*h2[0]*h1[1] + h1[0]*h2[0]*c[1] - c[0]*h1[0]*h2[1] +
            c_dot_h1*h2[1] + c_dot_h2*h1[1] - h1_dot_h2*c[1] +
            h2[0] * (c[2]*h1[3] - c[3]*h1[2]) +
            c[0] * (h1[2]*h2[3] - h1[3]*h2[2]) +
            h1[0] * (c[2]*h2[3] - c[3]*h2[2]))
    d[2] = (c[0]*h2[0]*h1[2] + h1[0]*h2[0]*c[2] - c[0]*h1[0]*h2[2] +
            c_dot_h1*h2[2] + c_dot_h2*h1[2] - h1_dot_h2*c[2] +
            h2[0] * (c[3]*h1[1] - c[1]*h1[3]) +
            c[0] * (h1[3]*h2[1] - h1[1]*h2[3]) +
            h1[0] * (c[3]*h2[1] - c[1]*h2[3]))
    d[3] = (c[0]*h2[0]*h1[3] + h1[0]*h2[0]*c[3] - c[0]*h1[0]*h2[3] +
            c_dot_h1*h2[3] + c_dot_h2*h1[3] - h1_dot_h2*c[3] +
            h2[0] * (c[1]*h1[2] - c[2] * h1[1]) +
            c[0] * (h1[1]*h2[2] - h1[2] * h2[1]) +
            h1[0] * (c[1]*h2[2] - c[2] * h2[1]))

    # Vectorised calculation of:
    # burg_trans * (d * burg_trans.inv)
    # labelled: beta_vars = b * (c * b.inv)
    b = burg_trans.quatCoef
    beta_vars = np.zeros_like(h1)

    b_dot_b = b[1]*b[1] + b[2]*b[2] + b[3]*b[3]
    b_dot_d = b[1]*d[1] + b[2]*d[2] + b[3]*d[3]

    beta_vars[0] = d[0] * (b[0]*b[0] + b_dot_b)
    beta_vars[1] = (d[1] * (b[0]*b[0] - b_dot_b) + 2*b_dot_d*b[1] +
                    2*b[0] * (b[2]*d[3] - b[3]*d[2]))
    beta_vars[2] = (d[2] * (b[0]*b[0] - b_dot_b) + 2*b_dot_d*b[2] +
                    2*b[0] * (b[3]*d[1] - b[1]*d[3]))
    beta_vars[3] = (d[3] * (b[0]*b[0] - b_dot_b) + 2*b_dot_d*b[3] +
                    2*b[0] * (b[1]*d[2] - b[2]*d[1]))

    # calculate misorientation to each of the cubic symmetries
    misoris = np.einsum("ij,ik->jk", beta_vars, unq_cub_sym_comps)
    misoris = np.abs(misoris)
    misoris[misoris > 1] = 1.
    misoris = 2 * np.arccos(misoris)

    # find the cubic symmetry with minimum misorientation for each of
    # the beta misorientation variants
    min_cub_sym_idx = np.argmin(misoris, axis=1)
    min_misoris = misoris[np.arange(144), min_cub_sym_idx]
    # reshape to 12 x 12 for each of the hex sym multiplications
    min_cub_sym_idx = min_cub_sym_idx.reshape((12, 12))
    min_misoris = min_misoris.reshape((12, 12))

    return min_misoris, min_cub_sym_idx


def calc_beta_oris_from_misori(
    alpha_ori: Quat,
    neighbour_oris: List[Quat],
    burg_tol: float = 5
) -> Tuple[List[List[Quat]], List[float]]:
    """Calculate the possible beta orientations for a given alpha
    orientation using the misorientation relation to neighbour orientations.

    Parameters
    ----------
    alpha_ori
        A quaternion representing the alpha orientation

    neighbour_oris
        Quaternions representing neighbour grain orientations

    burg_tol
        The threshold misorientation angle to determine neighbour relations

    Returns
    -------
    list of lists of defdap.Quat.quat
        Possible beta orientations, grouped by each neighbour. Any
        neighbour with deviation greater than the tolerance is excluded.
    list of float
        Deviations from perfect Burgers transformation

    """
    burg_tol *= np.pi / 180.
    # This needed to move further up calculation process
    unq_cub_sym_comps = Quat.extract_quat_comps(unq_cub_syms)

    alpha_ori_inv = alpha_ori.conjugate

    beta_oris = []
    beta_devs = []

    for neighbour_ori in neighbour_oris:

        min_misoris, min_cub_sym_idxs = calc_misori_of_variants(
            alpha_ori_inv, neighbour_ori, unq_cub_sym_comps
        )

        # find the hex symmetries (i, j) from give the minimum
        # deviation from the burgers relation for the minimum store:
        # the deviation, the hex symmetries (i, j) and the cubic
        # symmetry if the deviation is over a threshold then set
        # cubic symmetry to -1
        min_misori_idx = np.unravel_index(np.argmin(min_misoris),
                                          min_misoris.shape)
        burg_dev = min_misoris[min_misori_idx]

        if burg_dev < burg_tol:
            beta_oris.append(beta_oris_from_cub_sym(
                alpha_ori, min_cub_sym_idxs[min_misori_idx], int(min_misori_idx[0])
            ))
            beta_devs.append(burg_dev)

    return beta_oris, beta_devs


def calc_beta_oris_from_boundary_misori(
    grain: ebsd.Grain,
    neighbour_network: nx.Graph,
    quat_array: np.ndarray,
    alpha_phase_id : int,
    burg_tol: float = 5
) -> Tuple[List[List[Quat]], List[float], List[Quat]]:
    """Calculate the possible beta orientations for pairs of alpha and
    neighbour orientations using the misorientation relation to neighbour
    orientations.

    Parameters
    ----------
    grain
        The grain currently being reconstructed

    neighbour_network
        A neighbour network mapping grain boundary connectivity

    quat_array
        Array of quaternions, representing the orientations of the pixels of the EBSD map

    burg_tol :
        The threshold misorientation angle to determine neighbour relations

    Returns
    -------
    list of lists of defdap.Quat.quat
        Possible beta orientations, grouped by each neighbour. Any
        neighbour with deviation greater than the tolerance is excluded.
    list of float
        Deviations from perfect Burgers transformation
    list of Quat
        Alpha orientations

    """
    # This needed to move further up calculation process
    unq_cub_sym_comps = Quat.extract_quat_comps(unq_cub_syms)

    beta_oris = []
    beta_devs = []
    alpha_oris = []

    neighbour_grains = neighbour_network.neighbors(grain)
    neighbour_grains = [grain for grain in neighbour_grains
                        if grain.phaseID == alpha_phase_id]
    for neighbour_grain in neighbour_grains:

        bseg = neighbour_network[grain][neighbour_grain]['boundary']
        # check sense of bseg
        if grain is bseg.grain1:
            ipoint = 0
        else:
            ipoint = 1

        for boundary_point_pair in bseg.boundaryPointPairsX:
            point = boundary_point_pair[ipoint]
            alpha_ori = quat_array[point[1], point[0]]

            point = boundary_point_pair[ipoint - 1]
            neighbour_ori = quat_array[point[1], point[0]]

            min_misoris, min_cub_sym_idxs = calc_misori_of_variants(
                alpha_ori.conjugate, neighbour_ori, unq_cub_sym_comps
            )

            # find the hex symmetries (i, j) from give the minimum
            # deviation from the burgers relation for the minimum store:
            # the deviation, the hex symmetries (i, j) and the cubic
            # symmetry if the deviation is over a threshold then set
            # cubic symmetry to -1
            min_misori_idx = np.unravel_index(np.argmin(min_misoris),
                                              min_misoris.shape)
            burg_dev = min_misoris[min_misori_idx]

            if burg_dev < burg_tol / 180 * np.pi:
                beta_oris.append(beta_oris_from_cub_sym(
                    alpha_ori, min_cub_sym_idxs[min_misori_idx], int(min_misori_idx[0])
                ))
                beta_devs.append(burg_dev)
                alpha_oris.append(alpha_ori)

    return beta_oris, beta_devs, alpha_oris


def count_beta_variants(
    beta_oris: List[Quat],
    possible_beta_oris: List[List[Quat]],
    ori_tol: float
) -> np.ndarray:
    """

    Parameters
    ----------
    beta_oris
        Possible beta orientations from burgers relation, there are always 6
    possible_beta_oris
        Possible beta orientations from misorientations between neighbouring grains
    ori_tol
        Tolerance for binning of the orientations into the possible 6
    Returns
    -------
    list of int:
        The variant count for the grain

    """
    if not possible_beta_oris:
        return np.zeros(6, dtype=int)
    # divide 2 because of 2* in misorientation definition
    ori_tol = np.cos(ori_tol / 2 * np.pi / 180.)
    # flatten list of lists
    possible_beta_oris = [item for sublist in possible_beta_oris for item in sublist]

    misoris = np.empty((len(possible_beta_oris), 6))
    for ori_index, ori in enumerate(possible_beta_oris):
        for other_ori_index, other_ori in enumerate(beta_oris):
            misoris[ori_index, other_ori_index] = ori.misOri(other_ori, "cubic")

    # max is actually min because actual misorientation is arccos of this
    max_misoris_idx = np.nanargmax(misoris, axis=1)
    max_misoris = misoris[np.arange(len(possible_beta_oris)), max_misoris_idx]
    variant_count, _ = np.histogram(max_misoris_idx[max_misoris > ori_tol],
                                    range(0, 7))

    return variant_count


def load_map(
    ebsd_path: str,
    min_grain_size: int = 3,
    boundary_tolerance: int = 3,
    use_kuwahara: bool = False,
    kuwahara_tolerance: int = 5
) -> ebsd.Map:
    """Load in EBSD data and do the required prerequisite computations."""

    ebsd_path = pathlib.Path(ebsd_path)
    if ebsd_path.suffix == ".ctf":
        map_type = "OxfordText"
    elif ebsd_path.suffix == ".crc":
        map_type = "OxfordBinary"
    else:
        raise TypeError("Unknown EBSD map type. Can only read .ctf and .crc files.")

    ebsd_map = ebsd.Map(ebsd_path.with_suffix(""), dataType=map_type)
    ebsd_map.buildQuatArray()

    if use_kuwahara:
        ebsd_map.filterData(misOriTol=kuwahara_tolerance)

    ebsd_map.findBoundaries(boundDef=boundary_tolerance)
    ebsd_map.findGrains(minGrainSize=min_grain_size)

    ebsd_map.calcGrainAvOris()

    ebsd_map.buildNeighbourNetwork()

    return ebsd_map


def modal_variant(alpha_grains: List[ebsd.Grain]) -> np.ndarray:
    """Given a map of grains with variant counts, assign the prior beta
    orientation of the grains to the variant with the highest count.

    Parameters
    ----------
    ebsd_map
        EBSD map to assign the beta variants for.
    alpha_phase_id
        Index of the alpha phase in the EBSD map.

    """
    modal_variants = np.empty(len(alpha_grains), dtype=np.int8)
    for i, grain in enumerate(alpha_grains):
        variant_count = grain.variant_count
        mode_variant = np.where(variant_count == np.max(variant_count))[0]
        if len(mode_variant) == 1:
            mode_variant = mode_variant[0]
        else:
            #  multiple variants with same max
            mode_variant = -1
        modal_variants[i] = mode_variant

    # TODO: Change modeVariant to assigned_variant
    return modal_variants


def assign_beta_variants(
    ebsd_map: ebsd.Map,
    mode: str = "modal",
    alpha_phase_id: int = 0
):
    """Given a map of grains with variant counts, determine the prior
    beta orientation of the grains.

    Parameters
    ----------
    ebsd_map:
        EBSD map to assign the beta variants for.
    mode
        How to perform beta orientation assignment
            'modal': The beta orientation is assigned to the variant
                     with the highest count.
    alpha_phase_id
        Index of the alpha phase in the EBSD map.

    """
    alpha_grains = [grain for grain in ebsd_map
                    if grain.phaseID == alpha_phase_id]
    if mode == "modal":
        assigned_variants = modal_variant(alpha_grains)
    else:
        raise NotImplementedError(f"Mode '{mode}' is not a recognised "
                                  f"way to assign variants.")

    for grain, assigned_variant in zip(alpha_grains, assigned_variants):
        if assigned_variant >= 0:
            parent_beta_ori = grain.beta_oris[assigned_variant]
        else:
            # parent_beta_ori = Quat(1., 0., 0., 0.)
            parent_beta_ori = None

        grain.assigned_variant = assigned_variant
        grain.parent_beta_ori = parent_beta_ori
        # grain.modeVariant = mode_variant
        # grain.parentBetaOri = parent_beta_ori

    print("Assignment of beta variants complete.")


def construct_variant_map(
    ebsd_map: ebsd.Map,
    alpha_phase_id: int = 0
) -> np.ndarray:
    alpha_grains = (grain for grain in ebsd_map
                    if grain.phaseID == alpha_phase_id)
    all_lists = ((grain.grainID, grain.assigned_variant) for grain in alpha_grains)
    grain_ids, assigned_variants = zip(*all_lists)

    # points not part of a grain or other phases (-2) and
    # those that were not reconstructed (-1)
    return ebsd_map.grainDataToMapData(
        assigned_variants, grainIds=grain_ids, bg=-2
    )


def construct_beta_quat_array(
    ebsd_map: ebsd.Map,
    alpha_phase_id: int = 0,
    variant_map: np.ndarray = None,
) -> np.ndarray:
    """Construct

    Parameters
    ----------
    ebsd_map:
        EBSD map to assign the beta variants for.
    alpha_phase_id
        Index of the alpha phase in the EBSD map.

    """
    if variant_map is None:
        variant_map = construct_variant_map(ebsd_map, alpha_phase_id)

    transformations = []
    for sym in unq_hex_syms:
        transformations.append(burg_trans * sym.conjugate)
    trans_comps = Quat.extract_quat_comps(transformations)
    trans_comps = trans_comps[:, variant_map[variant_map >= 0]]

    quat_comps = Quat.extract_quat_comps(ebsd_map.quatArray[variant_map >= 0])
    quat_comps_beta = np.empty_like(quat_comps)

    # transformations[variant] * quat
    quat_comps_beta[0, :] = (trans_comps[0, :] * quat_comps[0, :]
                             - trans_comps[1, :] * quat_comps[1, :]
                             - trans_comps[2, :] * quat_comps[2, :]
                             - trans_comps[3, :] * quat_comps[3, :])
    quat_comps_beta[1, :] = (trans_comps[1, :] * quat_comps[0, :]
                             + trans_comps[0, :] * quat_comps[1, :]
                             - trans_comps[3, :] * quat_comps[2, :]
                             + trans_comps[2, :] * quat_comps[3, :])
    quat_comps_beta[2, :] = (trans_comps[2, :] * quat_comps[0, :]
                             + trans_comps[0, :] * quat_comps[2, :]
                             - trans_comps[1, :] * quat_comps[3, :]
                             + trans_comps[3, :] * quat_comps[1, :])
    quat_comps_beta[3, :] = (trans_comps[3, :] * quat_comps[0, :]
                             + trans_comps[0, :] * quat_comps[3, :]
                             - trans_comps[2, :] * quat_comps[1, :]
                             + trans_comps[1, :] * quat_comps[2, :])
    # swap into positive hemisphere if required
    quat_comps_beta[:, quat_comps_beta[0, :] < 0] *= -1

    beta_quat_array = np.empty_like(ebsd_map.quatArray)
    beta_quat_array[variant_map < 0] = Quat(1, 0, 0, 0)
    for i, idx in enumerate(zip(*np.where(variant_map >= 0))):
        beta_quat_array[idx] = Quat(quat_comps_beta[:, i])

    return beta_quat_array


def create_beta_ebsd_map(
    ebsd_map: ebsd.Map,
    mode: str = 'only_beta',
    beta_quat_array: np.ndarray = None,
    variant_map: np.array = None,
    alpha_phase_id: int = 0,
    beta_phase_id: int = 1,
) -> ebsd.Map:
    """

    Parameters
    ----------
    ebsd_map
    mode
        How to copy data from the input map
            'alone': Only include the reconstructed beta
            'append': Append reconstructed beta to present beta phase
            'add': Create a new phase for reconstructed beta
    beta_quat_array
    variant_map
    alpha_phase_id
    beta_phase_id
    """
    if variant_map is None:
        variant_map = construct_variant_map(
            ebsd_map, alpha_phase_id=alpha_phase_id
        )
    if beta_quat_array is None:
        beta_quat_array = construct_beta_quat_array(
            ebsd_map, variant_map=variant_map
        )

    if mode == 'alone':
        # Create map with only the reconstructed beta
        new_phase = copy.copy(ebsd_map.phases[beta_phase_id])
        new_phase.name += " (recon)"
        phases = [new_phase]

        out_phase_array = np.zeros_like(ebsd_map.phaseArray)
        out_phase_array[variant_map >= 0] = 1

        out_quat_array = beta_quat_array

    elif mode == 'append':
        # Append reconstructed beta to original beta phase
        phases = copy.copy(ebsd_map.phases)

        out_phase_array = np.copy(ebsd_map.phaseArray)
        out_phase_array[variant_map >= 0] = beta_phase_id + 1

        out_quat_array = np.where(variant_map >= 0, beta_quat_array,
                                  ebsd_map.quatArray)

    elif mode == 'add':
        # Create addition phase for the reconstructed beta
        phases = copy.copy(ebsd_map.phases)
        new_phase = copy.copy(ebsd_map.phases[beta_phase_id])
        new_phase.name += " (recon)"
        phases.append(new_phase)

        out_phase_array = np.copy(ebsd_map.phaseArray)
        out_phase_array[variant_map >= 0] = ebsd_map.numPhases + 1

        out_quat_array = np.where(variant_map >= 0, beta_quat_array,
                                  ebsd_map.quatArray)

    else:
        raise ValueError(f"Unknown beta map construction mode '{mode}'")

    out_euler_array = np.zeros((3,) + ebsd_map.shape)
    for i in range(ebsd_map.yDim):
        for j in range(ebsd_map.xDim):
            out_euler_array[:, i, j] = out_quat_array[i, j].eulerAngles()

    beta_ebsd_data = {
        'stepSize': ebsd_map.stepSize,
        'phases': phases,
        'phase': out_phase_array,
        'eulerAngle': out_euler_array,
        'bandContrast': ebsd_map.bandContrastArray
    }

    # TODO: Change so quats can be loaded instead of going via Euler angles
    beta_map = ebsd.Map(beta_ebsd_data, dataType="PythonDict")
    beta_map.quatArray = out_quat_array

    return beta_map


def do_reconstruction(
    ebsd_map: ebsd.Map,
    mode: str = 'average',
    burg_tol: float = 5,
    ori_tol: float = 3,
    alpha_phase_id: int = 0,
    beta_phase_id: int = 1
):
    """Apply beta reconstruction to a ebsd map object.

    The reconstructed beta map is stored directly in the ebsd map (this
    should probably change)

    Parameters
    ----------
    ebsd_map:
        EBSD map to apply reconstruction to
    mode
        How to perform reconstruction
            'average': grain average orientations
            'boundary': grain boundary orientations
            'beta': retained beta
    burg_tol
        Maximum deviation from the Burgers relation to allow (degrees)
    ori_tol: float
        Maximum deviation from a beta orientation (degrees)
    alpha_phase_id: int
        Index of the alpha phase in the EBSD map.
    beta_phase_id: int
        Index of the beta phase in the EBSD map.

    """
    # this is the only function that interacts with the ebsd map/grain objects
    alpha_grains = [grain for grain in ebsd_map
                    if grain.phaseID == alpha_phase_id]
    first = True
    for grain in tqdm(alpha_grains):

        beta_oris = calc_beta_oris(grain.refOri)
        variant_count = np.zeros(6, dtype=int)

        if mode == 'boundary':
            if first:
                print("Using boundary mode.")
                first = False
            possible_beta_oris, beta_deviations, alpha_oris = \
                calc_beta_oris_from_boundary_misori(
                    grain, ebsd_map.neighbourNetwork, ebsd_map.quatArray,
                    alpha_phase_id, burg_tol=burg_tol
                )

            for possible_beta_ori, beta_deviation, alpha_ori in zip(
                    possible_beta_oris, beta_deviations, alpha_oris):

                beta_oris_l = calc_beta_oris(alpha_ori)

                variant_count += count_beta_variants(
                    beta_oris_l, [possible_beta_ori], ori_tol
                )

        elif mode == 'beta':
            if first:
                print("Using beta mode.")
                first = False
            neighbour_grains = ebsd_map.neighbourNetwork.neighbors(grain)
            neighbour_oris = [[grain.refOri] for grain in neighbour_grains
                              if grain.phaseID == beta_phase_id]

            possible_beta_oris = neighbour_oris
            beta_deviations = [0.] * len(neighbour_oris)

            variant_count += count_beta_variants(
                beta_oris, possible_beta_oris, ori_tol
            )

        elif mode == 'average':
            if first:
                print("Using average mode.")
                first = False
            neighbour_grains = ebsd_map.neighbourNetwork.neighbors(grain)
            neighbour_oris = [grain.refOri for grain in neighbour_grains
                              if grain.phaseID == alpha_phase_id]

            # determine the possible beta orientations based on misorientation
            # between neighbouring alpha grains
            possible_beta_oris, beta_deviations = calc_beta_oris_from_misori(
                grain.refOri, neighbour_oris, burg_tol=burg_tol
            )

            variant_count += count_beta_variants(
                beta_oris, possible_beta_oris, ori_tol
            )

        else:
            raise ValueError(f"Unknown reconstruction mode '{mode}'")

        # save results in the grain objects
        grain.beta_oris = beta_oris
        grain.possible_beta_oris = possible_beta_oris
        grain.beta_deviations = beta_deviations
        grain.variant_count = variant_count
