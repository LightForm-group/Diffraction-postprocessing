import numpy as np
from matflow import utils
from defdap.quat import Quat
from defdap.crystal import CrystalStructure, Phase

def lattice_strain(workflow, phases, axis, tol=5):
    tensor_comp, unit_vector = utils.tensorcomps_fromaxis(axis)

    latticestrain = {}
    plane_intensity = {}
    for phase_name, phase in phases.items():
        print(f"Processing phase {phase_name}...")

        # define volume_element_response data from simulation...
        ve_response = workflow.tasks.simulate_volume_element_loading.elements[0].outputs.volume_element_response

        # get indicies for desired phase...
        phase_idx = ve_response['field_data']['phase']['meta']['phase_names'].index(phase_name)
        phase_mask = ve_response['field_data']['phase']['data'] == phase_idx

        # Using left Cauchy-Green defomation tensor for elastic strain...
        Ee = ve_response['field_data']['epsilon_V^2(F_e)']['data'][:, phase_mask, :, :]
        Ee_incs = ve_response['field_data']['epsilon_V^2(F_e)']['meta']['increments']

        # get number of increments...
        ori_in = ve_response['field_data']['O']['data']
        ori_incs = ve_response['field_data']['O']['meta']['increments']

        # ensure increments are correct
        assert Ee_incs == ori_incs
        incs = ori_incs
        assert ori_in['type'] == 'quat'
        assert ori_in['quat_component_ordering'] == 'scalar-vector'
        assert ori_in['unit_cell_alignment']['x'] == 'a' # SAME AS DAMASK

        # get orientation data for phase...
        quat_comps = ori_in['quaternions'][:, phase_mask, :]
        # convert to P = +1 convention
        if ori_in['P'] == -1:
            quat_comps[:, :, 1:] *= -1
        ori = np.empty(quat_comps.shape[:-1], dtype=object)
        for i, row in enumerate(quat_comps):
            for j, val in enumerate(row):
                ori[i, j] = Quat(*val) 

        latticestrain_phase = {}
        plane_intensity_phase = {}
        for plane_label, crystal_plane in phase.diffraction_planes.items():
            print(f"Processing plane {plane_label} in direction {unit_vector},[{tensor_comp-1},{tensor_comp-1}]...")
            latticestrain_plane = []
            plane_intensity_plane = []
            for inc in incs:
                try:
                    inc_idx = incs.index(inc)
                except ValueError:
                    print(f"Increment {inc} does not exist.")
                    continue
                # determine lit up material points
                lit_up = find_illuminated_points(ori[inc_idx], crystal_plane, phase, incs, unit_vector, tol)
                # filter data by illuminated points
                unit_vector /= np.sqrt(np.dot(unit_vector, unit_vector))
                # assert np.allclose(measure_dir, measure_dir), f'expect measure dir to be [1, 0, 0]'
                latticestrain_inc = Ee[inc_idx, lit_up, tensor_comp-1, tensor_comp-1]
            
                # append lattice strain for plane
                latticestrain_plane.append(latticestrain_inc)
                # account for multiplicity of planes here?
                plane_intensity_plane.append(np.count_nonzero(lit_up)/6)
            # append lattice strain for plane to phase
            latticestrain_phase[plane_label] = latticestrain_plane
            plane_intensity_phase[plane_label] = plane_intensity_plane
        # append lattice strain for phase to overall
        latticestrain[phase_name] = latticestrain_phase
        plane_intensity[phase_name] = plane_intensity_phase

    return latticestrain, plane_intensity


def find_illuminated_points(oris, crystal_plane, phase, incs, measure_dir, tol=5, batch_size=100000):
    """
    Find the crystal orientations that will contribute to the diffraction peak of
    a given crystal plane for a given load direction.
    
    Parameters
    ----------
    oris : numpy.ndarray(defdap.quat.Quat)
        Orientation of each point
    crystal_plane : numpy.ndarray
        Crystal plane expressed in Miller-Bravais indicies
    phase : defdap.crystal.Phase
        Defdap phase
    measure_dir : numpy.ndarray
        Direction strain is measured in. Currently only working for z direction (0,0,1)
    tolerance : float
        Tolerance for meeting the Bragg condition, in degrees
    batch_size : int, optional
        Number of orientations to process at one time
        
    """
    sym_group = phase.crystalStructure.name
    
    if sym_group == 'hexagonal':
        crystal_plane_ortho = transform_plane(crystal_plane, phase.cOverA) # only for HCP
    else:
        crystal_plane_ortho = crystal_plane.astype(float)
        crystal_plane_ortho /= np.sqrt(np.dot(crystal_plane_ortho, crystal_plane_ortho))
        
    measure_dir /= np.sqrt(np.dot(measure_dir, measure_dir))
    
    tol = np.cos(tol * np.pi / 180.)

    num_points = oris.shape[0]
    num_batchs = num_points // batch_size + 1

    lit_up = np.zeros(num_points, dtype=bool)

    for i in range(num_batchs):
        oris_batch = oris[i*batch_size:(i+1)*batch_size]
        
        # transform load_dir to crystal coords for all orientations and symmetries

        quat_comps_sym = Quat.calcSymEqvs(oris_batch, sym_group)
        
        # temp variables to use below
        quat_dot_vec = (quat_comps_sym[:, 1, :] * measure_dir[0] +
                        quat_comps_sym[:, 2, :] * measure_dir[1] +
                        quat_comps_sym[:, 3, :] * measure_dir[2])
        temp = (np.square(quat_comps_sym[:, 0, :]) - np.square(quat_comps_sym[:, 1, :]) -
                np.square(quat_comps_sym[:, 2, :]) - np.square(quat_comps_sym[:, 3, :]))
        
        # array to store crytal directions for all orientations and symmetries
        measure_dir_crystal = np.empty((3, quat_comps_sym.shape[0], quat_comps_sym.shape[2]))

        # (quat_comps_sym * measure_dir) * quat_comps_sym.conjugate
        measure_dir_crystal[0] = (2 * quat_dot_vec * quat_comps_sym[:, 1, :] +
                                  temp * measure_dir[0] +
                                  2 * quat_comps_sym[:, 0, :] * (quat_comps_sym[:, 2, :] * measure_dir[2] -
                                                                 quat_comps_sym[:, 3, :] * measure_dir[1]))
        measure_dir_crystal[1] = (2 * quat_dot_vec * quat_comps_sym[:, 2, :] +
                                  temp * measure_dir[1] +
                                  2 * quat_comps_sym[:, 0, :] * (quat_comps_sym[:, 3, :] * measure_dir[0] -
                                                                 quat_comps_sym[:, 1, :] * measure_dir[2]))
        measure_dir_crystal[2] = (2 * quat_dot_vec * quat_comps_sym[:, 3, :] +
                                  temp * measure_dir[2] +
                                  2 * quat_comps_sym[:, 0, :] * (quat_comps_sym[:, 1, :] * measure_dir[1] -
                                                                 quat_comps_sym[:, 2, :] * measure_dir[0]))

        # normalise vectors (just in case)
        measure_dir_crystal /= np.sqrt(np.einsum('ijk,ijk->jk', measure_dir_crystal, measure_dir_crystal))

        dot_prod = np.einsum('ijk,i->jk', measure_dir_crystal, crystal_plane_ortho)
        lit_up[i*batch_size:(i+1)*batch_size] = np.any(np.abs(dot_prod) > tol, axis=0)

#     print("\rDone. {:} points of {:} will contribute to diffraction peak.".format(np.count_nonzero(lit_up), num_points))
    
    return lit_up


def find_illuminated_points(oris, crystal_plane, phase, incs, measure_dir, tol=5, batch_size=100000):
    """
    Find the crystal orientations that will contribute to the diffraction peak of
    a given crystal plane for a given load direction.
    
    Parameters
    ----------
    oris : numpy.ndarray(defdap.quat.Quat)
        Orientation of each point
    crystal_plane : numpy.ndarray
        Crystal plane expressed in Miller-Bravais indicies
    phase : defdap.crystal.Phase
        Defdap phase
    measure_dir : numpy.ndarray
        Direction strain is measured in. Currently only working for z direction (0,0,1)
    tolerance : float
        Tolerance for meeting the Bragg condition, in degrees
    batch_size : int, optional
        Number of orientations to process at one time
        
    """
    sym_group = phase.crystalStructure.name
    
    if sym_group == 'hexagonal':
        crystal_plane_ortho = transform_plane(crystal_plane, phase.cOverA) # only for HCP
    else:
        crystal_plane_ortho = crystal_plane.astype(float)
        crystal_plane_ortho /= np.sqrt(np.dot(crystal_plane_ortho, crystal_plane_ortho))
        
    measure_dir /= np.sqrt(np.dot(measure_dir, measure_dir))
    
    tol = np.cos(tol * np.pi / 180.)

    num_points = oris.shape[0]
    num_batchs = num_points // batch_size + 1

    lit_up = np.zeros(num_points, dtype=bool)

    for i in range(num_batchs):
        oris_batch = oris[i*batch_size:(i+1)*batch_size]
        
        # transform load_dir to crystal coords for all orientations and symmetries

        quat_comps_sym = Quat.calcSymEqvs(oris_batch, sym_group)
        
        # temp variables to use below
        quat_dot_vec = (quat_comps_sym[:, 1, :] * measure_dir[0] +
                        quat_comps_sym[:, 2, :] * measure_dir[1] +
                        quat_comps_sym[:, 3, :] * measure_dir[2])
        temp = (np.square(quat_comps_sym[:, 0, :]) - np.square(quat_comps_sym[:, 1, :]) -
                np.square(quat_comps_sym[:, 2, :]) - np.square(quat_comps_sym[:, 3, :]))
        
        # array to store crytal directions for all orientations and symmetries
        measure_dir_crystal = np.empty((3, quat_comps_sym.shape[0], quat_comps_sym.shape[2]))

        # (quat_comps_sym * measure_dir) * quat_comps_sym.conjugate
        measure_dir_crystal[0] = (2 * quat_dot_vec * quat_comps_sym[:, 1, :] +
                                  temp * measure_dir[0] +
                                  2 * quat_comps_sym[:, 0, :] * (quat_comps_sym[:, 2, :] * measure_dir[2] -
                                                                 quat_comps_sym[:, 3, :] * measure_dir[1]))
        measure_dir_crystal[1] = (2 * quat_dot_vec * quat_comps_sym[:, 2, :] +
                                  temp * measure_dir[1] +
                                  2 * quat_comps_sym[:, 0, :] * (quat_comps_sym[:, 3, :] * measure_dir[0] -
                                                                 quat_comps_sym[:, 1, :] * measure_dir[2]))
        measure_dir_crystal[2] = (2 * quat_dot_vec * quat_comps_sym[:, 3, :] +
                                  temp * measure_dir[2] +
                                  2 * quat_comps_sym[:, 0, :] * (quat_comps_sym[:, 1, :] * measure_dir[1] -
                                                                 quat_comps_sym[:, 2, :] * measure_dir[0]))

        # normalise vectors (just in case)
        measure_dir_crystal /= np.sqrt(np.einsum('ijk,ijk->jk', measure_dir_crystal, measure_dir_crystal))

        dot_prod = np.einsum('ijk,i->jk', measure_dir_crystal, crystal_plane_ortho)
        lit_up[i*batch_size:(i+1)*batch_size] = np.any(np.abs(dot_prod) > tol, axis=0)

#     print("\rDone. {:} points of {:} will contribute to diffraction peak.".format(np.count_nonzero(lit_up), num_points))
    
    return lit_up


def transform_plane(crystal_plane, c_over_a):
    """Transform crystal plane expressed in Miller-Bravais indicies to a unit vector in orthonormal basis."""
    # Convert plane and dir from Miller-Bravais to Miller indicies
    crystal_plane_m = crystal_plane[[0, 1, 3]] # (M-B hkil -> M hkl)

    # Create L matrix. Transformation from crystal to orthonormal coords
    l_matrix = CrystalStructure.lMatrix(1, 1, c_over_a, np.pi / 2, np.pi / 2, np.pi * 2 / 3)

    # Create Q matrix for transforming planes to orthonormal coords
    q_matrix = CrystalStructure.qMatrix(l_matrix)

    # Transform into orthonormal basis and then normalise
    crystal_plane_ortho = np.matmul(q_matrix, crystal_plane_m)
    crystal_plane_ortho /= np.sqrt(np.dot(crystal_plane_ortho, crystal_plane_ortho))

    return crystal_plane_ortho


def process_strain(Ee, lit_up, measure_dir, comp):
    """
    Filter the strain data by the illuminated points.
    
    Parameters
    ----------
    Ee : numpy.ndarray
        Elastic strain tensor at each point, shape (num_points, 3, 3)
    lit_up : numpy.ndarray(bool)
        Boolean array indicating if each point is illuminated, shape (num_points,)
    measure_dir : numpy.ndarray
        Direction strain is to be measured in. Currently only working for z direction (0,0,1)
        
    """    
    measure_dir /= np.sqrt(np.dot(measure_dir, measure_dir))
#     assert np.allclose(measure_dir, measure_dir), f'expect measure dir to be [1, 0, 0]'
    
    measured_strain = Ee[lit_up, comp, comp]
    
    return measured_strain