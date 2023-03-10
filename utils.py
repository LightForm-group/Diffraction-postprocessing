import yaml
import numpy as np
from defdap.crystal import Phase


def tensorcomps_fromaxis(axis):
    """
    Given axis direction as string eg. "X"
    return corresponding tensor indecies eg. "11".
    return unit vector as numpy array eg. [1,0,0].
    Used for extracting tensor components in post-processing.
    """
    if axis=="X":
        tensor_comp = 1
        unit_vector = np.array([1.,0.,0.])
    if axis=="Y":
        tensor_comp = 2
        unit_vector = np.array([0.,1.,0.])
    if axis=="Z":
        tensor_comp = 3
        unit_vector = np.array([0.,0.,1.])
    return tensor_comp, unit_vector


def import_diff_params(*args):
    """
    Given phase_labels as any amount of strings eg. "Ti_alpha", "Ti_beta", ...
    Return phases dict with each phase as defdap phase object
    containing crystallographic parameters eg. diffraction planes.
    """
    yaml_file = open("./diffraction_parameters.yaml", 'r') # define phase parameters database path

    yaml_dict = yaml.safe_load(yaml_file) # use yaml to load database
    phases=dict();                        # define phases dict (function output)
    for phase in args:
        # write lattice parameters for chosen phase into defdap object
        name = yaml_dict[phase]['Phase']['name']
        luaeGroup = yaml_dict[phase]['Phase']['laueGroup']
        spaceGroup = yaml_dict[phase]['Phase']['spaceGroup']
        latticeParams = yaml_dict[phase]['Phase']['latticeParams']
        phases[phase] = Phase(name, luaeGroup, spaceGroup, latticeParams)

        # write diffraction planes for chosen phase into defdap object
        # planes in yaml in cryst direction convention[] should be planes{} in phases string keys
        phases[phase].diffraction_planes = {}
        planes = yaml_dict[phase]['diffraction_planes']
        for plane in planes:
            plane_str = '{'+str(plane).strip('[]').replace(',','').replace(' ','')+'}'
            plane_array = np.array(plane)
            phases[phase].diffraction_planes[plane_str] = plane_array
    yaml_file.close()
    return phases