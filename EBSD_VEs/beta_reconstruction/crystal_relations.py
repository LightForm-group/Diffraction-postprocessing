import numpy as np
from defdap.quat import Quat

hex_syms = Quat.symEqv("hexagonal")
# subset of hexagonal symmetries that give unique orientations when the
# Burgers transformation is applied
unq_hex_syms = [
    hex_syms[0],
    hex_syms[5],
    hex_syms[4],
    hex_syms[2],
    hex_syms[10],
    hex_syms[11]
]

cubic_syms = Quat.symEqv("cubic")
# subset of cubic symmetries that give unique orientations when the
# Burgers transformation is applied
unq_cub_syms = [
    cubic_syms[0],
    cubic_syms[7],
    cubic_syms[9],
    cubic_syms[1],
    cubic_syms[22],
    cubic_syms[16],
    cubic_syms[12],
    cubic_syms[15],
    cubic_syms[4],
    cubic_syms[8],
    cubic_syms[21],
    cubic_syms[20]
]

# HCP -> BCC
burg_eulers = np.array([135, 90, 354.74]) * np.pi / 180
burg_trans = Quat.fromEulerAngles(*burg_eulers).conjugate
