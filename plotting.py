
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import matplotlib.pyplot as plt


def plot_latticestrain(axs, true_stress, latticestrain, axis, incs="*", xlim=None, ylim=None):
    phase_names = [
             'Ti_alpha',
             'Ti_beta',
             ]

    custom_cycler = (cycler(color=[
                                   '#FFE800', '#9FDF00', '#0DCD52', '#00948D', # alpha colours
                                   '#67001F', '#DF2179', '#CDA0CD' # beta colours
                                  ]) +
                     cycler(marker=[
                                    's', 'h', '^', 'D', # alpha markers
                                    '+', 'x', 'P' # beta markers
                                   ]))
    axs.set_prop_cycle(custom_cycler)
    
    latticestrain_mean = {}
    for phase_name, latticestrain_phase in latticestrain[axis].items():
        latticestrain_phase_mean = {}
        for plane_label, strains in latticestrain_phase.items():
            latticestrain_phase_mean[plane_label] = np.array([strain.mean() for strain in strains])
        latticestrain_mean[phase_name] = latticestrain_phase_mean

    for phase in phase_names:
        for plane_label, mean_strain in latticestrain_mean[phase].items():
            total_microstrain = np.array([0] + mean_strain) * 1e6 # convert to microstrain (1e6)
            axs.plot(total_microstrain, true_stress, label=plane_label)

    axs.title.set_text(f"{axis}")
    axs.set_xlabel("Lattice Strain ($10^{-6}$)")
    axs.set_xlim([None, xlim])
    axs.set_ylabel("True Stress $\sigma$ (MPa)")
    axs.set_ylim([None, ylim])
    axs.legend()
    
    
def plot_truestrain_peakint(axs, plane_intensity, true_strain, axis, phase_labels, xlim=None, ylim=None):
    
    # maybe the colours are just wrong??
    custom_cycler = (cycler(color=[
                                   '#FFE800', '#9FDF00', '#0DCD52', '#00948D', # alpha
                                   '#67001F', '#DF2179', '#CDA0CD' # beta
                                  ]) +
                     cycler(marker=[
                                    's', 'h', '^', 'D', # alpha
                                    '+', 'x', 'P' # beta
                                   ]))
    
    custom_cycler = (cycler(color=[
                                   '#FFE800', '#9FDF00', '#0DCD52', '#00948D', # alpha
                                   '#67001F', '#DF2179', '#CDA0CD' # beta
                                  ]) +
                     cycler(marker=[
                                    's', 'h', '^', 'D', # alpha
                                    '+', 'x', 'P' # beta
                                   ]))
    
    axs.set_prop_cycle(custom_cycler)

    for phase in phase_labels:
        for plane_label, peakint in plane_intensity[axis][phase].items():
            
            # Elements of X-Ray Diffraction 2nd ed., B.D. Cullity, 1978 ISBN 0-201-01174-3
            # alpha
            if plane_label == "{0002}":
                peakint = [intensity / 2 for intensity in peakint]
            elif plane_label == "{10-10}":
                peakint = [intensity / 6 for intensity in peakint]
            elif plane_label == "{10-11}":
                peakint = [intensity / 12 for intensity in peakint]
            elif plane_label == "{11-20}":
                peakint = [intensity / 6 for intensity in peakint]
            # beta
            elif plane_label == "{200}":
                peakint = [intensity / 6 for intensity in peakint]
            elif plane_label == "{110}":
                peakint = [intensity / 12 for intensity in peakint]
                
#             ax.plot(true_strain, peakint, label=plane_label)
            axs.plot(true_strain, peakint, label=plane_label)
            
    axs.title.set_text(f"{axis}")
    axs.set_xlabel(f"True Strain $\epsilon_{axis}^{{VM}}$")
    axs.set_xlim([None, xlim])
    axs.set_ylabel("Material Point Count")
    axs.set_ylim([None, ylim])
    axs.legend()
    
    
def plot_lattice_strain_dist_inc(axs, latticestrain, axis, phase, inc, bins=20, xmin=None, xlim=None, ymin=None, ylim=None):
    
    if phase=="Ti_beta":
        colour=np.array(['#67001F', '#DF2179'], dtype='object')
        
    elif phase=="Ti_alpha":
        colour=np.array(['#FFE800', '#9FDF00', '#0DCD52', '#00948D'], dtype='object')

    print(f"phase: {phase} Direction: {axis}")
    for plane_num, plane in enumerate(latticestrain[axis][phase].keys()):
        print(f"\tplane: {plane}\tlattice strain ")
        
        # calculate how many numbers within stdev around mean
        lattstrain_dist = latticestrain[axis][phase][plane][inc]
        axs.hist(lattstrain_dist*1e6, bins=bins,
                alpha=0.4, color=colour[plane_num],
                label=plane)

    axs.legend()
    axs.title.set_text(f"{axis}")
    axs.set_xlabel("Lattice Strain")
    axs.set_xlim(xmin, xlim)
    axs.set_ylabel("Measurements")
    axs.set_ylim(ymin, ylim)
