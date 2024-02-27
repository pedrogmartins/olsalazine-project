import numpy as np
import matplotlib.pyplot as plt

# This file is adapted from script written by Alex Smith (graduate
#  student in Neaton Group at UC Berkeley). It extracts the NMR shifts
#  computed by VASP, and substracts the appropriate NMR shift references.
#  Shifts are saved in .csv files.

#References for PBE calculations (also work for PBE-D3 since
#  D3 doesn’t affect a static structure NMR calculations
#  from Forse, Alexander C., et al. Journal of the American Chemical
#  Society 140.51 (2018): 18016-18031.
ref_c = [-170.5]
ref_n = [-225.0]
ref_h = [-30.9]

# Define functions needed
def get_absolute_shifts(fn):
    ashift = np.loadtxt(fn, dtype=str)
    return [(a, float(b)) for (a,b) in ashift]

def get_shifts(shifts_list, elem = 'H', ref = ref_h[0]):
    s = []
    for (atom, abs_shift) in shifts_list:
        if atom == elem:
            s.append(-(ref - abs_shift))
    return s

def get_c_rel(fn): return get_shifts(get_absolute_shifts(fn), elem='C', ref=ref_c[0])
def get_h_rel(fn): return get_shifts(get_absolute_shifts(fn), elem='H', ref=ref_h[0])
def get_n_rel(fn): return get_shifts(get_absolute_shifts(fn), elem='N', ref=ref_n[0])


def plot_shifts_separate(shifts, colors, labels, elem):
    min_shift = 1000
    max_shift = -1000
    for s in shifts:
        for i in s:
            if i < min_shift:
                min_shift = i
            if i > max_shift:
                max_shift = i
    #print(min_shift, max_shift)
    fig, axs = plt.subplots(len(colors), figsize = [10.0, 6.4])
    for (i,(s, c, l)) in enumerate(zip(shifts, colors, labels)):
        #plt.eventplot(s, orientation=‘horizontal’, colors=c, label = l)
        axs[i].eventplot(s, orientation='horizontal', colors=c, label = l)
        axs[i].set_xlim([min_shift-0.3, max_shift+0.3])
    for i in range(len(colors)):
        axs[i].legend(loc=0)
    fig.suptitle('Relative ' + elem + ' Chemical Shift')
    fig.tight_layout()
    #plt.savefig(‘s2.png’, dpi=300)

#Extract NMR shifts for bare structure
bare_c = get_c_rel('bare.txt')
bare_h = get_h_rel('bare.txt')
bare_n = get_n_rel('bare.txt')

#Save numpy arrays as .csv files
np.asarray(bare_c).tofile('bare_c.csv', sep = ',')
np.asarray(bare_h).tofile('bare_h.csv', sep = ',')
np.asarray(bare_n).tofile('bare_n.csv', sep = ',')
