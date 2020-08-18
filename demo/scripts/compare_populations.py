#!t()usr/bin/env python
# -*- coding: utf-8 -*-
#
# -----------------------------------------------------------------------------|
#   Header
# -----------------------------------------------------------------------------|
from    matplotlib import rc
import  matplotlib.pyplot as plt
import  numpy as np

rc('text', usetex=True)
rc('font', size=10, family='serif')

# -----------------------------------------------------------------------------|
#   Load data
# -----------------------------------------------------------------------------|
dat_f = np.loadtxt('dat/hot_tails_fortran.dat')
dat_m = np.loadtxt('dat/hot_tails_matlab.dat')
dat_p = np.loadtxt('dat/hot_tails_python.dat')

# -----------------------------------------------------------------------------|
#   Prepare plots
# -----------------------------------------------------------------------------|
color = ['k', 'b', 'r']
label_sources = ['matlab', 'python', 'fortran']

ax_titles = [
        'Hot-tail population', 
        'Electron temperature', 
        'Critical velocity',
        'Electron density', 
        r'$\tau$', 
        'Electric field'
        ]

y_labels = [
        r'$n_{\rm hot}$~(10$^{16}$~ m$^{-3}$)',
        r'$T_{\rm e}$~(eV)',
        r'$v_{\rm c}/v_{T_0}$',
        r'$n_{\rm e}$~(10$^{19}$~m$^{-3}$)',
        r'$\tau$',
        r'$E$~(V/m)',
        ]

y_labels_rel = [
        r'$\Delta n_{\rm hot}$',
        r'$\Delta T_{\rm e}$',
        r'$\Delta v_{\rm c}/v_{T_0}$',
        r'$\Delta n_{\rm e}$',
        r'$\Delta \tau$',
        r'$\Delta E$',
        ]

# -----------------------------------------------------------------------------|
#   Plot data - absolute values
# -----------------------------------------------------------------------------|
fig_1, axs = plt.subplots(3, 2, figsize=(7,6))
axs = fig_1.axes

mult = [1, 1e-16, 1e-19, 1, 1, 1, 1]
for j, dat in enumerate([dat_m, dat_p, dat_f]):
    # ----- Plot traces -------------------------------------------------------|
    for i, ind in enumerate([1, 3, 5, 2, 6, 1]):
        axs[i].plot(1e3*dat[:,0], mult[ind]*dat[:,ind], c=color[j], 
                label=label_sources[j])

    # ----- Plot settings -----------------------------------------------------|
    axs[1].semilogy()
    axs[1].set_ylim(bottom=1)

for i, ax in enumerate(axs):
    ax.legend(loc=0)
    ax.set_title(ax_titles[i])
    ax.set_xlabel(r'$t$~(ms)')
    ax.set_xlim((1e3*dat[0,0], 1e3*dat[-1,0]))
    ax.set_ylabel(y_labels[i])

    if i != 1:
        ax.set_ylim(bottom=0)

fig_1.tight_layout()

# -----------------------------------------------------------------------------|
#   Plot data - relative difference
# -----------------------------------------------------------------------------|
fig_2, axs = plt.subplots(3, 1, figsize=(3.5,6))
axs = fig_2.axes

for j, dat in enumerate([dat_m, dat_p]):
    if np.all(dat[:,0] != dat_f[:,0]): continue

    c = ['k', 'b'][j]

    # ----- Plot traces -------------------------------------------------------|
    for i, ind in enumerate([1, 5, 6]):
        axs[i].plot(1e3*dat_f[:,0], 
                np.abs(dat_f[:,ind]/dat[:,ind]-1), 
                c=c, label=label_sources[j])

    # ----- Plot settings -----------------------------------------------------|
for i, ax in enumerate(axs):
    ax.legend(loc=1, title='fortran vs.')
    ax.semilogy()
    ax.set_title(ax_titles[2*i])
    ax.set_xlabel(r'$t$~(ms)')
    ax.set_xlim((1e3*dat_f[0,0], 1e3*dat_f[-1,0]))
    ax.set_ylabel(y_labels_rel[2*i])
    ax.set_ylim((1e-12, 1))

fig_2.tight_layout()

# -----------------------------------------------------------------------------|
#   Show figures
# -----------------------------------------------------------------------------|
plt.show()

fig_1.savefig('hot_tails_evolution.png', dpi=300)
fig_2.savefig('hot_tails_evolution_difference.png', dpi=300)

# -----------------------------------------------------------------------------|
