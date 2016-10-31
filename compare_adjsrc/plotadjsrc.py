#!/usr/bin/python
from __future__ import (absolute_import, division, print_function)

# import os
# import numpy as np
# import sys
# import subprocess
# import argparse
# import json
# from obspy import read
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def plot_adjsrc_onepair_withdata(ref, syn, obsd, synt, wins, adjsrc_type,
                                 filename):

    left_window = 60000
    right_window = 0
    for win in wins:
        left_window = min(left_window, win[0])
        right_window = max(right_window, win[1])

    buffwin = 600
    x_range = obsd.stats.endtime - obsd.stats.starttime
    left_window = max(0, left_window - buffwin)
    right_window = min(x_range, right_window + buffwin)

    print(left_window, right_window)

    plt.figure(figsize=(12, 6))
    plt.subplot(211)
    plt.plot(obsd.times(), obsd.data, color='0', label="Observed", lw=2)
    plt.plot(synt.times(), synt.data, color='r', label="Synthetic", lw=2)
    ylim = max(map(abs, plt.ylim()))

    for win in wins:
        re = patches.Rectangle((win[0], plt.ylim()[0]),
                               win[1] - win[0], plt.ylim()[1] - plt.ylim()[0],
                               color="blue", alpha=0.1)
        plt.gca().add_patch(re)

    plt.grid()
    plt.legend(fancybox=True, framealpha=0.5)
    plt.ylim(-ylim, ylim)
    plt.xlim(left_window, right_window)

    plt.subplot(212)
    plt.plot(synt.times(), ref.adjoint_source, color='b',
             label='Pyadjoint (misfit = %13.5le)' % ref.misfit, lw=3)
    plt.plot(synt.times(), syn.adjoint_source, color='r', ls='--',
             label='measure_adj (misfit = %13.5le)' % syn.misfit, lw=2)
    ylim = max(map(abs, plt.ylim()))

    for win in wins:
        re = patches.Rectangle((win[0], plt.ylim()[0]),
                               win[1] - win[0], plt.ylim()[1] - plt.ylim()[0],
                               color="blue", alpha=0.1)
        plt.gca().add_patch(re)

    plt.grid()
    plt.legend(fancybox=True, framealpha=0.5)
    plt.ylim(-ylim, ylim)
    plt.xlim(left_window, right_window)
    plt.title('%s adjoint source' % adjsrc_type)
    plt.xlabel('Time (sec)')
    plt.savefig(filename, bbox_inches='tight', dpi=300)


def plot_adjsrc_onepair(ref, syn, wins, adjsrc_type, filename):

    left_window = 60000
    right_window = 0
    for win in wins:
        left_window = min(left_window, win[0])
        right_window = max(right_window, win[1])

    x_range = ref[-1,0] - ref[0,0]
    buffwin = (right_window - left_window) * 0.3
    left_window = max(0, left_window - buffwin)
    right_window = min(x_range, right_window + buffwin)

    # print(left_window, right_window)

    plt.figure(figsize=(12, 3))
    plt.plot(ref[:,0], ref[:,1], color='b',
             label='Pyadjoint', lw=3)
             # label='Pyadjoint (misfit = %13.5le)' % ref.misfit, lw=3)
    plt.plot(syn[:,0], syn[:,1], color='r', ls='--',
             label='measure_adj', lw=2)
             # label='measure_adj (misfit = %13.5le)' % syn.misfit, lw=2)
    ylim = max(map(abs, plt.ylim()))

    for win in wins:
        re = patches.Rectangle((win[0], plt.ylim()[0]),
                               win[1] - win[0], plt.ylim()[1] - plt.ylim()[0],
                               color="gray", alpha=0.05)
        plt.gca().add_patch(re)

    plt.grid()
    plt.legend(fancybox=True, framealpha=0.5)
    plt.ylim(-ylim, ylim)
    plt.xlim(left_window, right_window)
    plt.title('%s adjoint source' % adjsrc_type)
    plt.xlabel('Time (sec)')
    plt.savefig(filename, bbox_inches='tight', dpi=300)


def plot_seismogram_onepair(obsd, synt, wins, adjsrc_type, filename):

    left_window = 60000
    right_window = 0
    for win in wins:
        left_window = min(left_window, win[0])
        right_window = max(right_window, win[1])

    x_range = synt.data[-1] - synt.data[0]
    buffwin = (right_window - left_window) * 0.3
    left_window = max(0, left_window - buffwin)
    right_window = min(x_range, right_window + buffwin)

    # print(left_window, right_window)

    plt.figure(figsize=(12, 3))
    plt.plot(obsd.times(), obsd.data, color='0', label="Observed", lw=2)
    plt.plot(synt.times(), synt.data, color='r', label="Synthetic", lw=2)
    ylim = max(map(abs, plt.ylim()))

    for win in wins:
        re = patches.Rectangle((win[0], plt.ylim()[0]),
                               win[1] - win[0], plt.ylim()[1] - plt.ylim()[0],
                               color="gray", alpha=0.05)
        plt.gca().add_patch(re)

    plt.grid()
    plt.legend(fancybox=True, framealpha=0.5)
    plt.ylim(-ylim, ylim)
    plt.xlim(left_window, right_window)
    plt.title("")
    plt.xlabel('Time (sec)')
    plt.savefig(filename, bbox_inches='tight', dpi=300)

def plot_adjsrc_txt(adjsrc, adjsrc_type, filename):

    x_range = adjsrc[-1,0] - adjsrc[0,0]
    left_window = 0
    right_window = x_range

    plt.figure(figsize=(12, 3))
    plt.plot(adjsrc[:,0], adjsrc[:,1], color='b',
             label='adjoint source', lw=3)
    ylim = max(map(abs, plt.ylim()))

    plt.grid()
    plt.legend(fancybox=True, framealpha=0.5)
    plt.ylim(-ylim, ylim)
    plt.xlim(left_window, right_window)
    plt.title('%s adjoint source' % adjsrc_type)
    plt.xlabel('Time (sec)')
    plt.show()
    plt.savefig(filename, bbox_inches='tight', dpi=300)
