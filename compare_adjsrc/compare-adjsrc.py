#!/bin/env python

import os  # sys
import argparse
import numpy as np
import json
from plotadjsrc import plot_adjsrc_onepair, plot_seismogram_onepair
from obspy import read

def compare_error(ref, syn):
    """
    calculate difference between two seismogram/adjoint source
    """

    corr_min = 1.
    err_max = 0.

    # correlation test
    corr_mat = np.corrcoef(ref, syn)
    corr = np.min(corr_mat)
    corr_min = min(corr, corr_min)

    # least square test
    norm = np.linalg.norm
    sqrt = np.sqrt
    err = norm(ref-syn)/sqrt(norm(ref)*norm(syn))
    err_max = max(err, err_max)

    return corr_min, err_max


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', action='store', dest='adjointsrc', required=True)
    parser.add_argument('-s', action='store', dest='seismogram', required=True)
    parser.add_argument('-w', action='store', dest='windowfile', required=True)
    parser.add_argument('-v', action='store_true', dest='verbose')
    args = parser.parse_args()

    figure = True

    tol_corr = 0.99
    tol_err = 3.E-02

    def _read_file_json(json_file):
        with open(json_file) as fh:
            return json.load(fh)

    adjointsrc = _read_file_json(args.adjointsrc)
    windows = _read_file_json(args.windowfile)
    seismogram = _read_file_json(args.seismogram)

    for key in adjointsrc:
        sta, ntwk, locid, dummy = key.split(".")
        sta_key = "%s.%s" % (sta, ntwk)
        wins = windows[sta_key][key]

        fobsd = seismogram[key]["obsd"]
        fsynt = seismogram[key]["synt"]

        obsd = read(fobsd)[0]
        synt = read(fsynt)[0]

        window = []
        for win in wins:
            wind = [win["relative_starttime"], win["relative_endtime"]]
            window.append(wind)

        for case in adjointsrc[key]:
            ref_file = case["ref"].encode('utf-8')
            syn_file = case["syn"].encode('utf-8')
            tag_plot = case["tag"].encode('utf-8')

            ref_tr = np.loadtxt(ref_file)
            syn_tr = np.loadtxt(syn_file)

            ref = ref_tr[:, 1]
            syn = syn_tr[:, 1]

            corr_min, err_max = compare_error(ref, syn)

            if corr_min >= tol_corr and err_max <= tol_err:
                ifpass = True
            else:
                ifpass = False

            fname = os.path.basename(syn_file)

            # print results to screen
            print("|%20s| %13.5le| %13.5le| %s" % (fname,
                                                   corr_min,
                                                   err_max,
                                                   str(ifpass)))

            if figure:
                filename = "%s.eps" % fname
                adjsrc_type = tag_plot
                plot_adjsrc_onepair(ref_tr, syn_tr, window, adjsrc_type,
                                    filename)

        if figure:
            filename = "%s.eps" % key
            plot_seismogram_onepair(obsd, synt, window, adjsrc_type, filename)
