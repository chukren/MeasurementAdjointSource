import os
from obspy import read, read_inventory, readEvents
# from obspy.core.event import readEvents
import pyadjoint
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.transforms as trf
from pytomo3d.adjoint.adjsrc import postprocess_adjsrc
import logging
logger = logging.getLogger("pyadjoint")
logger.setLevel(logging.DEBUG)

obsd_dir = "../adjsrc_fortran/sac/C201009031635A.proc_obsd_60_100.sac"
synt_dir = "../adjsrc_fortran/sac/C201009031635A.proc_synt_60_100.sac"
quakeml = "./C201009031635A.xml"
stationxml = "./IU.KBL.xml"

src_type = "multitaper_misfit"
tag_type = "mt"

config = pyadjoint.ConfigMultiTaper(min_period=60.0, max_period=100.0,
                          lnpt=15,
                          transfunc_waterlevel=1.0E-10,
                          water_threshold=0.02,
                          ipower_costaper=10,
                          min_cycle_in_window=0.5,
                          taper_percentage=1.0,
                          taper_type='cos_p10',
                          mt_nw=4,
                          num_taper=5,
                          dt_fac=0.2,
                          phase_step=1.5,
                          err_fac=2.5,
                          dt_max_scale=3.5,
                          dt_sigma_min=1.0,
                          dlna_sigma_min=0.5,
                          use_cc_error=False,
                          use_mt_error=False)

event = readEvents(quakeml)[0]
inventory = read_inventory(stationxml)

for catlog in ("dt", "am"):

    print("========= %s ========" % catlog)
    config.measure_type = catlog
    fwindow_chi = "%s_%s_window_chi.txt" % (tag_type, catlog)
    fchi = open(fwindow_chi, 'w')

    adjsrc_all = []
    # R
    fobsd = os.path.join(obsd_dir, "IU.KBL..BHR.proc_obsd_60_100.sac")
    fsynt = os.path.join(synt_dir, "IU.KBL.S3.MXR.proc_synt_60_100.sac")

    obs = read(fobsd)
    syn = read(fsynt)

    fadjsrc_r = "IU.KBL..BHR.%s.%s.adj" % (tag_type, config.measure_type)
    window_r = [[4032., 4304.6]]
    adjsrc_r =\
        pyadjoint.calculate_adjoint_source(adj_src_type=src_type,
                                           observed=obs[0], synthetic=syn[0],
                                           config=config, window=window_r,
                                           adjoint_src=True, plot=False)

    adjsrc_r.write(filename=fadjsrc_r, format="SPECFEM", time_offset=0)
    fchi.write("%s %f\n" % (fadjsrc_r, adjsrc_r.misfit))
    print("%s %s tr_chi %f am_chi %f" %
          (fadjsrc_r, adjsrc_r.measurement[0]["type"],
           adjsrc_r.measurement[0]["misfit_dt"],
           adjsrc_r.measurement[0]["misfit_dlna"]))
    adjsrc_all.append(adjsrc_r)

    # T
    fobsd = os.path.join(obsd_dir, "IU.KBL..BHT.proc_obsd_60_100.sac")
    fsynt = os.path.join(synt_dir, "IU.KBL.S3.MXT.proc_synt_60_100.sac")

    obs = read(fobsd)
    syn = read(fsynt)

    fadjsrc_t = "IU.KBL..BHT.%s.%s.adj" % (tag_type, config.measure_type)
    window = [[3008.2, 3345.8]]

    adjsrc_t =\
        pyadjoint.calculate_adjoint_source(adj_src_type=src_type,
                                           observed=obs[0], synthetic=syn[0],
                                           config=config, window=window,
                                           adjoint_src=True, plot=False)
    adjsrc_t.write(filename=fadjsrc_t, format="SPECFEM", time_offset=0)
    fchi.write("%s %f\n" % (fadjsrc_t, adjsrc_t.misfit))
    print("%s %s tr_chi %f am_chi %f" %
          (fadjsrc_t, adjsrc_t.measurement[0]["type"],
           adjsrc_t.measurement[0]["misfit_dt"],
           adjsrc_t.measurement[0]["misfit_dlna"]))
    adjsrc_all.append(adjsrc_t)

    # Z
    fobsd = os.path.join(obsd_dir, "IU.KBL..BHZ.proc_obsd_60_100.sac")
    fsynt = os.path.join(synt_dir, "IU.KBL.S3.MXZ.proc_synt_60_100.sac")

    obs = read(fobsd)
    syn = read(fsynt)

    fadjsrc_z = "IU.KBL..BHZ.%s.%s.adj" % (tag_type, config.measure_type)
    window = [[3313.6, 3756.0]]

    adjsrc_z =\
        pyadjoint.calculate_adjoint_source(adj_src_type=src_type,
                                           observed=obs[0], synthetic=syn[0],
                                           config=config, window=window,
                                           adjoint_src=True, plot=False)
    adjsrc_z.write(filename=fadjsrc_z, format="SPECFEM", time_offset=0)
    fchi.write("%s %f\n" % (fadjsrc_z, adjsrc_z.misfit))
    print("%s %s tr_chi %f am_chi %f" %
          (fadjsrc_z, adjsrc_z.measurement[0]["type"],
           adjsrc_z.measurement[0]["misfit_dt"],
           adjsrc_z.measurement[0]["misfit_dlna"]))
    adjsrc_all.append(adjsrc_z)

    # print(adjsrc_t.starttime,adjsrc_r.starttime,adjsrc_z.starttime)
    # print(adjsrc_t.dt,adjsrc_r.dt,adjsrc_z.dt)
    # print(len(adjsrc_t.adjoint_source),len(adjsrc_r.adjoint_source),
    #    len(adjsrc_z.adjoint_source))

    interp_starttime = adjsrc_t.starttime
    interp_delta = adjsrc_t.dt
    interp_npts = len(adjsrc_t.adjoint_source)

    adjsrc_final =\
        postprocess_adjsrc(adjsrc_all, interp_starttime, interp_delta,
                           interp_npts, rotate_flag=True, inventory=inventory,
                           event=event, sum_over_comp_flag=False,
                           weight_flag=False, weight_dict=None,
                           filter_flag=False, pre_filt=None)

    for adj in adjsrc_final:
        fadjsrc = "%s.%s.%s.%s.%s.%s.adj" %\
            (adj.network, adj.station, adj.location, adj.component, tag_type,
             config.measure_type)

        adj.write(filename=fadjsrc, format="SPECFEM", time_offset=0)

fchi.close()
