import multiprocessing
from functools import partial

import h5py
import matplotlib.pyplot as plt
import itertools
from scipy.optimize import minimize, curve_fit
import numpy as np
from matplotlib.colors import LogNorm
import os

from OffAnaClass import *


def fit_position_time(xyz,
                      times,
                      c,
                      tresid_pdf,
                      seed=None):
    '''
    Runs a position time fit on a collection of hits given a speed'o'light and
    callable object to evaluate negative log likelihood values for hit time
    residuals.

    Generates a seed using the average hit position and emission time most
    consistent with this position if seed is None.

    Runs simplex to find minima region, then default quasi Newton BFGS minimizer

    returns a MinimizerResult for parameters (x,y,z,t)
    '''

    def sqresid(par):  # unused
        x, y, z, t = par
        pos = np.asarray([x, y, z])
        P = xyz - pos
        D = np.sqrt(np.sum(np.square(P), axis=1))
        T = times - t
        tresid = T - D / c
        return np.sum(np.square(tresid))

    def nll_postime(par, verbose=True):
        x, y, z, t = par
        pos = np.asarray([x, y, z])
        P = xyz - pos
        D = np.sqrt(np.sum(np.square(P), axis=1))
        T = times - t
        tresid = T - D / c
        # this may not be correct

        return np.sum(tresid_pdf(tresid))

    # guess pos/time seed by minimizing time residuals
    if seed is None:
        pos_guess = np.mean(xyz, axis=0)
        P = xyz - pos_guess
        D = np.sqrt(np.sum(np.square(P), axis=1))
        mean = tresid_pdf.mean()
        t_guess = np.mean(times - D / c) - mean
        guess = np.concatenate([pos_guess, [t_guess]])
    else:
        guess = seed[:4]

    # find best minima coarsely with simplex
    m = minimize(nll_postime, x0=guess, method='Nelder-Mead')
    # compute proper minima
    m = minimize(nll_postime, x0=m.x)
    # print('stage1:',guess,m.x)

    return m


def fit_direction(prompt_hits, pos, costheta_pdf, seed=None):
    '''
    Runs a direction fit using a fixed position and a callable object to evaluate
    negative log likelihoods for cos(theta) values.

    Generates a seed using the average hit direction if seed is None

    returns MinimizerResult for parameters (theta,phi)
    '''

    def nll_dir(par, verbose=False):
        theta, phi = par
        P = prompt_hits - pos
        D = np.sqrt(np.sum(np.square(P), axis=1))
        dvec = np.asarray([np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)])
        costheta = np.sum(P * dvec, axis=1) / D
        return np.sum(costheta_pdf(costheta))

    # fit direction with pos/time fixed
    if seed is None:
        P = prompt_hits - pos
        D = np.sqrt(np.sum(np.square(P), axis=1))
        avg_dir = np.mean(P.T / D, axis=1)
        avg_dir = avg_dir / np.sqrt(np.sum(np.square(avg_dir)))
        costheta = avg_dir[2]
        cosphi = avg_dir[0] / np.sqrt(1 - costheta ** 2)
        guess = (np.arccos(costheta), np.arccos(cosphi))
    else:
        guess = seed[-2:]

    # find best minima coarsely with simplex
    m_dir = minimize(nll_dir, x0=guess, method='Nelder-Mead')
    # compute proper minima
    m_dir = minimize(nll_dir, x0=m_dir.x)

    return m_dir


def position_time_direction_staged(xyz, times, rindex, MCID,
                                   prompt_cut=None, seed=None, coordinate=None, pdfs=None,
                                   verbose=False, return_nll=False):
    '''Staged position+time -> direction fit

    Inspired by the SNO FTP algorithm
        - Uses 1D pdfs for hit time residuals and cos(theta)
        - All hits used for position+time component
        - Prompt hit time residuals used for direction component (pos,time fixed)
        - Works well for most materials with proper choice of prompt_cut

    xyz is a numpy array of [x,y,z] photon hit positions
    times is a numpy array of the times of the hits in the xyz array
    rindex is the characteristic refractive index for computing hit time residuals
    prompt_cut is either None for no cut, or a time in nanoseconds used to select
        prompt hit-time residuals for a direction fit
    seed is either None which uses a fixed seed or an (x,y,z,t,theta,phi) tuple
    coordinate changes the behavior of the function to return hit time residuals
        and cos(theta) values to be used in PDF generation. In this mode coordinate
        should be a tuple (x,y,z,t,theta,phi) defining the TRUE event information
    pdfs should be set when not coordinating to provide callable objects for to
        evaluate negative log likelihood values in the fit.
    return_nll will return the negative log likelihood function instead of
        minimizing it.

    result will either be (tresids,costhetas) if coordinating or (success,res)
        where res is a MinimizerResult containing the fit information for
        parameters (x,y,z,theta,phi)'''

    c = 299792458 * 10 ** 3 * 10 ** -9 / rindex  # mm/ns

    if coordinate is not None:  # coordinate variable should be the truth information
        x, y, z, t, theta, phi = coordinate
        pos = np.asarray([x, y, z])
        P = xyz - pos
        D = np.sqrt(np.sum(np.square(P), axis=1))
        T = times - t
        tresid = T - D / c
        dvec = np.asarray([np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)])
        costheta = np.sum(P * dvec, axis=1) / D
        if prompt_cut is None:
            return tresid, costheta
        else:
            mask = tresid < prompt_cut
            return tresid, costheta[mask]
    else:  # pdfs should contain callable objects to evaluate -log(pdf) for each element in numpy array
        if prompt_cut is None:
            tresid_pdf, costheta_pdf = pdfs
        else:
            tresid_pdf, costheta_pdf, angular_norm = pdfs

        m = fit_position_time(xyz, times, c, tresid_pdf, seed=seed)

        ID = MCID[0]

        x, y, z, t = m.x
        pos = np.asarray([x, y, z])
        P = xyz - pos
        D = np.sqrt(np.sum(np.square(P), axis=1))
        T = times - t
        tresid = T - D / c

        if prompt_cut is None:
            prompt_hits = xyz
        else:
            prompt_hits = xyz[tresid < prompt_cut]

        m_dir = fit_direction(prompt_hits, pos, costheta_pdf, seed=seed)

        return (m.success and m_dir.success), np.append(np.concatenate([m.x, m_dir.x]), ID)


invrt2pi = 1 / np.sqrt(np.pi * 2)


def gauss(x, mean, sigma, norm):
    '''gaussian distribution in a convenient for for scipy.optimize.curve_fit'''
    return norm * invrt2pi / sigma * np.exp(-np.square((x - mean) / sigma) / 2)


def CollectEvtsHitInfo(vEvts):
    xyz = []
    times = []
    source = []
    MCID = []
    for evt in vEvts:
        xyz.append(vEvts[evt][:, :3])
        times.append(vEvts[evt][:, 3])
        source.append(vEvts[evt][:, 4])
        MCID.append(vEvts[evt][:, 5])
    return xyz, times, source, MCID


def GetRIndex(xyz, times, sources):
    all_pos = np.concatenate(xyz)[sources == 0]  # Cerenkov Only
    all_time = np.concatenate(times)[sources == 0]  # Cerenkov Only
    all_dist = np.sqrt(np.sum(np.square(all_pos - x0_true[:3]), axis=1))

    ''' Extract rindex of material
    Creates histogram of distance from true origin pos to PMTs, divided by time of hits
    Gaussian fit and mean value gives you average speed of light in specific target material
    '''

    speed_o_light = all_dist / all_time
    fastest_speed = 200
    counts, edges = np.histogram(speed_o_light, bins=np.linspace(fastest_speed - 80, fastest_speed + 80, 40))
    centers = (edges[1:] + edges[:-1]) / 2

    errs = np.sqrt(counts)
    errs[errs == 0] = 1
    p, pcov = curve_fit(gauss,
                        centers,
                        counts,
                        p0=(centers[np.argmax(counts)], 1, np.sum(counts)),
                        sigma=errs)
    c = p[0]

    rindex = 299792458 * 10 ** 3 * 10 ** -9 / c
    return rindex


def reconstruct_from_ratpac(tag, evts,
                            x0_true=None,
                            cores=8,
                            max_fits=None,
                            prompt_cut=None,
                            seed=None,
                            pdfs=None,
                            rindex=None,
                            t_max=800):
    if x0_true is None:
        x0_true = [0, 0, 0, 0, 0, 0]
    X_TRUE, Y_TRUE, Z_TRUE, T_TRUE, THETA_TRUE, PHI_TRUE = x0_true

    xyz, times, source, MCID = CollectEvtsHitInfo(evts)
    sources = np.concatenate(source)
    if rindex is None:
        rindex = GetRIndex(xyz, times, sources)
    else:
        rindex = rindex

    minTRes = -5
    maxTRes = t_max
    nbBinsTRes = int(2 * (t_max + 5))
    minCT = -1
    maxCT = 1
    nbBinsCT = 18

    if pdfs is None:
        ''' Create PDFs '''

        # first coordinate
        try:
            pool = multiprocessing.Pool(cores)
            coordinate = partial(position_time_direction_staged, coordinate=x0_true, prompt_cut=prompt_cut)
            pdf_vals = pool.starmap(coordinate, zip(xyz, times, itertools.repeat(rindex)))
            tresids = np.concatenate([t for t, _ in pdf_vals])
            costhetas = np.concatenate([cth for _, cth in pdf_vals])
            angular_norm = len(costhetas) / len(tresids)
        except:
            raise
        finally:
            pool.close()

    else:
        tresids, costhetas = pdfs

    counts, edges = np.histogram(tresids, bins=np.linspace(minTRes, maxTRes, nbBinsTRes))
    tresid_pdf = NLLPDF(counts, edges, pull=1)
    counts, edges = np.histogram(costhetas, bins=np.linspace(minCT, maxCT, nbBinsCT))
    costhetas_pdf = CosThetaPDF(counts, edges)

    if prompt_cut is not None:
        pdfs = tresid_pdf, costhetas_pdf, angular_norm
    else:
        pdfs = tresid_pdf, costhetas_pdf

    ''' Now FIT '''

    # then fit
    try:
        pool = multiprocessing.Pool(cores)
        if seed is not None:
            fitter = partial(position_time_direction_staged, pdfs=pdfs, seed=seed, prompt_cut=prompt_cut)
        else:
            fitter = partial(position_time_direction_staged, pdfs=pdfs, prompt_cut=prompt_cut)
            fits = pool.starmap(fitter,
                                zip(xyz if max_fits is None else xyz[:max_fits], times, itertools.repeat(rindex), MCID))
            success = np.asarray([z for z, _ in fits])
            fits = np.asarray([f for _, f in fits])
    except:
        raise
    finally:
        pool.close()

    cache_file = tag + '_CACHE' + '.h5'
    with h5py.File(cache_file, 'w') as hf:
        hf['rindex'] = rindex
        hf['fits'] = fits
        hf['tresids'] = tresids
        hf['sources'] = sources

    return tag, rindex, fits, tresids, sources


import uproot


def Export2TTree(sOutputName, fits):
    with uproot.recreate(sOutputName) as f:
        f["Recon"] = uproot.newtree({"X": uproot.newbranch(float, title="X"),
                                     "Y": uproot.newbranch(float, title="Y"),
                                     "Z": uproot.newbranch(float, title="Z"),
                                     "T": uproot.newbranch(float, title="T"),
                                     "Theta": uproot.newbranch(float, title="Theta"),
                                     "Phi": uproot.newbranch(float, title="Phi"),
                                     "MCID": uproot.newbranch(int, title="MCID")})
        f["Recon"].extend({"X": fits[:, 0],
                           "Y": fits[:, 1],
                           "Z": fits[:, 2],
                           "T": fits[:, 3],
                           "Theta": fits[:, 4],
                           "Phi": fits[:, 5],
                           "MCID": fits[:, 6]})


''' Create TRes hist for an evt, 
    from a user input speed-of-light (SoL) and true origin TrueOrigin and time TrueTime. 
    SoL should be in mm/ns, TrueOrigin a list [x, y, z] in mm and TrueTime in ns (or any consistent S.U)'''


def CreateEvtTRes(evt, SoL, TrueOrigin, TrueTime):
    TrueOrigin = np.asarray(TrueOrigin)
    x, y, z, t = evt[:, :4].T
    pos = np.asarray([x, y, z])
    P = pos - TrueOrigin[:, np.newaxis]
    D = np.sqrt(np.sum(np.square(P.T), axis=1))
    T = t - np.full((len(t),), TrueTime)
    tresid = T - D / SoL
    return tresid


''' Create costheta hist for an evt, 
    from a user input true origin TrueOrigin, and dir TrueDir. 
    TrueOrigin a list [x, y, z], TrueDir a list [x, y, z] NORMED'''


def CreateEvtCTheta(evt, TrueOrigin, TrueDir):
    TrueOrigin = np.asarray(TrueOrigin)
    TrueDir = np.asarray(TrueDir)
    x, y, z = evt[:, :3].T
    pos = np.asarray([x, y, z])
    P = pos - TrueOrigin[:, np.newaxis]
    D = np.sqrt(np.sum(np.square(P.T), axis=1))
    costheta = np.dot(P.T, TrueDir) / D
    return costheta


''' Create TRes and CTheta hist for an evt '''


def CreateEvtTResAndCTheta(evt, SoL, TrueOrigin, TrueTime, TrueDir):
    TrueOrigin = np.asarray(TrueOrigin)
    TrueDir = np.asarray(TrueDir)
    x, y, z, t = evt[:, :4].T
    pos = np.asarray([x, y, z])
    P = pos - TrueOrigin[:, np.newaxis]
    D = np.sqrt(np.sum(np.square(P.T), axis=1))
    T = t - np.full((len(t),), TrueTime)
    tresid = T - D / SoL
    costheta = np.dot(P.T, TrueDir) / D
    return tresid, costheta


def GetAverageRindex(evts, x0_true):
    xyz, times, source, MCID = CollectEvtsHitInfo(evts)
    sources = np.concatenate(source)
    all_pos = np.concatenate(xyz)[sources == 0]  # Cerenkov Only
    all_time = np.concatenate(times)[sources == 0]  # Cerenkov Only
    all_dist = np.sqrt(np.sum(np.square(all_pos - x0_true[:3]), axis=1))

    ''' Extract rindex of material
    Creates histogram of distance from true origin pos to PMTs, divided by time of hits
    Gaussian fit and mean value gives you average speed of light in specific target material
    '''

    speed_o_light = all_dist / all_time
    fastest_speed = 200
    counts, edges = np.histogram(speed_o_light, bins=np.linspace(fastest_speed - 80, fastest_speed + 80, 40))
    centers = (edges[1:] + edges[:-1]) / 2

    errs = np.sqrt(counts)
    errs[errs == 0] = 1
    p, pcov = curve_fit(gauss,
                        centers,
                        counts,
                        p0=(centers[np.argmax(counts)], 1, np.sum(counts)),
                        sigma=errs)
    c = p[0]

    rindex = 299792458 * 10 ** 3 * 10 ** -9 / c
    return rindex


def PlotPDFs(tresids, costhetas, minTRes=-5, maxTRes=200, nbBinsTRes=205, minCT=-1, maxCT=1, nbBinsCT=18):
    plt.figure()
    plt.hist2d(tresids,
               costhetas,
               bins=(np.linspace(minTRes, maxTRes, nbBinsTRes), np.linspace(minCT, maxCT, nbBinsCT)),
               norm=LogNorm())
    plt.colorbar()
    #     plt.savefig(tag+"_tresidsVScostheta.png")
    plt.figure()
    plt.hist(tresids,
             bins=np.linspace(minTRes,
                              maxTRes,
                              nbBinsTRes))
    #     plt.savefig(tag+"_tresids.png")
    plt.figure()
    plt.hist(costhetas, np.linspace(minCT, maxCT, nbBinsCT))
#     plt.savefig(tag+"_costheta.png")
