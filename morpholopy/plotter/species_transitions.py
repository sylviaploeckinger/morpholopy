import numpy as np
from matplotlib import gridspec

import matplotlib.pyplot as plt
from plotter.helpers import (
    get_Schruba_data,
    get_Schruba_upperlimits,
    get_Bigiel2008_data,
)
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
import scipy.ndimage
from mpl_toolkits.axes_grid1 import make_axes_locatable

from scipy.stats import gaussian_kde

cmRdBu = plt.cm.get_cmap("RdYlBu_r")

Sigma_min = 0.0
Sigma_max = 99.0
H2_over_HI_min = 10.0 ** (-1.9)
H2_over_HI_max = 10.0 ** (+2.4)

twelve_plus_logOH_solar = 8.69
twelve_plus_logOH_min = twelve_plus_logOH_solar - 0.5
twelve_plus_logOH_max = twelve_plus_logOH_solar + 0.5

logSFRmin = -12.2
logSFRmax = -5.5
logmolmin = -2.2
logmolmax = 3.2
logstarmin = -1.5
logstarmax = 4.25
nbins = 100
logstarbins = np.linspace(logstarmin, logstarmax, nbins)
logmolbins = np.linspace(logmolmin, logmolmax, nbins)
nbins_hist = 50


def species_transitions_combined(sim_name, output_path, halo_min_stellar_mass):

    # Plots using azimuthally averaged quantities
    filename_text_f1 = output_path + "/Azimuth_averaged_" + sim_name + ".dat"
    (
        redshift,
        binsize,
        SigmaHIH2,
        SigmaH2_over_HI,
        twelveplusOH,
        Sigma_star,
    ) = np.genfromtxt(
        filename_text_f1, comments="#", unpack=True, usecols=(2, 4, 5, 6, 7, 8)
    )

    text = sim_name + ", z = %.1f" % (redshift[0]) + "\n"
    text += "Limited to galaxies with:\n"
    text += r"log$_{\mathrm{10}}$ M$_{\star}$ [M$_{\odot}$] > %.2f" % (
        np.log10(halo_min_stellar_mass)
    )

    #######################
    # Sigma HI / Sigma H2  (data: Schruba+2011, color: metallicity)
    #######################

    fig = plt.figure(figsize=(3.5, 3.0))
    fig.subplots_adjust(left=0.1, right=0.8, top=0.73, bottom=0.17)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.0, 0.05], wspace=0.0)

    ax = plt.subplot(gs[0])
    ax.text(
        -0.1, 1.01, text, ha="left", va="bottom", transform=ax.transAxes, fontsize=8
    )
    ax.set_aspect(
        (Sigma_max - Sigma_min) / (np.log10(H2_over_HI_max) - np.log10(H2_over_HI_min))
    )
    ax.set_xlim(Sigma_min, Sigma_max)
    ax.set_ylim(H2_over_HI_min, H2_over_HI_max)
    ax.set_yscale("log")
    ax.set_yticks([0.1, 1.0, 10.0, 100.0])
    ax.set_yticklabels(["0.1", "1", "10", "100"])
    ax.set_xlabel(
        r"$\Sigma_{\mathrm{HI}}$ + $\Sigma_{\mathrm{H2}}$ [M$_{\odot}$ pc$^{-2}$]"
    )
    ax.set_ylabel(r"$\Sigma_{\mathrm{H2}}$ / $\Sigma_{\mathrm{HI}}$")

    # Schruba+2011 data
    S_HI_obs, S_H2_obs, flag_obs = get_Schruba_upperlimits()
    ax.scatter(
        S_HI_obs[flag_obs == 0] + S_H2_obs[flag_obs == 0],
        S_H2_obs[flag_obs == 0] / S_HI_obs[flag_obs == 0],
        color="darkgrey",
        marker="v",
        facecolor="none",
        label="Schruba+2011 (u.limit)",
    )
    S_HI_obs, S_H2_obs, flag_obs = get_Schruba_data()
    ax.scatter(
        S_HI_obs[flag_obs == 2] + S_H2_obs[flag_obs == 2],
        S_H2_obs[flag_obs == 2] / S_HI_obs[flag_obs == 2],
        color="lightgrey",
        marker="+",
        label="Schruba+2011 (detect.)",
    )

    radialbinsizes = [0.22, 0.8, 1.8]
    markers = ["o", "s", "D"]

    for irad in range(len(radialbinsizes)):
        radialbin_kpc = radialbinsizes[irad]
        x = SigmaHIH2[binsize == radialbin_kpc]
        y = SigmaH2_over_HI[binsize == radialbin_kpc]
        c = twelveplusOH[binsize == radialbin_kpc]

        sc = ax.scatter(
            x,
            y,
            c=c,
            vmin=twelve_plus_logOH_min,
            vmax=twelve_plus_logOH_max,
            cmap=cmRdBu,
            edgecolors="black",
            label="Azim.avg. %.2f kpc" % (radialbin_kpc),
            marker=markers[irad],
        )

    ax.legend(bbox_to_anchor=(0.85, 1.03), loc="lower left", ncol=1, fontsize=6)

    cax = plt.subplot(gs[1])
    cb = fig.colorbar(sc, cax=cax, orientation="vertical")
    cb.set_ticks(
        [twelve_plus_logOH_min, twelve_plus_logOH_solar, twelve_plus_logOH_max]
    )
    cb.set_label(r"12 + log$_{\mathrm{10}}$(O/H)$_{\mathrm{diffuse}}$")

    fig.savefig(
        output_path + "/Species_transition_" + sim_name + "_Schruba2011.png", dpi=150
    )
    plt.close()

    #######################
    # Sigma HI / Sigma H2  (data: Schruba+2011, color: stellar surface density)
    #######################
    fig = plt.figure(figsize=(3.5, 3.0))
    fig.subplots_adjust(left=0.1, right=0.8, top=0.73, bottom=0.17)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.0, 0.05], wspace=0.0)

    ax = plt.subplot(gs[0])
    ax.text(
        -0.1, 1.01, text, ha="left", va="bottom", transform=ax.transAxes, fontsize=8
    )
    ax.set_aspect(
        (Sigma_max - Sigma_min) / (np.log10(H2_over_HI_max) - np.log10(H2_over_HI_min))
    )
    ax.set_xlim(Sigma_min, Sigma_max)
    ax.set_ylim(H2_over_HI_min, H2_over_HI_max)
    ax.set_yscale("log")
    ax.set_yticks([0.1, 1.0, 10.0, 100.0])
    ax.set_yticklabels(["0.1", "1", "10", "100"])
    ax.set_xlabel(
        r"$\Sigma_{\mathrm{HI}}$ + $\Sigma_{\mathrm{H2}}$ [M$_{\odot}$ pc$^{-2}$]"
    )
    ax.set_ylabel(r"$\Sigma_{\mathrm{H2}}$ / $\Sigma_{\mathrm{HI}}$")

    # Schruba+2011 data
    S_HI_obs, S_H2_obs, flag_obs = get_Schruba_upperlimits()
    ax.scatter(
        S_HI_obs[flag_obs == 0] + S_H2_obs[flag_obs == 0],
        S_H2_obs[flag_obs == 0] / S_HI_obs[flag_obs == 0],
        color="darkgrey",
        marker="v",
        facecolor="none",
        label="Schruba+2011 (u.limit)",
    )
    S_HI_obs, S_H2_obs, flag_obs = get_Schruba_data()
    ax.scatter(
        S_HI_obs[flag_obs == 2] + S_H2_obs[flag_obs == 2],
        S_H2_obs[flag_obs == 2] / S_HI_obs[flag_obs == 2],
        color="lightgrey",
        marker="+",
        label="Schruba+2011 (detect.)",
    )

    radialbinsizes = [0.22, 0.8, 1.8]
    markers = ["o", "s", "D"]

    for irad in range(len(radialbinsizes)):
        radialbin_kpc = radialbinsizes[irad]
        x = SigmaHIH2[binsize == radialbin_kpc]
        y = SigmaH2_over_HI[binsize == radialbin_kpc]
        c = np.log10(Sigma_star[binsize == radialbin_kpc])

        sc = ax.scatter(
            x,
            y,
            c=c,
            vmin=logstarmin,
            vmax=logstarmax,
            cmap=cmRdBu,
            edgecolors="black",
            label="Azim.avg. %.2f kpc" % (radialbin_kpc),
            marker=markers[irad],
        )

    ax.legend(bbox_to_anchor=(0.85, 1.03), loc="lower left", ncol=1, fontsize=6)

    cax = plt.subplot(gs[1])
    cb = fig.colorbar(sc, cax=cax, orientation="vertical")
    cb.set_label(r"log$_{\mathrm{10}} \Sigma_{\star}$ [M$_{\odot}$ pc$^{-2}$]")

    fig.savefig(
        output_path + "/Species_transition_" + sim_name + "_Schruba2011_Sigma_star.png",
        dpi=150,
    )
    plt.close()

    # Plots using grid averaged quantities
    filename_text_f2 = output_path + "/Grid_averaged_" + sim_name + ".dat"
    (
        redshift_fp,
        binsize_fp,
        log10sfr_fp,
        log10H2_fp,
        log10HI_fp,
        log10star_fp,
    ) = np.genfromtxt(
        filename_text_f2, comments="#", unpack=True, usecols=(2, 4, 5, 6, 7, 8)
    )

    #######################
    # Sigma HI / Sigma H2  (data: Bigiel+2008, scatter)
    #######################

    fig = plt.figure(figsize=(3.5, 3.0))
    fig.subplots_adjust(left=0.1, right=0.8, top=0.73, bottom=0.17)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.0, 0.05], wspace=0.0)

    ax = plt.subplot(gs[0])
    text += "\n" + "grid =%.2f kpc" % (binsize_fp[0])
    ax.text(
        -0.1, 1.01, text, ha="left", va="bottom", transform=ax.transAxes, fontsize=8
    )

    ax.set_aspect(
        (Sigma_max - Sigma_min) / (np.log10(H2_over_HI_max) - np.log10(H2_over_HI_min))
    )
    ax.set_xlim(Sigma_min, Sigma_max)
    ax.set_ylim(np.log10(H2_over_HI_min), np.log10(H2_over_HI_max))
    ax.set_yticks([-1.0, 0.0, 1.0, 2.0])
    ax.set_yticklabels(["0.1", "1", "10", "100"])
    ax.set_xlabel(
        r"$\Sigma_{\mathrm{HI}}$ + $\Sigma_{\mathrm{H2}}$ [M$_{\odot}$ pc$^{-2}$]"
    )
    ax.set_ylabel(r"$\Sigma_{\mathrm{H2}}$ / $\Sigma_{\mathrm{HI}}$")

    S_HI_obs, S_H2_obs = get_Bigiel2008_data()
    ax.scatter(
        (S_HI_obs + S_H2_obs),
        np.log10(S_H2_obs / S_HI_obs),
        color="lightgrey",
        marker="+",
        label="Bigiel+2008",
    )

    x = np.power(10.0, log10H2_fp) + np.power(10.0, log10HI_fp)
    y = log10H2_fp - log10HI_fp

    sc = ax.scatter(
        np.power(10.0, log10H2_fp) + np.power(10.0, log10HI_fp),
        log10H2_fp - log10HI_fp,
        c="black",
        marker="o",
        s=10,
        zorder=3,
    )

    ax.legend(loc="lower right", fontsize=6)
    fig.savefig(
        output_path + "/Species_transition_" + sim_name + "_Bigiel2008_scatter.png",
        dpi=150,
    )
    plt.close()

    #######################
    # Sigma HI / Sigma H2  (data: Bigiel+2008, contour)
    #######################

    fig = plt.figure(figsize=(3.5, 3.0))
    fig.subplots_adjust(left=0.1, right=0.8, top=0.73, bottom=0.17)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.0, 0.05], wspace=0.0)

    ax = plt.subplot(gs[0])
    text += "\n" + "grid =%.2f kpc" % (binsize_fp[0])
    ax.text(
        -0.1, 1.01, text, ha="left", va="bottom", transform=ax.transAxes, fontsize=8
    )

    ax.set_aspect(
        (Sigma_max - Sigma_min) / (np.log10(H2_over_HI_max) - np.log10(H2_over_HI_min))
    )
    ax.set_xlim(Sigma_min, Sigma_max)
    ax.set_ylim(np.log10(H2_over_HI_min), np.log10(H2_over_HI_max))
    ax.set_yticks([-1.0, 0.0, 1.0, 2.0])
    ax.set_yticklabels(["0.1", "1", "10", "100"])
    ax.set_xlabel(
        r"$\Sigma_{\mathrm{HI}}$ + $\Sigma_{\mathrm{H2}}$ [M$_{\odot}$ pc$^{-2}$]"
    )
    ax.set_ylabel(r"$\Sigma_{\mathrm{H2}}$ / $\Sigma_{\mathrm{HI}}$")

    S_HI_obs, S_H2_obs = get_Bigiel2008_data()
    ax.scatter(
        (S_HI_obs + S_H2_obs),
        np.log10(S_H2_obs / S_HI_obs),
        color="lightgrey",
        marker="+",
        label="Bigiel+2008",
    )

    x = np.power(10.0, log10H2_fp) + np.power(10.0, log10HI_fp)
    y = log10H2_fp - log10HI_fp

    data = np.vstack([x, y])
    kde = gaussian_kde(data, bw_method=0.1)
    xgrid = np.linspace(Sigma_min, Sigma_max, nbins_hist)
    ygrid = np.linspace(np.log10(H2_over_HI_min), np.log10(H2_over_HI_max), nbins_hist)

    Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)
    Z = kde.evaluate(np.vstack([Xgrid.ravel(), Ygrid.ravel()]))

    levels = (
        np.nanmax(Z) * 0.01,
        np.nanmax(Z) * 0.05,
        np.nanmax(Z) * 0.1,
        np.nanmax(Z) * 0.25,
        np.nanmax(Z) * 0.5,
        np.nanmax(Z) * 0.75,
        np.nanmax(Z) * 0.9,
        np.nanmax(Z),
    )
    ax.contour(
        xgrid,
        ygrid,
        Z.reshape(Xgrid.shape),
        origin="lower",
        extent=[xgrid[0], xgrid[-1], ygrid[0], ygrid[-1]],
        zorder=5,
        colors="black",
        levels=levels,
        linestyles="solid",
    )

    ax.legend(loc="lower right", fontsize=6)
    fig.savefig(
        output_path + "/Species_transition_" + sim_name + "_Bigiel2008_contour.png",
        dpi=150,
    )
    plt.close()

    #######################
    # Sigma_SFR vs. Sigma_star
    #######################
    fig = plt.figure(figsize=(3.0, 3.0))
    fig.subplots_adjust(left=0.10, right=0.95, top=0.75, bottom=0.17)
    gs = gridspec.GridSpec(1, 1)

    ax = plt.subplot(gs[0])
    ax.text(0.0, 1.01, text, ha="left", va="bottom", transform=ax.transAxes, fontsize=8)
    ax.set_xlim(logstarmin, logstarmax)
    ymax = -7.0
    ymin = -13.0
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(r"log $\Sigma_{\star}$ [M$_{\odot}$ pc$^{-2}$]")
    ax.set_ylabel(r"log $\Sigma_{\mathrm{sfr}}$ / $\Sigma_{\star}$ [yr$^{-1}$]")
    ax.set_xticks([0.0, 2.0, 4.0])
    ax.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax.set_yticks([-12.0, -10.0, -8.0])
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax.set_aspect((logstarmax - logstarmin) / (ymax - ymin))

    x = log10star_fp
    y = log10sfr_fp - log10star_fp

    hist, xe, ye = np.histogram2d(
        x, y, bins=nbins_hist, range=[[logstarmin, logstarmax], [ymin, ymax]]
    )
    hist_resampled = scipy.ndimage.zoom(hist, 3)
    hist = np.log10(hist.T)
    hist_resampled = np.log10(hist_resampled.T)
    xe_resampled = np.linspace(xe[0], xe[-1], hist_resampled.shape[0], endpoint=True)
    ye_resampled = np.linspace(ye[0], ye[-1], hist_resampled.shape[1], endpoint=True)
    # plot individual points separately
    hist[hist <= np.log10(2.0)] = np.nan

    ax.scatter(x, y, facecolors="none", edgecolors="grey", zorder=0, s=4)
    ax.imshow(
        hist,
        interpolation="none",
        origin="lower",
        extent=[xe[0], xe[-1], ye[0], ye[-1]],
        vmin=0.0,
        zorder=1,
        aspect=(logstarmax - logstarmin) / (ymax - ymin),
        cmap="Greys",
    )

    ax.plot(
        logstarbins,
        -10.10 + 0.02 * logstarbins,
        color="#CC0066",
        linestyle="dashed",
        zorder=10,
    )
    ax.plot(
        logstarbins[logstarbins > 1.2],
        -10.10 + 0.02 * logstarbins[logstarbins > 1.2],
        color="#CC0066",
        linestyle="solid",
        zorder=10,
        label="Sanchez+2021(EDGE fit)",
    )
    ax.legend(loc="upper right", fontsize=6)

    fig.savefig(output_path + "/Grid_averaged_" + sim_name + "_SFR_Mstar.png", dpi=150)
    plt.close()

    #######################
    # Sigma_SFR vs. Sigma_mol
    #######################
    fig = plt.figure(figsize=(3.0, 3.0))
    fig.subplots_adjust(left=0.10, right=0.95, top=0.75, bottom=0.17)
    gs = gridspec.GridSpec(1, 1)

    ax = plt.subplot(gs[0])
    ax.text(0.0, 1.01, text, ha="left", va="bottom", transform=ax.transAxes, fontsize=8)
    ax.set_xlim(logmolmin, logmolmax)
    ymax = -7.0
    ymin = -11.0
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(r"log $\Sigma_{\mathrm{mol}}$ [M$_{\odot}$ pc$^{-2}$]")
    ax.set_ylabel(r"log $\Sigma_{\mathrm{sfr}}$ / $\Sigma_{\mathrm{mol}}$ [yr$^{-1}$]")
    ax.set_xticks([-2.0, 0.0, 2.0])
    ax.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax.set_yticks([-10.0, -8.0])
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax.set_aspect((logmolmax - logmolmin) / (ymax - ymin))

    x = log10H2_fp
    y = log10sfr_fp - log10H2_fp

    hist, xe, ye = np.histogram2d(
        x, y, bins=nbins_hist, range=[[logmolmin, logmolmax], [ymin, ymax]]
    )
    hist_resampled = scipy.ndimage.zoom(hist, 3)
    hist = np.log10(hist.T)
    hist_resampled = np.log10(hist_resampled.T)
    xe_resampled = np.linspace(xe[0], xe[-1], hist_resampled.shape[0], endpoint=True)
    ye_resampled = np.linspace(ye[0], ye[-1], hist_resampled.shape[1], endpoint=True)
    # plot individual points separately
    hist[hist <= np.log10(2.0)] = np.nan

    ax.scatter(x, y, facecolors="none", edgecolors="grey", zorder=0, s=4)
    ax.imshow(
        hist,
        interpolation="none",
        origin="lower",
        extent=[xe[0], xe[-1], ye[0], ye[-1]],
        vmin=0.0,
        zorder=1,
        aspect=(logmolmax - logmolmin) / (ymax - ymin),
        cmap="Greys",
    )

    ax.plot(
        logmolbins,
        -9.01 + (0.98 - 1.0) * logmolbins,
        color="#CC0066",
        linestyle="dashed",
        zorder=10,
    )
    ax.plot(
        logmolbins[logmolbins > 0.3],
        -9.01 + (0.98 - 1.0) * logmolbins[logmolbins > 0.3],
        color="#CC0066",
        linestyle="solid",
        zorder=10,
        label="Sanchez+2021(EDGE fit)",
    )

    ax.legend(loc="upper right", fontsize=6)

    fig.savefig(output_path + "/Grid_averaged_" + sim_name + "_SFR_Mmol.png", dpi=150)
    plt.close()

    #######################
    # Sigma_mol vs. Sigma_star
    #######################
    fig = plt.figure(figsize=(3.0, 3.0))
    fig.subplots_adjust(left=0.10, right=0.95, top=0.75, bottom=0.17)
    gs = gridspec.GridSpec(1, 1)

    ax = plt.subplot(gs[0])
    ax.text(0.0, 1.01, text, ha="left", va="bottom", transform=ax.transAxes, fontsize=8)

    ax.set_xlim(logstarmin, logstarmax)
    ymin = -3.0
    ymax = +1.0
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(r"log $\Sigma_{\star}$ [M$_{\odot}$ pc$^{-2}$]")
    ax.set_ylabel(r"log $\Sigma_{\mathrm{mol}}$ / $\Sigma_{\star}$")
    ax.set_xticks([0.0, 2.0, 4.0])
    ax.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax.set_yticks([-2.0, 0.0])
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax.set_aspect((logstarmax - logstarmin) / (ymax - ymin))

    x = log10star_fp
    y = log10H2_fp - log10star_fp

    hist, xe, ye = np.histogram2d(
        x, y, bins=nbins_hist, range=[[logstarmin, logstarmax], [ymin, ymax]]
    )
    hist_resampled = scipy.ndimage.zoom(hist, 3)
    hist = np.log10(hist.T)
    hist_resampled = np.log10(hist_resampled.T)
    xe_resampled = np.linspace(xe[0], xe[-1], hist_resampled.shape[0], endpoint=True)
    ye_resampled = np.linspace(ye[0], ye[-1], hist_resampled.shape[1], endpoint=True)
    # plot individual points separately
    hist[hist <= np.log10(2.0)] = np.nan

    ax.scatter(x, y, facecolors="none", edgecolors="grey", zorder=0, s=4)
    ax.imshow(
        hist,
        interpolation="none",
        origin="lower",
        extent=[xe[0], xe[-1], ye[0], ye[-1]],
        vmin=0.0,
        zorder=1,
        aspect=(logstarmax - logstarmin) / (ymax - ymin),
        cmap="Greys",
    )

    ax.plot(
        logstarbins[logstarbins > 1.2],
        -0.91 + (0.93 - 1.0) * logstarbins[logstarbins > 1.2],
        color="#CC0066",
        linestyle="solid",
        label="Sanchez+2021(EDGE fit)",
        zorder=10,
    )
    ax.plot(
        logstarbins,
        -0.91 + (0.93 - 1.0) * logstarbins,
        color="#CC0066",
        linestyle="dashed",
        zorder=10,
    )

    ax.legend(loc="upper right", fontsize=6)
    fig.savefig(output_path + "/Grid_averaged_" + sim_name + "_Mmol_Mstar.png", dpi=150)
    plt.close()

    return
