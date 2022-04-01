import numpy as np
from matplotlib import gridspec
import matplotlib.pyplot as plt

cols_to_use = np.asarray(2)
cols_to_use = np.append(cols_to_use, np.arange(4,17,1))

def scale_height_length_plots(sim_name, output_path, halo_min_stellar_mass):
    filename_text_f3 = output_path + "/Scale_height_length_" + sim_name + ".dat"
    f3 = open(filename_text_f3, 'r')
    (
        redshift, 
        Mstar,
        MHI,
        MH2,
        H_kpc_g,
        H_kpc_r,
        H_kpc_i,
        H_kpc_HI,
        H_kpc_H2,
        L_kpc_g,
        L_kpc_r,
        L_kpc_i,
        L_kpc_HI,
        L_kpc_H2,
    ) = np.genfromtxt (
        filename_text_f3, comments="#", unpack=True, usecols=cols_to_use
    )
    f3.close()

    text = sim_name + ", z = %.1f" % (redshift[0]) + "\n"
    text += "Limited to galaxies with:\n"
    text += r"log$_{\mathrm{10}}$ M$_{\star}$ [M$_{\odot}$] > %.2f" % (
        np.log10(halo_min_stellar_mass)
    )

    x = np.log10(Mstar)
    xmin = 6.
    xmax = 11.

    #######################
    # Scale height stars vs. stellar mass
    #######################

    ymin = 0.
    ymax = 1.9

    fig = plt.figure(figsize=(3.5, 3.0))
    fig.subplots_adjust(left=0.1, right=0.8, top=0.73, bottom=0.17)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.0, 0.05], wspace=0.0)
    ax = plt.subplot(gs[0])
    ax.text(
        -0.1, 1.01, text, ha="left", va="bottom", transform=ax.transAxes, fontsize=8
    )
    ax.set_aspect((xmax - xmin)/(ymax - ymin))
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xticks(np.arange(xmin, xmax + 1., 1.))
    ax.set_yticks(np.arange(ymin, ymax , 0.2))
    ax.set_xlabel(r"log$_{\mathrm{10}}$ M$_{\star\mathrm{,30kpc}}$ [M$_{\odot}$]")
    ax.set_ylabel(r"Scale height [kpc]")
    ax.scatter(x, H_kpc_g, edgecolor = 'black', facecolor = 'none', marker="o", label = "g")
    ax.scatter(x, H_kpc_r, edgecolor = 'black', facecolor = 'none', marker="s", label = "r")
    ax.scatter(x, H_kpc_i, edgecolor = 'black', facecolor = 'none', marker="D", label = "i")
    ax.legend(title = 'GAMA filter',bbox_to_anchor=(1.0, 1.0), loc = "upper left", fontsize = 8)

    fig.savefig(
        output_path + "/Scale_height_stars_" + sim_name + ".png", dpi = 150
    )
    plt.close()

    #######################
    # Scale height gas vs. stellar mass
    #######################

    ymin = 0.
    ymax = 1.9

    fig = plt.figure(figsize=(3.5, 3.0))
    fig.subplots_adjust(left=0.1, right=0.8, top=0.73, bottom=0.17)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.0, 0.05], wspace=0.0)
    ax = plt.subplot(gs[0])
    ax.text(
        -0.1, 1.01, text, ha="left", va="bottom", transform=ax.transAxes, fontsize=8
    )
    ax.set_aspect((xmax - xmin)/(ymax - ymin))
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xticks(np.arange(xmin, xmax + 1., 1.))
    ax.set_yticks(np.arange(ymin, ymax , 0.2))
    ax.set_xlabel(r"log$_{\mathrm{10}}$ M$_{\star\mathrm{,30kpc}}$ [M$_{\odot}$]")
    ax.set_ylabel(r"Scale height [kpc]")
    ax.scatter(x, H_kpc_HI, edgecolor = 'black', facecolor = 'none', marker="o", label = "HI")
    ax.scatter(x, H_kpc_H2, edgecolor = 'black', facecolor = 'none', marker="s", label = "H2")
    ax.legend(bbox_to_anchor=(1.0, 1.0), loc = "upper left", fontsize = 8)

    fig.savefig(
        output_path + "/Scale_height_gas_" + sim_name + ".png", dpi = 150
    )
    plt.close()

    #######################
    # Scale height / scale_length stars vs. stellar mass
    #######################

    ymin = 0.
    ymax = 1.2

    fig = plt.figure(figsize=(3.5, 3.0))
    fig.subplots_adjust(left=0.1, right=0.8, top=0.73, bottom=0.17)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.0, 0.05], wspace=0.0)
    ax = plt.subplot(gs[0])
    ax.text(
        -0.1, 1.01, text, ha="left", va="bottom", transform=ax.transAxes, fontsize=8
    )
    ax.set_aspect((xmax - xmin)/(ymax - ymin))
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xticks(np.arange(xmin, xmax + 1., 1.))
    ax.set_yticks(np.arange(ymin, ymax , 0.2))
    ax.set_xlabel(r"log$_{\mathrm{10}}$ M$_{\star\mathrm{,30kpc}}$ [M$_{\odot}$]")
    ax.set_ylabel(r"Scale height / scale length")
    ax.scatter(
        x[H_kpc_g>0.], 
        H_kpc_g[H_kpc_g>0.] / L_kpc_g[H_kpc_g>0.], 
        edgecolor = 'black', 
        facecolor = 'none', 
        marker="o", 
        label = "g"
    )
    ax.scatter(
        x[H_kpc_r>0.], 
        H_kpc_r[H_kpc_r>0.] / L_kpc_r[H_kpc_r>0.], 
        edgecolor = 'black', 
        facecolor = 'none', 
        marker="s", 
        label = "r"
    )
    ax.scatter(
        x[H_kpc_i>0.],
        H_kpc_i[H_kpc_i>0.] / L_kpc_i[H_kpc_i>0.],
        edgecolor = 'black',
        facecolor = 'none',
        marker="D",
        label = "i"
    )
    ax.legend(title = 'GAMA filter',bbox_to_anchor=(1.0, 1.0), loc = "upper left", fontsize = 8)

    fig.savefig(
        output_path + "/Scale_height_length_stars_" + sim_name + ".png", dpi = 150
    )
    plt.close()

    #######################
    # Scale height /scale length gas vs. stellar mass
    #######################

    ymin = 0.
    ymax = 1.2

    fig = plt.figure(figsize=(3.5, 3.0))
    fig.subplots_adjust(left=0.1, right=0.8, top=0.73, bottom=0.17)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.0, 0.05], wspace=0.0)
    ax = plt.subplot(gs[0])
    ax.text(
        -0.1, 1.01, text, ha="left", va="bottom", transform=ax.transAxes, fontsize=8
    )
    ax.set_aspect((xmax - xmin)/(ymax - ymin))
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xticks(np.arange(xmin, xmax + 1., 1.))
    ax.set_yticks(np.arange(ymin, ymax , 0.2))
    ax.set_xlabel(r"log$_{\mathrm{10}}$ M$_{\star\mathrm{,30kpc}}$ [M$_{\odot}$]")
    ax.set_ylabel(r"Scale height / scale length")
    ax.scatter(
        x[H_kpc_HI>0.], 
        H_kpc_HI[H_kpc_HI>0.] / L_kpc_HI[H_kpc_HI>0.], 
        edgecolor = 'black', 
        facecolor = 'none', 
        marker="o", 
        label = "HI",
    )
    ax.scatter(
        x[H_kpc_H2>0.], 
        H_kpc_H2[H_kpc_H2>0.] / L_kpc_H2[H_kpc_H2>0.],
        edgecolor = 'black', 
        facecolor = 'none', 
        marker="s", 
        label = "H2",
    )
    ax.legend(bbox_to_anchor=(1.0, 1.0), loc = "upper left", fontsize = 8)

    fig.savefig(
        output_path + "/Scale_height_length_gas_" + sim_name + ".png", dpi = 150
    )
    plt.close()

    return


