import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import matplotlib.ticker as ticker
from matplotlib.ticker import FuncFormatter
from matplotlib import gridspec
from sphviewer.tools import QuickView
from matplotlib.colors import LogNorm

def plot_galaxy(pstars,kappa,ihalo):
    # partsDATA contains particles data and is structured as follow
    # [ (:3)Position[Mpc]: (0)X | (1)Y | (2)Z ]
    cmap_stars = plt.cm.magma
    
    r_img = 5
    xmin = -r_img
    ymin = -r_img
    xmax = r_img
    ymax = r_img
    
    # Plot parameters
    params = {
        "font.size": 11,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (7, 4),
        "figure.subplot.left": 0.12,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.2,
        "figure.subplot.top": 0.9,
        "figure.subplot.wspace": 0.4,
        "figure.subplot.hspace": 0.4,
        "lines.markersize": 0.5,
        "lines.linewidth": 0.2,
    }
    rcParams.update(params)
    fig = plt.figure()
    ax = plt.subplot(1,2,1)
    ax.set_title(r"$\kappa_{\mathrm{co}} = $%0.2f / face on" % (kappa))

    ###### plot one side ########################
    qv = QuickView(pstars, mass=np.ones(len(pstars)), plot=False,
                   r='infinity', p=0, t=0, extent=[xmin, xmax, ymin, ymax],
                   x=0, y=0, z=0)
    img = qv.get_image()
    ext = qv.get_extent()
    ax.tick_params(labelleft=True, labelbottom=True, length=0)
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    im = ax.imshow(img, cmap=cmap_stars, norm=LogNorm(), extent=ext)
    ax.autoscale(False)
    
    ###### plot another side ########################
    ax = plt.subplot(1,2,2)
    ax.set_title("Edge on")
    
    qv = QuickView(pstars, mass=np.ones(len(pstars)), plot=False,
                   r='infinity', p=90, t=0, extent=[xmin, xmax, ymin, ymax],
                   x=0, y=0, z=0)
    img = qv.get_image()
    ext = qv.get_extent()
    ax.tick_params(labelleft=True, labelbottom=True, length=0)
    plt.xlabel('x [kpc]')
    plt.ylabel('z [kpc]')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    im = ax.imshow(img, cmap=cmap_stars, norm=LogNorm(), extent=ext)
    ax.autoscale(False)
   
    outfile = "galaxy_%i.png" % (ihalo)
    fig.savefig(outfile, dpi=150)
    print("Saved: %s" % (outfile))
    plt.close("all")

def plot_momentum(stellar_mass,momentum,parttype):
    
    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (5, 4),
        "figure.subplot.left": 0.15,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.18,
        "figure.subplot.top": 0.9,
        "lines.markersize": 6,
        "lines.linewidth": 2.0,
    }
    rcParams.update(params)

    figure()
    ax = plt.subplot(1,1,1)
    if parttype==4:
        ax.set_title("Stellar component")
        color = 'tab:blue'
        ylabel = "$j_{\mathrm{stars}}$ [kpc km/s]"

    if parttype==0:
        ax.set_title("HI+H2 gas")
        color = 'tab:green'
        ylabel = "$j_{\mathrm{gas}}$ [kpc km/s]"

    plt.grid("True")

    plt.plot(stellar_mass, momentum, 'o', color=color)
    
    plt.xlabel("Stellar Mass [M$_{\odot}$]")
    plt.ylabel(ylabel)
    plt.savefig("./momentum_parttype_%i.png"%parttype, dpi=200)
    plt.close()

def plot_kappa(stellar_mass,kappa,parttype):
    
    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (5, 4),
        "figure.subplot.left": 0.15,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.18,
        "figure.subplot.top": 0.9,
        "lines.markersize": 6,
        "lines.linewidth": 2.0,
    }
    rcParams.update(params)

    figure()
    ax = plt.subplot(1,1,1)
    if parttype==4:
        ax.set_title("Stellar component")
        color = 'tab:blue'
    if parttype==0:
        ax.set_title("HI+H2 gas")
        color = 'tab:green'

    plt.grid("True")
    
    plt.plot(stellar_mass, kappa, 'o', color=color)
    
    plt.xlabel("Stellar Mass [M$_{\odot}$]")
    plt.ylabel(r"$\kappa_{\mathrm{co}}$")
    plt.savefig("./Kappa_co_parttype_%i.png"%parttype, dpi=200)
    plt.close()

def plot_axis_ratios(stellar_mass,axis_ratios,parttype):
    
    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (9, 3.5),
        "figure.subplot.left": 0.1,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.15,
        "figure.subplot.top": 0.9,
        "figure.subplot.wspace": 0.4,
        "figure.subplot.hspace": 0.4,
        "lines.markersize": 6,
        "lines.linewidth": 2.0,
    }
    rcParams.update(params)

    if parttype==4:
        title = "Stellar component"
        color = 'tab:blue'
    if parttype==0:
        title = "HI+H2 gas"
        color = 'tab:green'

    ########
    figure()
    ax = plt.subplot(1,3,1)
    ax.set_title(title)
    plt.grid("True")

    plt.plot(stellar_mass, axis_ratios[:,0], 'o', color=color)

    plt.xlabel("Stellar Mass [M$_{\odot}$]")
    plt.ylabel("c/a")

    ########
    ax = plt.subplot(1,3,2)
    plt.grid("True")

    plt.plot(stellar_mass, axis_ratios[:,1], 'o', color=color)
    plt.xlabel("Stellar Mass [M$_{\odot}$]")
    plt.ylabel("c/b")

    ########
    ax = plt.subplot(1,3,3)
    plt.grid("True")
    
    plt.plot(stellar_mass, axis_ratios[:,2], 'o', color=color)
    
    plt.xlabel("Stellar Mass [M$_{\odot}$]")
    plt.ylabel("b/a")

    plt.savefig("./Axis_ratios_parttype_%i.png"%parttype, dpi=200)
    plt.close()


def plot_morphology(galaxy_data):

    # plot kappa for stars and gas :
    plot_kappa(galaxy_data.stellar_mass,galaxy_data.kappa_co,4)
    plot_kappa(galaxy_data.stellar_mass,galaxy_data.gas_kappa_co,0)

    # plot specific angular momentum  for stars and gas :
    plot_momentum(galaxy_data.stellar_mass,galaxy_data.momentum,4)
    plot_momentum(galaxy_data.stellar_mass,galaxy_data.gas_momentum,0)

    # plot axis ratios
    axis_ratios = np.zeros((len(galaxy_data.axis_ca),3))
    axis_ratios[:,0] = galaxy_data.axis_ca
    axis_ratios[:,1] = galaxy_data.axis_cb
    axis_ratios[:,2] = galaxy_data.axis_ba
    plot_axis_ratios(galaxy_data.stellar_mass,axis_ratios,4)

    axis_ratios = np.zeros((len(galaxy_data.gas_axis_ca),3))
    axis_ratios[:,0] = galaxy_data.gas_axis_ca
    axis_ratios[:,1] = galaxy_data.gas_axis_cb
    axis_ratios[:,2] = galaxy_data.gas_axis_ba
    plot_axis_ratios(galaxy_data.stellar_mass,axis_ratios,0)

def plot_galaxy_sparts(partsDATA,kappa,ihalo):
    # partsDATA contains particles data and is structured as follow
    # [ (:3)Position[kpc]: (0)X | (1)Y | (2)Z ]
    
    # Plot ranges
    r_img = 5
    xmin = -r_img
    ymin = -r_img
    xmax = r_img
    ymax = r_img
    
    pstars = partsDATA[:,0:3].copy()
    
    # Plot parameters
    params = {
        "font.size": 11,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (7, 4),
        "figure.subplot.left": 0.1,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.2,
        "figure.subplot.top": 0.9,
        "figure.subplot.wspace": 0.3,
        "figure.subplot.hspace": 0.3,
        "lines.markersize": 0.5,
        "lines.linewidth": 0.2,
    }
    rcParams.update(params)

    fig = plt.figure()
    ax = plt.subplot(1,2,1)
    ax.set_title(r"Stellar component - $\kappa_{\mathrm{co}} = $%0.2f" % (kappa))
    ax.tick_params(labelleft=True, labelbottom=True, length=0)
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    plt.plot(pstars[:,0],pstars[:,1],'o')
    ax.autoscale(False)

    ax = plt.subplot(1,2,2)
    ax.set_title(r"Stellar component - $\kappa_{\mathrm{co}} = $%0.2f" % (kappa))
    ax.tick_params(labelleft=True, labelbottom=True, length=0)
    plt.xlabel('x [kpc]')
    plt.ylabel('z [kpc]')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    plt.plot(pstars[:,0],pstars[:,2],'o')
    ax.autoscale(False)
    outfile = "galaxy_stars_%i.png" % (ihalo)
    fig.savefig(outfile, dpi=150)
    print("Saved: %s" % (outfile))
    plt.close("all")

def plot_galaxy_gas_parts(partsDATA,kappa,ihalo):
    # partsDATA contains particles data and is structured as follow
    # [ (:3)Position[kpc]: (0)X | (1)Y | (2)Z ]
    
    # Plot ranges
    r_img = 5
    xmin = -r_img
    ymin = -r_img
    xmax = r_img
    ymax = r_img
    
    pstars = partsDATA[:,0:3].copy()
    
    # Plot parameters
    params = {
        "font.size": 11,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (7, 4),
        "figure.subplot.left": 0.1,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.2,
        "figure.subplot.top": 0.9,
        "figure.subplot.wspace": 0.3,
        "figure.subplot.hspace": 0.3,
        "lines.markersize": 1.0,
        "lines.linewidth": 0.2,
    }
    rcParams.update(params)

    fig = plt.figure()
    ax = plt.subplot(1,2,1)
    ax.set_title(r"HI+H2 gas - $\kappa_{\mathrm{co}} = $%0.2f" % (kappa))
    ax.tick_params(labelleft=True, labelbottom=True, length=0)
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    plt.plot(pstars[:,0],pstars[:,1],'o',color='tab:green')
    ax.autoscale(False)

    ax = plt.subplot(1,2,2)
    ax.set_title(r"HI+H2 gas - $\kappa_{\mathrm{co}} = $%0.2f" % (kappa))
    ax.tick_params(labelleft=True, labelbottom=True, length=0)
    plt.xlabel('x [kpc]')
    plt.ylabel('z [kpc]')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    plt.plot(pstars[:,0],pstars[:,2],'o',color='tab:green')
    ax.autoscale(False)
    outfile = "galaxy_gas_%i.png" % (ihalo)
    fig.savefig(outfile, dpi=150)
    print("Saved: %s" % (outfile))
    plt.close("all")

