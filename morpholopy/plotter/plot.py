from matplotlib.colors import LogNorm
from pylab import *
from sphviewer.tools import QuickView
from plotter.html import add_web_section

def get_normalized_image(image,vmin=None,vmax=None):
    if(vmin==None):
        vmin = np.min(image[image>0])
        image[image==0] = vmin
    if(vmax==None):
        vmax = np.max(image)

    image=np.clip(image,vmin,vmax)
    image=(image-vmin)/(vmax-vmin)
    return image


def plot_galaxy(parts_data, kappa, mass, ihalo, parttype, GalPlotsInWeb):
    # partsDATA contains particles data and is structured as follow
    # [ (:3)Position[Mpc]: (0)X | (1)Y | (2)Z ]
    if parttype == 4: cmap = plt.cm.magma
    if parttype == 0: cmap = plt.cm.viridis

    pstars = parts_data[:, :3]
    hsml_parts = parts_data[:, 7]
    mass = parts_data[:, 3]

    r_img = 8
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
    ax = plt.subplot(1, 2, 1)
    if parttype == 4: title = r"Stellar component, $\kappa_{\mathrm{co}} = $%0.2f" % (kappa)
    if parttype == 0: title = r"HI+H2 gas, $\kappa_{\mathrm{co}} = $%0.2f" % (kappa)
    ax.set_title(title)

    ###### plot one side ########################
    qv = QuickView(pstars, mass=mass, logscale=True, hsml=hsml_parts, plot=False,
                   r='infinity', p=0, t=0, extent=[xmin, xmax, ymin, ymax],
                   x=0, y=0, z=0)
    img = qv.get_image()
    ext = qv.get_extent()
    ax.tick_params(labelleft=True, labelbottom=True, length=0)
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    img = get_normalized_image(img)
    im = ax.imshow(img, cmap=cmap, extent=ext)
    ax.autoscale(False)

    ###### plot another side ########################
    ax = plt.subplot(1, 2, 2)
    if parttype == 4: title = r"Stellar component, $\log_{10} M_{*}/M_{\odot} = $%0.2f" % (mass)
    if parttype == 0: title = r"HI+H2 gas, $\log_{10} M_{gas}/M_{\odot} = $%0.2f" % (mass)
    ax.set_title(title)

    qv = QuickView(pstars, hsml=hsml_parts, logscale=True, plot=False,
                   r='infinity', p=90, t=0, extent=[xmin, xmax, ymin, ymax],
                   x=0, y=0, z=0)
    img = qv.get_image()
    ext = qv.get_extent()
    ax.tick_params(labelleft=True, labelbottom=True, length=0)
    plt.xlabel('x [kpc]')
    plt.ylabel('z [kpc]')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    img = get_normalized_image(img)
    im = ax.imshow(img, cmap=cmap, extent=ext)
    ax.autoscale(False)

    if parttype == 0: outfile = "galaxy_gas_%i.png" % (ihalo)
    if parttype == 4: outfile = "galaxy_stars_%i.png" % (ihalo)
    fig.savefig(outfile, dpi=150)
    print("Saved: %s" % (outfile))
    plt.close("all")

    if parttype == 0: title = "Galaxy %i / Gas component" % (ihalo)
    if parttype == 4: title = "Galaxy %i / Stellar component" % (ihalo)
    caption = "Face-on (left) and edge-on (right)."
    if parttype == 0: id = abs(hash("galaxy gas"))
    if parttype == 4: id = abs(hash("galaxy stars"))
    GalPlotsInWeb.load_plots(title, caption, outfile, id)

def plot_momentum(stellar_mass,momentum,parttype,MorphologyPlotsInWeb):
    
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
    plt.yscale('log')
    plt.xscale('log')
    plt.axis([1e6,1e10,1e-1,1e3])
    plt.savefig("momentum_parttype_%i.png"%parttype, dpi=200)
    plt.close()

    if parttype==4: title = "Specific angular momentum / Stars"
    if parttype==0: title = "Specific angular momentum / HI+H2 gas"

    caption = "Ratio between the total angular momentum of stars (or gas) within 30 kpc of "
    caption += "aperture divided by the total mass in stars (or gas)."
    filename = "momentum_parttype_%i.png"%parttype
    id = abs(hash("momentum %i"%parttype))

    MorphologyPlotsInWeb.load_plots(title, caption, filename, id)

def plot_kappa(stellar_mass,kappa,parttype,MorphologyPlotsInWeb):
    
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
    
    plt.xscale('log')
    plt.xlabel("Stellar Mass [M$_{\odot}$]")
    plt.ylabel(r"$\kappa_{\mathrm{co}}$")
    plt.savefig("Kappa_co_parttype_%i.png"%parttype, dpi=200)
    plt.close()

    if parttype==4: title = "Kappa corotation / Stars"
    if parttype==0: title = "Kappa corotation / HI+H2 gas"

    caption = "Kappa corotation is defined as the fraction of kinetic energy in a galaxy "
    caption += "that is in ordered rotation. Note that the rotating contribution is calculated "
    caption += "only for prograde rotation."
    filename = "Kappa_co_parttype_%i.png"%parttype
    id = abs(hash("kappa co %i"%parttype))
    MorphologyPlotsInWeb.load_plots(title, caption, filename, id)


def plot_axis_ratios(stellar_mass,axis_ratios,parttype,MorphologyPlotsInWeb):
    
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

    plt.xscale('log')
    plt.xlabel("Stellar Mass [M$_{\odot}$]")
    plt.ylabel("c/a")

    ########
    ax = plt.subplot(1,3,2)
    plt.grid("True")

    plt.plot(stellar_mass, axis_ratios[:,1], 'o', color=color)

    plt.xscale('log')
    plt.xlabel("Stellar Mass [M$_{\odot}$]")
    plt.ylabel("c/b")

    ########
    ax = plt.subplot(1,3,3)
    plt.grid("True")
    
    plt.plot(stellar_mass, axis_ratios[:,2], 'o', color=color)
    
    plt.xscale('log')
    plt.xlabel("Stellar Mass [M$_{\odot}$]")
    plt.ylabel("b/a")

    plt.savefig("Axis_ratios_parttype_%i.png"%parttype, dpi=200)
    plt.close()

    if parttype==4: title = "Axis ratios / Stars"
    if parttype==0: title = "Axis ratios / HI+H2 gas"
    caption = "Axial ratios of galaxies more massive than 1e6 Msun in stellar mass. "
    caption += "a, b and c (a >= b >= c) represent the lengths of the primary axes. "
    caption += "Ratios have been calculated following eqs. (1) and (2) from Trayford+2018."
    filename = "Axis_ratios_parttype_%i.png"%parttype
    id = abs(hash("galaxy axis %i"%parttype))
    MorphologyPlotsInWeb.load_plots(title, caption, filename, id)


def plot_morphology(galaxy_data,web,MorphologyPlotsInWeb):

    # plot kappa for stars and gas :
    plot_kappa(10**galaxy_data.stellar_mass,galaxy_data.kappa_co,4,MorphologyPlotsInWeb)
    plot_kappa(10**galaxy_data.stellar_mass,galaxy_data.gas_kappa_co,0,MorphologyPlotsInWeb)

    title = 'Kappa corotation'
    id = abs(hash("Kappa corotation"))
    plots = MorphologyPlotsInWeb.plots_details
    add_web_section(web, title, id, plots)

    MorphologyPlotsInWeb.reset_plots_list()

    # plot specific angular momentum  for stars and gas :
    plot_momentum(10**galaxy_data.stellar_mass,galaxy_data.momentum,4,MorphologyPlotsInWeb)
    plot_momentum(10**galaxy_data.stellar_mass,galaxy_data.gas_momentum,0,MorphologyPlotsInWeb)

    title = 'Specific angular momentum'
    id = abs(hash("angular momentum"))
    plots = MorphologyPlotsInWeb.plots_details
    add_web_section(web, title, id, plots)

    MorphologyPlotsInWeb.reset_plots_list()

    # plot axis ratios
    axis_ratios = np.zeros((len(galaxy_data.axis_ca),3))
    axis_ratios[:,0] = galaxy_data.axis_ca
    axis_ratios[:,1] = galaxy_data.axis_cb
    axis_ratios[:,2] = galaxy_data.axis_ba
    plot_axis_ratios(10**galaxy_data.stellar_mass,axis_ratios,4,MorphologyPlotsInWeb)

    axis_ratios = np.zeros((len(galaxy_data.gas_axis_ca),3))
    axis_ratios[:,0] = galaxy_data.gas_axis_ca
    axis_ratios[:,1] = galaxy_data.gas_axis_cb
    axis_ratios[:,2] = galaxy_data.gas_axis_ba
    plot_axis_ratios(10**galaxy_data.stellar_mass,axis_ratios,0,MorphologyPlotsInWeb)

    title = 'Axis ratios'
    id = abs(hash("axis ratios"))
    plots = MorphologyPlotsInWeb.plots_details
    add_web_section(web, title, id, plots)


def plot_galaxy_sparts(partsDATA,kappa,mass,ihalo,PlotsInWeb):
    # partsDATA contains particles data and is structured as follow
    # [ (:3)Position[kpc]: (0)X | (1)Y | (2)Z ]
    
    # Plot ranges
    r_img = 15
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
        "lines.markersize": 0.2,
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
    ax.set_title(r"Stellar component - $\log_{10}$ $M_{*}/M_{\odot} = $%0.2f" % (mass))
    ax.tick_params(labelleft=True, labelbottom=True, length=0)
    plt.xlabel('x [kpc]')
    plt.ylabel('z [kpc]')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    plt.plot(pstars[:,0],pstars[:,2],'o')
    ax.autoscale(False)
    outfile = "galaxy_sparts_%i.png" % (ihalo)
    fig.savefig(outfile, dpi=150)
    print("Saved: %s" % (outfile))
    plt.close("all")

    title = "Galaxy %i / Star particles" % (ihalo)
    caption = "Face-on (left) and edge-on (right)."
    id = abs(hash("galaxy stars parts"))
    PlotsInWeb.load_plots(title, caption, outfile, id)

def plot_galaxy_gas_parts(partsDATA,kappa,mass,ihalo,PlotsInWeb):
    # partsDATA contains particles data and is structured as follow
    # [ (:3)Position[kpc]: (0)X | (1)Y | (2)Z ]
    
    # Plot ranges
    r_img = 15
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
        "lines.markersize": 0.2,
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
    ax.set_title(r"HI+H2 gas - $\log_{10}$ $M_{gas}/M_{\odot} = $%0.2f" % (mass))
    ax.tick_params(labelleft=True, labelbottom=True, length=0)
    plt.xlabel('x [kpc]')
    plt.ylabel('z [kpc]')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    plt.plot(pstars[:,0],pstars[:,2],'o',color='tab:green')
    ax.autoscale(False)
    outfile = "galaxy_parts_%i.png" % (ihalo)
    fig.savefig(outfile, dpi=150)
    print("Saved: %s" % (outfile))
    plt.close("all")

    title = "Galaxy %i / Gas particles" % (ihalo)
    caption = "Face-on (left) and edge-on (right)."
    id = abs(hash("galaxy gas parts"))
    PlotsInWeb.load_plots(title, caption, outfile, id)