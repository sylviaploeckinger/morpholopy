from pylab import *
from sphviewer.tools import QuickView
from swiftsimio.visualisation.rotation import rotation_matrix_from_vector
from numpy import matmul

import unyt
from swiftsimio import load, mask
from swiftsimio.visualisation.projection import project_pixel_grid
from matplotlib.colors import LogNorm

def plot_galaxy_parts(partsDATA, parttype, ang_momentum, halo_data, index, PlotsInWeb, output_path):

    # Plot ranges
    r_img = 5 * halo_data.halfmass_radius_star[index]
    xmin = -r_img
    ymin = -r_img
    xmax = r_img
    ymax = r_img

    pos_parts = partsDATA[:, 0:3].copy()

    face_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum)
    edge_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum,axis="y")

    pos_face_on = matmul(face_on_rotation_matrix, pos_parts.T)
    pos_face_on = pos_face_on.T
    pos_edge_on = matmul(edge_on_rotation_matrix, pos_parts.T)
    pos_edge_on = pos_edge_on.T

    # Plot parameters
    params = {
        "font.size": 11,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (7, 4),
        "figure.subplot.left": 0.1,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.1,
        "figure.subplot.top": 0.8,
        "figure.subplot.wspace": 0.3,
        "figure.subplot.hspace": 0.3,
        "lines.markersize": 0.2,
        "lines.linewidth": 0.2,
    }
    rcParams.update(params)

    fig = plt.figure()
    ax = plt.subplot(1, 2, 1)
    if parttype ==4:
        title = "Stellar component"
        color = 'tab:blue'
    if parttype ==0:
        title = "Gas component"
        color = 'tab:green'
    ax.set_title(title)
    ax.tick_params(labelleft=True, labelbottom=True, length=0)
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    plt.plot(pos_face_on[:, 0], pos_face_on[:, 1], 'o', color=color)
    ax.autoscale(False)

    ax = plt.subplot(1, 2, 2)
    if parttype ==4:
        kappa = halo_data.kappa_co[index]
        mass = halo_data.stellar_mass[index]
        ac = halo_data.axis_ca[index]
        cb = halo_data.axis_cb[index]
        ba = halo_data.axis_ba[index]
        title = r" $\kappa_{\mathrm{co}} = $%0.2f" % (kappa)
        title += " - $\log_{10}$ $M_{*}/M_{\odot} = $%0.2f" % (mass)
        title += " \n c/a = %0.2f," % (ac)
        title += " c/b = %0.2f," % (cb)
        title += " b/a = %0.2f" % (ba)
        color = 'tab:blue'
    if parttype ==0:
        kappa = halo_data.gas_kappa_co[index]
        mass = halo_data.gas_mass[index]
        ac = halo_data.gas_axis_ca[index]
        cb = halo_data.gas_axis_cb[index]
        ba = halo_data.gas_axis_ba[index]
        title = r" $\kappa_{\mathrm{co}} = $%0.2f" % (kappa)
        title += " - $\log_{10}$ $M_{gas}/M_{\odot} = $%0.2f" % (mass)
        title += " \n c/a = %0.2f," % (ac)
        title += " c/b = %0.2f," % (cb)
        title += " b/a = %0.2f" % (ba)
        color = 'tab:green'
    ax.set_title(title)

    ax.tick_params(labelleft=True, labelbottom=True, length=0)
    plt.xlabel('x [kpc]')
    plt.ylabel('z [kpc]')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    plt.plot(pos_edge_on[:, 0], pos_edge_on[:, 1], 'o', color=color)
    ax.autoscale(False)
    if parttype ==4: outfile = f"{output_path}/galaxy_sparts_%i.png" % (index)
    if parttype ==0: outfile = f"{output_path}/galaxy_parts_%i.png" % (index)
    fig.savefig(outfile, dpi=150)
    plt.close("all")

    if parttype ==4:
        title = "Star particles"
        caption = "Projection of stars within 5 times %0.1f kpc" % (halo_data.halfmass_radius_star[index])
        caption += " (the galaxy's  stellar half mass radius). Face-on (left) and edge-on (right)."
        id = abs(hash("galaxy stars parts %i" % (index)))
        outfile = "galaxy_sparts_%i.png" % (index)

    if parttype ==0:
        title = "Gas particles"
        caption = "Projection of gas within 5 times %0.1f kpc" % (halo_data.halfmass_radius_star[index])
        caption += " (the galaxy's  stellar half mass radius). Face-on (left) and edge-on (right)."
        id = abs(hash("galaxy gas parts %i" % (index)))
        outfile = "galaxy_parts_%i.png" % (index)

    PlotsInWeb.load_plots(title, caption, outfile, id)


def get_normalized_image(image,vmin=None,vmax=None):
    if(vmin==None):
        vmin = np.min(image[image>0])
        image[image==0] = vmin
    if(vmax==None):
        vmax = np.max(image)

    image=np.clip(image,vmin,vmax)
    image=(image-vmin)/(vmax-vmin)
    return image

def plot_galaxy(parts_data, parttype, ang_momentum, halo_data, index, GalPlotsInWeb, output_path):
    # partsDATA contains particles data and is structured as follow
    # [ (:3)Position[Mpc]: (0)X | (1)Y | (2)Z ]
    if parttype == 4: cmap = plt.cm.magma
    if parttype == 0: cmap = plt.cm.viridis

    pos_parts = parts_data[:, 0:3].copy()

    face_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum)
    edge_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum,axis="y")

    pos_face_on = matmul(face_on_rotation_matrix, pos_parts.T)
    pos_face_on = pos_face_on.T
    pos_edge_on = matmul(edge_on_rotation_matrix, pos_parts.T)
    pos_edge_on = pos_edge_on.T

    hsml_parts = parts_data[:, 7]
    mass = parts_data[:, 3]

    r_img = 5 * halo_data.halfmass_radius_star[index]
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
        "figure.subplot.bottom": 0.1,
        "figure.subplot.top": 0.8,
        "figure.subplot.wspace": 0.4,
        "figure.subplot.hspace": 0.4,
        "lines.markersize": 0.5,
        "lines.linewidth": 0.2,
    }
    rcParams.update(params)
    fig = plt.figure()
    ax = plt.subplot(1, 2, 1)
    if parttype == 4: title = "Stellar component"
    if parttype == 0: title = "HI+H2 gas"
    ax.set_title(title)

    ###### plot one side ########################
    qv = QuickView(pos_face_on, mass=mass, hsml=hsml_parts, logscale=True, plot=False,
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
    ax.imshow(img, cmap=cmap, extent=ext)
    ax.autoscale(False)

    ###### plot another side ########################
    ax = plt.subplot(1, 2, 2)
    if parttype ==4:
        kappa = halo_data.kappa_co[index]
        mass_galaxy = halo_data.stellar_mass[index]
        ac = halo_data.axis_ca[index]
        cb = halo_data.axis_cb[index]
        ba = halo_data.axis_ba[index]
        title = r" $\kappa_{\mathrm{co}} = $%0.2f" % (kappa)
        title += " - $\log_{10}$ $M_{*}/M_{\odot} = $%0.2f" % (mass_galaxy)
        title += " \n c/a = %0.2f," % (ac)
        title += " c/b = %0.2f," % (cb)
        title += " b/a = %0.2f" % (ba)
    if parttype ==0:
        kappa = halo_data.gas_kappa_co[index]
        mass_galaxy = halo_data.gas_mass[index]
        ac = halo_data.gas_axis_ca[index]
        cb = halo_data.gas_axis_cb[index]
        ba = halo_data.gas_axis_ba[index]
        title = r" $\kappa_{\mathrm{co}} = $%0.2f" % (kappa)
        title += " - $\log_{10}$ $M_{gas}/M_{\odot} = $%0.2f" % (mass_galaxy)
        title += " \n c/a = %0.2f," % (ac)
        title += " c/b = %0.2f," % (cb)
        title += " b/a = %0.2f" % (ba)
    ax.set_title(title)

    qv = QuickView(pos_edge_on, mass=mass, hsml=hsml_parts, logscale=True, plot=False,
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
    ax.imshow(img, cmap=cmap, extent=ext)
    ax.autoscale(False)

    if parttype == 0: outfile = f"{output_path}/galaxy_gas_%i.png" % (index)
    if parttype == 4: outfile = f"{output_path}/galaxy_stars_%i.png" % (index)
    fig.savefig(outfile, dpi=150)
    plt.close("all")

    if parttype == 0:
        title = "Gas component"
        caption = "Projection of gas within 5 times %0.1f kpc" % (halo_data.halfmass_radius_star[index])
        caption += " (the galaxy's  stellar half mass radius). Face-on (left) and edge-on (right)."
        id = abs(hash("galaxy gas %i" % (index)))
        outfile = "galaxy_gas_%i.png" % (index)

    if parttype == 4:
        title = "Stellar component"
        caption = "Projection of stars within 5 times %0.1f kpc" % (halo_data.halfmass_radius_star[index])
        caption += " (the galaxy's  stellar half mass radius). Face-on (left) and edge-on (right)."
        id = abs(hash("stars galaxy %i" % (index)))
        outfile = "galaxy_stars_%i.png" % (index)

    GalPlotsInWeb.load_plots(title, caption, outfile, id)

def plot_galaxy_swiftsimio(stars_ang_momentum, gas_ang_momentum, halo_data, index, GalPlotsInWeb, output_path, siminfo):

    # Radius around which to load data, we will visualise half of this
    size = 0.06 * unyt.Mpc

    x = halo_data.xminpot[index] * unyt.Mpc
    y = halo_data.yminpot[index] * unyt.Mpc
    z = halo_data.zminpot[index] * unyt.Mpc

    region = [
        [x - size, x + size],
        [y - size, y + size],
        [z - size, z + size],
    ]

    visualise_region = [
        x - 0.5 * size, x + 0.5 * size,
        y - 0.5 * size, y + 0.5 * size,
    ]

    extent = [-0.5 * size.value, 0.5 * size.value, -0.5 * size.value, 0.5 * size.value]

    size_pixel = 1e-3 * unyt.Mpc
    number_of_pixels = int(size / size_pixel + 1)

    data_mask = mask(siminfo.snapshot)
    data_mask.constrain_spatial(region)
    data = load(siminfo.snapshot, mask=data_mask)

    stars_face_on_rotation_matrix = rotation_matrix_from_vector(stars_ang_momentum)
    stars_edge_on_rotation_matrix = rotation_matrix_from_vector(stars_ang_momentum,axis="y")
    gas_face_on_rotation_matrix = rotation_matrix_from_vector(gas_ang_momentum)
    gas_edge_on_rotation_matrix = rotation_matrix_from_vector(gas_ang_momentum,axis="y")


    common_arguments = dict(data=data.stars,
                            boxsize=data.metadata.boxsize,
                            resolution=number_of_pixels,
                            parallel=False,
                            project="masses",
                            region=visualise_region)

    stars_face_on = project_pixel_grid(
        **common_arguments,
        rotation_center=unyt.unyt_array([x, y, z]),
        rotation_matrix=stars_face_on_rotation_matrix,
    )
    stars_face_on = stars_face_on.T

    stars_edge_on = project_pixel_grid(
        **common_arguments,
        rotation_center=unyt.unyt_array([x, y, z]),
        rotation_matrix=stars_edge_on_rotation_matrix,
    )
    stars_edge_on = stars_edge_on.T

    common_arguments = dict(data=data.gas,
                            boxsize=data.metadata.boxsize,
                            resolution=number_of_pixels,
                            parallel=False,
                            region=visualise_region)

    gas_face_on = project_pixel_grid(
        **common_arguments,
        rotation_center=unyt.unyt_array([x, y, z]),
        rotation_matrix=gas_face_on_rotation_matrix,
    )
    gas_face_on = gas_face_on.T

    gas_edge_on = project_pixel_grid(
        **common_arguments,
        rotation_center=unyt.unyt_array([x, y, z]),
        rotation_matrix=gas_edge_on_rotation_matrix,
    )
    gas_edge_on = gas_edge_on.T

    r_img = 0.03
    xmin = -r_img
    ymin = -r_img
    xmax = r_img
    ymax = r_img

    cmap_stars = plt.cm.magma
    cmap_gas = plt.cm.viridis


    # Plot parameters
    params = {
        "font.size": 11,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (7, 4),
        "figure.subplot.left": 0.12,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.1,
        "figure.subplot.top": 0.8,
        "figure.subplot.wspace": 0.4,
        "figure.subplot.hspace": 0.4,
        "lines.markersize": 0.5,
        "lines.linewidth": 0.2,
    }
    rcParams.update(params)

    fig = plt.figure()
    ax = plt.subplot(1, 2, 1)
    title = "Stellar component"
    ax.set_title(title)

    ###### plot one side ########################
    ax.tick_params(labelleft=True, labelbottom=True, length=0)
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    im = ax.imshow(LogNorm()(stars_face_on), cmap=cmap_stars, extent=extent)
    ax.autoscale(False)

    ###### plot the other side ###################
    ax = plt.subplot(1, 2, 2)
    kappa = halo_data.kappa_co[index]
    mass_galaxy = halo_data.stellar_mass[index]
    ac = halo_data.axis_ca[index]
    cb = halo_data.axis_cb[index]
    ba = halo_data.axis_ba[index]
    title = r" $\kappa_{\mathrm{co}} = $%0.2f" % (kappa)
    title += " - $\log_{10}$ $M_{*}/M_{\odot} = $%0.2f" % (mass_galaxy)
    title += " \n c/a = %0.2f," % (ac)
    title += " c/b = %0.2f," % (cb)
    title += " b/a = %0.2f" % (ba)
    ax.set_title(title)

    ax.tick_params(labelleft=True, labelbottom=True, length=0)
    plt.xlabel('x [kpc]')
    plt.ylabel('z [kpc]')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    im = ax.imshow(LogNorm()(stars_edge_on), cmap=cmap_stars, extent=extent)
    ax.autoscale(False)

    outfile = f"{output_path}/galaxy_stars_%i.png" % (index)
    fig.savefig(outfile, dpi=150)
    print("Saved: %s" % (outfile))
    plt.close("all")

    title = "Stellar component" % (index)
    caption = "Face-on (left) and edge-on (right)."
    id = abs(hash("stars galaxy %i" % (index)))
    outfile = "galaxy_stars_%i.png" % (index)
    GalPlotsInWeb.load_plots(title, caption, outfile, id)

    """
    fig = plt.figure()
    ax = plt.subplot(1, 2, 1)
    if parttype == 4: title = "Stellar component"
    if parttype == 0: title = "HI+H2 gas"
    ax.set_title(title)

    cmap_stars = plt.cm.magma
    cmap_gas = plt.cm.viridis

    ###### plot one side ########################
    qv = QuickView(pos_face_on, mass=mass, hsml=hsml_parts, logscale=True, plot=False,
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
    if parttype ==4:
        kappa = halo_data.kappa_co[index]
        mass_galaxy = halo_data.stellar_mass[index]
        ac = halo_data.axis_ca[index]
        cb = halo_data.axis_cb[index]
        ba = halo_data.axis_ba[index]
        title = r" $\kappa_{\mathrm{co}} = $%0.2f" % (kappa)
        title += " - $\log_{10}$ $M_{*}/M_{\odot} = $%0.2f" % (mass_galaxy)
        title += " \n c/a = %0.2f," % (ac)
        title += " c/b = %0.2f," % (cb)
        title += " b/a = %0.2f" % (ba)
    if parttype ==0:
        kappa = halo_data.gas_kappa_co[index]
        mass_galaxy = halo_data.gas_mass[index]
        ac = halo_data.gas_axis_ca[index]
        cb = halo_data.gas_axis_cb[index]
        ba = halo_data.gas_axis_ba[index]
        title = r" $\kappa_{\mathrm{co}} = $%0.2f" % (kappa)
        title += " - $\log_{10}$ $M_{gas}/M_{\odot} = $%0.2f" % (mass_galaxy)
        title += " \n c/a = %0.2f," % (ac)
        title += " c/b = %0.2f," % (cb)
        title += " b/a = %0.2f" % (ba)
    ax.set_title(title)

    qv = QuickView(pos_edge_on, mass=mass, hsml=hsml_parts, logscale=True, plot=False,
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

    if parttype == 0: outfile = f"{output_path}/galaxy_gas_%i.png" % (index)
    if parttype == 4: outfile = f"{output_path}/galaxy_stars_%i.png" % (index)
    fig.savefig(outfile, dpi=150)
    print("Saved: %s" % (outfile))
    plt.close("all")

    if parttype == 0: title = "Galaxy %i / Gas component" % (index)
    if parttype == 4: title = "Galaxy %i / Stellar component" % (index)
    caption = "Face-on (left) and edge-on (right)."
    if parttype == 0: id = abs(hash("galaxy gas %i" % (index)))
    if parttype == 4: id = abs(hash("stars galaxy %i" % (index)))
    if parttype == 0: outfile = "galaxy_gas_%i.png" % (index)
    if parttype == 4: outfile = "galaxy_stars_%i.png" % (index)
    GalPlotsInWeb.load_plots(title, caption, outfile, id)
    """

def visualize_galaxy(stars_data, gas_data, stars_ang_momentum, gas_ang_momentum,
                             halo_data, i, GalPlotsInWeb, output_path):

    #plot_galaxy(stars_data, 4, stars_ang_momentum, halo_data, i, GalPlotsInWeb, output_path)
    #plot_galaxy(gas_data, 0, gas_ang_momentum, halo_data, i, GalPlotsInWeb, output_path)

    plot_galaxy_parts(stars_data, 4, stars_ang_momentum, halo_data, i, GalPlotsInWeb, output_path)
    plot_galaxy_parts(gas_data, 0, gas_ang_momentum, halo_data, i, GalPlotsInWeb, output_path)
