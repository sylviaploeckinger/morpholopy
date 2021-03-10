from pylab import *
from sphviewer.tools import QuickView
from swiftsimio.visualisation.rotation import rotation_matrix_from_vector
from numpy import matmul

def rotation_matrix(vector: float64, pos_parts: float64, axis: str = "z"):
    normed_vector = vector / norm(vector)

    # Directional vector describing the axis we wish to look 'down'
    original_direction = array([0.0, 0.0, 0.0], dtype=float64)
    switch = {"x": 0, "y": 1, "z": 2}

    try:
        original_direction[switch[axis]] = 1.0
    except KeyError:
        raise ValueError(
            f"Parameter axis must be one of x, y, or z. You supplied {axis}."
        )

    dot_product = dot(original_direction, normed_vector)
    cross_product = cross(original_direction, normed_vector)
    mod_cross_product = norm(cross_product)
    cross_product /= mod_cross_product
    theta = arccos(dot_product)

    q0 = cos(theta/2)
    q1 = sin(theta/2) * cross_product[0]
    q2 = sin(theta/2) * cross_product[1]
    q3 = sin(theta/2) * cross_product[2]

    # Skew symmetric matrix for cross product
    Q = array(
        [
            [q0**2+q1**2-q2**2-q3**2, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2)],
            [2*(q2*q1+q0*q3), q0**2-q1**2+q2**2-q3**2, 2*(q3*q2-q0*q1)],
            [2*(q1*q3-q0*q2), 2*(q3*q2+q0*q1), q0**2-q1**2-q2**2+q3**2],
        ]
    )

    i = np.array([1,0,0])
    j = np.array([0,1,0])
    k = np.array([0,0,1])

    u = matmul(Q, i)
    v = matmul(Q, j)
    w = matmul(Q, k)

    pos_face_on = pos_parts.copy()
    pos_face_on[:,0] = dot(pos_parts, u)
    pos_face_on[:,1] = dot(pos_parts, v)
    pos_face_on[:,2] = dot(pos_parts, w)

    return pos_face_on

def plot_galaxy_parts(partsDATA, parttype, ang_momentum, halo_data, index, PlotsInWeb, output_path):

    # Plot ranges
    r_limit = 5 * halo_data.halfmass_radius_star[index]
    r_img = 30.
    if r_limit < r_img: r_img = r_limit
    xmin = -r_img
    ymin = -r_img
    xmax = r_img
    ymax = r_img

    pos_parts = partsDATA[:, 0:3].copy()
    check_axis = np.where(np.abs(ang_momentum) == np.max(np.abs(ang_momentum)))[0]
    main_axis = ["x","y","z"]
    secondary_axis = ["y","z","x"]

    pos_face_on = rotation_matrix(ang_momentum,pos_parts,axis=main_axis[check_axis[0]])
    pos_edge_on = rotation_matrix(ang_momentum,pos_parts,axis=secondary_axis[check_axis[0]])
    #face_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum,axis=main_axis[check_axis[0]])
    #edge_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum,axis=secondary_axis[check_axis[0]])

    #pos_face_on = matmul(face_on_rotation_matrix, pos_parts.T)
    #pos_face_on = pos_face_on.T
    #pos_edge_on = matmul(edge_on_rotation_matrix, pos_parts.T)
    #pos_edge_on = pos_edge_on.T


    if parttype == 0: density = np.log10(partsDATA[:,11])
    if parttype == 4: density = np.log10(partsDATA[:, 8])

    # Sort particles for better viewing
    arg_sort = np.argsort(density)
    density = density[arg_sort]
    pos_face_on = pos_face_on[arg_sort,:]
    pos_edge_on = pos_edge_on[arg_sort,:]

    # Plot parameters
    params = {
        "font.size": 11,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (7.5, 4),
        "figure.subplot.left": 0.1,
        "figure.subplot.right": 0.85,
        "figure.subplot.bottom": 0.15,
        "figure.subplot.top": 0.8,
        "figure.subplot.wspace": 0.3,
        "figure.subplot.hspace": 0.3,
        "lines.markersize": 0.2,
        "lines.linewidth": 0.2,
    }
    rcParams.update(params)

    fig = plt.figure()
    ax = plt.subplot(1, 2, 1)
    if parttype ==4: title = "Stellar component"
    if parttype ==0: title = "Gas component"
    ax.set_title(title)
    ax.tick_params(labelleft=True, labelbottom=True, length=0)
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    plt.scatter(pos_face_on[:, 0], pos_face_on[:, 1], c=density, alpha=1, s=10,
                vmin=6, vmax=10, cmap='magma', edgecolors='none')
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
    ax.set_title(title)

    ax.tick_params(labelleft=True, labelbottom=True, length=0)
    plt.xlabel('x [kpc]')
    plt.ylabel('z [kpc]')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    plt.scatter(pos_edge_on[:, 0], pos_edge_on[:, 1], c=density, alpha=1, s=10,
                vmin=6, vmax=10, cmap='magma', edgecolors='none')

    ax.autoscale(False)

    cbar_ax = fig.add_axes([0.86, 0.22, 0.018, 0.5])
    cbar_ax.tick_params(labelsize=15)
    cb = plt.colorbar(ticks=[4, 6, 8, 10], cax=cbar_ax)
    cb.set_label(label=r'$\log_{10}$ $\rho$ [M$_{\odot}$/kpc$^{3}$]', labelpad=0.5)

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
        #vmin = np.min(image[image>0])
        #image[image==0] = vmin
        image[image<4] = 4
    #if(vmax==None):
        #vmax = np.max(image)
        #image[image>8] = 8
    #image=np.clip(image,vmin,vmax)
    #image=(image-vmin)/(vmax-vmin)
    return image

def plot_galaxy(parts_data, parttype, ang_momentum, halo_data, index, GalPlotsInWeb, output_path):
    # partsDATA contains particles data and is structured as follow
    # [ (:3)Position[Mpc]: (0)X | (1)Y | (2)Z ]
    if parttype == 4: cmap = plt.cm.magma
    if parttype == 0: cmap = plt.cm.viridis

    pos_parts = parts_data[:, 0:3].copy()
    check_axis = np.where(np.abs(ang_momentum) == np.max(np.abs(ang_momentum)))[0]
    main_axis = ["x","y","z"]
    secondary_axis = ["y","z","x"]

    pos_face_on = rotation_matrix(ang_momentum,pos_parts,axis=main_axis[check_axis[0]])
    pos_edge_on = rotation_matrix(ang_momentum,pos_parts,axis=secondary_axis[check_axis[0]])

    #face_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum,axis=main_axis[check_axis[0]])
    #edge_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum,axis=secondary_axis[check_axis[0]])

    #pos_face_on = matmul(face_on_rotation_matrix, pos_parts.T)
    #pos_face_on = pos_face_on.T
    #pos_edge_on = matmul(edge_on_rotation_matrix, pos_parts.T)
    #pos_edge_on = pos_edge_on.T

    hsml_parts = parts_data[:, 7]
    mass = parts_data[:, 3]

    r_limit = 5 * halo_data.halfmass_radius_star[index]
    r_img = 30.
    if r_limit < r_img: r_img = r_limit
    xmin = -r_img
    ymin = -r_img
    xmax = r_img
    ymax = r_img

    # Plot parameters
    params = {
        "font.size": 11,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (7.5, 4),
        "figure.subplot.left": 0.1,
        "figure.subplot.right": 0.85,
        "figure.subplot.bottom": 0.15,
        "figure.subplot.top": 0.8,
        "figure.subplot.wspace": 0.3,
        "figure.subplot.hspace": 0.3,
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
                   r='infinity', p=0, t=0, extent=[xmin, xmax, ymin, ymax],
                   x=0, y=0, z=0)
    img = qv.get_image()
    ext = qv.get_extent()
    ax.tick_params(labelleft=True, labelbottom=True, length=0)
    plt.xlabel('x [kpc]')
    plt.ylabel('z [kpc]')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    img = get_normalized_image(img)
    ims = ax.imshow(img, cmap=cmap, extent=ext)
    ax.autoscale(False)

    cbar_ax = fig.add_axes([0.86, 0.22, 0.018, 0.5])
    cbar_ax.tick_params(labelsize=15)
    cb = plt.colorbar(ims, ticks=[4, 6, 8, 10, 12], cax=cbar_ax)
    cb.set_label(label=r'$\log_{10}$ $\Sigma$ [M$_{\odot}$/kpc$^{2}$]', labelpad=0.5)


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


def visualize_galaxy(stars_data, gas_data, stars_ang_momentum, gas_ang_momentum,
                             halo_data, i, GalPlotsInWeb, output_path):

    plot_galaxy(stars_data, 4, stars_ang_momentum, halo_data, i, GalPlotsInWeb, output_path)
    plot_galaxy(gas_data, 0, stars_ang_momentum, halo_data, i, GalPlotsInWeb, output_path)

    plot_galaxy_parts(stars_data, 4, stars_ang_momentum, halo_data, i, GalPlotsInWeb, output_path)
    plot_galaxy_parts(gas_data, 0, stars_ang_momentum, halo_data, i, GalPlotsInWeb, output_path)
