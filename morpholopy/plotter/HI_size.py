import numpy as np
from swiftsimio.visualisation.rotation import rotation_matrix_from_vector
from swiftsimio.visualisation.projection_backends import backends, backends_parallel
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl
import scipy.optimize as opt


def gauss_curve(x, A, a, b, c, x0, y0):
    """
    General 2D elliptical Gaussian function, based on Rajohnson et al. (2022),
    equation (2).

    Parameters
    ----------
    x: (N,2) array containing positions (in kpc)

    A: central surface density (in Msun pc^-2)
    a, b, c: rotated major and minor axis of the ellipse (in kpc^-2)
    x0, y0: centre of the ellipse (in kpc)

    Returns
    -------
    (N,) array containing the predicted surface density for each position
    """
    return A * np.exp(
        -(
            a * (x[:, 0] - x0) ** 2
            + 2.0 * b * (x[:, 0] - x0) * (x[:, 1] - y0)
            + c * (x[:, 1] - y0) ** 2
        )
    )


def calculate_HI_size(
    data, ang_momentum, galaxy_data, output_path, index, resolution=64
):
    """
    Calculate the HI size of the given galaxy by fitting a 2D elliptical function
    to the surface density profile and measuring the diameter along the major axis
    of the 1 Msun pc^-2 contour (see Rajohnson et al. (2022) section 4).
    The result is stored in the given HaloCatalogue (in units of kpc).

    Parameters
    ----------
    data: gas data (as returned by simulation_data.make_particle_data() - note
          that the positions have been recentred onto the centre of the halo)
    ang_momentum: angular momentum vector of the galaxy
    galaxy_data: HaloCatalogue object to store the HI size in
    output_path: directory in which to store images
    index: halo index, used to index images
    resolution: resolution of the image that is used to fit the elliptical profile
                the default value (128) provides a good mix of accuracy and speed
    """

    # get a local copy of the particle coordinates (that we can rotate)
    pos_parts = data[:, 0:3].copy()

    # determine the size of the galaxy from the maximum of the coordinates
    R = np.sqrt(
        pos_parts[:, 0] ** 2 + pos_parts[:, 1] ** 2 + pos_parts[:, 2] ** 2
    ).max()

    # rotate the coordinates to a face-on projection
    face_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum, axis="z")
    pos_face_on = np.matmul(face_on_rotation_matrix, pos_parts.T).T

    # get the smoothing lengths and HI mass
    hsml = data[:, 7]
    mass = data[:, 8]

    # rescale the coordinates to a square box [-R,R] -> [0, 1]
    x = 0.5 * (pos_face_on[:, 0] + R) / R
    y = 0.5 * (pos_face_on[:, 1] + R) / R
    h = 0.5 * hsml / R

    # create the projected image using the swiftsimio backend
    image = backends["subsampled"](x=x, y=y, m=mass, h=h, res=resolution)
    # units are Msun / (2*R)^2
    # convert to Msun / pc^2
    image *= 2.5e-7 / R ** 2

    # get the limits (for plotting)
    imax = image.max()
    imin = max(1.0e-4 * imax, image.min())

    # set up a grid of positions for the fit
    xg, yg = np.meshgrid(np.linspace(-R, R, resolution), np.linspace(-R, R, resolution))
    xs = np.zeros((resolution ** 2, 2))
    xs[:, 0] = xg.flatten()
    xs[:, 1] = yg.flatten()
    # initial values for the fit (based on Rajohnson et al.)
    A = imax
    sigX = 40.0 * R / resolution
    sigY = 40.0 * R / resolution
    a = 0.5 / sigX ** 2
    b = 0.0
    c = 0.5 / sigY ** 2
    # perform the fit
    params, _ = opt.curve_fit(
        gauss_curve, xs, image.flatten(), p0=(A, a, b, c, 0.0, 0.0)
    )

    # extract the fitted parameters
    A, a, b, c, x0, y0 = params

    # the HI size is determined from the 1 Msun pc^-2 contour level. If the fitted
    # profile has a central surface density that is lower, the HI size is undefined
    if A <= 1.0:
        print(
            "Central surface density below 1 Msun pc^-2 limit, no HI size measurement possible!"
        )
    else:
        # convert from the general elliptical coordinates to a coordinate frame where
        # the major axis is aligned with the x axis (and the minor axis with the y
        # axis). Determine the rotation angle theta of the general ellipse.
        # The relations below have been derived from Rajohnson et al., eq (3)-(5),
        # using some basic trigonometry.
        # Note that the derivation depends on the condition 1/sigY2 > 1/sigX2, so
        # that sigX2 is guaranteed to correspond to the major axis
        sigX2 = 1.0 / (a + c - np.sqrt((a - c) ** 2 + 4.0 * b ** 2))
        sigY2 = 1.0 / (np.sqrt((a - c) ** 2 + 4.0 * b ** 2) + a + c)
        theta = 0.5 * np.arctan2(-2.0 * b, c - a)
        Dx = np.sqrt(2.0 * sigX2 * np.log(A))
        Dy = np.sqrt(2.0 * sigY2 * np.log(A))

        if Dx < Dy:
            print("Error: major axis smaller than minor axis!")

        # Compute the HI size
        HIsize = 2.0 * Dx
        # Compute the HI mass as the integrated surface density of the image
        HImass = image.sum() * R ** 2 * 4.0e6 / resolution ** 2
        # Set the values for this galaxy
        galaxy_data.add_HI_size_mass(HIsize, HImass, index)

    rcparams = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (5.5, 4),
        "figure.subplot.left": 0.05,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.15,
        "figure.subplot.top": 0.9,
        "figure.subplot.wspace": 0.3,
        "figure.subplot.hspace": 0.3,
        "lines.markersize": 2,
        "lines.linewidth": 1.0,
        "xtick.top": True,
        "ytick.right": True,
    }
    pl.rcParams.update(rcparams)

    # Create a figure to show how good the fit was
    fig, ax = pl.subplots(1, 1)

    # First, plot the surface density map
    levels = np.logspace(np.log10(imin), np.log10(imax), 100)
    cs = ax.contourf(
        xg,
        yg,
        image,
        norm=matplotlib.colors.LogNorm(vmin=imin, vmax=imax),
        levels=levels,
    )

    if A > 1.0:
        # Now overplot the 1 Msun pc^-2 contour, using the parametric equation for
        # the general ellipse
        tpar = np.linspace(0.0, 2.0 * np.pi, 1000)
        xpar = Dx * np.cos(theta) * np.cos(tpar) - Dy * np.sin(theta) * np.sin(tpar)
        ypar = Dx * np.sin(theta) * np.cos(tpar) + Dy * np.cos(theta) * np.sin(tpar)
        ax.plot(x0 + xpar, y0 + ypar, color="w", linestyle="--")

    # Overplot 9 contour levels of the Gaussian function
    image = gauss_curve(xs, *params).reshape((resolution, resolution))
    levels = np.logspace(np.log10(imin), np.log10(imax), 10)
    ax.contour(
        xg,
        yg,
        image,
        norm=matplotlib.colors.LogNorm(vmin=imin, vmax=imax),
        levels=levels,
        colors="w",
    )

    # Finalise and save the plot
    ax.set_aspect("equal")
    ax.set_xlabel("x (kpc)")
    ax.set_ylabel("y (kpc)")

    fig.colorbar(
        cs,
        label="Surface density (M$_\\odot{}$ pc$^{-2}$)",
        ticks=matplotlib.ticker.LogLocator(),
        format=matplotlib.ticker.LogFormatterMathtext(),
    )

    pl.tight_layout()
    pl.savefig(f"{output_path}/HI_size_{index}.png", dpi=300)
    pl.close()


def plot_HI_size_mass(output_path, name_list):

    fit_slope = 0.501
    fit_intercept = -3.252
    fit_intercept_sigma_low = 0.074
    fit_intercept_sigma_high = 0.073

    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (5.5, 4),
        "figure.subplot.left": 0.05,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.15,
        "figure.subplot.top": 0.9,
        "figure.subplot.wspace": 0.3,
        "figure.subplot.hspace": 0.3,
        "lines.markersize": 4,
        "lines.linewidth": 2.0,
        "xtick.top": True,
        "xtick.labeltop": True,
        "ytick.right": True,
        "ytick.labelright": True,
    }
    pl.rcParams.update(params)

    logMfit = np.array([6.5, 12.0])
    pl.fill_between(
        10.0 ** logMfit,
        10.0 ** (fit_slope * logMfit + fit_intercept - fit_intercept_sigma_low),
        10.0 ** (fit_slope * logMfit + fit_intercept + fit_intercept_sigma_high),
        alpha=0.1,
        color="k",
    )
    pl.plot(
        10.0 ** logMfit,
        10.0 ** (fit_slope * logMfit + fit_intercept),
        "--",
        color="k",
        label="Rajohnson et al. (2022)",
    )

    for i, name in enumerate(name_list):
        data = np.loadtxt(
            f"{output_path}/HI_size_mass_{name}.txt",
            dtype=[("size", "f8"), ("mass", "f8")],
        )
        pl.loglog(data["mass"], data["size"], "o", color=f"C{i}", label=name)

    pl.xlabel("M$_{\\rm{}HI}$ (M$_\\odot$)")
    pl.ylabel("D$_{\\rm{}HI}$ (kpc)")

    pl.legend(loc="best")

    pl.tight_layout()
    pl.savefig(f"{output_path}/HI_size_mass.png", dpi=200)
    pl.close()
