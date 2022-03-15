from swiftsimio import load, mask
from velociraptor import load as load_catalogue
from swiftsimio.visualisation.rotation import rotation_matrix_from_vector
from swiftsimio.visualisation.projection import project_gas
from swiftsimio.visualisation.projection import project_pixel_grid
from swiftsimio.visualisation.slice import kernel_gamma
from swiftsimio.visualisation.smoothing_length_generation import (
    generate_smoothing_lengths,
)
from swiftsimio import swift_cosmology_to_astropy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import gridspec

import unyt
from unyt import pc, kpc, msun, cm, yr
from unyt import proton_mass_cgs as mH
from unyt import boltzmann_constant_cgs as kB
from unyt import gravitational_constant_cgs as G
from unyt import unyt_array

from astropy.visualization import make_lupton_rgb

from scipy.optimize import curve_fit

def exponential(x, Sigma0, H, offset):
        return Sigma0 * np.exp(-np.abs(x+offset)/H)

def calculate_scaleheight_fit(mass_map, r_img_kpc):
        r_abs_max_kpc = 7.5
        xx = np.linspace(-r_img_kpc.value, r_img_kpc.value, len(mass_map[:,0]),endpoint=True)
        z  = (np.tile(xx, (len(xx),1))).T
        z_1D = np.ravel(z       [:, (np.abs(xx) < r_abs_max_kpc)])
        S_1D = np.ravel(mass_map[:, (np.abs(xx) < r_abs_max_kpc)])

        popt, pcov = curve_fit(exponential, z_1D[np.isfinite(S_1D)], S_1D[np.isfinite(S_1D)])

        return popt[1]

def get_radial_profile(mass_map, radialbin_kpc, pixsize_kpc, r_img_kpc):
        nbins_radial = int(r_img_kpc.value / radialbin_kpc)

        y, x = np.indices((mass_map.shape))
        npix = len(mass_map[0,:])
        pcenter = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
        r = np.hypot(x - pcenter[0], y - pcenter[1])
        r_1D = r.ravel() * pixsize_kpc

        # from mass per pc-2 to mass per pixel
        mass_map_per_pixel = mass_map * pixsize_kpc * pixsize_kpc * kpc * kpc
        mass_map_per_pixel.convert_to_units(msun)
        mass_map_1D = (mass_map_per_pixel.ravel()).value

        rmin = 0.
        rmax = nbins_radial * radialbin_kpc

        M_hist, r_hist = np.histogram(r_1D, weights = mass_map_1D, bins = nbins_radial, range = (rmin, rmax), density = False)

        # calculate the area of each radial bin
        A = np.zeros(len(r_hist)-1)
        for i in range(len(r_hist)-1):
                A[i] = np.pi * (r_hist[i+1] * r_hist[i+1] - r_hist[i] * r_hist[i])

        A_unyt = A * kpc * kpc
        Sigma = M_hist * msun / A_unyt
        Sigma.convert_to_units(msun / pc**2)

        r_central = np.zeros(len(r_hist) - 1)
        for i in range(len(r_hist)-1):
                r_central[i] = 0.5 * (r_hist[i+1] + r_hist[i])

        r_central = r_central * kpc

        return r_central, Sigma

def get_Schruba_data():
        S_HI, S_H2, flag = np.loadtxt("plotter/obs_data/Schruba_2011_original_data.txt", skiprows = 14, usecols = (3,5,8), unpack = True)
        return S_HI, S_H2, flag

def get_Schruba_upperlimits():
        S_HI, S_H2, flag = np.loadtxt("plotter/obs_data/Schruba_2011_original_data.txt", skiprows = 14, usecols = (3,7,8), unpack = True)
        return S_HI, S_H2, flag

def get_angular_momentum_vector(rparticles, vparticles, rgalaxy, vgalaxy, mparticles):
        # make sure everything is in the same unit system
        rparticles = rparticles.to_physical()
        vparticles = vparticles.to_physical() 

        r = rparticles - rgalaxy
        v = vparticles - vgalaxy
        m = mparticles

        d = np.linalg.norm(r, axis = 1)

        # mask out the innermost x% of star particles
        dmin = np.quantile(d, 0.3, axis = 0)
        dmax = np.quantile(d, 0.5, axis = 0)
        m[d<dmin] = 0.
        m[d>dmax] = 0.

        L = np.cross(r,v)
        L[:,0] *= m
        L[:,1] *= m 
        L[:,2] *= m

        Ltotal = np.sum(L, axis = 0)
        
        Ltotal = Ltotal / np.linalg.norm(Ltotal)

        face_on_rotation_matrix = rotation_matrix_from_vector(
           Ltotal
        )
        edge_on_rotation_matrix = rotation_matrix_from_vector(
           Ltotal,
           axis="y"
        )

        return face_on_rotation_matrix, edge_on_rotation_matrix


def get_gas_surface_density_map(catalogue, halo_id, plottype, snapshot_filename, size, npixlocal):
        # center of the halo in physical coordinates
        x = catalogue.positions.xcmbp[halo_id]
        y = catalogue.positions.ycmbp[halo_id]
        z = catalogue.positions.zcmbp[halo_id]

        # angular momentum of the stars (for projection)
        lx = catalogue.angular_momentum.lx_star[halo_id]
        ly = catalogue.angular_momentum.ly_star[halo_id]
        lz = catalogue.angular_momentum.lz_star[halo_id]

        angular_momentum_vector = np.array([lx.value, ly.value, lz.value])
        angular_momentum_vector /= np.linalg.norm(angular_momentum_vector)

        #face_on_rotation_matrix = rotation_matrix_from_vector(
        #   angular_momentum_vector
        #)
        #edge_on_rotation_matrix = rotation_matrix_from_vector(
        #   angular_momentum_vector,
        #   axis="y"
        #)
        
        # needs to be in comoving coordinates for the mask
        region = [
           [x / catalogue.scale_factor - size / catalogue.scale_factor, x / catalogue.scale_factor + size / catalogue.scale_factor],
           [y / catalogue.scale_factor - size / catalogue.scale_factor, y / catalogue.scale_factor + size / catalogue.scale_factor],
           [z / catalogue.scale_factor - size / catalogue.scale_factor, z / catalogue.scale_factor + size / catalogue.scale_factor],
        ]

        visualise_region = [
           x - 0.5 * size, x + 0.5 * size,
           y - 0.5 * size, y + 0.5 * size,
        ]

        data_mask = mask(snapshot_filename)
        data_mask.constrain_spatial(region)
        data = load(snapshot_filename, mask=data_mask)
        data.gas.coordinates = data.gas.coordinates.to_physical()


        face_on_rotation_matrix, edge_on_rotation_matrix = get_angular_momentum_vector(data.stars.coordinates, \
                                                                                       data.stars.velocities, \
                                                                                       [x,y,z], \
                                                                                       [catalogue.velocities.vxcmbp[halo_id], catalogue.velocities.vycmbp[halo_id], catalogue.velocities.vzcmbp[halo_id]], 
                                                                                       data.stars.masses)

        if plottype == 'HI':
                data.gas.usermass = data.gas.masses * data.gas.species_fractions.HI * data.gas.element_mass_fractions.hydrogen
        elif plottype == 'H2':
                data.gas.usermass = data.gas.masses * 2. * data.gas.species_fractions.H2 * data.gas.element_mass_fractions.hydrogen
        elif 'GK11' in  plottype:
                data.gas.usermass = get_GK11fractions(data, plottype) * data.gas.masses * data.gas.element_mass_fractions.hydrogen
        elif 'hydrogen' in plottype:
                data.gas.usermass = data.gas.masses * data.gas.element_mass_fractions.hydrogen
        elif 'helium' in plottype:
                data.gas.usermass = data.gas.masses * data.gas.element_mass_fractions.helium
        elif 'totaloxygen' in plottype:
                data.gas.usermass = data.gas.masses * data.gas.element_mass_fractions.oxygen
        elif 'diffuseoxygen' in plottype:
                data.gas.usermass = data.gas.diffuse_oxygen_masses_from_table
        elif 'sfr' in plottype:
                data.gas.star_formation_rates.convert_to_units(msun / yr)
                data.gas.star_formation_rates[data.gas.star_formation_rates<0.] = 0.
                data.gas.usermass = data.gas.star_formation_rates
        else:
                print ('Unknown plottype: ', plottype)
                import sys
                sys.exit()

        if not 'sfr' in plottype:
                data.gas.usermass.convert_to_units("Msun")

        # Face on projection
        mass_map_face = project_gas(data, resolution=int(npixlocal), project="usermass", parallel=True, region = visualise_region,
                               rotation_center=unyt.unyt_array([x, y, z]), rotation_matrix=face_on_rotation_matrix, backend = "subsampled")

        if 'sfr' in plottype:
                mass_map_face.convert_to_units(msun / yr / pc**2)
        else:
                mass_map_face.convert_to_units(msun / pc**2)

        pixelsize = (visualise_region[1] - visualise_region[0]) / float(npixlocal)
        totalmass = np.sum(mass_map_face) * pixelsize * pixelsize
        if 'sfr' in plottype:
                totalmass.convert_to_units('Msun/yr')
        else:
                totalmass.convert_to_units('Msun')

        # Edge on projection
        mass_map_edge = project_gas(data, resolution=int(npixlocal), project="usermass", parallel=True, region = visualise_region,
                               rotation_center=unyt.unyt_array([x, y, z]), rotation_matrix=edge_on_rotation_matrix, backend = "subsampled")

        if 'sfr' in plottype:
                mass_map_edge.convert_to_units(msun / yr / pc**2)
        else:
                mass_map_edge.convert_to_units(msun / pc**2)

        # mask with circle
        lx, ly = mass_map_edge.shape
        X, Y = np.ogrid[0:lx, 0:ly]
        mask_circle = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 > lx * ly / 4
        mass_map_edge[mask_circle] = np.nan

        # mask with circle
        lx, ly = mass_map_face.shape
        X, Y = np.ogrid[0:lx, 0:ly]
        mask_circle = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 > lx * ly / 4
        mass_map_face[mask_circle] = np.nan
 
        return mass_map_face.T, mass_map_edge.T, visualise_region, x, y, totalmass
         
def get_stars_surface_brightness_map(catalogue, halo_id, snapshot_filename, size, npix, r_img_kpc):
        # center of the halo
        x = catalogue.positions.xcmbp[halo_id]
        y = catalogue.positions.ycmbp[halo_id]
        z = catalogue.positions.zcmbp[halo_id]

        # angular momentum of the stars (for projection)
        lx = catalogue.angular_momentum.lx_star[halo_id]
        ly = catalogue.angular_momentum.ly_star[halo_id]
        lz = catalogue.angular_momentum.lz_star[halo_id]

        angular_momentum_vector = np.array([lx.value, ly.value, lz.value])
        angular_momentum_vector /= np.linalg.norm(angular_momentum_vector)

        #face_on_rotation_matrix = rotation_matrix_from_vector(
        #   angular_momentum_vector
        #)
        #edge_on_rotation_matrix = rotation_matrix_from_vector(
        #   angular_momentum_vector,
        #   axis="y"
        #)

        # needs to be in comoving coordinates for the mask
        region = [
           [x / catalogue.scale_factor - size / catalogue.scale_factor, x / catalogue.scale_factor + size / catalogue.scale_factor],
           [y / catalogue.scale_factor - size / catalogue.scale_factor, y / catalogue.scale_factor + size / catalogue.scale_factor],
           [z / catalogue.scale_factor - size / catalogue.scale_factor, z / catalogue.scale_factor + size / catalogue.scale_factor],
        ]

        visualise_region = [
           x - 0.5 * size, x + 0.5 * size,
           y - 0.5 * size, y + 0.5 * size,
        ]

        data_mask = mask(snapshot_filename)
        data_mask.constrain_spatial(region)
        data = load(snapshot_filename, mask=data_mask)
        data.stars.coordinates = data.stars.coordinates.to_physical()

        data.stars.smoothing_lengths = generate_smoothing_lengths(
                coordinates=data.stars.coordinates,
                boxsize=data.metadata.boxsize,
                kernel_gamma=kernel_gamma,
                neighbours=11,
                speedup_fac=1,
                dimension=3,
            )

        face_on_rotation_matrix, edge_on_rotation_matrix = get_angular_momentum_vector(data.stars.coordinates, \
                                                                                       data.stars.velocities, \
                                                                                       [x,y,z], \
                                                                                       [catalogue.velocities.vxcmbp[halo_id], catalogue.velocities.vycmbp[halo_id], catalogue.velocities.vzcmbp[halo_id]], 
                                                                                       data.stars.masses)
        
        luminosities = [data.stars.luminosities.GAMA_i, data.stars.luminosities.GAMA_r, data.stars.luminosities.GAMA_g]
        rgb_image_face = np.zeros((npix, npix, len(luminosities)))

        for ilum in range(len(luminosities)):
                # Face on projection
                data.stars.usermass = luminosities[ilum]
                pixel_grid = project_pixel_grid(data.stars, resolution=int(npix), project='usermass', parallel=True, region = visualise_region,
                                       rotation_center=unyt.unyt_array([x, y, z]), rotation_matrix=face_on_rotation_matrix, boxsize=data.metadata.boxsize, backend = "subsampled")

                x_range = visualise_region[1] - visualise_region[0]
                y_range = visualise_region[3] - visualise_region[2]
                units = 1.0 / (x_range * y_range)
                # Unfortunately this is required to prevent us from {over,under}flowing
                # the units...
                units.convert_to_units(1.0 / (x_range.units * y_range.units))

                mass_map_face = unyt_array(pixel_grid, units=units)
                mass_map_face.convert_to_units(1. / pc**2)
                try:
                        mass_map_face[mass_map_face == 0.] = mass_map_face[mass_map_face > 0.].min()
                except:
                        mass_map_face[mass_map_face == 0.] = 1.e-10


                rgb_image_face[:,:,ilum] = mass_map_face.T

        image_face = make_lupton_rgb(rgb_image_face[:,:,0], rgb_image_face[:,:,1], rgb_image_face[:,:,2], Q=10, stretch=0.5)
        # mask with circle
        lx, ly = mass_map_face.shape
        X, Y = np.ogrid[0:lx, 0:ly]
        mask_circle = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 > lx * ly / 4
        image_face[mask_circle,:] = 255

        H_kpc_gri = np.zeros(len(luminosities))
        rgb_image_edge = np.zeros((npix, npix, len(luminosities)))
        for ilum in range(len(luminosities)):
                # Face on projection
                data.stars.usermass = luminosities[ilum]
                pixel_grid = project_pixel_grid(data.stars, resolution=int(npix), project='usermass', parallel=True, region = visualise_region,
                                       rotation_center=unyt.unyt_array([x, y, z]), rotation_matrix=edge_on_rotation_matrix, boxsize=data.metadata.boxsize, backend = "subsampled")

                x_range = visualise_region[1] - visualise_region[0]
                y_range = visualise_region[3] - visualise_region[2]
                units = 1.0 / (x_range * y_range)
                # Unfortunately this is required to prevent us from {over,under}flowing
                # the units...
                units.convert_to_units(1.0 / (x_range.units * y_range.units))

                mass_map_edge = unyt_array(pixel_grid, units=units)
                mass_map_edge.convert_to_units(1. / pc**2)
                try:
                        mass_map_edge[mass_map_edge == 0.] = mass_map_edge[mass_map_edge > 0.].min()
                except:
                        mass_map_edge[mass_map_edge == 0.] = 1.e-10

                try:
                        H_kpc_gri[ilum] = calculate_scaleheight_fit ( mass_map_edge.T, r_img_kpc )
                except:
                        H_kpc_gri[ilum] = -1.
                rgb_image_edge[:,:,ilum] = mass_map_edge.T

        print('H (gri): ', H_kpc_gri)

        image_edge = make_lupton_rgb(rgb_image_edge[:,:,0], rgb_image_edge[:,:,1], rgb_image_edge[:,:,2], Q=10, stretch=0.5)
        # mask with circle
        lx, ly = mass_map_edge.shape
        X, Y = np.ogrid[0:lx, 0:ly]
        mask_circle = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 > lx * ly / 4
        image_edge[mask_circle,:] = 255

        return image_face, image_edge, visualise_region, x, y, -1., H_kpc_gri

def get_stars_surface_density_map(catalogue, halo_id, plottype, snapshot_filename, size, npixloc):
        # center of the halo
        x = catalogue.positions.xcmbp[halo_id]
        y = catalogue.positions.ycmbp[halo_id]
        z = catalogue.positions.zcmbp[halo_id]

        # angular momentum of the stars (for projection)
        lx = catalogue.angular_momentum.lx_star[halo_id]
        ly = catalogue.angular_momentum.ly_star[halo_id]
        lz = catalogue.angular_momentum.lz_star[halo_id]

        angular_momentum_vector = np.array([lx.value, ly.value, lz.value])
        angular_momentum_vector /= np.linalg.norm(angular_momentum_vector)

        #face_on_rotation_matrix = rotation_matrix_from_vector(
        #   angular_momentum_vector
        #)
        #edge_on_rotation_matrix = rotation_matrix_from_vector(
        #   angular_momentum_vector,
        #   axis="y"
        #)

        # needs to be in comoving coordinates for the mask
        region = [
           [x / catalogue.scale_factor - size / catalogue.scale_factor, x / catalogue.scale_factor + size / catalogue.scale_factor],
           [y / catalogue.scale_factor - size / catalogue.scale_factor, y / catalogue.scale_factor + size / catalogue.scale_factor],
           [z / catalogue.scale_factor - size / catalogue.scale_factor, z / catalogue.scale_factor + size / catalogue.scale_factor],
        ]

        visualise_region = [
           x - 0.5 * size, x + 0.5 * size,
           y - 0.5 * size, y + 0.5 * size,
        ]

        data_mask = mask(snapshot_filename)
        data_mask.constrain_spatial(region)
        data = load(snapshot_filename, mask=data_mask)
        data.stars.coordinates = data.stars.coordinates.to_physical()

        face_on_rotation_matrix, edge_on_rotation_matrix = get_angular_momentum_vector(data.stars.coordinates, \
                                                                                       data.stars.velocities, \
                                                                                       [x,y,z], \
                                                                                       [catalogue.velocities.vxcmbp[halo_id], catalogue.velocities.vycmbp[halo_id], catalogue.velocities.vzcmbp[halo_id]], 
                                                                                       data.stars.masses)
        data.stars.smoothing_lengths = generate_smoothing_lengths(
                coordinates=data.stars.coordinates,
                boxsize=data.metadata.boxsize,
                kernel_gamma=kernel_gamma,
                neighbours=11,
                speedup_fac=1,
                dimension=3,
            )        

        # Face on projection
        pixel_grid = project_pixel_grid(data.stars, resolution=int(npixloc), project="masses", parallel=True, region = visualise_region,
                               rotation_center=unyt.unyt_array([x, y, z]), rotation_matrix=face_on_rotation_matrix, boxsize=data.metadata.boxsize, backend = "subsampled")

        x_range = visualise_region[1] - visualise_region[0]
        y_range = visualise_region[3] - visualise_region[2]
        units = 1.0 / (x_range * y_range)
        # Unfortunately this is required to prevent us from {over,under}flowing
        # the units...
        units.convert_to_units(1.0 / (x_range.units * y_range.units))
        units *= getattr(data.stars, 'masses').units

        mass_map_face = unyt_array(pixel_grid, units=units)
        mass_map_face.convert_to_units(msun / pc**2)
        mass_map_face[mass_map_face == 0.] = 1.e-20

        pixelsize = (visualise_region[1] - visualise_region[0]) / float(npixloc)
        totalmass = np.sum(mass_map_face) * pixelsize * pixelsize
        totalmass.convert_to_units('Msun')

        # Edge on projection
        pixel_grid = project_pixel_grid(data.stars, resolution=int(npixloc), project="masses", parallel=True, region = visualise_region,
                               rotation_center=unyt.unyt_array([x, y, z]), rotation_matrix=edge_on_rotation_matrix, boxsize=data.metadata.boxsize, backend = "subsampled")

        x_range = visualise_region[1] - visualise_region[0]
        y_range = visualise_region[3] - visualise_region[2]
        units = 1.0 / (x_range * y_range)
        # Unfortunately this is required to prevent us from {over,under}flowing
        # the units...
        units.convert_to_units(1.0 / (x_range.units * y_range.units))
        units *= getattr(data.stars, 'masses').units

        mass_map_edge = unyt_array(pixel_grid, units=units)
        mass_map_edge.convert_to_units(msun / pc**2)
        try:
                mass_map_edge[mass_map_edge == 0.] = mass_map_edge[mass_map_edge > 0.].min()
        except:
                mass_map_edge[mass_map_edge == 0.] = 1.e-10


        # mask with circle
        lx, ly = mass_map_edge.shape
        X, Y = np.ogrid[0:lx, 0:ly]
        mask_circle = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 > lx * ly / 4
        mass_map_edge[mask_circle] = np.nan

        # mask with circle
        lx, ly = mass_map_face.shape
        X, Y = np.ogrid[0:lx, 0:ly]
        mask_circle = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 > lx * ly / 4
        mass_map_face[mask_circle] = np.nan

        return mass_map_face.T, mass_map_edge.T, visualise_region, x, y, totalmass


def get_Bigiel2008_data():
        # from morpholopy
        with open("KS_relation_Bigiel2010.txt") as f:
            lines = f.readlines()

        size_array = len(lines) - 49

        sigma_HI = -5 * np.ones(size_array)
        sigma_HI_err = -5 * np.ones(size_array)
        sigma_H2 = -5 * np.ones(size_array)
        sigma_H2_err = -5 * np.ones(size_array)
        sigma_SFR = -5 * np.ones(size_array)
        sigma_SFR_err = -5 * np.ones(size_array)

        for i in range(49, len(lines)):
            k = i - 49
            word1 = lines[i][25:29]
            word2 = lines[i][30:34]
            word3 = lines[i][35:39]
            word4 = lines[i][40:44]
            word5 = lines[i][45:50]
            word6 = lines[i][52:56]
            if word1 != "    ":
                sigma_HI[k] = float(word1)
            if word2 != "    ":
                sigma_HI_err[k] = float(word2)
            if word3 != "    ":
                sigma_H2[k] = float(word3)
            if word4 != "    ":
                sigma_H2_err[k] = float(word4)
            if word5 != "     ":
                sigma_SFR[k] = float(word5)
            if word6 != "    ":
                sigma_SFR_err[k] = float(word6)

        return np.power(10., sigma_HI)/1.36, np.power(10., sigma_H2)/1.36

from math import floor, log10

# Define function for string formatting of scientific notation
def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.
    """
    try:
        if exponent is None:
            exponent = int(floor(log10(abs(num))))
        coeff = round(num / float(10**exponent), decimal_digits)
        if precision is None:
            precision = decimal_digits

        return r"${0:.{2}f}\times10^{{{1:d}}}$".format(coeff, exponent, precision)

    except:
        return "0"
