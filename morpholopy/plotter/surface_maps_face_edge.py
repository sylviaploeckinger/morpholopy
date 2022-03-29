from velociraptor import load as load_catalogue
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
from matplotlib import gridspec
from unyt import pc, kpc, msun, yr, Myr

from plotter.helpers import (
    get_stars_surface_density_map,
    get_gas_surface_density_map,
    calculate_scaleheight_fit,
    get_stars_surface_brightness_map,
    get_radial_profile,
    get_Schruba_data,
    get_Schruba_upperlimits,
    sci_notation,
)


cmRdBu = plt.cm.get_cmap("RdYlBu_r")
Sigma_min = 0.0
Sigma_max = 99.0
H2_over_HI_min = 10.0 ** (-1.9)
H2_over_HI_max = 10.0 ** (+2.4)

twelve_plus_logOH_solar = 8.69
twelve_plus_logOH_min = twelve_plus_logOH_solar - 0.5
twelve_plus_logOH_max = twelve_plus_logOH_solar + 0.5

Zsun = 0.0134
twelve_plus_log_OH_solar = 8.69

#########################
# nr of map plots
#########################

nr_text_boxes = 1
nr_map_plots = 6
nr_scatter_plots = 1
nr_total_plots = nr_text_boxes + nr_map_plots + nr_scatter_plots

#################################
# Miscellaneous
#################################

npix = int(512 / 2)

vmin = -0.98
vmax = 2.98

r_img_kpc = 30.0 * kpc
lbar_kpc = 15.0 * kpc
ypos_bar = 20.0 * kpc

size = 2.0 * r_img_kpc

radialbinsizes = [0.22, 0.8, 1.8]  # in kpc
pixsize_kpc = r_img_kpc.value / npix

bingrid = 0.8 * kpc
npix_coarse = int(2.0 * r_img_kpc / bingrid)

cmap = matplotlib.cm.get_cmap("Greens")
usergreen = cmap(0.7)
cmap = matplotlib.cm.get_cmap("Blues")
userblue = cmap(0.7)


def surface_densities_overview(
    sim_name,
    directory,
    snapshot,
    catalogue_file,
    output_path,
    nhalos,
    halo_min_stellar_mass,
    halo_ids_sample,
):

    filename_text_f1 = output_path + "/Azimuth_averaged_" + sim_name + ".dat"
    f1 = open(filename_text_f1, "w")
    f1.write("# column 0: Simulations name\n")
    f1.write("# column 1: Snapshot number\n")
    f1.write("# column 2: Redshift\n")
    f1.write("# column 3: Halo id\n")
    f1.write("# column 4: Radial bin size [kpc]\n")
    f1.write("# column 5: Sigma (HI+H2) [Msol pc-2]\n")
    f1.write("# column 6: Sigma H2 / Sigma HI\n")
    f1.write("# column 7: Metallicity [12 + log10(O/H)] (diffuse)\n")
    f1.write("# column 8: Stellar mass surface density [Msol pc-2]\n")
    f1.close()

    filename_text_f2 = output_path + "/Grid_averaged_" + sim_name + ".dat"
    f2 = open(filename_text_f2, "w")
    f2.write("# column 0: Simulations name\n")
    f2.write("# column 1: Snapshot number\n")
    f2.write("# column 2: Redshift\n")
    f2.write("# column 3: Halo id\n")
    f2.write("# column 4: Grid pixel size [kpc]\n")
    f2.write("# column 5: log10 Star formation rate surface density [Msol yr-1 pc-2]\n")
    f2.write("# column 6: log10 H2 surface density [Msol pc-2]\n")
    f2.write("# column 7: log10 HI surface density [Msol pc-2]\n")
    f2.write("# column 8: log10 Stellar mass surface density [Msol pc-2]\n")
    f2.close()

    snapshot_filename = directory + "/" + snapshot

    catalogue_filename = glob.glob(directory + catalogue_file + "*")[0]
    catalogue = load_catalogue(catalogue_filename)
    catalogue.masses.mass_star_30kpc.convert_to_units("Msun")

    catalogue.masses.mass_200crit.convert_to_units("Msun")
    catalogue.masses.mass_gas_30kpc.convert_to_units("Msun")
    catalogue.gas_hydrogen_species_masses.HI_mass_30_kpc.convert_to_units("Msun")
    catalogue.gas_hydrogen_species_masses.H2_mass_30_kpc.convert_to_units("Msun")
    catalogue.apertures.mass_star_100_kpc.convert_to_units("Msun")

    Zdiffuse = catalogue.cold_dense_gas_properties.cold_dense_diffuse_metal_mass_100_kpc
    ColdGas = catalogue.cold_dense_gas_properties.cold_dense_gas_mass_100_kpc
    twelve_plus_logOH = np.log10(Zdiffuse / (Zsun * ColdGas)) + twelve_plus_log_OH_solar

    snapnum = "".join([s for s in snapshot if s.isdigit()])

    for ihalo, halo_id in enumerate(halo_ids_sample):

        fig = plt.figure(figsize=(15.0, 3.5))
        fig.subplots_adjust(left=0.01, right=0.95, top=0.85, bottom=0.12)
        gs = gridspec.GridSpec(
            2, nr_total_plots, wspace=0.0, hspace=0.15, height_ratios=[0.05, 1.0]
        )

        # General information
        text = sim_name + ", z = %.1f" % (catalogue.redshift) + "\n\n"
        text += r"${\bf" + "VR\ halo\ id:\ \ \ %3.3i" % (halo_id) + r"}$" + "\n"
        text += (
            r"M$_{\mathrm{200,crit}}$ = "
            + sci_notation(catalogue.masses.mass_200crit[halo_id].value)
            + r" M$_{\odot}$"
            + "\n"
        )
        text += (
            r"M$_{\mathrm{*,30kpc}}$ = "
            + sci_notation(catalogue.masses.mass_star_30kpc[halo_id].value)
            + r" M$_{\odot}$"
            + "\n"
        )
        text += (
            r"M$_{\mathrm{gas,30kpc}}$ = "
            + sci_notation(catalogue.masses.mass_gas_30kpc[halo_id].value)
            + r" M$_{\odot}$"
            + "\n"
        )
        text += (
            r"M$_{\mathrm{HI,30kpc}}$ = "
            + sci_notation(
                catalogue.gas_hydrogen_species_masses.HI_mass_30_kpc[halo_id].value
            )
            + r" M$_{\odot}$"
            + "\n"
        )
        text += (
            r"M$_{\mathrm{H2,30kpc}}$ = "
            + sci_notation(
                catalogue.gas_hydrogen_species_masses.H2_mass_30_kpc[halo_id].value
            )
            + r" M$_{\odot}$"
            + "\n"
        )

        # catalogue.apertures.sfr_gas_100_kpc is in units of '97.78*Msun/yr' and somehow the unit conversion does not work
        # if I use catalogue.apertures.sfr_gas_100_kpc.convert_to_units('Msun/yr') it multiplies by 97.78 instead of dividing it by 97.78
        # as a workaround do the unit conversion manually (checked that it agrees with sSFR from the pipeline)
        sSFR = (
            catalogue.apertures.sfr_gas_100_kpc[halo_id].value
            / 97.78
            / catalogue.apertures.mass_star_100_kpc[halo_id].value
            * 1.0e9
        )
        print(halo_id, sSFR)
        if np.isfinite(sSFR):
            text += r"sSFR$_{\mathrm{100}}$ = " + sci_notation(sSFR) + r" Gyr$^{-1}$"

        ax = plt.subplot(gs[nr_total_plots])
        ax.set_aspect("equal")
        ax.tick_params(labelleft=False, labelbottom=False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(
            0.05, 0.95, text, ha="left", va="top", transform=ax.transAxes, fontsize=8
        )
        ax.text(
            0.05,
            1.20,
            "%i Galaxy" % (ihalo + 1),
            ha="left",
            va="bottom",
            transform=ax.transAxes,
            fontsize=14,
        )

        # Stars gri face-on
        ax = plt.subplot(gs[nr_total_plots + 1])
        ax.set_title("Stars (gri) - face")
        (
            mass_map_face,
            mass_map_edge,
            visualise_region,
            x,
            y,
            totalmass,
            H_kpc_gri,
        ) = get_stars_surface_brightness_map(
            catalogue, halo_id, snapshot_filename, size, npix, r_img_kpc
        )
        mass_map_face_plot = mass_map_face
        mass_map_edge_plot = mass_map_edge
        ax.tick_params(labelleft=False, labelbottom=False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_axis_off()
        im = ax.imshow(mass_map_face_plot, extent=visualise_region)
        circle = plt.Circle(
            (x, y),
            (0.99 * r_img_kpc.value) / 1000.0,
            color="black",
            fill=False,
            linewidth=2,
        )
        ax.add_artist(circle)
        ax.plot(
            [x - lbar_kpc / 2.0, x + lbar_kpc / 2.0],
            [y + ypos_bar, y + ypos_bar],
            color="white",
            linewidth=2,
            linestyle="solid",
        )
        ax.text(
            x,
            y + ypos_bar,
            "%i kpc" % (int(lbar_kpc.value)),
            color="white",
            verticalalignment="bottom",
            horizontalalignment="center",
        )
        # Stars gri edge-on
        ax = plt.subplot(gs[nr_total_plots + 2])
        ax.set_title("Stars (gri) - edge")
        ax.tick_params(labelleft=False, labelbottom=False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_axis_off()
        im = ax.imshow(mass_map_edge_plot, extent=visualise_region)
        circle = plt.Circle(
            (x, y),
            (0.99 * r_img_kpc.value) / 1000.0,
            color="black",
            fill=False,
            linewidth=2,
        )
        ax.add_artist(circle)
        ax.text(
            0.5,
            0.2,
            r"H$_{r}$ = %.2f kpc" % (H_kpc_gri[1]),
            ha="center",
            va="top",
            color="white",
            transform=ax.transAxes,
            fontsize=8,
        )

        # HI face-on
        ax = plt.subplot(gs[nr_total_plots + 3])
        ax.set_title("Gas (HI) - face")
        (
            mass_map_face,
            mass_map_edge,
            visualise_region,
            x,
            y,
            totalmass,
        ) = get_gas_surface_density_map(
            catalogue, halo_id, "HI", snapshot_filename, size, npix
        )
        mass_map_face.convert_to_units("Msun / pc**2")
        mass_map_edge.convert_to_units("Msun / pc**2")
        mass_map_face_plot_HI = mass_map_face  # save for H2 ratio plot
        totalmass_H2 = totalmass  # save for H2 ratio plot
        mass_map_face_plot = np.log10(mass_map_face.value)
        mass_map_edge_plot = np.log10(mass_map_edge.value)
        ax.tick_params(labelleft=False, labelbottom=False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_axis_off()
        colormap = "Blues"
        cmap_loc = matplotlib.cm.get_cmap(colormap)
        ccolor = cmap_loc(0.5)
        im = ax.imshow(
            mass_map_face_plot,
            cmap=colormap,
            extent=visualise_region,
            vmin=vmin,
            vmax=vmax,
        )
        circle = plt.Circle(
            (x, y),
            (0.99 * r_img_kpc.value) / 1000.0,
            color=ccolor,
            fill=False,
            linewidth=2,
        )
        ax.add_artist(circle)

        # HI edge-on
        ax = plt.subplot(gs[nr_total_plots + 4])
        ax.set_title("Gas (HI) - edge")
        ax.tick_params(labelleft=False, labelbottom=False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_axis_off()
        im = ax.imshow(
            mass_map_edge_plot,
            cmap=colormap,
            extent=visualise_region,
            vmin=vmin,
            vmax=vmax,
        )
        circle = plt.Circle(
            (x, y),
            (0.99 * r_img_kpc.value) / 1000.0,
            color=ccolor,
            fill=False,
            linewidth=2,
        )
        ax.add_artist(circle)
        try:
            H_kpc = calculate_scaleheight_fit(mass_map_edge.value)
            ax.text(
                0.5,
                0.2,
                r"H$_{\mathrm{HI}}$ = %.2f kpc" % (H_kpc),
                ha="center",
                va="top",
                color="black",
                transform=ax.transAxes,
                fontsize=8,
            )
        except:
            H_kpc = -1.0

        # HI colorbar
        cax = plt.subplot(gs[3:5])
        cb = fig.colorbar(im, cax=cax, orientation="horizontal")
        cb.set_label("log $\Sigma_{\mathrm{%s}}$ [M$_{\odot}$ pc$^{-2}$]" % ("HI"))
        cb.ax.xaxis.set_ticks_position("top")
        cb.ax.xaxis.set_label_position("top")
        cb.set_ticks(np.arange(round(vmin), vmax + 0.5, 0.5))

        # H2 face-on
        ax = plt.subplot(gs[nr_total_plots + 5])
        ax.set_title("Gas (H2) - face")
        (
            mass_map_face,
            mass_map_edge,
            visualise_region,
            x,
            y,
            totalmass,
        ) = get_gas_surface_density_map(
            catalogue, halo_id, "H2", snapshot_filename, size, npix
        )
        mass_map_face.convert_to_units("Msun / pc**2")
        mass_map_edge.convert_to_units("Msun / pc**2")
        mass_map_face_plot_H2 = mass_map_face  # save for H2 ratio plot
        mass_map_face_plot = np.log10(mass_map_face.value)
        mass_map_edge_plot = np.log10(mass_map_edge.value)
        ax.tick_params(labelleft=False, labelbottom=False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_axis_off()
        colormap = "Greens"
        cmap_loc = matplotlib.cm.get_cmap(colormap)
        ccolor = cmap_loc(0.5)
        im = ax.imshow(
            mass_map_face_plot,
            cmap=colormap,
            extent=visualise_region,
            vmin=vmin,
            vmax=vmax,
        )
        circle = plt.Circle(
            (x, y),
            (0.99 * r_img_kpc.value) / 1000.0,
            color=ccolor,
            fill=False,
            linewidth=2,
        )
        ax.add_artist(circle)

        # H2 edge-on
        ax = plt.subplot(gs[nr_total_plots + 6])
        ax.set_title("Gas (H2) - edge")
        ax.tick_params(labelleft=False, labelbottom=False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_axis_off()
        im = ax.imshow(
            mass_map_edge_plot,
            cmap=colormap,
            extent=visualise_region,
            vmin=vmin,
            vmax=vmax,
        )
        circle = plt.Circle(
            (x, y),
            (0.99 * r_img_kpc.value) / 1000.0,
            color=ccolor,
            fill=False,
            linewidth=2,
        )
        ax.add_artist(circle)
        try:
            H_kpc = calculate_scaleheight_fit(mass_map_edge.value)
            ax.text(
                0.5,
                0.2,
                r"H$_{\mathrm{H2}}$ = %.2f kpc" % (H_kpc),
                ha="center",
                va="top",
                color="black",
                transform=ax.transAxes,
                fontsize=8,
            )
        except:
            H_kpc = -1.0

        # H2 colorbar
        cax = plt.subplot(gs[5:7])
        cb = fig.colorbar(im, cax=cax, orientation="horizontal")
        cb.set_label("log $\Sigma_{\mathrm{%s}}$ [M$_{\odot}$ pc$^{-2}$]" % ("H2"))
        cb.ax.xaxis.set_ticks_position("top")
        cb.ax.xaxis.set_label_position("top")
        cb.set_ticks(np.arange(round(vmin), vmax + 0.5, 0.5))

        # Sigma H2 / Sigma HI vs. Sigma HI+H2
        ax = plt.subplot(gs[nr_total_plots + 7])
        ax.set_title("HI - H2 transition")
        ax.set_aspect(
            (Sigma_max - Sigma_min)
            / (np.log10(H2_over_HI_max) - np.log10(H2_over_HI_min))
        )
        ax.set_xlim(Sigma_min, Sigma_max)
        ax.set_ylim(H2_over_HI_min, H2_over_HI_max)
        ax.set_yscale("log")
        ax.set_yticks([0.1, 1.0, 10.0, 100.0])
        ax.set_yticklabels(["0.1", "1", "10", "100"])
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
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

        if totalmass_H2 > 1.0e7 * msun:
            markers = ["o", "s", "D"]
            f1 = open(filename_text_f1, "a")
            for irad in range(len(radialbinsizes)):
                radialbin_kpc = radialbinsizes[irad]
                # datapoint position
                _, Sigma_HIH2_sub = get_radial_profile(
                    mass_map_face_plot_HI + mass_map_face_plot_H2,
                    radialbin_kpc,
                    pixsize_kpc,
                    r_img_kpc,
                )
                _, Sigma_HI_sub = get_radial_profile(
                    mass_map_face_plot_HI, radialbin_kpc, pixsize_kpc, r_img_kpc
                )
                _, Sigma_H2_sub = get_radial_profile(
                    mass_map_face_plot_H2, radialbin_kpc, pixsize_kpc, r_img_kpc
                )
                # datapoint color (metallicity)
                (
                    mass_map_face_plot_H,
                    _,
                    _,
                    _,
                    _,
                    totalmass_H,
                ) = get_gas_surface_density_map(
                    catalogue, halo_id, "hydrogen", snapshot_filename, size, npix
                )
                (
                    mass_map_face_plot_O_diffuse,
                    _,
                    _,
                    _,
                    _,
                    totalmass_O_diffuse,
                ) = get_gas_surface_density_map(
                    catalogue, halo_id, "diffuseoxygen", snapshot_filename, size, npix
                )
                (
                    mass_map_face_plot_star,
                    _,
                    _,
                    _,
                    _,
                    totalmass_star,
                ) = get_stars_surface_density_map(
                    catalogue, halo_id, "stars", snapshot_filename, size, npix
                )

                _, Sigma_H_sub = get_radial_profile(
                    mass_map_face_plot_H, radialbin_kpc, pixsize_kpc, r_img_kpc
                )
                _, Sigma_O_diffuse_sub = get_radial_profile(
                    mass_map_face_plot_O_diffuse, radialbin_kpc, pixsize_kpc, r_img_kpc
                )
                _, Sigma_star = get_radial_profile(
                    mass_map_face_plot_star, radialbin_kpc, pixsize_kpc, r_img_kpc
                )

                sc = ax.scatter(
                    Sigma_HIH2_sub,
                    Sigma_H2_sub / Sigma_HI_sub,
                    c=12.0 + np.log10(Sigma_O_diffuse_sub / Sigma_H_sub / 16.0),
                    vmin=twelve_plus_logOH_min,
                    vmax=twelve_plus_logOH_max,
                    cmap=cmRdBu,
                    edgecolors="black",
                    label="Azim.avg. %.2f kpc" % (radialbin_kpc),
                    marker=markers[irad],
                )

                for i in range(len(Sigma_HIH2_sub)):
                    if Sigma_H2_sub[i] > 0.0 and Sigma_HI_sub[i] > 0.0:
                        f1.write(
                            sim_name
                            + "\t"
                            + snapnum
                            + "\t%.2f\t%i\t%.2f\t%.2e\t%.2e\t%.2e\t%.2e\n"
                            % (
                                catalogue.redshift,
                                halo_id,
                                radialbin_kpc,
                                Sigma_HIH2_sub[i],
                                Sigma_H2_sub[i] / Sigma_HI_sub[i],
                                12.0
                                + np.log10(
                                    Sigma_O_diffuse_sub[i] / Sigma_H_sub[i] / 16.0
                                ),
                                Sigma_star[i],
                            )
                        )
            f1.close()

        ax.legend(bbox_to_anchor=(-0.1, -0.1), ncol=3, loc="upper right", fontsize=8)

        # Metallicity colorbar
        if totalmass_H2 > 1.0e7 * msun:
            cax = plt.subplot(gs[7])
            cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
            cb.set_label(r"12 + log$_{\mathrm{10}}$(O/H)$_{\mathrm{diffuse}}$")
            cb.ax.xaxis.set_ticks_position("top")
            cb.ax.xaxis.set_label_position("top")
            cb.set_ticks(
                [twelve_plus_logOH_min, twelve_plus_logOH_solar, twelve_plus_logOH_max]
            )

        ax.text(
            0.95,
            0.05,
            "12 + log$_{\mathrm{10}}$(O/H) = %.2f" % (twelve_plus_logOH[halo_id]),
            va="bottom",
            ha="right",
            transform=ax.transAxes,
            fontsize=8,
        )

        fig.savefig(
            output_path + "/surface_overview_halo%3.3i_" % (ihalo) + sim_name + ".png",
            dpi=150,
        )
        plt.close()

        # produce the grid averaged properties

        mass_map_face_plot_sfr, _, _, _, _, _ = get_gas_surface_density_map(
            catalogue, halo_id, "sfr", snapshot_filename, size, npix_coarse
        )
        mass_map_face_plot_H2, _, _, _, _, _ = get_gas_surface_density_map(
            catalogue, halo_id, "H2", snapshot_filename, size, npix_coarse
        )
        mass_map_face_plot_HI, _, _, _, _, _ = get_gas_surface_density_map(
            catalogue, halo_id, "HI", snapshot_filename, size, npix_coarse
        )
        mass_map_face_plot_star, _, _, _, _, _ = get_stars_surface_density_map(
            catalogue, halo_id, "Stars", snapshot_filename, size, npix_coarse
        )

        mass_map_face_plot_sfr.convert_to_units(msun / yr / pc ** 2)
        mass_map_face_plot_H2.convert_to_units(msun / pc ** 2)
        mass_map_face_plot_HI.convert_to_units(msun / pc ** 2)
        mass_map_face_plot_star.convert_to_units(msun / pc ** 2)

        f2 = open(filename_text_f2, "a")

        for i in range(len(mass_map_face_plot_sfr.ravel())):
            logsfr = np.log10(mass_map_face_plot_sfr.ravel()[i])
            logH2 = np.log10(mass_map_face_plot_H2.ravel()[i])
            logHI = np.log10(mass_map_face_plot_HI.ravel()[i])
            logstar = np.log10(mass_map_face_plot_star.ravel()[i])
            if np.isfinite(logsfr) and np.isfinite(logH2):
                f2.write(
                    sim_name
                    + "\t"
                    + snapnum
                    + "\t%.2f\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n"
                    % (
                        catalogue.redshift,
                        halo_id,
                        bingrid.value,
                        logsfr,
                        logH2,
                        logHI,
                        logstar,
                    )
                )

        f2.close()

    return
