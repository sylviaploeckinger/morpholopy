import numpy as np
from .html import add_web_section, PlotsInPipeline
from typing import List


def loadGalaxyPlots(
    web,
    output_path: str,
    num_galaxies_to_show: int,
    name_list: List[int],
    min_stellar_mass,
):
    """
    @TODO Create separate .yaml config containing all necessary information about the plots
    """

    PlotsInWeb = PlotsInPipeline()

    title = "Specific angular momentum / Stars"
    caption = "Ratio between the total angular momentum of stars within 30 kpc of "
    caption += "aperture divided by the total mass in stars."
    filename = "momentum_parttype_%i.png" % 4
    id = abs(hash("momentum %i" % 4))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Specific angular momentum / HI+H2 gas"
    caption = "Ratio between the total angular momentum of gas within 30 kpc of "
    caption += "aperture divided by the total mass in gas."
    filename = "momentum_parttype_%i.png" % 0
    id = abs(hash("momentum %i" % 0))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Specific angular momentum"
    id = abs(hash("angular momentum"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    title = "Kappa corotation / Stars"
    caption = (
        "Kappa corotation is defined as the fraction of kinetic energy in a galaxy "
    )
    caption += "that is in ordered rotation. Note that the rotating contribution is calculated "
    caption += "only for prograde rotation."
    filename = "Kappa_co_parttype_%i.png" % 4
    id = abs(hash("kappa co %i" % 4))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Kappa corotation / HI+H2 gas"
    caption = (
        "Kappa corotation is defined as the fraction of kinetic energy in a galaxy "
    )
    caption += "that is in ordered rotation. Note that the rotating contribution is calculated "
    caption += "only for prograde rotation."
    filename = "Kappa_co_parttype_%i.png" % 0
    id = abs(hash("kappa co %i" % 0))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Kappa corotation"
    id = abs(hash("Kappa corotation"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    title = "Axis ratios / Stars"
    caption = "Axial ratios of galaxies more massive than 1e6 Msun in stellar mass. "
    caption += "a, b and c (a >= b >= c) represent the lengths of the primary axes. "
    caption += (
        "Ratios have been calculated following eqs. (1) and (2) from Trayford+2018."
    )
    filename = "Axis_ratios_parttype_%i.png" % 4
    id = abs(hash("galaxy axis %i" % 4))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Axis ratios / HI+H2 gas"
    caption = "Axial ratios of galaxies more massive than 1e6 Msun in stellar mass. "
    caption += "a, b and c (a >= b >= c) represent the lengths of the primary axes. "
    caption += (
        "Ratios have been calculated following eqs. (1) and (2) from Trayford+2018."
    )
    filename = "Axis_ratios_parttype_%i.png" % 0
    id = abs(hash("galaxy axis %i" % 0))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Axis ratios"
    id = abs(hash("axis ratios"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    title = "HI size vs HI mass"
    caption = "HI size as a function of stellar mass"
    filename = "HI_size_mass.png"
    id = abs(hash("HI size vs mass"))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "HI size"
    id = abs(hash("HI size"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for name in name_list:
        title = "Integrated surface densities (" + name + ")"
        caption = "Integrated surface densities of H2+HI gas and star-forming gas for each individual galaxy. "
        caption += "Quantities are calculated summing up all gas (and SFR) within the galaxies' stellar half mass radius."
        filename = "surface_density_H2_" + name + ".png"
        id = abs(hash("surface_density_H2_" + name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for name in name_list:
        title = "Integrated surface densities (" + name + ")"
        caption = "Integrated surface densities of H2 gas and star-forming gas for each individual galaxy. "
        caption += "Quantities are calculated summing up all gas (and SFR) within the galaxies' stellar half mass radius."
        filename = "surface_density_gas_" + name + ".png"
        id = abs(hash("surface_density_gas_" + name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Integrated surface densities"
    id = abs(hash("Integrated Surface density"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for name in name_list:
        title = (
            "Ratio between molecular and atomic neutral surface density vs. neutral gas surface density ("
            + name
            + ")"
        )

        caption = "Combined galaxy sample: central galaxies with stellar masses "
        caption += r"log M$_{\star}$ [M$_{\odot}$] > %.2f." % (
            np.log10(min_stellar_mass)
        )
        caption += "\nAzimuthally averaged measurements, color coded by the mean oxygen abundance of the diffuse component"
        caption += (
            " (i.e. excluding the oxygen in dust) for each bin. Each galaxy is binned with three"
            " different radial bin sizes (see legend). The grey crosses (triangles) show the detections"
            " (upper limits) from Schruba et al. (2011)"
        )
        filename = "Species_transition_" + name + "_Schruba2011.png"
        id = abs(hash("species_transition_Schruba2011" + name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for name in name_list:
        title = (
            "Ratio between molecular and atomic neutral surface density vs. neutral gas surface density ("
            + name
            + ")"
        )
        caption = "Combined galaxy sample: central galaxies with stellar masses "
        caption += r"log M$_{\star}$ [M$_{\odot}$] > %.2f." % (
            np.log10(min_stellar_mass)
        )
        caption += "\nAzimuthally averaged measurements, color coded by the stellar surface density for each bin."
        caption += (
            " Each galaxy is binned with three"
            " different radial bin sizes (see legend). The grey crosses (triangles) show the detections"
            " (upper limits) from Schruba et al. (2011)"
        )
        filename = "Species_transition_" + name + "_Schruba2011_Sigma_star.png"
        id = abs(hash("species_transition_Schruba2011_Sigma_star" + name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for name in name_list:
        title = (
            "Ratio between molecular and atomic neutral surface density vs. neutral gas surface density ("
            + name
            + ")"
        )
        caption = "Combined galaxy sample: central galaxies with stellar masses "
        caption += r"log M$_{\star}$ [M$_{\odot}$] > %.2f." % (
            np.log10(min_stellar_mass)
        )
        caption += "\nSpatially resolved (grid averaged) measurements of galaxies."
        caption += (
            " The grey crosses show the detections"
            " from Bigiel et al. (2008). The simulation data is displayed with points."
        )
        filename = "Species_transition_" + name + "_Bigiel2008_scatter.png"
        id = abs(hash("species_transition_Bigiel2008 scatter" + name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for name in name_list:
        title = (
            "Ratio between molecular and atomic neutral surface density vs. neutral gas surface density ("
            + name
            + ")"
        )
        caption = "Combined galaxy sample: central galaxies with stellar masses "
        caption += r"log M$_{\star}$ [M$_{\odot}$] > %.2f." % (
            np.log10(min_stellar_mass)
        )
        caption += "\nSpatially resolved (grid averaged) measurements of galaxies."
        caption += (
            " The grey crosses show the detections"
            " from Bigiel et al. (2008). The simulation data is displayed with"
            " contours of a kernel-density estimate using Gaussian kernels (scipy.stats.gaussian_kde)."
        )
        filename = "Species_transition_" + name + "_Bigiel2008_contour.png"
        id = abs(hash("species_transition_Bigiel2008 contour" + name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for name in name_list:
        title = "Molecular surface density vs. stellar surface density (" + name + ")"
        caption = "Combined galaxy sample: central galaxies with stellar masses "
        caption += r"log M$_{\star}$ [M$_{\odot}$] > %.2f." % (
            np.log10(min_stellar_mass)
        )
        caption += "\nSpatially resolved (grid averaged) measurements."
        caption += (
            " Simulation data is shown in grey circles and/or a greyscale 2d histogram,"
            " if points saturate."
        )
        filename = "Grid_averaged_" + name + "_Mmol_Mstar.png"
        id = abs(hash("Grid_averaged_" + name + "_Mmol_Mstar"))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for name in name_list:
        title = (
            " Star formation rate surface density vs. molecular surface density ("
            + name
            + ")"
        )
        caption = "Combined galaxy sample: central galaxies with stellar masses "
        caption += r"log M$_{\star}$ [M$_{\odot}$] > %.2f." % (
            np.log10(min_stellar_mass)
        )
        caption += "\nSpatially resolved (grid averaged) measurements."
        caption += (
            " Simulation data is shown in grey circles and/or a greyscale 2d histogram,"
            " if enough points are present."
        )
        filename = "Grid_averaged_" + name + "_SFR_Mmol.png"
        id = abs(hash("Grid_averaged_" + name + "_SFR_Mmol"))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for name in name_list:
        title = (
            " Star formation rate surface density vs. stellar surface density ("
            + name
            + ")"
        )
        caption = "Combined galaxy sample: central galaxies with stellar masses "
        caption += r"log M$_{\star}$ [M$_{\odot}$] > %.2f." % (
            np.log10(min_stellar_mass)
        )
        caption += "\nSpatially resolved (grid averaged) measurements."
        caption += (
            " Simulation data is shown in grey circles and/or a greyscale 2d histogram,"
            " if enough points are present."
        )
        filename = "Grid_averaged_" + name + "_SFR_Mstar.png"
        id = abs(hash("Grid_averaged_" + name + "_SFR_Mstar"))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for name in name_list:
        title = (
            "Molecular-to-neutral surface density vs. neutral gas surface density ("
            + name
            + ")"
        )
        caption = "Combined spatially resolved measurements from N most massive individual galaxies,"
        caption += (
            " coloured by the mean metallicity of the resolved pixel. The surface densities were calculated"
            " using the grid method with a pixel size of 750pc. Coloured solid lines show the median relations"
            " considering only cells with fixed metallicity (as indicated in the legends). The grey solid line"
            " shows the median relation for all pixels."
        )
        filename = "combined_surface_density_ratios_" + name + ".png"
        id = abs(hash("combined_surface_density_ratios_" + name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for name in name_list:
        title = "SFR surface density vs. neutral gas surface density (" + name + ")"
        caption = "Combined spatially resolved measurements from N most massive individual galaxies,"
        caption += (
            " coloured by the mean metallicity of the resolved pixel. The X axis shows the surface density of neutral"
            " gas and the Y axis shows the star formation rate surface density. The surface densities were calculated"
            " using the grid method with a pixel size of 750pc. Coloured lines show in the median relations"
            " considering only cells with fixed metallicity (as indicated in the legends). The grey solid line"
            " shows the median relation for all pixels, whereas the black solid line shows the relation"
            " only for pixels that have SFR surface density >0."
        )
        filename = "combined_surface_density_gas_" + name + ".png"
        id = abs(hash("combined_surface_density_gas_" + name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for name in name_list:
        title = "SFR surface density vs. molecular gas surface density (" + name + ")"
        caption = "Combined spatially resolved measurements from N most massive individual galaxies,"
        caption += (
            " coloured by the mean metallicity of the resolved pixel. The X axis shows the surface density of molecular"
            " gas and the Y axis shows the star formation rate surface density. The surface densities were calculated"
            " using the grid method with a pixel size of 750pc. Coloured lines show in the median relations"
            " considering only cells with fixed metallicity (as indicated in the legends). The grey solid line"
            " shows the median relation for all pixels, whereas the black solid line shows the relation"
            " only for pixels that have SFR surface density >0."
        )
        filename = "combined_surface_density_H2_" + name + ".png"
        id = abs(hash("combined_surface_density_H2_" + name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for name in name_list:
        title = "SFR surface density vs. atomic gas surface density (" + name + ")"
        caption = "Combined spatially resolved measurements from N most massive individual galaxies,"
        caption += (
            " coloured by the mean metallicity of the resolved pixel. The X axis shows the surface density of atomic"
            " gas and the Y axis shows the star formation rate surface density. The surface densities were calculated"
            " using the grid method with a pixel size of 750pc. Coloured lines show in the median relations"
            " considering only cells with fixed metallicity (as indicated in the legends). The grey solid line"
            " shows the median relation for all pixels, whereas the black solid line shows the relation"
            " only for pixels that have SFR surface density >0."
        )
        filename = "combined_surface_density_HI_" + name + ".png"
        id = abs(hash("combined_surface_density_HI_" + name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for name in name_list:
        title = "Neutral KS relation (azimuthally-averaged) (" + name + ")"
        caption = "Combined spatially resolved measurements from N most massive individual galaxies,"
        caption += (
            " The X axis shows the surface density of neutral"
            " gas and the Y axis shows the star formation rate surface density. The surface densities were calculated"
            " using the azimuthally-averaging method with a radial width of 750pc. Coloured lines show in the median relations"
            " considering only cells with fixed metallicity (as indicated in the legends). The grey solid line"
            " shows the median relation for all pixels, whereas the black solid line shows the relation"
            " only for pixels that have SFR surface density >0."
        )
        filename = "radii_combined_surface_density_gas_" + name + ".png"
        id = abs(hash("radii_combined_surface_density_gas_" + name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for name in name_list:
        title = "Molecular KS relation (azimuthally-averaged) (" + name + ")"
        caption = "Combined spatially resolved measurements from N most massive individual galaxies,"
        caption += (
            " The X axis shows the surface density of molecular"
            " gas and the Y axis shows the star formation rate surface density. The surface densities were calculated"
            " using the azimuthally-averaging method with a radial width of 750pc. Coloured lines show in the median relations"
            " considering only cells with fixed metallicity (as indicated in the legends). The grey solid line"
            " shows the median relation for all pixels, whereas the black solid line shows the relation"
            " only for pixels that have SFR surface density >0."
        )
        filename = "radii_combined_surface_density_H2_" + name + ".png"
        id = abs(hash("radii_combined_surface_density_H2_" + name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for name in name_list:
        title = "Atomic KS relation (azimuthally-averaged) (" + name + ")"
        caption = "Combined spatially resolved measurements from N most massive individual galaxies,"
        caption += (
            " The X axis shows the surface density of atomic"
            " gas and the Y axis shows the star formation rate surface density. The surface densities were calculated"
            " using the azimuthally-averaging method with a radial width of 750pc. Coloured lines show in the median relations"
            " considering only cells with fixed metallicity (as indicated in the legends). The grey solid line"
            " shows the median relation for all pixels, whereas the black solid line shows the relation"
            " only for pixels that have SFR surface density >0."
        )
        filename = "radii_combined_surface_density_HI_" + name + ".png"
        id = abs(hash("radii_combined_surface_density_HI_" + name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for name in name_list:
        title = "Depletion time vs. surface density, neutral gas (" + name + ")"
        caption = "Depletion time of neutral gas vs. neutral gas surface density from N most massive individual galaxies,"
        caption += (
            " coloured by the mean metallicity of the resolved pixel. The surface densities"
            " were calculated using a grid with pixel size of 750 pc. Coloured lines show in the median relations"
            " considering only cells with fixed metallicity (as indicated in the legends). The grey solid line"
            " shows the median relation for all pixels, whereas the black solid line shows the relation"
            " only for pixels that have SFR surface density >0."
        )
        filename = "depletion_time_combined_surface_density_gas_" + name + ".png"
        id = abs(hash("depletion_time_combined_surface_density_gas_" + name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for name in name_list:
        title = "Depletion time vs. surface density, molecular gas (" + name + ")"
        caption = "Depletion time of molecular gas vs. molecular gas surface density from N most massive individual galaxies,"
        caption += (
            " coloured by the mean metallicity of the resolved pixel. The surface densities"
            " were calculated using a grid with pixel size of 750 pc. Coloured lines show in the median relations"
            " considering only cells with fixed metallicity (as indicated in the legends). The grey solid line"
            " shows the median relation for all pixels, whereas the black solid line shows the relation"
            " only for pixels that have SFR surface density >0."
        )
        filename = "depletion_time_combined_surface_density_H2_" + name + ".png"
        id = abs(hash("depletion_time_combined_surface_density_H2_" + name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for name in name_list:
        title = "Depletion time vs. surface density, atomic gas (" + name + ")"
        caption = "Depletion time of atomic gas vs. atomic gas surface density from N most massive individual galaxies,"
        caption += (
            " coloured by the mean metallicity of the resolved pixel. The surface densities"
            " were calculated using a grid with pixel size of 750 pc. Coloured lines show in the median relations"
            " considering only cells with fixed metallicity (as indicated in the legends). The grey solid line"
            " shows the median relation for all pixels, whereas the black solid line shows the relation"
            " only for pixels that have SFR surface density >0."
        )
        filename = "depletion_time_combined_surface_density_HI_" + name + ".png"
        id = abs(hash("depletion_time_combined_surface_density_HI_" + name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Combined surface densities"
    id = abs(hash("Surface density"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    # Individual galaxy gallery
    for index in range(num_galaxies_to_show):

        for name in name_list:
            title = "Overview plot (" + name + ")"
            caption = "Surface brightness (gri-composite) and surface density maps"
            id = abs(hash("galaxy overview %i " % (index) + name))

            outfile = "surface_overview_halo%3.3i_" % (index) + name + ".png"
            PlotsInWeb.load_plots(title, caption, outfile, id)

    title = "Gallery"
    caption = " "
    id = abs(hash("galaxy overview gallery %i" % index))
    plots = PlotsInWeb.plots_details
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    # Individual galaxy plots
    for index in range(num_galaxies_to_show):

        for name in name_list:
            title = "Gas component (" + name + ")"
            caption = (
                "Projection of gas within 5 times the galaxy's  stellar "
                "half mass radius. Face-on (left) and edge-on (right)."
            )
            id = abs(hash("galaxy gas %i " % (index) + name))
            outfile = "galaxy_gas_%i_" % (index) + name + ".png"
            # Don't show these plots
            # PlotsInWeb.load_plots(title, caption, outfile, id)

        for name in name_list:
            title = "Stellar component (" + name + ")"
            caption = (
                "Projection of stars within 5 times the galaxy's  stellar "
                "half mass radius. Face-on (left) and edge-on (right)."
            )
            id = abs(hash("stars galaxy %i " % (index) + name))
            outfile = "galaxy_stars_%i_" % (index) + name + ".png"
            # Don't show these plots
            # PlotsInWeb.load_plots(title, caption, outfile, id)

        for name in name_list:
            title = "Star particles (" + name + ")"
            caption = "Projection of stars within 5 times the galaxy's stellar"
            caption += " half mass radius. Face-on (left) and edge-on (right)."
            id = abs(hash("galaxy stars parts %i " % (index) + name))
            outfile = "galaxy_sparts_%i_" % (index) + name + ".png"
            # Don't show these plots
            # PlotsInWeb.load_plots(title, caption, outfile, id)

        for name in name_list:
            title = "Gas particles (" + name + ")"
            caption = "Projection of gas within 5 times the galaxy's stellar"
            caption += " half mass radius. Face-on (left) and edge-on (right)."
            id = abs(hash("galaxy gas parts %i " % (index) + name))
            outfile = "galaxy_parts_%i_" % (index) + name + ".png"
            # Don't show these plots
            # PlotsInWeb.load_plots(title, caption, outfile, id)

        for name in name_list:
            for filtname in ["u", "r", "K"]:
                title = f"Rest frame {filtname}-band light "
                caption = "Projection %s-band emission within 5 times " % (filtname)
                caption += "the galaxy's stellar half mass radius. Face-on (left) and edge-on (right)."
                id = abs(hash("%s-band emission galaxy %i" % (filtname, index) + name))
                outfile = "galaxy_%s_map_%i_" % (filtname, index) + name + ".png"
                # Don't show these plots
                # PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "KS relation (data: H2+HI mass, method: grid)"
        id = abs(hash("galaxy KS relation H2+HI grid %i" % (index)))
        outfile = "KS_relation_best_grid_%i.png" % (index)
        caption = "KS relation. Surface densities were calculated using a grid with pixel size of 750 pc."
        caption += " Each blue dot shows the total SFR and H2 mass in the pixel divided by the pixel area."
        caption += " Black solid line indicates the median relation and shaded area the 84-16th percentiles."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "KS relation (data: H2 mass, method: grid)"
        id = abs(hash("galaxy KS relation H2 grid %i" % (index)))
        outfile = "KS_molecular_relation_grid_%i.png" % (index)
        caption = "KS relation. Surface densities were calculated using a grid with pixel size of 750 pc."
        caption += " Each blue dot shows the total SFR and H2 mass in the pixel divided by the pixel area."
        caption += " Black solid line indicates the median relation and shaded area the 84-16th percentiles."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "KS relation (data: HI mass, method: grid)"
        id = abs(hash("galaxy KS relation HI grid %i" % (index)))
        outfile = "KS_atomic_relation_grid_%i.png" % (index)
        caption = "KS relation. Surface densities were calculated using a grid with pixel size of 750 pc."
        caption += " Each blue dot shows the total SFR and H2 mass in the pixel divided by the pixel area."
        caption += " Black solid line indicates the median relation and shaded area the 84-16th percentiles."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "KS relation (data: H2+HI mass, method: Azimuthal average)"
        id = abs(hash("galaxy KS relation H2+HI radii %i" % (index)))
        outfile = "KS_relation_best_radii_%i.png" % (index)
        caption = "KS relation. Surface densities were calculated by azimuthally averaging radial concentric shells"
        caption += " of 750 pc of width. The shells are centered in the minimum of the dark matter potential."
        caption += " Each blue dot shows the total SFR and H2 mass in the shell divided by the shell area."
        caption += " Black solid line indicates the median relation and shaded area the 84-16th percentiles."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "KS relation (data: H2 mass, method: Azimuthal average)"
        id = abs(hash("galaxy KS relation H2 radii %i" % (index)))
        outfile = "KS_molecular_relation_radii_%i.png" % (index)
        caption = "KS relation. Surface densities were calculated by azimuthally averaging radial concentric shells"
        caption += " of 750 pc of width. The shells are centered in the minimum of the dark matter potential."
        caption += " Each blue dot shows the total SFR and H2 mass in the shell divided by the shell area."
        caption += " Black solid line indicates the median relation and shaded area the 84-16th percentiles."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "KS relation (data: HI mass, method: Azimuthal average)"
        id = abs(hash("galaxy KS relation HI radii %i" % (index)))
        outfile = "KS_atomic_relation_radii_%i.png" % (index)
        caption = "KS relation. Surface densities were calculated by azimuthally averaging radial concentric shells"
        caption += " of 750 pc of width. The shells are centered in the minimum of the dark matter potential."
        caption += " Each blue dot shows the total SFR and H2 mass in the shell divided by the shell area."
        caption += " Black solid line indicates the median relation and shaded area the 84-16th percentiles."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "Depletion time (data: H2+HI mass, method: grid)"
        id = abs(hash("galaxy depletion H2+HI grid %i" % (index)))
        outfile = "gas_depletion_timescale_best_grid_%i.png" % (index)
        caption = "Gas depletion times. The surface densities were calculated using a grid with pixel size of 750 pc."
        caption += " Black solid line indicates the median relation, shaded area the 84-16th percentiles, "
        caption += "and the observational data-points correspond to Bigiel et al. (2008, 2010) inner, same as in KS relation (H2+HI mass) figure."
        # Don't shown dipletion time plots
        # PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "Depletion time (data: H2 mass, method: grid)"
        id = abs(hash("galaxy depletion H2 grid %i" % (index)))
        outfile = "molecular_gas_depletion_timescale_grid_%i.png" % (index)
        caption = "Gas depletion times. The surface densities were calculated using a grid with pixel size of 750 pc."
        caption += " Black solid line indicates the median relation, shaded area the 84-16th percentiles, "
        caption += "and the observational data-points correspond to Bigiel et al. (2008) inner, same as in KS relation (H2 mass) figure."
        # Don't shown dipletion time plots
        # PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "Depletion time (data: H2 mass, method: Azimuthal average)"
        id = abs(hash("galaxy depletion H2 radii %i" % (index)))
        outfile = "molecular_gas_depletion_timescale_radii_%i.png" % (index)
        caption = "Gas depletion times. The surface densities were calculated by azimuthally averaging radial concentric shells"
        caption += " of 800 pc of width. The shells were centered in the minimum of the dark matter potential."
        caption += " Black solid line indicates the median relation, shaded area the 84-16th percentiles, "
        caption += "and the observational data-points correspond to Bigiel et al. (2008) inner, same as in KS relation (H2 mass) figure."
        # Don't shown dipletion time plots
        # PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "Depletion time (data: H2+HI mass, method: Azimuthal average)"
        id = abs(hash("galaxy depletion H2+HI radii %i" % (index)))
        outfile = "gas_depletion_timescale_best_radii_%i.png" % (index)
        caption = "Gas depletion times. The surface densities were calculated by azimuthally averaging radial concentric shells"
        caption += " of 800 pc of width. The shells were centered in the minimum of the dark matter potential."
        caption += " Black solid line indicates the median relation, shaded area the 84-16th percentiles, "
        caption += "and the observational data-points correspond to Bigiel et al. (2008, 2010) inner, same as in KS relation (H2+HI mass) figure."
        # Don't shown dipletion time plots
        # PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "Surface density ratios (method: grid)"
        id = abs(hash("density ratio H2+HI grid %i" % (index)))
        outfile = "Surface_density_ratio_grid_%i.png" % (index)
        caption = "Surface density ratios. The y-axis shows the ratio between surface densities calculated using a grid"
        caption += " with pixel size of 750 pc. Red dashed line corresponds to Krumholz+ (2009) semi-analytic model, the"
        caption += " black solid line indicates the median relation and the shaded area the 84-16th percentiles, "
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "Surface density ratios (method: Azimuthal average)"
        id = abs(hash("density ratio H2+HI radii %i" % (index)))
        outfile = "Surface_density_ratio_radii_%i.png" % (index)
        caption = "Surface density ratios. The y-axis shows the ratio between surface densities calculated"
        caption += " by azimuthally averaging radial concentric shells of 800 pc of width. The shells were centered"
        caption += " in the minimum of the dark matter potential."
        caption += " The red dashed and dotted lines correspond to Krumholz+ (2009) semi-analytic model,"
        caption += " black solid line indicates the median relation and the shaded area the 84-16th percentiles, "
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "%i Galaxy " % (index + 1)
        caption = " "
        for name in name_list:
            data = np.loadtxt(f"{output_path}/galaxy_data_" + name + ".txt")

            try:
                sfr_galaxy = data[:, 0]
                mass_galaxy = data[:, 1]
                gas_mass_galaxy = data[:, 2]
                mass_halo = data[:, 3]
                galaxy_metallicity_gas_sfr = data[:, 4]
                galaxy_metallicity_gas = data[:, 5]

                caption += "<strong>Simulation: " + name + "</strong>. Galaxy details: "
                caption += r"$\log_{10}$ M$_{200}$/M$_{\odot} = $%0.2f," % (
                    mass_halo[index]
                )
                caption += " SFR = %0.1f M$_{\odot}$/yr," % (sfr_galaxy[index])
                caption += " Z$_{\mathrm{SFR}>0}$ = %0.3f," % (
                    galaxy_metallicity_gas_sfr[index]
                )
                caption += " Z$_{\mathrm{gas}}$ = %0.3f," % (
                    galaxy_metallicity_gas[index]
                )
                caption += " $\log_{10}$ M$_{*}$/M$_{\odot} = $%0.2f" % (
                    mass_galaxy[index]
                )
                caption += (
                    ' $\&$ $\log_{10}$ M$_{\mathrm{gas}}$/M$_{\odot} = $%0.2f.</p><p style="font-size:18px;">'
                    % (gas_mass_galaxy[index])
                )

            # If only one galaxy is available, the data array is one dimensional
            except IndexError:

                sfr_galaxy = data[0]
                mass_galaxy = data[1]
                gas_mass_galaxy = data[2]
                mass_halo = data[3]
                galaxy_metallicity_gas_sfr = data[4]
                galaxy_metallicity_gas = data[5]
                caption += "<strong>Simulation: " + name + "</strong>. Galaxy details: "
                caption += r"$\log_{10}$ M$_{200}$/M$_{\odot} = $%0.2f," % (mass_halo)
                caption += " SFR = %0.1f M$_{\odot}$/yr," % (sfr_galaxy)
                caption += " Z$_{\mathrm{SFR}>0}$ = %0.3f," % (
                    galaxy_metallicity_gas_sfr
                )
                caption += " Z$_{\mathrm{gas}}$ = %0.3f," % (galaxy_metallicity_gas)
                caption += " $\log_{10}$ M$_{*}$/M$_{\odot} = $%0.2f" % (mass_galaxy)
                caption += (
                    ' $\&$ $\log_{10}$ M$_{\mathrm{gas}}$/M$_{\odot} = $%0.2f.</p><p style="font-size:18px;">'
                    % (gas_mass_galaxy)
                )

        id = abs(hash("galaxy and ks relation %i" % index))
        plots = PlotsInWeb.plots_details
        add_web_section(web, title, caption, id, plots)
        PlotsInWeb.reset_plots_list()
