import numpy as np
from .html import add_web_section, PlotsInPipeline
from typing import List


def loadGalaxyPlots(
    web, output_path: str, num_galaxies_to_show: int, name_list: List[int]
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
            "Molecular-to-neutral surface density vs. neutral gas surface density ("
            + name
            + ")"
        )
        caption = "Combined spatially resolved measurements from N most massive individual galaxies,"
        caption += (
            " coloured by the mean metallicity of the resolved pixel. The surface densities were calculated"
            " using the grid method with a pixel size of 250pc. Coloured solid lines show the median relations"
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
            " using the grid method with a pixel size of 250pc. Coloured lines show in the median relations"
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
            " using the grid method with a pixel size of 250pc. Coloured lines show in the median relations"
            " considering only cells with fixed metallicity (as indicated in the legends). The grey solid line"
            " shows the median relation for all pixels, whereas the black solid line shows the relation"
            " only for pixels that have SFR surface density >0."
        )
        filename = "combined_surface_density_H2_" + name + ".png"
        id = abs(hash("combined_surface_density_H2_" + name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for name in name_list:
        title = "Depletion time vs. surface density, neutral gas (" + name + ")"
        caption = "Depletion time of neutral gas vs. neutral gas surface density from N most massive individual galaxies,"
        caption += (
            " coloured by the mean metallicity of the resolved pixel. The surface densities"
            " were calculated using a grid with pixel size of 250 pc. Coloured lines show in the median relations"
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
            " were calculated using a grid with pixel size of 250 pc. Coloured lines show in the median relations"
            " considering only cells with fixed metallicity (as indicated in the legends). The grey solid line"
            " shows the median relation for all pixels, whereas the black solid line shows the relation"
            " only for pixels that have SFR surface density >0."
        )
        filename = "depletion_time_combined_surface_density_H2_" + name + ".png"
        id = abs(hash("depletion_time_combined_surface_density_H2_" + name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Combined surface densities"
    id = abs(hash("Surface density"))
    plots = PlotsInWeb.plots_details
    caption = " "
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
            PlotsInWeb.load_plots(title, caption, outfile, id)

        for name in name_list:
            title = "Stellar component (" + name + ")"
            caption = (
                "Projection of stars within 5 times the galaxy's  stellar "
                "half mass radius. Face-on (left) and edge-on (right)."
            )
            id = abs(hash("stars galaxy %i " % (index) + name))
            outfile = "galaxy_stars_%i_" % (index) + name + ".png"
            PlotsInWeb.load_plots(title, caption, outfile, id)

        for name in name_list:
            title = "Star particles (" + name + ")"
            caption = "Projection of stars within 5 times the galaxy's stellar"
            caption += " half mass radius. Face-on (left) and edge-on (right)."
            id = abs(hash("galaxy stars parts %i " % (index) + name))
            outfile = "galaxy_sparts_%i_" % (index) + name + ".png"
            PlotsInWeb.load_plots(title, caption, outfile, id)

        for name in name_list:
            title = "Gas particles (" + name + ")"
            caption = "Projection of gas within 5 times the galaxy's stellar"
            caption += " half mass radius. Face-on (left) and edge-on (right)."
            id = abs(hash("galaxy gas parts %i " % (index) + name))
            outfile = "galaxy_parts_%i_" % (index) + name + ".png"
            PlotsInWeb.load_plots(title, caption, outfile, id)

        for name in name_list:
            for filtname in ["u", "r", "K"]:
                title = f"Rest frame {filtname}-band light "
                caption = "Projection %s-band emission within 5 times " % (filtname)
                caption += "the galaxy's stellar half mass radius. Face-on (left) and edge-on (right)."
                id = abs(hash("%s-band emission galaxy %i" % (filtname, index) + name))
                outfile = "galaxy_%s_map_%i_" % (filtname, index) + name + ".png"
                # Don't show these plots
                # PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "KS relation (data: H2 mass, method: grid)"
        id = abs(hash("galaxy KS relation H2 grid %i" % (index)))
        outfile = "KS_molecular_relation_grid_%i.png" % (index)
        caption = "KS relation. Surface densities were calculated using a grid with pixel size of 250 pc."
        caption += " Each blue dot shows the total SFR and H2 mass in the pixel divided by the pixel area."
        caption += " Black solid line indicates the median relation and shaded area the 84-16th percentiles."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "KS relation (data: H2+HI mass, method: grid)"
        id = abs(hash("galaxy KS relation H2+HI grid %i" % (index)))
        outfile = "KS_relation_best_grid_%i.png" % (index)
        caption = "KS relation. Surface densities were calculated using a grid with pixel size of 250 pc."
        caption += " Each blue dot shows the total SFR and H2 mass in the pixel divided by the pixel area."
        caption += " Black solid line indicates the median relation and shaded area the 84-16th percentiles."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "KS relation (data: H2 mass, method: Azimuthal average)"
        id = abs(hash("galaxy KS relation H2 radii %i" % (index)))
        outfile = "KS_molecular_relation_radii_%i.png" % (index)
        caption = "KS relation. Surface densities were calculated by azimuthally averaging radial concentric shells"
        caption += " of 800 pc of width. The shells are centered in the minimum of the dark matter potential."
        caption += " Each blue dot shows the total SFR and H2 mass in the shell divided by the shell area."
        caption += " Black solid line indicates the median relation and shaded area the 84-16th percentiles."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "KS relation (data: H2+HI mass, method: Azimuthal average)"
        id = abs(hash("galaxy KS relation H2+HI radii %i" % (index)))
        outfile = "KS_relation_best_radii_%i.png" % (index)
        caption = "KS relation. Surface densities were calculated by azimuthally averaging radial concentric shells"
        caption += " of 800 pc of width. The shells are centered in the minimum of the dark matter potential."
        caption += " Each blue dot shows the total SFR and H2 mass in the shell divided by the shell area."
        caption += " Black solid line indicates the median relation and shaded area the 84-16th percentiles."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "Depletion time (data: H2 mass, method: grid)"
        id = abs(hash("galaxy depletion H2 grid %i" % (index)))
        outfile = "molecular_gas_depletion_timescale_grid_%i.png" % (index)
        caption = "Gas depletion times. The surface densities were calculated using a grid with pixel size of 250 pc."
        caption += " Black solid line indicates the median relation, shaded area the 84-16th percentiles, "
        caption += "and the observational data-points correspond to Bigiel et al. (2008) inner, same as in KS relation (H2 mass) figure."
        # Don't shown dipletion time plots
        # PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "Depletion time (data: H2+HI mass, method: grid)"
        id = abs(hash("galaxy depletion H2+HI grid %i" % (index)))
        outfile = "gas_depletion_timescale_best_grid_%i.png" % (index)
        caption = "Gas depletion times. The surface densities were calculated using a grid with pixel size of 250 pc."
        caption += " Black solid line indicates the median relation, shaded area the 84-16th percentiles, "
        caption += "and the observational data-points correspond to Bigiel et al. (2008, 2010) inner, same as in KS relation (H2+HI mass) figure."
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
        caption += " with pixel size of 250 pc. Red dashed line corresponds to Krumholz+ (2009) semi-analytic model, the"
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


def loadAbundancePlots(
    web, output_path: str, name_list: List[int]
):
    """
    @TODO Create separate .yaml config containing all necessary information about the plots
    """

    PlotsInWeb = PlotsInPipeline()

    num_sims = len(name_list)

    for i in range(num_sims):
        title = name_list[i]
        caption = "Carbon abundance [C/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Grey dots correspond to individual stars in MW-type haloes, solid line shows the median relation. "
        caption += "Contours corresponds to GALAH DR3 data (Buder+21) and are constructed from histogram of star counts "
        caption += "using a log scale with a minimum star count of 10."
        filename = "C_Fe_"+name_list[i]+".png"
        id = abs(hash("Carbon %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    if num_sims > 1:
        title = "Comparison"
        caption = "Comparison between the [C/Fe]-[Fe/H] median relations from each simulation listed in this catalogue."
        filename = "C_Fe_comparison.png"
        id = abs(hash("Carbon comparison %i" % i))
        PlotsInWeb.load_plots(title, caption, filename, id)


    title = "Carbon"
    id = abs(hash("Carbon section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for i in range(num_sims):
        title = name_list[i]
        caption = "Nitrogen abundance [N/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Grey dots correspond to individual stars within MW-type haloes, the solid line shows the median relation. "
        filename = "N_Fe_"+name_list[i]+".png"
        id = abs(hash("Nitrogen %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    if num_sims > 1:
        title = "Comparison"
        caption = "Comparison between the [N/Fe]-[Fe/H] median relations from each simulation listed in this catalogue."
        filename = "N_Fe_comparison.png"
        id = abs(hash("Nitrogen_comparison %i" % i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Nitrogen"
    id = abs(hash("Nitrogen section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for i in range(num_sims):
        title = name_list[i]
        caption = "Oxygen abundance [O/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Left panel shows the stellar abundance of MW-type haloes, whereas right panel shows the abundance of satellites. "
        caption += "Grey dots correspond to individual stars, solid line shows the median relation. "
        caption += "Contours corresponds to GALAH DR3 data (Buder+21) and are constructed from histogram of star counts "
        caption += "using a log scale with a minimum star count of 10. "
        caption += "The observational data for MW, Carina, Fornax, Sculptor and Sagittarious, corresponds to a data compilation presented by Tolstoy, Hill & Tosi (2009)"
        filename = "O_Fe_"+name_list[i]+".png"
        id = abs(hash("Oxygen %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    if num_sims > 1:
        title = "Comparison"
        caption = "Comparison between the [O/Fe]-[Fe/H] median relations from each simulation listed in this catalogue."
        filename = "O_Fe_comparison.png"
        id = abs(hash("Oxygen comparison %i" % i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Oxygen"
    id = abs(hash("Oxygen section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for i in range(num_sims):
        title = name_list[i]
        caption = "Neon abundance [Ne/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Grey dots correspond to individual stars within MW-type haloes, the solid line shows the median relation. "
        filename = "Ne_Fe_"+name_list[i]+".png"
        id = abs(hash("Neon %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    if num_sims > 1:
        title = "Comparison"
        caption = "Comparison between the [Ne/Fe]-[Fe/H] median relations from each simulation listed in this catalogue."
        filename = "Ne_Fe_comparison.png"
        id = abs(hash("Neon_comparison %i" % i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Neon"
    id = abs(hash("Neon section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for i in range(num_sims):
        title = name_list[i]
        caption = "Magensium abundance [Mg/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Left panel shows the stellar abundance of MW-type haloes, whereas right panel shows the abundance of satellites. "
        caption += "Grey dots correspond to individual stars, solid line shows the median relation. "
        caption += "Contours corresponds to GALAH DR3 data (Buder+21) and are constructed from histogram of star counts "
        caption += "using a log scale with a minimum star count of 10. "
        caption += "The observational data for MW, Carina, Fornax, Sculptor and Sagittarious, corresponds to a data compilation presented by Tolstoy, Hill & Tosi (2009)"
        filename = "Mg_Fe_"+name_list[i]+".png"
        id = abs(hash("Magnesium %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    if num_sims > 1:
        title = "Comparison"
        caption = "Comparison between the [Mg/Fe]-[Fe/H] median relations from each simulation listed in this catalogue."
        filename = "Mg_Fe_comparison.png"
        id = abs(hash("Magnesium comparison %i" % i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Magnesium"
    id = abs(hash("Magnesium section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for i in range(num_sims):
        title = name_list[i]
        caption = "Silicon abundance [Si/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Grey dots correspond to individual stars within MW-type haloes, the solid line shows the median relation. "
        caption += "Contours corresponds to GALAH DR3 data (Buder+21) and are constructed from histogram of star counts "
        caption += "using a log scale with a minimum star count of 10. "
        filename = "Si_Fe_"+name_list[i]+".png"
        id = abs(hash("Silicon %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    if num_sims > 1:
        title = "Comparison"
        caption = "Comparison between the [Si/Fe]-[Fe/H] median relations from each simulation listed in this catalogue."
        filename = "Si_Fe_comparison.png"
        id = abs(hash("Silicon comparison %i" % i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Silicon"
    id = abs(hash("Silicon section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for i in range(num_sims):
        title = name_list[i]
        caption = "Strontium abundance [Sr/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Grey dots correspond to individual stars within MW-type haloes, the solid line shows the median relation. "
        caption += "Coloured symbols correspond to observational data from Roederer et al. (2014) and Spite et al. (2018)."
        filename = "Sr_Fe_"+name_list[i]+".png"
        id = abs(hash("Strontium %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    if num_sims > 1:
        title = "Comparison"
        caption = "Comparison between the [Sr/Fe]-[Fe/H] median relations from each simulation listed in this catalogue."
        filename = "Sr_Fe_comparison.png"
        id = abs(hash("Strontium_comparison %i" % i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Strontium"
    id = abs(hash("Strontium section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for i in range(num_sims):
        title = name_list[i]
        caption = "Barium abundance [Ba/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Grey dots correspond to individual stars within MW-type haloes, the solid line shows the median relation. "
        caption += "Contours corresponds to GALAH DR3 data (Buder+21) and are constructed from histogram of star counts "
        caption += "using a log scale with a minimum star count of 10. "
        filename = "Ba_Fe_"+name_list[i]+".png"
        id = abs(hash("Barium %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    if num_sims > 1:
        title = "Comparison"
        caption = "Comparison between the [Ba/Fe]-[Fe/H] median relations from each simulation listed in this catalogue."
        filename = "Ba_Fe_comparison.png"
        id = abs(hash("Barium comparison %i" % i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Barium"
    id = abs(hash("Barium section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for i in range(num_sims):
        title = name_list[i]
        caption = "Europium abundance [Eu/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Grey dots correspond to individual stars within MW-type haloes, the solid line shows the median relation. "
        caption += "Contours corresponds to GALAH DR3 data (Buder+21) and are constructed from histogram of star counts "
        caption += "using a log scale with a minimum star count of 10. "
        filename = "Eu_Fe_"+name_list[i]+".png"
        id = abs(hash("Europium %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    if num_sims > 1:
        title = "Comparison"
        caption = "Comparison between the [Eu/Fe]-[Fe/H] median relations from each simulation listed in this catalogue."
        filename = "Eu_Fe_comparison.png"
        id = abs(hash("Europium comparison %i" % i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Europium"
    id = abs(hash("Europium section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    title = "Stellar Mass - Z/Zsun relation (light-weighted r-band, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-Z/Zsun median relations from each simulation listed in this catalogue. "
    caption += "The values of Z/Zsun are obtained by calculated the light-weighted r-band average of the stars metallicity."
    filename = "Mstellar_Z_light_weighted_r_band_comparison.png"
    id = abs(hash("Mstellar Z light weighted r band comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - Z/Zsun relation (mass-weighted, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-Z/Zsun median relations from each simulation listed in this catalogue. "
    caption += "The values of Z/Zsun are obtained by calculated the mass-weighted average of the stars total metallicity."
    filename = "Mstellar_Z_mass_weighted_comparison.png"
    id = abs(hash("Mstellar Z mass weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Fe/H] relation (light-weight r-band of mass, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Fe/H] relations from each simulation listed in this catalogue. "
    caption += "The values of [Fe/H] are obtained by calculated the light-weighted r-band of total Fe and H, "
    caption += "then the log of the ratio is calculated and normalised by the corresponding solar abundances."
    filename = "Mstellar_Fe_H_light_weighted_r_band_with_mass_comparison.png"
    id = abs(hash("Mstellar FeH light weighted r band mass comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Fe/H] relation (light-weight r-band of ratio, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Fe/H] relations from each simulation listed in this catalogue. "
    caption += "The values of [Fe/H] are obtained by calculated the light-weighted r-band average of the ratio Fe/H. "
    caption += "Then the log is calculated and normalised by the corresponding solar abundances."
    filename = "Mstellar_Fe_H_light_weighted_r_band_comparison.png"
    id = abs(hash("Mstellar FeH light weighted r band ratio comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Fe/H] relation (mass-weighted, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Fe/H] relations from each simulation listed in this catalogue. "
    caption += "The values of [Fe/H] are obtained by calculating the total mass of Fe and H. "
    caption += "Then the log of the ratio is calculated and normalised by the corresponding solar abundances."
    filename = "Mstellar_Fe_H_mass_weighted_comparison.png"
    id = abs(hash("Mstellar FeH mass weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Fe/H] relation (mass-weight of ratio, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Fe/H] relations from each simulation listed in this catalogue. "
    caption += "The values of [Fe/H] are obtained by calculated the mass-weighted average of the ratio Fe/H. "
    caption += "Then the log is calculated and normalised by the corresponding solar abundances."
    filename = "Mstellar_Fe_H_ratio_weighted_comparison.png"
    id = abs(hash("Mstellar FeH ratio weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Fe/H] relation (median-of-log, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Fe/H] relations from each simulation listed in this catalogue. "
    caption += "The values of [Fe/H] are obtained by selecting star particles with (log of) ratio > -3 and calculating the median."
    filename = "Mstellar_Fe_H_median_comparison.png"
    id = abs(hash("Mstellar FeH comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Metallicity"
    id = abs(hash("Stellar Metallicity section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    ##############################

    ## Oxygen ##

    title = "Stellar Mass - [O/Fe] relation (light-weight r-band of mass)"
    caption = "Comparison between the Stellar mass-[O/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [O/Fe] are obtained by calculated the light-weighted r-band average of the total Iron and Oxygen in each galaxy, "
    caption += "then the log of the ratio, O/Fe, is calculated and normalised by the corresponding solar abundances."
    filename = "Mstellar_O_Fe_light_weighted_r_band_with_mass_comparison.png"
    id = abs(hash("Mstellar OFe light-weighted mass comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [O/Fe] relation (light-weighted r-band of ratio)"
    caption = "Comparison between the Stellar mass-[O/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [O/Fe] are obtained by calculated the light-weighted r-band average of the ratio O/Fe in each galaxy, "
    caption += "then the log is calculated and normalised by the corresponding solar abundances."
    filename = "Mstellar_O_Fe_light_weighted_r_band_comparison.png"
    id = abs(hash("Mstellar OFe light-weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [O/Fe] relation (mass-weighted, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[O/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [O/Fe] are obtained by calculating the total mass of Iron and Oxygen in each galaxy. "
    caption += "Then the log of the ratio, O/Fe, is calculated and normalised by the corresponding solar abundances."
    filename = "Mstellar_O_Fe_mass_weighted_comparison.png"
    id = abs(hash("Mstellar OFe mass weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [O/Fe] relation (mass-weighted of ratio, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[O/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [O/Fe] are obtained by calculated the mass-weighted average of the ratio O/Fe."
    caption += "Then the log is calculated and normalised by the corresponding solar abundances."
    filename = "Mstellar_O_Fe_ratio_weighted_comparison.png"
    id = abs(hash("Mstellar OFe ratio weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [O/Fe] relation (median-of-log, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[O/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [O/Fe] are obtained by selecting star particles with (log of) ratio of [Fe/H] > -3 and calculating the median."
    filename = "Mstellar_O_Fe_median_comparison.png"
    id = abs(hash("Mstellar OFe comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass-[O/Fe]"
    id = abs(hash("Stellar Mass-[O/Fe] section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    ## Magnesium ##

    title = "Stellar Mass - [Mg/Fe] relation (light-weight r-band of mass)"
    caption = "Comparison between the Stellar mass-[Mg/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [Mg/Fe] are obtained by calculated the light-weighted r-band average of the total Iron and Magnesium in each galaxy, "
    caption += "then the log of the ratio, Mg/Fe, is calculated and normalised by the corresponding solar abundances."
    filename = "Mstellar_Mg_Fe_light_weighted_r_band_with_mass_comparison.png"
    id = abs(hash("Mstellar MgFe light-weighted mass comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Mg/Fe] relation (light-weighted r-band of ratio)"
    caption = "Comparison between the Stellar mass-[Mg/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [Mg/Fe] are obtained by calculated the light-weighted r-band average of the ratio Mg/Fe in each galaxy, "
    caption += "then the log is calculated and normalised by the corresponding solar abundances."
    filename = "Mstellar_Mg_Fe_light_weighted_r_band_comparison.png"
    id = abs(hash("Mstellar MgFe light-weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Mg/Fe] relation (mass-weighted)"
    caption = "Comparison between the Stellar mass-[Mg/Fe] median relations from each simulation listed in this catalogue. "
    caption += "The values of [Mg/Fe] are obtained by computing the (log of) ratio between the total Magnesium mass in stars and "
    caption += "the total Iron mass in stars, and then normalising it by the corresponding solar abundances."
    filename = "Mstellar_Mg_Fe_mass_weighted_comparison.png"
    id = abs(hash("Mstellar MgFe total comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Mg/Fe] relation (mass-weighted of ratio, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Mg/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [Mg/Fe] are obtained by calculated the mass-weighted average of the ratio Mg/Fe."
    caption += "Then the log is calculated and normalised by the corresponding solar abundances."
    filename = "Mstellar_Mg_Fe_ratio_weighted_comparison.png"
    id = abs(hash("Mstellar MgFe ratio weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Mg/Fe] relation (median-of-log, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Mg/Fe] median relations from each simulation listed in this catalogue. "
    caption += "The values of [Mg/Fe] are obtained by selecting star particles with (log of) ratio of [Fe/H] > -3 and calculating the median."
    filename = "Mstellar_Mg_Fe_median_comparison.png"
    id = abs(hash("Mstellar MgFe comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass-[Mg/Fe]"
    id = abs(hash("Stellar Mass-[Mg/Fe] section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    ##############################

    for i in range(num_sims):
        title = name_list[i]
        caption = "Fe SNIa Fractions as a function of Iron abundance [Fe/H]. "
        caption += "Grey dots correspond to individual stars within MW-type haloes, the solid line shows the median relation. "
        filename = "FeSNIa_Fe_"+name_list[i]+".png"
        id = abs(hash("FeSNIa Fe %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    if num_sims > 1:
        title = "Fe SNIa Fractions"
        caption = "Comparison between the fraction of Fe from SNIa from each simulation listed in this catalogue. "
        caption += "The curves show the median relations for the mass fraction of Fe from SNIa vs [Fe/H]."
        filename = "FeSNIa_Fe_comparison.png"
        id = abs(hash("FeSNIa Fe comparison %i" % i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Fe(SNIa)/Fe"
    id = abs(hash("Fe(SNIa)/Fe section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    title = "SNIa Rates"
    caption = "Comparison between the SNIa rates from each simulation listed in this catalogue. "
    filename = "SNIa_rates_comparison.png"
    id = abs(hash("SNIa rates comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "SFH"
    caption = "Comparison between the SFH from each simulation listed in this catalogue. "
    filename = "SFH_comparison.png"
    id = abs(hash("SFH comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "SNIa Rates"
    id = abs(hash("SNIa rates section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    ##############################

    # title = "Stellar Mass - Z/Zsun relation (light-weighted i-band, 100 kpc aperture)"
    # caption = "Comparison between the Stellar mass-Z/Zsun median relations from each simulation listed in this catalogue. "
    # caption += "The values of Z/Zsun are obtained by calculated the light-weighted i-band average of the stars metallicity."
    # filename = "Mstellar_Z_light_weighted_i_band_comparison.png"
    # id = abs(hash("Mstellar Z light weighted i band comparison %i" % i))
    # PlotsInWeb.load_plots(title, caption, filename, id)
    #
    # title = "Stellar Mass - Z/Zsun relation (light-weighted z-band, 100 kpc aperture)"
    # caption = "Comparison between the Stellar mass-Z/Zsun median relations from each simulation listed in this catalogue. "
    # caption += "The values of Z/Zsun are obtained by calculated the light-weighted z-band average of the stars metallicity."
    # filename = "Mstellar_Z_light_weighted_z_band_comparison.png"
    # id = abs(hash("Mstellar Z light weighted z band comparison %i" % i))
    # PlotsInWeb.load_plots(title, caption, filename, id)
    #
    # title = "Stellar Mass - [Fe/H] relation (light-weighted i-band, 100 kpc aperture)"
    # caption = "Comparison between the Stellar mass-[Fe/H] median relations from each simulation listed in this catalogue. "
    # caption += "The values of [Fe/H] are obtained by calculated the light-weighted i-band average of [Fe/H]."
    # filename = "Mstellar_Fe_H_light_weighted_i_band_comparison.png"
    # id = abs(hash("Mstellar FeH light weighted i band comparison %i" % i))
    # PlotsInWeb.load_plots(title, caption, filename, id)
    #
    # title = "Stellar Mass - [Fe/H] relation (light-weighted z-band, 100 kpc aperture)"
    # caption = "Comparison between the Stellar mass-[Fe/H] median relations from each simulation listed in this catalogue. "
    # caption += "The values of [Fe/H] are obtained by calculated the light-weighted z-band average of [Fe/H]."
    # filename = "Mstellar_Fe_H_light_weighted_z_band_comparison.png"
    # id = abs(hash("Mstellar FeH light weighted z band comparison %i" % i))
    # PlotsInWeb.load_plots(title, caption, filename, id)
    #
    # title = "Extra plots"
    # id = abs(hash("Extra plots section"))
    # plots = PlotsInWeb.plots_details
    # caption = " "
    # add_web_section(web, title, caption, id, plots)
    # PlotsInWeb.reset_plots_list()
