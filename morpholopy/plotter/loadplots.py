import numpy as np
from plotter.html import add_web_section, PlotsInPipeline

def loadGalaxyPlots(web,name_list,output_path):

    PlotsInWeb = PlotsInPipeline()

    for index in range(10):

        for name in name_list:
            title = "Gas component ("+name+")"
            caption = "Projection of gas within 5 times the galaxy's  stellar " \
                      "half mass radius. Face-on (left) and edge-on (right)."
            id = abs(hash("galaxy gas %i " %(index) + name))
            outfile = "galaxy_gas_%i_"%(index)+name+".png"
            PlotsInWeb.load_plots(title, caption, outfile, id)

        for name in name_list:
            title = "Stellar component ("+name+")"
            caption = "Projection of stars within 5 times the galaxy's  stellar " \
                      "half mass radius. Face-on (left) and edge-on (right)."
            id = abs(hash("stars galaxy %i " %(index) + name))
            outfile = "galaxy_stars_%i_"%(index)+name+".png"
            PlotsInWeb.load_plots(title, caption, outfile, id)

        for name in name_list:
            title = "Star particles ("+name+")"
            caption = "Projection of stars within 5 times the galaxy's stellar"
            caption += " half mass radius. Face-on (left) and edge-on (right)."
            id = abs(hash("galaxy stars parts %i " % (index) + name))
            outfile = "galaxy_sparts_%i_"%(index)+name+".png"
            PlotsInWeb.load_plots(title, caption, outfile, id)

        for name in name_list:
            title = "Gas particles ("+name+")"
            caption = "Projection of gas within 5 times the galaxy's stellar"
            caption += " half mass radius. Face-on (left) and edge-on (right)."
            id = abs(hash("galaxy gas parts %i " % (index) + name))
            outfile = "galaxy_parts_%i_"%(index)+name+".png"
            PlotsInWeb.load_plots(title, caption, outfile, id)

        for name in name_list:
            for filtname in ['u', 'r', 'K']:
                title = f"Rest frame {filtname}-band light "
                caption = "Projection %s-band emission within 5 times " % (filtname)
                caption += "the galaxy's stellar half mass radius. Face-on (left) and edge-on (right)."
                id = abs(hash("%s-band emission galaxy %i" % (filtname, index)+name))
                outfile = "galaxy_%s_map_%i_" % (filtname, index) + name+".png"
                PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "KS relation (data: H2 mass, method: grid)"
        id = abs(hash("galaxy KS relation H2 grid %i" % (index)))
        outfile = "KS_molecular_relation_grid_%i.png"%(index)
        caption = "KS relation. Surface densities were calculated using a grid with pixel size of 250 pc."
        caption += " Each blue dot shows the total SFR and H2 mass in the pixel divided by the pixel area."
        caption += " Black solid line indicates the median relation and shaded area the 84-16th percentiles."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "KS relation (data: H2+HI mass, method: grid)"
        id = abs(hash("galaxy KS relation H2+HI grid %i" % (index)))
        outfile = "KS_relation_best_grid_%i.png"%(index)
        caption = "KS relation. Surface densities were calculated using a grid with pixel size of 250 pc."
        caption += " Each blue dot shows the total SFR and H2 mass in the pixel divided by the pixel area."
        caption += " Black solid line indicates the median relation and shaded area the 84-16th percentiles."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "KS relation (data: H2 mass, method: Azimuthal average)"
        id = abs(hash("galaxy KS relation H2 radii %i" % (index)))
        outfile = "KS_molecular_relation_radii_%i.png"%(index)
        caption = "KS relation. Surface densities were calculated by azimuthally averaging radial concentric shells"
        caption += " of 800 pc of width. The shells are centered in the minimum of the dark matter potential."
        caption += " Each blue dot shows the total SFR and H2 mass in the shell divided by the shell area."
        caption += " Black solid line indicates the median relation and shaded area the 84-16th percentiles."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "KS relation (data: H2+HI mass, method: Azimuthal average)"
        id = abs(hash("galaxy KS relation H2+HI radii %i" % (index)))
        outfile = "KS_relation_best_radii_%i.png"%(index)
        caption = "KS relation. Surface densities were calculated by azimuthally averaging radial concentric shells"
        caption += " of 800 pc of width. The shells are centered in the minimum of the dark matter potential."
        caption += " Each blue dot shows the total SFR and H2 mass in the shell divided by the shell area."
        caption += " Black solid line indicates the median relation and shaded area the 84-16th percentiles."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "Depletion time (data: H2 mass, method: grid)"
        id = abs(hash("galaxy depletion H2 grid %i" % (index)))
        outfile = "molecular_gas_depletion_timescale_grid_%i.png"%(index)
        caption = "Gas depletion times. The surface densities were calculated using a grid with pixel size of 250 pc."
        caption += " Black solid line indicates the median relation, shaded area the 84-16th percentiles, "
        caption += "and the observational data-points correspond to Bigiel et al. (2008) inner, same as in KS relation (H2 mass) figure."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "Depletion time (data: H2+HI mass, method: grid)"
        id = abs(hash("galaxy depletion H2+HI grid %i" % (index)))
        outfile = "gas_depletion_timescale_best_grid_%i.png"%(index)
        caption = "Gas depletion times. The surface densities were calculated using a grid with pixel size of 250 pc."
        caption += " Black solid line indicates the median relation, shaded area the 84-16th percentiles, "
        caption += "and the observational data-points correspond to Bigiel et al. (2008, 2010) inner, same as in KS relation (H2+HI mass) figure."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "Depletion time (data: H2 mass, method: Azimuthal average)"
        id = abs(hash("galaxy depletion H2 radii %i" % (index)))
        outfile = "molecular_gas_depletion_timescale_radii_%i.png"%(index)
        caption = "Gas depletion times. The surface densities were calculated by azimuthally averaging radial concentric shells"
        caption += " of 800 pc of width. The shells were centered in the minimum of the dark matter potential."
        caption += " Black solid line indicates the median relation, shaded area the 84-16th percentiles, "
        caption += "and the observational data-points correspond to Bigiel et al. (2008) inner, same as in KS relation (H2 mass) figure."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "Depletion time (data: H2+HI mass, method: Azimuthal average)"
        id = abs(hash("galaxy depletion H2+HI radii %i" % (index)))
        outfile = "gas_depletion_timescale_best_radii_%i.png"%(index)
        caption = "Gas depletion times. The surface densities were calculated by azimuthally averaging radial concentric shells"
        caption += " of 800 pc of width. The shells were centered in the minimum of the dark matter potential."
        caption += " Black solid line indicates the median relation, shaded area the 84-16th percentiles, "
        caption += "and the observational data-points correspond to Bigiel et al. (2008, 2010) inner, same as in KS relation (H2+HI mass) figure."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "Surface density ratios (method: grid)"
        id = abs(hash("density ratio H2+HI grid %i" % (index)))
        outfile = "Surface_density_ratio_grid_%i.png"%(index)
        caption = "Surface density ratios. The y-axis shows the ratio between surface densities calculated using a grid"
        caption += " with pixel size of 250 pc. Red dashed line corresponds to Krumholz+ (2009) semi-analytic model, the"
        caption += " black solid line indicates the median relation and the shaded area the 84-16th percentiles, "
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "Surface density ratios (method: Azimuthal average)"
        id = abs(hash("density ratio H2+HI radii %i" % (index)))
        outfile = "Surface_density_ratio_radii_%i.png"%(index)
        caption = "Surface density ratios. The y-axis shows the ratio between surface densities calculated"
        caption += " by azimuthally averaging radial concentric shells of 800 pc of width. The shells were centered"
        caption += " in the minimum of the dark matter potential."
        caption += " The red dashed and dotted lines correspond to Krumholz+ (2009) semi-analytic model,"
        caption += " black solid line indicates the median relation and the shaded area the 84-16th percentiles, "
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = '%i Galaxy ' % (index + 1)
        caption = " "
        for name in name_list:
            data = np.loadtxt(f"{output_path}/galaxy_data_"+name+".txt")
            sfr_galaxy = data[:,0]
            mass_galaxy = data[:,1]
            gas_mass_galaxy = data[:,2]
            mass_halo = data[:,3]
            galaxy_metallicity_gas_sfr = data[:,4]
            galaxy_metallicity_gas = data[:,5]

            caption += "<strong>Simulation: "+name+"</strong>. Galaxy details: "
            caption += r"$\log_{10}$ M$_{200}$/M$_{\odot} = $%0.2f," % (mass_halo[index])
            caption += " SFR = %0.1f M$_{\odot}$/yr," % (sfr_galaxy[index])
            caption += " Z$_{\mathrm{SFR}>0}$ = %0.3f," % (galaxy_metallicity_gas_sfr[index])
            caption += " Z$_{\mathrm{gas}}$ = %0.3f," % (galaxy_metallicity_gas[index])
            caption += " $\log_{10}$ M$_{*}$/M$_{\odot} = $%0.2f" % (mass_galaxy[index])
            caption += " $\&$ $\log_{10}$ M$_{\mathrm{gas}}$/M$_{\odot} = $%0.2f.</p><p style=\"font-size:18px;\">" % (gas_mass_galaxy[index])

        id = abs(hash("galaxy and ks relation %i" % index))
        plots = PlotsInWeb.plots_details
        add_web_section(web, title, caption, id, plots)
        PlotsInWeb.reset_plots_list()


    title = "Specific angular momentum / Stars"
    caption = "Ratio between the total angular momentum of stars within 30 kpc of "
    caption += "aperture divided by the total mass in stars."
    filename = "momentum_parttype_%i.png" %4
    id = abs(hash("momentum %i" %4))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Specific angular momentum / HI+H2 gas"
    caption = "Ratio between the total angular momentum of gas within 30 kpc of "
    caption += "aperture divided by the total mass in gas."
    filename = "momentum_parttype_%i.png" %0
    id = abs(hash("momentum %i" %0))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = 'Specific angular momentum'
    id = abs(hash("angular momentum"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    title = "Kappa corotation / Stars"
    caption = "Kappa corotation is defined as the fraction of kinetic energy in a galaxy "
    caption += "that is in ordered rotation. Note that the rotating contribution is calculated "
    caption += "only for prograde rotation."
    filename = "Kappa_co_parttype_%i.png" %4
    id = abs(hash("kappa co %i" %4))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Kappa corotation / HI+H2 gas"
    caption = "Kappa corotation is defined as the fraction of kinetic energy in a galaxy "
    caption += "that is in ordered rotation. Note that the rotating contribution is calculated "
    caption += "only for prograde rotation."
    filename = "Kappa_co_parttype_%i.png" %0
    id = abs(hash("kappa co %i" %0))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = 'Kappa corotation'
    id = abs(hash("Kappa corotation"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    title = "Axis ratios / Stars"
    caption = "Axial ratios of galaxies more massive than 1e6 Msun in stellar mass. "
    caption += "a, b and c (a >= b >= c) represent the lengths of the primary axes. "
    caption += "Ratios have been calculated following eqs. (1) and (2) from Trayford+2018."
    filename = "Axis_ratios_parttype_%i.png" %4
    id = abs(hash("galaxy axis %i" %4))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Axis ratios / HI+H2 gas"
    caption = "Axial ratios of galaxies more massive than 1e6 Msun in stellar mass. "
    caption += "a, b and c (a >= b >= c) represent the lengths of the primary axes. "
    caption += "Ratios have been calculated following eqs. (1) and (2) from Trayford+2018."
    filename = "Axis_ratios_parttype_%i.png" %0
    id = abs(hash("galaxy axis %i" %0))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = 'Axis ratios'
    id = abs(hash("axis ratios"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for name in name_list:
        title = "Integrated surface densities (" + name + ")"
        caption = "Integrated surface densities of H2+HI gas and star-forming gas for each individual galaxy. "
        caption += "Quantities are calculated summing up all gas (and SFR) within the galaxies' stellar half mass radius."
        filename = "surface_density_H2_"+name+".png"
        id = abs(hash("surface_density_H2_"+name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for name in name_list:
        title = "Integrated surface densities (" + name + ")"
        caption = "Integrated surface densities of H2 gas and star-forming gas for each individual galaxy. "
        caption += "Quantities are calculated summing up all gas (and SFR) within the galaxies' stellar half mass radius."
        filename = "surface_density_gas_"+name+".png"
        id = abs(hash("surface_density_gas_"+name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = 'Integrated surface densities'
    id = abs(hash("Integrated Surface density"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for name in name_list:
        title = "Combined spatially resolved surface density ratios (" + name + ")"
        caption = "Combined spatially resolved measurements from the ten most massive individual galaxies,"
        caption += " coloured by the mean metallicity of the resolved pixel. The surface densities were calculated" \
                   " using the grid method with a pixel size of 250pc. Coloured solid lines show the median relations" \
                   " considering only cells with fixed metallicity (as indicated in the legends). The grey solid line" \
                   " shows the median relation for all pixels."
        filename = "combined_surface_density_ratios_"+name+".png"
        id = abs(hash("combined_surface_density_ratios_"+name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for name in name_list:
        title = "Combined spatially resolved surface density (" + name + ")"
        caption = "Combined spatially resolved measurements from the ten most massive individual galaxies,"
        caption += " coloured by the mean metallicity of the resolved pixel. The surface densities were calculated" \
                   " using the grid method with a pixel size of 250pc. Coloured lines show in the median relations" \
                   " considering only cells with fixed metallicity (as indicated in the legends). The grey solid line" \
                   " shows the median relation for all pixels, whereas the dashed black line shows the relation" \
                   " only for pixels that have SFR surface density >0."
        filename = "combined_surface_density_gas_"+name+".png"
        id = abs(hash("combined_surface_density_gas_"+name))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = 'Combined surface densities'
    id = abs(hash("Surface density"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()