
from plotter.html import add_web_section, PlotsInPipeline

def loadGalaxyPlots(web,name_list):

    PlotsInWeb = PlotsInPipeline()

    for index in range(10):

        for name in name_list:
            title = "Star particles ("+name+")"
            caption = "Projection of stars within 5 times the galaxy's stellar"
            caption += " half mass radius. Face-on (left) and edge-on (right)."
            id = abs(hash("galaxy stars parts %i" % (index)))
            outfile = "galaxy_sparts_%i_"%(index)+name+".png"
            PlotsInWeb.load_plots(title, caption, outfile, id)

        for name in name_list:
            title = "Gas particles ("+name+")"
            caption = "Projection of gas within 5 times the galaxy's stellar"
            caption += " half mass radius. Face-on (left) and edge-on (right)."
            id = abs(hash("galaxy gas parts %i" % (index)))
            outfile = "galaxy_parts_%i_"%(index)+name+".png"
            PlotsInWeb.load_plots(title, caption, outfile, id)

        name = name_list[0]

        title = "KS relation (data: H2 mass, method: grid)"
        id = abs(hash("galaxy KS relation H2 grid %i" % (index)))
        outfile = "KS_molecular_relation_grid_%i_"%(index)+name+".png"
        caption = "KS relation. Surface densities were calculated using a grid with pixel size of 250 pc."
        caption += " Each blue dot shows the total SFR and H2 mass in the pixel divided by the pixel area."
        caption += " Black solid line indicates the median relation and shaded area the 84-16th percentiles."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "KS relation (data: H2+HI mass, method: grid)"
        id = abs(hash("galaxy KS relation H2+HI grid %i" % (index)))
        outfile = "KS_relation_best_grid_%i_"%(index)+name+".png"
        caption = "KS relation. Surface densities were calculated using a grid with pixel size of 250 pc."
        caption += " Each blue dot shows the total SFR and H2 mass in the pixel divided by the pixel area."
        caption += " Black solid line indicates the median relation and shaded area the 84-16th percentiles."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "KS relation (data: H2 mass, method: Azimuthal average)"
        id = abs(hash("galaxy KS relation H2 radii %i" % (index)))
        outfile = "KS_molecular_relation_radii_%i_"%(index)+name+".png"
        caption = "KS relation. Surface densities were calculated by azimuthally averaging radial concentric shells"
        caption += " of 800 pc of width. The shells are centered in the minimum of the dark matter potential."
        caption += " Each blue dot shows the total SFR and H2 mass in the shell divided by the shell area."
        caption += " Black solid line indicates the median relation and shaded area the 84-16th percentiles."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "KS relation (data: H2+HI mass, method: Azimuthal average)"
        id = abs(hash("galaxy KS relation H2+HI radii %i" % (index)))
        outfile = "KS_relation_best_radii_%i_"%(index)+name+".png"
        caption = "KS relation. Surface densities were calculated by azimuthally averaging radial concentric shells"
        caption += " of 800 pc of width. The shells are centered in the minimum of the dark matter potential."
        caption += " Each blue dot shows the total SFR and H2 mass in the shell divided by the shell area."
        caption += " Black solid line indicates the median relation and shaded area the 84-16th percentiles."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "Depletion time (data: H2 mass, method: grid)"
        id = abs(hash("galaxy depletion H2 grid %i" % (index)))
        outfile = "molecular_gas_depletion_timescale_grid_%i_"%(index)+name+".png"
        caption = "Gas depletion times. The surface densities were calculated using a grid with pixel size of 250 pc."
        caption += " Black solid line indicates the median relation, shaded area the 84-16th percentiles, "
        caption += "and the observational data-points correspond to Bigiel et al. (2008) inner, same as in KS relation (H2 mass) figure."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "Depletion time (data: H2+HI mass, method: grid)"
        id = abs(hash("galaxy depletion H2+HI grid %i" % (index)))
        outfile = "gas_depletion_timescale_best_grid_%i_"%(index)+name+".png"
        caption = "Gas depletion times. The surface densities were calculated using a grid with pixel size of 250 pc."
        caption += " Black solid line indicates the median relation, shaded area the 84-16th percentiles, "
        caption += "and the observational data-points correspond to Bigiel et al. (2008, 2010) inner, same as in KS relation (H2+HI mass) figure."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "Depletion time (data: H2 mass, method: Azimuthal average)"
        id = abs(hash("galaxy depletion H2 radii %i" % (index)))
        outfile = "molecular_gas_depletion_timescale_radii_%i_"%(index)+name+".png"
        caption = "Gas depletion times. The surface densities were calculated by azimuthally averaging radial concentric shells"
        caption += " of 800 pc of width. The shells were centered in the minimum of the dark matter potential."
        caption += " Black solid line indicates the median relation, shaded area the 84-16th percentiles, "
        caption += "and the observational data-points correspond to Bigiel et al. (2008) inner, same as in KS relation (H2 mass) figure."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "Depletion time (data: H2+HI mass, method: Azimuthal average)"
        id = abs(hash("galaxy depletion H2+HI radii %i" % (index)))
        outfile = "gas_depletion_timescale_best_radii_%i_"%(index)+name+".png"
        caption = "Gas depletion times. The surface densities were calculated by azimuthally averaging radial concentric shells"
        caption += " of 800 pc of width. The shells were centered in the minimum of the dark matter potential."
        caption += " Black solid line indicates the median relation, shaded area the 84-16th percentiles, "
        caption += "and the observational data-points correspond to Bigiel et al. (2008, 2010) inner, same as in KS relation (H2+HI mass) figure."
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "Surface density ratios (method: grid)"
        id = abs(hash("density ratio H2+HI grid %i" % (index)))
        outfile = "Surface_density_ratio_grid_%i_"%(index)+name+".png"
        caption = "Surface density ratios. The y-axis shows the ratio between surface densities calculated using a grid"
        caption += " with pixel size of 250 pc. Red dashed line corresponds to Krumholz+ (2009) semi-analytic model, the"
        caption += " black solid line indicates the median relation and the shaded area the 84-16th percentiles, "
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = "Surface density ratios (method: Azimuthal average)"
        id = abs(hash("density ratio H2+HI radii %i" % (index)))
        outfile = "Surface_density_ratio_radii_%i_"%(index)+name+".png"
        caption = "Surface density ratios. The y-axis shows the ratio between surface densities calculated"
        caption += " by azimuthally averaging radial concentric shells of 800 pc of width. The shells were centered"
        caption += " in the minimum of the dark matter potential."
        caption += " The red dashed and dotted lines correspond to Krumholz+ (2009) semi-analytic model,"
        caption += " black solid line indicates the median relation and the shaded area the 84-16th percentiles, "
        PlotsInWeb.load_plots(title, caption, outfile, id)

        title = '%i Galaxy ' % (index + 1)
        id = abs(hash("galaxy and ks relation %i" % index))
        plots = PlotsInWeb.plots_details
        add_web_section(web, title, id, plots)
        PlotsInWeb.reset_plots_list()
