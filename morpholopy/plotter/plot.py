from pylab import *
from plotter.html import add_web_section
from plotter.KS_relation import median_relations
import numpy as np

def plot_momentum(stellar_mass,momentum,parttype,MorphologyPlotsInWeb,output_path):
    
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
    plt.xlim(1e6, 1e12)
    plt.ylim(1e-1, 1e4)
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.savefig(f"{output_path}/momentum_parttype_%i.png"%parttype, dpi=200)
    plt.close()

    if parttype==4: title = "Specific angular momentum / Stars"
    if parttype==0: title = "Specific angular momentum / HI+H2 gas"

    caption = "Ratio between the total angular momentum of stars (or gas) within 30 kpc of "
    caption += "aperture divided by the total mass in stars (or gas)."
    filename = "momentum_parttype_%i.png"%parttype
    id = abs(hash("momentum %i"%parttype))

    MorphologyPlotsInWeb.load_plots(title, caption, filename, id)

def plot_kappa(stellar_mass,kappa,parttype,MorphologyPlotsInWeb,output_path):
    
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
    plt.xlim(1e6, 1e12)
    plt.ylim(0, 1)
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)

    plt.savefig(f"{output_path}/Kappa_co_parttype_%i.png"%parttype, dpi=200)
    plt.close()

    if parttype==4: title = "Kappa corotation / Stars"
    if parttype==0: title = "Kappa corotation / HI+H2 gas"

    caption = "Kappa corotation is defined as the fraction of kinetic energy in a galaxy "
    caption += "that is in ordered rotation. Note that the rotating contribution is calculated "
    caption += "only for prograde rotation."
    filename = "Kappa_co_parttype_%i.png"%parttype
    id = abs(hash("kappa co %i"%parttype))
    MorphologyPlotsInWeb.load_plots(title, caption, filename, id)


def plot_axis_ratios(stellar_mass,axis_ratios,parttype,MorphologyPlotsInWeb,output_path):
    
    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (11, 3.5),
        "figure.subplot.left": 0.05,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.15,
        "figure.subplot.top": 0.9,
        "figure.subplot.wspace": 0.3,
        "figure.subplot.hspace": 0.3,
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
    plt.xlim(1e6, 1e12)
    plt.ylim(0.0, 1.0)
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)

    ########
    ax = plt.subplot(1,3,2)
    plt.grid("True")

    plt.plot(stellar_mass, axis_ratios[:,1], 'o', color=color)

    plt.xscale('log')
    plt.xlabel("Stellar Mass [M$_{\odot}$]")
    plt.ylabel("c/b")
    plt.xlim(1e6, 1e12)
    plt.ylim(0.0, 1.0)
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)

    ########
    ax = plt.subplot(1,3,3)
    plt.grid("True")
    
    plt.plot(stellar_mass, axis_ratios[:,2], 'o', color=color)
    
    plt.xscale('log')
    plt.xlabel("Stellar Mass [M$_{\odot}$]")
    plt.ylabel("b/a")
    plt.xlim(1e6, 1e12)
    plt.ylim(0.0, 1.0)
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.savefig(f"{output_path}/Axis_ratios_parttype_%i.png"%parttype, dpi=200)
    plt.close()

    if parttype==4: title = "Axis ratios / Stars"
    if parttype==0: title = "Axis ratios / HI+H2 gas"
    caption = "Axial ratios of galaxies more massive than 1e6 Msun in stellar mass. "
    caption += "a, b and c (a >= b >= c) represent the lengths of the primary axes. "
    caption += "Ratios have been calculated following eqs. (1) and (2) from Trayford+2018."
    filename = "Axis_ratios_parttype_%i.png"%parttype
    id = abs(hash("galaxy axis %i"%parttype))
    MorphologyPlotsInWeb.load_plots(title, caption, filename, id)


def plot_surface_densities(sigma_SFR, sigma_gas, sigma_H2, stellar_mass, MorphologyPlotsInWeb, output_path):

    # Get the default KS relation for correct IMF
    def KS(sigma_g, n, A):
        return A * sigma_g ** n

    Sigma_g = np.logspace(-2, 3, 1000)
    Sigma_star = KS(Sigma_g, 1.4, 1.515e-4)

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (5.5, 4),
        "figure.subplot.left": 0.15,
        "figure.subplot.right": 0.85,
        "figure.subplot.bottom": 0.18,
        "figure.subplot.top": 0.9,
        "lines.markersize": 4,
        "lines.linewidth": 2.0,
    }
    rcParams.update(params)

    fig = figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(np.log10(Sigma_g), np.log10(Sigma_star), '--',color='grey',label=r"1.51e-4 $\times$ $\Sigma_{g}^{1.4}$")
    plt.scatter(sigma_H2, sigma_SFR, c=stellar_mass, alpha=.9, s=25,
                vmin=6, vmax=10, cmap='CMRmap_r', edgecolors='none', zorder=2)

    plt.xlabel("log $\\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
    plt.ylabel("log $\\Sigma_{\\rm SFR}$ $[{\\rm M_\\odot \\cdot yr^{-1} \\cdot kpc^{-2}}]$")
    plt.xlim(-1.0, 3.0)
    plt.ylim(-6.0, 0.0)

    plt.legend()
    cbar_ax = fig.add_axes([0.87, 0.18, 0.018, 0.5])
    cbar_ax.tick_params(labelsize=15)
    cb = plt.colorbar(ticks=[6,7,8,9,10,11,12], cax=cbar_ax)
    cb.set_label(label='$\log_{10}$ M$_{*}$/M$_{\odot}$', labelpad=0.5)
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.savefig(f"{output_path}/surface_density_gas.png", dpi=200)
    plt.close()

    caption = "Integrated surface densities of H2 gas and star-forming gas for each individual galaxy. "
    caption += "Quantities are calculated summing up all gas (and SFR) within the galaxies' stellar half mass radius."
    filename = "surface_density_gas.png"
    id = abs(hash("surface_density_gas"))

    MorphologyPlotsInWeb.load_plots(title, caption, filename, id)

    #######
    fig = figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(np.log10(Sigma_g), np.log10(Sigma_star), '--',color='grey',label=r"1.51e-4 $\times$ $\Sigma_{g}^{1.4}$")
    plt.scatter(sigma_gas, sigma_SFR, c=stellar_mass, alpha=.9, s=25,
                vmin=6, vmax=10, cmap='CMRmap_r', edgecolors='none', zorder=2)


    plt.xlabel("log $\\Sigma_{HI}+ \\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
    plt.ylabel("log $\\Sigma_{\\rm SFR}$ $[{\\rm M_\\odot \\cdot yr^{-1} \\cdot kpc^{-2}}]$")
    plt.xlim(-1.0, 3.0)
    plt.ylim(-6.0, 0.0)
    plt.legend()
    cbar_ax = fig.add_axes([0.87, 0.18, 0.018, 0.5])
    cbar_ax.tick_params(labelsize=15)
    cb = plt.colorbar(ticks=[6,7,8,9,10,11,12], cax=cbar_ax)
    cb.set_label(label='$\log_{10}$ M$_{*}$/M$_{\odot}$', labelpad=0.5)
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.savefig(f"{output_path}/surface_density_H2.png", dpi=200)
    plt.close()

    caption = "Integrated surface densities of H2+HI gas and star-forming gas for each individual galaxy. "
    caption += "Quantities are calculated summing up all gas (and SFR) within the galaxies' stellar half mass radius."
    filename = "surface_density_H2.png"
    id = abs(hash("surface_density_H2"))

    MorphologyPlotsInWeb.load_plots(title, caption, filename, id)


def plot_accumulative_densities(galaxy_data, MorphologyPlotsInWeb, output_path):

    # Get the default KS relation for correct IMF
    def KS(sigma_g, n, A):
        return A * sigma_g ** n

    Sigma_g = np.logspace(-2, 3, 1000)
    Sigma_star = KS(Sigma_g, 1.4, 1.515e-4)

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (5.5, 4),
        "figure.subplot.left": 0.15,
        "figure.subplot.right": 0.85,
        "figure.subplot.bottom": 0.18,
        "figure.subplot.top": 0.9,
        "lines.markersize": 2,
        "lines.linewidth": 2.0,
    }
    rcParams.update(params)

    fig = figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(np.log10(Sigma_g), np.log10(Sigma_star), '--',color='grey',label=r"1.51e-4 $\times$ $\Sigma_{g}^{1.4}$")

    sigma_gas = galaxy_data.surface_density[1:]
    sigma_SFR = galaxy_data.SFR_density[1:]
    metals = galaxy_data.metallicity[1:]
    metals[metals==0] = 1e-6
    metals = np.log10(np.array(metals, dtype='float'))
    arg_sort = np.argsort(metals)
    metals = metals[arg_sort[::-1]]
    sigma_gas = sigma_gas[arg_sort[::-1]]
    sigma_SFR = sigma_SFR[arg_sort[::-1]]


    x, y, y_down, y_up = median_relations(sigma_gas, sigma_SFR)

    plt.scatter(sigma_gas, sigma_SFR, c=metals, alpha=.9, s=10,
                vmin=-3, vmax=1, cmap='CMRmap_r', edgecolors='none', zorder=2)
    plt.plot(x, y, '-', color='grey')

    plt.xlabel("log $\\Sigma_{HI} + \\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
    plt.ylabel("log $\\Sigma_{\\rm SFR}$ $[{\\rm M_\\odot \\cdot yr^{-1} \\cdot kpc^{-2}}]$")
    plt.xlim(-1.0, 3.0)
    plt.ylim(-6.0, 0.0)
    plt.legend()

    cbar_ax = fig.add_axes([0.87, 0.18, 0.018, 0.5])
    cbar_ax.tick_params(labelsize=15)
    cb = plt.colorbar(ticks=[-3,-2,-1,0,1], cax=cbar_ax)
    cb.set_label(label='Z$_{gas}$/Z$_{\odot}$', labelpad=0.5)
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.savefig(f"{output_path}/accumulative_surface_density_gas.png", dpi=200)
    plt.close()

    caption = "Combined spatially resolved measurements from the ten most massive individual galaxies,"
    caption += " coloured by the mean metallicity of the resolved pixel. The surface densities were calculated" \
               "using the grid method."
    filename = "accumulative_surface_density_gas.png"
    id = abs(hash("accumulative_surface_density_gas"))

    MorphologyPlotsInWeb.load_plots(title, caption, filename, id)

    #######
    fig = figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    sigma_gas = galaxy_data.surface_density[1:]
    sigma_ratio = galaxy_data.ratio_densities[1:]
    metals = galaxy_data.metallicity[1:]
    metals[metals==0] = 1e-6
    metals = np.log10(np.array(metals, dtype='float'))
    arg_sort = np.argsort(metals)
    metals = metals[arg_sort[::-1]]
    sigma_gas = sigma_gas[arg_sort[::-1]]
    sigma_ratio = sigma_ratio[arg_sort[::-1]]

    x, y, y_down, y_up = median_relations(sigma_gas, sigma_ratio)

    plt.scatter(sigma_gas, sigma_ratio, c=metals, alpha=.9, s=10,
                vmin=-3, vmax=1, cmap='CMRmap_r', edgecolors='none', zorder=2)
    plt.plot(x, y, '-', color='grey')

    plt.xlabel("log $\\Sigma_{HI} + \\Sigma_{H_2}$  $[{\\rm M_\\odot\\cdot pc^{-2}}]$")
    plt.ylabel(r"log $\Sigma_{\mathrm{H2}} / (\Sigma_{\mathrm{HI}}+\Sigma_{\mathrm{H2}})$")
    plt.xlim(-1.0, 3.0)
    plt.ylim(-8.0, 0.5)

    cbar_ax = fig.add_axes([0.87, 0.18, 0.018, 0.5])
    cbar_ax.tick_params(labelsize=15)
    cb = plt.colorbar(ticks=[-3,-2,-1,0,1], cax=cbar_ax)
    cb.set_label(label='Z$_{gas}$/Z$_{\odot}$', labelpad=0.5)
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.savefig(f"{output_path}/accumulative_surface_density_ratios.png", dpi=200)
    plt.close()

    caption = "Combined spatially resolved measurements from the ten most massive individual galaxies,"
    caption += " coloured by the mean metallicity of the resolved pixel. The surface densities were calculated" \
               "using the grid method."
    filename = "accumulative_surface_density_ratios.png"
    id = abs(hash("accumulative_surface_density_ratios"))

    MorphologyPlotsInWeb.load_plots(title, caption, filename, id)


def plot_morphology(galaxy_data,web,MorphologyPlotsInWeb,output_path):

    # plot kappa for stars and gas :
    plot_kappa(10**galaxy_data.stellar_mass,galaxy_data.kappa_co,4,MorphologyPlotsInWeb,output_path)
    plot_kappa(10**galaxy_data.stellar_mass,galaxy_data.gas_kappa_co,0,MorphologyPlotsInWeb,output_path)

    title = 'Kappa corotation'
    id = abs(hash("Kappa corotation"))
    plots = MorphologyPlotsInWeb.plots_details
    add_web_section(web, title, id, plots)

    MorphologyPlotsInWeb.reset_plots_list()

    # plot specific angular momentum  for stars and gas :
    plot_momentum(10**galaxy_data.stellar_mass,galaxy_data.momentum,4,MorphologyPlotsInWeb,output_path)
    plot_momentum(10**galaxy_data.stellar_mass,galaxy_data.gas_momentum,0,MorphologyPlotsInWeb,output_path)

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
    plot_axis_ratios(10**galaxy_data.stellar_mass,axis_ratios,4,MorphologyPlotsInWeb,output_path)

    axis_ratios = np.zeros((len(galaxy_data.gas_axis_ca),3))
    axis_ratios[:,0] = galaxy_data.gas_axis_ca
    axis_ratios[:,1] = galaxy_data.gas_axis_cb
    axis_ratios[:,2] = galaxy_data.gas_axis_ba
    plot_axis_ratios(10**galaxy_data.stellar_mass,axis_ratios,0,MorphologyPlotsInWeb,output_path)

    title = 'Axis ratios'
    id = abs(hash("axis ratios"))
    plots = MorphologyPlotsInWeb.plots_details
    add_web_section(web, title, id, plots)

    MorphologyPlotsInWeb.reset_plots_list()

    # plot surface densities
    plot_surface_densities(galaxy_data.sigma_SFR,galaxy_data.sigma_gas,galaxy_data.sigma_H2,
                           galaxy_data.stellar_mass,MorphologyPlotsInWeb,output_path)

    title = 'Integrated surface densities'
    id = abs(hash("Integrated Surface density"))
    plots = MorphologyPlotsInWeb.plots_details
    add_web_section(web, title, id, plots)

    MorphologyPlotsInWeb.reset_plots_list()

    # plot surface densities
    plot_accumulative_densities(galaxy_data, MorphologyPlotsInWeb, output_path)

    title = 'Combined surface densities'
    id = abs(hash("Surface density"))
    plots = MorphologyPlotsInWeb.plots_details
    add_web_section(web, title, id, plots)
