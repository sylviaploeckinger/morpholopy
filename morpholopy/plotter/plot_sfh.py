import matplotlib.pylab as plt
from matplotlib.pylab import rcParams
import numpy as np
from astropy.cosmology import FlatLambdaCDM, z_at_value
import astropy.units as u


def read_SFH(sim_info, data):
    kpctoMpc = 1e-3
    file = sim_info.SFH_output
    data_file = np.loadtxt(file, usecols=[3,7])
    z = data_file[:,0]
    SFR = data_file[:,1] * 1.022690e-02  #Msun/yr
    SFH = SFR / (sim_info.boxSize * kpctoMpc) ** 3 #Msun/yr/Mpc^3
    counter = np.array([len(z)])

    if data == None:
        data = {'z': z, 'SFH': SFH, 'counter': counter}
    else:
        z_data = data['z']
        z_data = np.append(z_data, z)
        SFH_data = data['SFH']
        SFH_data = np.append(SFH_data, SFH)
        counter_data = data['counter']
        counter_data = np.append(counter_data, counter)
        data = {'z': z_data, 'SFH': SFH_data, 'counter':counter_data}
    return data

def plot_SFH_obs():
    h = 0.7

    # SFR Observational data from Hopkins 2004
    hcorr = np.log10(h) - np.log10(0.7)  # h^-2 for SFR, h^-3 for volume
    (z, z_up, z_down, lgrho, lgrho_up, lgrho_down) = np.loadtxt("./plotter/obs_data/sfr_hopkins2004_cor.dat",
                                                                unpack=True)
    lgrho = lgrho + np.log10(0.6) + hcorr
    obs1_a = 1. / (1. + z)
    obs1_rho = 10 ** lgrho
    obs1_rho_err = np.array([obs1_rho - 10 ** lgrho_down, 10 ** lgrho_up - obs1_rho])

    # SFR Observational data from Karim 2011
    (z, rho, err_up, err_down) = np.loadtxt("./plotter/obs_data/sfr_karim2011.dat", unpack=True)
    obs2_a = 1. / (1. + z)
    obs2_rho = rho * 0.6777 / 0.7
    obs2_rho_err = np.array([-err_down, err_up])

    # SFR Observational data from Bouwens 2012
    z, rhostar = np.loadtxt("./plotter/obs_data/bouwens_2012_sfrd_no_dust.txt", unpack=True)
    z, rhostar_dust = np.loadtxt("./plotter/obs_data/bouwens_2012_sfrd_dustcorr.txt", unpack=True)
    rhostar = (rhostar / 1.8) * 0.6777 / 0.7  # convert to Chabrier IMF from Salpeter and adjust for change in cosmology
    rhostar_dust = (rhostar_dust / 1.8) * 0.6777 / 0.7  # convert to Chabrier IMF from Salpeter and adjust for change in cosmology
    obs3_a = 1. / (1. + z)
    obs3_rho = rhostar
    obs3_rho_dust = rhostar_dust

    # SFR Observational data from Rodighierio 2012
    z, rhostar, err_m, err_p = np.loadtxt("./plotter/obs_data/sfr_rodighiero2010.dat", unpack=True)
    rhostar = (rhostar / 1.65) * 0.6777 / 0.7  # convert to Chabrier IMF from Salpeter and adjust for change in cosmology
    obs4_a = 1. / (1. + z)
    obs4_rho = rhostar
    obs4_rho_err = np.array([-err_m / 1.65, err_p / 1.65])

    # SFR Observational data from Cucciati 2012
    z, rhostar, err_m, err_p = np.loadtxt("./plotter/obs_data/sfr_cucciati2011.dat", unpack=True)
    rhostar = (rhostar - 2. * np.log10(0.7) + 2. * np.log10(0.6777))
    obs5_a = 1. / (1 + z)
    obs5_rho = 10 ** rhostar / 1.65
    obs5_rho_err = 10 ** np.array([rhostar + (err_m), rhostar + err_p]) / 1.65
    obs5_rho_err[0] = -obs5_rho_err[0] + obs5_rho
    obs5_rho_err[1] = -obs5_rho + obs5_rho_err[1]

    # SFR Observational data from Magnelli 2013
    z, rhostar, err_p, err_m = np.loadtxt("./plotter/obs_data/sfr_magnelli2013.dat", unpack=True)
    obs6_a = 1. / (1 + z)
    obs6_rho = 10 ** rhostar / 1.65  # convert to Chabrier IMF from Salpeter
    obs6_rho_err = 10 ** np.array([rhostar + (err_m), rhostar + err_p]) / 1.65
    obs6_rho_err[0] = -obs6_rho_err[0] + obs6_rho
    obs6_rho_err[1] = -obs6_rho + obs6_rho_err[1]

    # SFR Observational data from Gruppioni 2013
    z, rhostar, err_p, err_m = np.loadtxt("./plotter/obs_data/sfr_gruppioni2013.dat", unpack=True)
    obs7_a = 1. / (1 + z)
    obs7_rho = 10 ** rhostar / 1.65  # convert to Chabrier IMF from Salpeter
    obs7_rho_err = 10 ** np.array([rhostar + (err_m), rhostar + err_p]) / 1.65
    obs7_rho_err[0] = -obs7_rho_err[0] + obs7_rho
    obs7_rho_err[1] = -obs7_rho + obs7_rho_err[1]
    ##################

    # Observational data
    plt.errorbar(obs4_a, obs4_rho, yerr=obs4_rho_err, fmt='s', mec='0.3', color='0.3', markersize=4,
                 markeredgewidth=0.5, linewidth=0.5, mfc='w')#, label="${\\rm Rodighiero~et~al.~(2010)~(24\\mu m)}$")
    plt.errorbar(obs2_a, obs2_rho, yerr=obs2_rho_err, fmt='.', mec='0.3', color='0.3', markersize=7,
                 markeredgewidth=0.5, linewidth=0.5, mfc='w')#, label="${\\rm Karim~et~al.~(2011)~(radio)}$")
    plt.errorbar(obs5_a, obs5_rho, yerr=obs5_rho_err, fmt='^', mec='0.3', color='0.3', markersize=4,
                 markeredgewidth=0.5, linewidth=0.5, mfc='w')#, label="${\\rm Cucciati~et~al.~(2012)~(FUV)}$")
    plt.plot(obs3_a, obs3_rho_dust, 'd', mec='0.3', color='0.3', markersize=4, markeredgewidth=0.5, linewidth=1,
             mfc='w')#, label="${\\rm Bouwens~et~al.~(2012)~(UV,~no~dust)}$")
    plt.errorbar(obs6_a, obs6_rho, yerr=obs6_rho_err, fmt='>', mec='0.3', color='0.3', markersize=4,
                 markeredgewidth=0.5, linewidth=0.5, mfc='w')#, label="${\\rm Magnelli~et~al.~(2013)~(IR)}$")
    plt.errorbar(obs7_a, obs7_rho, yerr=obs7_rho_err, fmt='<', mec='0.3', color='0.3', markersize=4,
                 markeredgewidth=0.5, linewidth=0.5, mfc='w')#, label="${\\rm Gruppioni~et~al.~(2013)~(IR)}$")

def cosmic_sfr(z):
    phi = 0.01 * (1.+z)**2.6 / (1.0+((1.+z)/3.2)**(6.2))  #Msun Mpc^-3 yr^-1
    return phi


def Madau2014(z):
    return 0.015 * (1. + z)**2.7/(1. + ((1.+z)/2.9)**5.6)

def plot_SFH_models():
    redshift = np.arange(0,10,0.1)
    IMF = 1.65
    plt.plot(1. / (1+redshift), Madau2014(redshift) / IMF,'--',color='darkblue')
    plt.plot(1. / (1 + redshift), cosmic_sfr(redshift) / IMF,'--',color='darkgreen')


def plot_SFH(sim_data, output_name_list, output_path):

    z = sim_data['z']
    SFH = sim_data['SFH']
    counter = sim_data['counter']

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (4, 3),
        "figure.subplot.left": 0.18,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.15,
        "figure.subplot.top": 0.95,
        "figure.subplot.wspace": 0.3,
        "figure.subplot.hspace": 0.3,
        "lines.markersize": 0.5,
        "lines.linewidth": 0.2,
    }
    rcParams.update(params)

    plt.figure()

    # Box stellar abundance --------------------------------
    ax = plt.subplot(1, 1, 1)

    plot_SFH_obs()
    plot_SFH_models()

    count = 0
    color = ['tab:blue', 'tab:green', 'tab:orange', 'crimson', 'tab:purple']
    for i in range(len(output_name_list)):
        xm = z[count:count + counter[i]]
        xm = 1./ (1. + xm)
        ym = SFH[count:count + counter[i]]
        count += counter[i]
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

    plt.ylabel("SFR Density [M$_{\odot}$ yr$^{-1}$ Mpc$^{-3}$]", labelpad=2)
    plt.xlabel("Redshift ($z$)", labelpad=2)
    plt.yscale('log')

    redshift_ticks = np.array([0., 0.2, 0.5, 1., 2., 3., 5., 10., 20., 50., 100.])
    redshift_labels = ["$0$", "$0.2$", "$0.5$", "$1$", "$2$", "$3$", "$5$", "$10$", "$20$", "$50$", "$100$"]
    a_ticks = 1. / (redshift_ticks + 1.)
    plt.xticks(a_ticks, redshift_labels)
    plt.gca().tick_params(axis='x', which='minor', bottom=False)
    ax.tick_params(direction='in', axis='both', which='both')

    plt.xlim(1.02, 0.07)
    plt.ylim(1e-4, 1)

    plt.legend(loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               fontsize=8, columnspacing=0.02)

    plt.savefig(f"{output_path}/SFH_comparison.png", dpi=200)
