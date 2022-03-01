import matplotlib.pylab as plt
from matplotlib.pylab import rcParams
import numpy as np
from astropy.cosmology import FlatLambdaCDM, z_at_value
import astropy.units as u
import scipy.integrate as sci
from tqdm import tqdm


def read_SNIa_rates(sim_info, data):

    file = sim_info.SNIa_output
    data_file = np.loadtxt(file, usecols=[6,11])
    z = data_file[:,0]
    NumSNIa = data_file[:,1] * 1.022690e-12 #yr^-1 Mpc^-3
    counter = np.array([len(z)])

    if data == None:
        data = {'z': z, 'NumSNIa': NumSNIa, 'counter': counter}
    else:
        z_data = data['z']
        z_data = np.append(z_data, z)
        NumSNIa_data = data['NumSNIa']
        NumSNIa_data = np.append(NumSNIa_data, NumSNIa)
        counter_data = data['counter']
        counter_data = np.append(counter_data, counter)
        data = {'z': z_data, 'NumSNIa': NumSNIa_data, 'counter':counter_data}
    return data

def plot_SNIa_obs():

    h = 0.6777

    z, ratenu, sys_err_p, sys_err_m, stat = np.loadtxt("./plotter/obs_data/SNIa_rate_frohmaier.dat", unpack=True)
    SNI_obs1_z = z
    SNI_obs1_a = 1. / (1. + z)
    err_p = (sys_err_p ** 2 + stat ** 2) ** .5 * (h / 0.7) ** 3
    err_m = (sys_err_m ** 2 + stat ** 2) ** .5 * (h / 0.7) ** 3
    SNI_obs1_rnu = ratenu * 1e-5 * (h / 0.7) ** 3
    SNI_obs1_err = np.zeros((2, 1))
    SNI_obs1_err[0] = err_m * 1e-5
    SNI_obs1_err[1] = err_p * 1e-5

    z, ratenu, sys_err_p, sys_err_m, stat_p, stat_m = np.loadtxt("./plotter/obs_data/SNIa_rate_dilday.dat",
                                                                 unpack=True)
    SNI_obs2_z = z
    SNI_obs2_a = 1. / (1. + z)
    err_p = (sys_err_p ** 2 + stat_p ** 2) ** .5 * (h / 0.7) ** 3
    err_m = (sys_err_m ** 2 + stat_m ** 2) ** .5 * (h / 0.7) ** 3
    SNI_obs2_rnu = ratenu * 1e-5 * (h / 0.7) ** 3
    SNI_obs2_err = np.zeros((2, len(SNI_obs2_rnu)))
    SNI_obs2_err[0] = err_m * 1e-5
    SNI_obs2_err[1] = err_p * 1e-5

    z, ratenu, sys_err_p, sys_err_m, stat_p, stat_m = np.loadtxt("./plotter/obs_data/SNIa_rate_perrett.dat",
                                                                 unpack=True)
    SNI_obs3_z = z
    SNI_obs3_a = 1. / (1. + z)
    err_p = (sys_err_p ** 2 + stat_p ** 2) ** .5 * (h / 0.7) ** 3
    err_m = (sys_err_m ** 2 + stat_m ** 2) ** .5 * (h / 0.7) ** 3
    SNI_obs3_rnu = ratenu * 1e-5 * (h / 0.7) ** 3
    SNI_obs3_err = np.zeros((2, len(SNI_obs3_rnu)))
    SNI_obs3_err[0] = err_m * 1e-5
    SNI_obs3_err[1] = err_p * 1e-5

    z, ratenu, sys_err_p, sys_err_m = np.loadtxt("./plotter/obs_data/SNIa_rate_graur.dat", unpack=True)
    SNI_obs4_z = z
    SNI_obs4_a = 1. / (1. + z)
    err_p = sys_err_p * (h / 0.7) ** 3
    err_m = -sys_err_m * (h / 0.7) ** 3
    SNI_obs4_rnu = ratenu * 1e-5 * (h / 0.7) ** 3
    SNI_obs4_err = np.zeros((2, len(SNI_obs4_rnu)))
    SNI_obs4_err[0] = err_m * 1e-5
    SNI_obs4_err[1] = err_p * 1e-5

    z, ratenu, sys_err_p, sys_err_m, stat_p, stat_m = np.loadtxt("./plotter/obs_data/SNIa_rate_rodney.dat",
                                                                 unpack=True)
    SNI_obs5_z = z
    SNI_obs5_a = 1. / (1. + z)
    err_p = (sys_err_p ** 2 + stat_p ** 2) ** .5 * (h / 0.7) ** 3
    err_m = (sys_err_m ** 2 + stat_m ** 2) ** .5 * (h / 0.7) ** 3
    SNI_obs5_rnu = ratenu * 1e-5 * (h / 0.7) ** 3
    SNI_obs5_err = np.zeros((2, len(SNI_obs5_rnu)))
    SNI_obs5_err[0] = err_m * 1e-5
    SNI_obs5_err[1] = err_p * 1e-5

    plt.errorbar(1 / (1. + SNI_obs1_z), SNI_obs1_rnu, yerr=SNI_obs1_err, fmt='s', mec='0.3', color='0.3', markersize=4,
                 markeredgewidth=0.5, linewidth=0.5, mfc='w')  # , label="${\\rm Frohmaier~et~al.~(2019)~(PTF)}$")
    plt.errorbar(1 / (1. + SNI_obs2_z), SNI_obs2_rnu, yerr=SNI_obs2_err, fmt='.', mec='0.3', color='0.3', markersize=4,
                 markeredgewidth=0.5, linewidth=0.5, mfc='w')  # , label="${\\rm Dilday~et~al.~(2010)~(SDSS)}$")
    plt.errorbar(1 / (1. + SNI_obs3_z), SNI_obs3_rnu, yerr=SNI_obs3_err, fmt='^', mec='0.3', color='0.3', markersize=4,
                 markeredgewidth=0.5, linewidth=0.5, mfc='w')  # , label="${\\rm Perrett~et~al.~(2013)~(SNLS)}$")
    plt.errorbar(1 / (1. + SNI_obs4_z), SNI_obs4_rnu, yerr=SNI_obs4_err, fmt='<', mec='0.3', color='0.3', markersize=4,
                 markeredgewidth=0.5, linewidth=0.5, mfc='w')  # , label="${\\rm Graur~et~al.~(2011)~(SDF)}$")
    plt.errorbar(1 / (1. + SNI_obs5_z), SNI_obs5_rnu, yerr=SNI_obs5_err, fmt='>', mec='0.3', color='0.3', markersize=4,
                 markeredgewidth=0.5, linewidth=0.5, mfc='w')  # , label="${\\rm Rodney~et~al.~(2014)~(CANDELS)}$")

def prefacpow(beta,tau,tstart,tend):
    firstpar = 1 - beta
    secondpar = (tend/tau)**(1.-beta) - (tstart/tau)**(1.-beta)
    thirdpar = tau**beta
    return firstpar / secondpar * thirdpar

def Madau2014(z):
    return 0.015 * (1. + z)**2.7/(1. + ((1.+z)/2.9)**5.6)

def DTD(t,beta):
    return t**(-beta)

def DTDexp(t,beta):
    return np.exp(-t/beta)

def function4(tau,t,beta):
    cosmo = FlatLambdaCDM(H0=67.77 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.307)
    return Madau2014(z_at_value(cosmo.age,(t-tau)*u.Gyr))*DTDexp(tau,beta)*np.heaviside(tau-0.04,1.)

def function2(tau,t,beta):
    cosmo = FlatLambdaCDM(H0=67.77 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.307)
    return Madau2014(z_at_value(cosmo.age,(t-tau)*u.Gyr))*DTD(tau,beta)*np.heaviside(tau-0.04,1.)

def plot_SNIa_models():

    cosmo = FlatLambdaCDM(H0=67.77 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.307)
    N = 40
    time = np.linspace(0.25, 13.4617, N)
    values = np.zeros((len(time), 15))
    zvalues = np.zeros(len(time))
    t_init = 0.040
    t_H = 13.6
    IMF = 1.65 # IMF correction factor

    # Calculate the functional form of the SNIa rate
    for i in tqdm(range(0, len(values[:, 0]))):
        # Let's do the integral for the power law function
        #integral = sci.quad(function2, t_init, time[i], args=(time[i], 1.0), epsrel=1e-4)[0]

        # tau = 2.0
        # beta = 1.0394117647058825
        # nu = 0.0011818181818181819
        # values[i, 0] = integral * prefacpow(beta, tau, t_init, t_H) * nu / tau / IMF

        beta = 0.5
        nu = 0.001
        tau = 1.0
        integral = sci.quad(function2, t_init, time[i], args=(time[i], beta), epsrel=1e-4)[0]
        values[i, 1] = integral * prefacpow(beta, tau, t_init, t_H) * nu / tau / IMF

        tau = 2.0
        nu = 2e-3
        # Let's do the integral for the exponential function
        values[i, 2] = sci.quad(function4, t_init, time[i], args=(time[i], tau), epsrel=1e-4)[0] * nu / tau / IMF

        # beta = 2
        # nu = 0.006
        # tau = 1.0
        # values[i, 3] = integral * prefacpow(beta, tau, t_init, t_H) * nu / tau / IMF

        beta = 0.8
        nu = 0.001
        tau = 1.0
        integral = sci.quad(function2, t_init, time[i], args=(time[i], beta), epsrel=1e-4)[0]
        values[i, 4] = integral * prefacpow(beta, tau, t_init, t_H) * nu / tau / IMF
        zvalues[i] = z_at_value(cosmo.age, time[i] * u.Gyr)

    plt.plot(1. / (1. + zvalues), values[:, 2], label='EAGLE model ($\\tau = 2, \\nu=2\cdot 10^{-3}$)', color='black',
             linewidth=0.6, linestyle='-.')
    # plt.plot(1. / (1. + zvalues), values[:, 3], label='Power law ($\\beta=2, \\tau=1, \\nu=6\cdot 10^{-3}$)',
    #          color='tab:red', linewidth=0.6, linestyle='-.')
    # plt.plot(1. / (1. + zvalues), values[:, 0], label='Power law ($\\beta=1.03, \\tau=2, \\nu=1.3\cdot 10^{-3}$)',
    #          color='tab:orange', linewidth=0.6, linestyle='-.')
    plt.plot(1. / (1. + zvalues), values[:, 4], label='Power law ($\\beta=0.8, \\tau=1, \\nu=1\cdot 10^{-3}$)',
             color='darkgreen', linewidth=0.6, linestyle='-.')
    plt.plot(1. / (1. + zvalues), values[:, 1], label='Power law ($\\beta=0.5, \\tau=1, \\nu=1\cdot 10^{-3}$)',
             color='darkblue', linewidth=0.6, linestyle='-.')


def plot_SNIa_rates(sim_data, output_name_list, output_path):


    z = sim_data['z']
    NumSNIa = sim_data['NumSNIa']
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

    plot_SNIa_obs()
    plot_SNIa_models()

    count = 0
    color = ['tab:blue', 'tab:orange', 'crimson', 'tab:green']
    for i in range(len(output_name_list)):
        xm = z[count:count + counter[i]]
        xm = 1./ (1. + xm)
        ym = NumSNIa[count:count + counter[i]]
        count += counter[i]
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

    plt.ylabel("SNIa Rate [yr$^{-1}$ Mpc$^{-3}$]", labelpad=2)
    plt.xlabel("Redshift ($z$)", labelpad=2)
    plt.yscale('log')

    redshift_ticks = np.array([0., 0.2, 0.5, 1., 2., 3., 5., 10., 20., 50., 100.])
    redshift_labels = ["$0$", "$0.2$", "$0.5$", "$1$", "$2$", "$3$", "$5$", "$10$", "$20$", "$50$", "$100$"]
    a_ticks = 1. / (redshift_ticks + 1.)
    plt.xticks(a_ticks, redshift_labels)
    plt.gca().tick_params(axis='x', which='minor', bottom=False)
    ax.tick_params(direction='in', axis='both', which='both')

    plt.xlim(1.02, 0.07)
    plt.ylim(1e-5, 1e-3)

    # handles, labels = plt.gca().get_legend_handles_labels()
    # order = np.arange(len(handles)-2)+2
    # order = np.append(order,0)
    # order = np.append(order,1)

    # plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order],
    #            loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
    #            fontsize=9, columnspacing=0.02)
    plt.legend(loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               fontsize=8, columnspacing=0.02)

    plt.savefig(f"{output_path}/SNIa_rates_comparison.png", dpi=200)



