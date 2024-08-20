import numpy
import mpmath
import matplotlib.pyplot as plt

# Número de espécies
Nspecies = 2

# Definição dos valores de kappa que você deseja comparar
k = [500, 50, 25, 6, 2]

npara = numpy.zeros(Nspecies, dtype='i4')
nperp = numpy.zeros(Nspecies, dtype='i4')

vparamin = numpy.zeros(Nspecies)
vparamax = numpy.zeros(Nspecies)

vperpmin = numpy.zeros(Nspecies)
vperpmax = numpy.zeros(Nspecies)

dens = numpy.zeros(Nspecies)
mu = numpy.zeros(Nspecies)

beta_para = numpy.zeros(Nspecies)
beta_perp = numpy.zeros(Nspecies)

theta_para = numpy.zeros(Nspecies)
theta_perp = numpy.zeros(Nspecies)
vdrift = numpy.zeros(Nspecies)

# Configurações para as duas espécies
npara[0], nperp[0], dens[0], mu[0], beta_para[0], beta_perp[0] = 255, 64, 1.0, 1.0, 4.0, 2.0
npara[1], nperp[1], dens[1], mu[1], beta_para[1], beta_perp[1] = 255, 64, 1.0, 1836.0, 1.0, 1.0

vparamin[0], vparamax[0], vperpmin[0], vperpmax[0] = -5.0, 5.0, 0.0, 10.0
vparamin[1], vparamax[1], vperpmin[1], vperpmax[1] = -500.0, 500.0, 0.0, 260.0

vdrift[:] = 0.0

limit = 10.0 ** (-300)

def dist_bimax(vpar, vper, n, mu, beta_par, beta_per, drift):
    bimax = numpy.exp(-n*(vpar-drift)**2/beta_par/mu -n*vper**2/beta_per/mu)* n**1.5 /(mu**1.5 *numpy.pi**1.5 *beta_per*numpy.sqrt(beta_par))
    return bimax.ravel()

def dist_bikappa(vpar, vper, n, mu, kappa,theta_par, theta_per, drift):
    theta_par_per2 = numpy.outer(theta_par, theta_per**2)
    bikappa = (n/(kappa**1.5*numpy.pi**1.5*theta_par_per2))*(mpmath.gamma(kappa+1)/(mpmath.gamma(kappa-0.5)))*(1+(vpar-drift)**2/(kappa*theta_par**2)+vper**2/(kappa*theta_per**2))**(-kappa-1)
    return bikappa.ravel()

for ispecies in range(0, Nspecies):
    # Calcula a distribuição bi-maxwelliana
    vpara = numpy.linspace(vparamin[ispecies], vparamax[ispecies], npara[ispecies])
    vperp = numpy.linspace(vperpmin[ispecies], vperpmax[ispecies], nperp[ispecies])

    vpara2, vperp2 = numpy.meshgrid(vpara, vperp)

    databimax = dist_bimax(vpara2, vperp2, dens[ispecies], mu[ispecies], beta_para[ispecies], beta_perp[ispecies], vdrift[ispecies])
    databimax = databimax.reshape(nperp[ispecies], npara[ispecies])

    data_new_bimax = numpy.zeros((nperp[ispecies], npara[ispecies]))

    for i in range(0, npara[ispecies]):
        for j in range(0, nperp[ispecies]):
            if databimax[j, i] > limit:
                data_new_bimax[j, i] = databimax[j, i] * dens[ispecies]
            else:
                data_new_bimax[j, i] = 0.0

    # Plot da distribuição bi-maxwelliana
    plt.plot(vpara, numpy.log10(data_new_bimax[:, :].sum(axis=0)), color='r', label='Maxwellian')

    for kappa in k:
        # Cálculo da distribuição bi-kappa
        theta_para[ispecies] = (((kappa - 1.5) / kappa) * beta_para[ispecies] * mu[ispecies] / dens[ispecies]) ** 0.5
        theta_perp[ispecies] = (((kappa - 1.5) / kappa) * beta_perp[ispecies] * mu[ispecies] / dens[ispecies]) ** 0.5
        
        databikappa = dist_bikappa(vpara2, vperp2, dens[ispecies], mu[ispecies], kappa, theta_para[ispecies], theta_perp[ispecies], vdrift[ispecies])
        databikappa = databikappa.reshape(nperp[ispecies], npara[ispecies])

        # Cálculo do valor médio de v_perp para cada v_para
        vperp_mean = numpy.zeros(npara[ispecies])
        for i in range(npara[ispecies]):
            # Calcula o valor médio de v_perp para cada valor de v_para
            vperp_mean[i] = numpy.sum(vperp * databikappa[:, i]) / numpy.sum(databikappa[:, i])
        
        # Plot do valor médio de v_perp para cada v_para
        plt.plot(vpara, vperp_mean, label=f'Kappa={kappa}')
        
    plt.xlabel('vpara')
    plt.ylabel('Average vperp')
    plt.title(f'Species {ispecies + 1} - Average vperp')
    plt.legend()
    plt.show()