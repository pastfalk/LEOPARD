import numpy
import mpmath
import matplotlib.pyplot as plt

# Número de espécies
Nspecies = 2

# Valores de kappa para comparação
k = [50, 10, 6, 4]

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

# Configurações das espécies
npara[0] = 255
nperp[0] = 64
dens[0] = 1.0
mu[0] = 1.0
beta_para[0] = 4.0
beta_perp[0] = 2.0
theta_para[0] = (((k[0] - 1.5) / k[0]) * beta_para[0] * mu[0] / dens[0]) ** 0.5
theta_perp[0] = (((k[0] - 1.5) / k[0]) * beta_perp[0] * mu[0] / dens[0]) ** 0.5
vdrift[0] = -0.0
vparamin[0] = -5.0
vparamax[0] = 5.0
vperpmin[0] = 0.0
vperpmax[0] = 10.0

npara[1] = 255
nperp[1] = 64
dens[1] = 1.0
mu[1] = 1836.0
beta_para[1] = 1.0
beta_perp[1] = 1.0
theta_para[1] = (((k[1] - 1.5) / k[1]) * beta_para[1] * mu[1] / dens[1]) ** 0.5
theta_perp[1] = (((k[1] - 1.5) / k[1]) * beta_perp[1] * mu[1] / dens[1]) ** 0.5
vdrift[1] = -0.0
vparamin[1] = -260.0
vparamax[1] = 260.0
vperpmin[1] = 0.0
vperpmax[1] = 260.0

limit = 10.0 ** (-300)

def dist_bimax(vpar, vper, n, mu, beta_par, beta_per, drift):
    bimax = numpy.exp(-n * (vpar - drift) ** 2 / beta_par / mu - n * vper ** 2 / beta_per / mu) * n ** 1.5 / (mu ** 1.5 * numpy.pi ** 1.5 * beta_per * numpy.sqrt(beta_par))
    return bimax.ravel()

def dist_bikappa(vpar, vper, n, mu, kappa, theta_par, theta_per, drift):
    theta_par_per2 = numpy.outer(theta_par, theta_per ** 2)
    bikappa = (n / (kappa ** 1.5 * numpy.pi ** 1.5 * theta_par_per2)) * (mpmath.gamma(kappa + 1) / (mpmath.gamma(kappa - 0.5))) * (1 + (vpar - drift) ** 2 / (kappa * theta_par ** 2) + vper ** 2 / (kappa * theta_per ** 2)) ** (-kappa - 1)
    return bikappa.ravel()

index = 0  # Índice inicial para os arquivos

for ispecies in range(0, Nspecies):
    # Identificação da espécie: iarb
    iarb = str(ispecies + 1)  # Número da espécie de 1 a 9
    index = 0
    for kappa in k:
        # Ajuste o theta para o kappa atual
        theta_para[ispecies] = (((kappa - 1.5) / kappa) * beta_para[ispecies] * mu[ispecies] / dens[ispecies]) ** 0.5
        theta_perp[ispecies] = (((kappa - 1.5) / kappa) * beta_perp[ispecies] * mu[ispecies] / dens[ispecies]) ** 0.5
        
        # Atualize o índice do arquivo
        index += 1
        set_num = str(index).zfill(3)  # Preenche com zeros à esquerda para 3 dígitos
        
        # Nomes dos arquivos ajustados para incluir apenas o índice e o número da espécie
        file_name_bimax = f'distribution_{set_num}-{iarb}.dat'
        file_name_bikappa = f'distribution_{set_num}-{iarb}.dat'

        # Calcula a distribuição bi-maxwelliana
        vpara = numpy.linspace(vparamin[ispecies], vparamax[ispecies], npara[ispecies])
        vperp = numpy.linspace(vperpmin[ispecies], vperpmax[ispecies], nperp[ispecies])
        vpara2, vperp2 = numpy.meshgrid(vpara, vperp)

        databimax = dist_bimax(vpara2, vperp2, dens[ispecies], mu[ispecies], beta_para[ispecies], beta_perp[ispecies], vdrift[ispecies])

        databikappa = dist_bikappa(vpara2, vperp2, dens[ispecies], mu[ispecies], kappa, theta_para[ispecies], theta_perp[ispecies], vdrift[ispecies])
        databikappa = databikappa.reshape(nperp[ispecies], npara[ispecies])

        data_new_bikappa = numpy.zeros((nperp[ispecies], npara[ispecies]))

        for i in range(0, npara[ispecies]):
            for j in range(0, nperp[ispecies]):
                if databikappa[j, i] > limit:
                    data_new_bikappa[j, i] = databikappa[j, i] * dens[ispecies]
                else:
                    data_new_bikappa[j, i] = 0.0
        # Cálculo da integral da distribuição bi-kappa para verificar a normalização
        integral_bikappa = numpy.sum(data_new_bikappa) * (vparamax[ispecies] - vparamin[ispecies]) / npara[ispecies] * (vperpmax[ispecies] - vperpmin[ispecies]) / nperp[ispecies]

        # Imprimir o valor da integral e da densidade esperada para verificação
        print(f'Integral da distribuição bi-kappa para espécie {ispecies + 1} e kappa = {kappa}: {integral_bikappa:.6e}')
        print(f'Densidade esperada: {dens[ispecies]:.6e}')
        dat_fin_bimax = open(file_name_bimax, 'w')
        dat_fin_bikappa = open(file_name_bikappa, 'w')
        for i in range(0, npara[ispecies]):
            for j in range(0, nperp[ispecies]):
                dat_fin_bimax.write(f"{vpara[i]:.6e} {vperp[j]:.6e} {data_new_bikappa[j, i]:.6e}\n")
                dat_fin_bikappa.write(f"{vpara[i]:.6e} {vperp[j]:.6e} {data_new_bikappa[j, i]:.6e}\n")


        dat_fin_bimax.close()
        dat_fin_bikappa.close()

        # Plot da distribuição bi-kappa para o valor atual de kappa
        plt.plot(vpara, numpy.log10(data_new_bikappa[:, :].sum(axis=0)), label=rf'$\kappa={kappa}$')

    plt.xlabel(r'$v_\parallel$')
    plt.ylabel(r'$\log{f(v)}$')
    plt.title(f'Species {ispecies + 1} Distribution')
    plt.legend()
    plt.savefig(f'distribution_{set_num}-{iarb}.png')  # Salva o gráfico como imagem
    plt.show()
