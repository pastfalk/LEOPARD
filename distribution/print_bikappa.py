import numpy as np
import mpmath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Número de espécies
Nspecies = 2

# Valores de kappa para comparação
k = [2, 5, 10, 50]

npara = np.zeros(Nspecies, dtype='i4')
nperp = np.zeros(Nspecies, dtype='i4')

vparamin = np.zeros(Nspecies)
vparamax = np.zeros(Nspecies)

vperpmin = np.zeros(Nspecies)
vperpmax = np.zeros(Nspecies)

dens = np.zeros(Nspecies)
mu = np.zeros(Nspecies)

beta_para = np.zeros(Nspecies)
beta_perp = np.zeros(Nspecies)

theta_para = np.zeros(Nspecies)
theta_perp = np.zeros(Nspecies)
vdrift = np.zeros(Nspecies)

# Configurações das espécies
npara[0] = 255
nperp[0] = 64
dens[0] = 1.0
mu[0] = 1.0
beta_para[0] = 1.0
beta_perp[0] = 1.0
theta_para[0] = (((k[0] - 1.5) / k[0]) * beta_para[0] * mu[0] / dens[0]) ** 0.5
theta_perp[0] = (((k[0] - 1.5) / k[0]) * beta_perp[0] * mu[0] / dens[0]) ** 0.5
vdrift[0] = -0.0
vparamin[0] = -12.0
vparamax[0] = 12.0
vperpmin[0] = 0.0
vperpmax[0] = 10.0

npara[1] = 255
nperp[1] = 64
dens[1] = 1.0
mu[1] = 1836.0
beta_para[1] = 20.0
beta_perp[1] = 1.0
theta_para[1] = (((k[0] - 1.5) / k[0]) * beta_para[1] * mu[1] / dens[1]) ** 0.5
theta_perp[1] = (((k[0] - 1.5) / k[0]) * beta_perp[1] * mu[1] / dens[1]) ** 0.5
vdrift[1] = -0.0
vparamin[1] = -260.0
vparamax[1] = 260.0
vperpmin[1] = 0.0
vperpmax[1] = 260.0

# npara[2] = 255
# nperp[2] = 64
# dens[2] = 1.0
# mu[2] = 1836.0
# beta_para[2] = 16.0
# beta_perp[2] = 1.0
# theta_para[2] = (((k[0] - 1.5) / k[0]) * beta_para[2] * mu[2] / dens[2]) ** 0.5
# theta_perp[2] = (((k[0] - 1.5) / k[0]) * beta_perp[2] * mu[2] / dens[2]) ** 0.5
# vdrift[2] = -0.0
# vparamin[2] = -260.0
# vparamax[2] = 260.0
# vperpmin[2] = 0.0
# vperpmax[2] = 260.0

# npara[3] = 255
# nperp[3] = 64
# dens[3] = 1.0
# mu[3] = 1836.0
# beta_para[3] = 18.0
# beta_perp[3] = 1.0
# theta_para[3] = (((k[0] - 1.5) / k[0]) * beta_para[3] * mu[3] / dens[3]) ** 0.5
# theta_perp[3] = (((k[0] - 1.5) / k[0]) * beta_perp[3] * mu[3] / dens[3]) ** 0.5
# vdrift[3] = -0.0
# vparamin[3] = -260.0
# vparamax[3] = 260.0
# vperpmin[3] = 0.0
# vperpmax[3] = 260.0

# npara[4] = 255
# nperp[4] = 64
# dens[4] = 1.0
# mu[4] = 1836.0
# beta_para[4] = 20.0
# beta_perp[4] = 1.0
# theta_para[4] = (((k[0] - 1.5) / k[0]) * beta_para[4] * mu[4] / dens[4]) ** 0.5
# theta_perp[4] = (((k[0] - 1.5) / k[0]) * beta_perp[4] * mu[4] / dens[4]) ** 0.5
# vdrift[4] = -0.0
# vparamin[4] = -260.0
# vparamax[4] = 260.0
# vperpmin[4] = 0.0
# vperpmax[4] = 260.0

# npara[5] = 255
# nperp[5] = 64
# dens[5] = 1.0
# mu[5] = 1836.0
# beta_para[5] = 30.0
# beta_perp[5] = 1.0
# theta_para[5] = (((k[0] - 1.5) / k[0]) * beta_para[5] * mu[5] / dens[5]) ** 0.5
# theta_perp[5] = (((k[0] - 1.5) / k[0]) * beta_perp[5] * mu[5] / dens[5]) ** 0.5
# vdrift[5] = -0.0
# vparamin[5] = -260.0
# vparamax[5] = 260.0
# vperpmin[5] = 0.0
# vperpmax[5] = 260.0

limit = 10.0 ** (-20)

def dist_bimax(vpar, vper, n, mu, beta_par, beta_per, drift):
    bimax = np.exp(-n * (vpar - drift) ** 2 / beta_par / mu - n * vper ** 2 / beta_per / mu) * n ** 1.5 / (mu ** 1.5 * np.pi ** 1.5 * beta_per * np.sqrt(beta_par))
    return bimax.ravel()

def dist_bikappa(vpar, vper, n, mu, kappa, theta_par, theta_per, drift):
    theta_par_per2 = np.outer(theta_par, theta_per ** 2)
    bikappa = (n / (kappa ** 1.5 * np.pi ** 1.5 * theta_par_per2)) * (mpmath.gamma(kappa + 1) / (mpmath.gamma(kappa - 0.5))) * (1 + (vpar - drift) ** 2 / (kappa * theta_par ** 2) + vper ** 2 / (kappa * theta_per ** 2)) ** (-kappa - 1)
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
        file_name_bikappa = f'distribution{set_num}-{iarb}.dat'

        # Calcula a distribuição bi-maxwelliana
        vpara = np.linspace(vparamin[ispecies], vparamax[ispecies], npara[ispecies])
        vperp = np.linspace(vperpmin[ispecies], vperpmax[ispecies], nperp[ispecies])
        vpara2, vperp2 = np.meshgrid(vpara, vperp)

        databimax = dist_bimax(vpara2, vperp2, dens[ispecies], mu[ispecies], beta_para[ispecies], beta_perp[ispecies], vdrift[ispecies])

        databikappa = dist_bikappa(vpara2, vperp2, dens[ispecies], mu[ispecies], kappa, theta_para[ispecies], theta_perp[ispecies], vdrift[ispecies])
        databikappa = databikappa.reshape(nperp[ispecies], npara[ispecies])

        data_new_bikappa = np.zeros((nperp[ispecies], npara[ispecies]))

        for i in range(0, npara[ispecies]):
            for j in range(0, nperp[ispecies]):
                if databikappa[j, i] > limit:
                    data_new_bikappa[j, i] = databikappa[j, i] * dens[ispecies]
                else:
                    data_new_bikappa[j, i] = 0.0
        
        dat_fin_bikappa = open(file_name_bikappa, 'w')
        for i in range(0, npara[ispecies]):
            for j in range(0, nperp[ispecies]):
                dat_fin_bikappa.write(f"{vpara[i]:.6e} {vperp[j]:.6e} {data_new_bikappa[j, i]:.6e}\n")

        dat_fin_bikappa.close()

        # Plot da distribuição bi-kappa para o valor atual de kappa
        vperp_fix=0
        idx_vperp = (np.abs(vperp - vperp_fix)).argmin()
        plt.plot(vpara, np.log10(data_new_bikappa[idx_vperp, :]), label=rf'$\kappa={kappa}$, $\beta_\parallel={beta_para[ispecies]}$')

        #plt.plot(vpara, np.log10(data_new_bikappa[:, :].sum(axis=0)), label=rf'$\kappa={kappa}$, $\beta_\parallel={beta_para[1]}$')

    plt.xlabel(r'$v_\parallel$')
    plt.ylabel(r'$\log{f(v)}$')
    if ispecies == 0:
        plt.title(f'Distribuição de íons para $v_\perp = {vperp_fix}$')
        plt.legend()
        plt.savefig(f'distribution_íons_vperp{vperp_fix}.png')  # Salva o gráfico como imagem
         # Configuração do gráfico 3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Criação da superfície 3D
        vpara2, vperp2 = np.meshgrid(vpara, vperp)
        ax.plot_surface(vpara2, vperp2, np.log10(data_new_bikappa), cmap='viridis')

        # Rótulos dos eixos
        ax.set_xlabel(r'$v_\parallel$')
        ax.set_ylabel(r'$v_\perp$')
        ax.set_zlabel(r'$\log{f(v_\parallel, v_\perp)}$')
        # Definir os limites dos eixos
        # plt.xlim(vparamin[ispecies], vparamax[ispecies])
        # plt.ylim(np.min(np.log10(data_new_bikappa)), np.max(np.log10(data_new_bikappa)))
        plt.title(f'Distribuição de íons para $v_\perp = {vperp_fix}$ e $\kappa={kappa}$')
    else:
        # Definir os limites dos eixos
        # plt.xlim(vparamin[ispecies], vparamax[ispecies])
        # plt.ylim(np.min(np.log10(data_new_bikappa)), np.max(np.log10(data_new_bikappa)))
        plt.title(f'Distribuição de elétrons para $v_\perp = {vperp_fix}$')
        plt.legend()
        plt.savefig(f'distribution_elétrons_beta{beta_para[1]}_vperp{vperp_fix}.png')  # Salva o gráfico como imagem
         # Configuração do gráfico 3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Criação da superfície 3D
        vpara2, vperp2 = np.meshgrid(vpara, vperp)
        ax.plot_surface(vpara2, vperp2, np.log10(data_new_bikappa), cmap='viridis')

        # Rótulos dos eixos
        ax.set_xlabel(r'$v_\parallel$')
        ax.set_ylabel(r'$v_\perp$')
        ax.set_zlabel(r'$\log{f(v_\parallel, v_\perp)}$')

        plt.title(f'Distribuição bi-Kappa para $\kappa={kappa}$')
    plt.show()
