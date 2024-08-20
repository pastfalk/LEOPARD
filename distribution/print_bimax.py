import numpy as np

#for python 2.7


#######################

Nspecies=2

#######################


npara=np.zeros(Nspecies,dtype='i4')
nperp=np.zeros(Nspecies,dtype='i4')

vparamin=np.zeros(Nspecies)
vparamax=np.zeros(Nspecies)

vperpmin=np.zeros(Nspecies)
vperpmax=np.zeros(Nspecies)

dens=np.zeros(Nspecies)
mu=np.zeros(Nspecies)
beta_para=np.zeros(Nspecies)
beta_perp=np.zeros(Nspecies)
vdrift=np.zeros(Nspecies)


########################

#species 1
ciarb=1
npara[0]=255
nperp[0]=64

dens[0]=1.0
mu[0]=1.0
beta_para[0]=1.0
beta_perp[0]=1.0
vdrift[0]=-0.0

vparamin[0]=-12.0
vparamax[0]=12.0

vperpmin[0]=0.0
vperpmax[0]=10.0

#species 2

npara[1]=255
nperp[1]=64

dens[1]=1.0
mu[1]=1836.0
beta_para[1]=20.0
beta_perp[1]=1.0
vdrift[1]=-0.0

vparamin[1]=-260.0
vparamax[1]=260.0

vperpmin[1]=0.0
vperpmax[1]=260.0 

#############################

limit=10.0**(-300)


def dist_bimax(vpar,vper, n,m,beta_par,beta_per,drift):
	bimax=np.exp(-n*(vpar-drift)**2/beta_par/m -n*vper**2/beta_per/m)* n**1.5 /(m**1.5 *np.pi**1.5 *beta_per*np.sqrt(beta_par))
	return bimax.ravel()

def dist_bikappa(vpar,vper, n,m,beta_par,beta_per,drift):
	bimax=np.exp(-n*(vpar-drift)**2/beta_par/m -n*vper**2/beta_per/m)* n**1.5 /(m**1.5 *np.pi**1.5 *beta_per*np.sqrt(beta_par))
	return bimax.ravel()


for ispecies in range(0,Nspecies):

	file_name='distribution00'+str(ispecies+1)+'-'+str(ciarb)+'.dat'

	vpara = np.linspace(vparamin[ispecies],vparamax[ispecies], npara[ispecies])
	vperp = np.linspace(vperpmin[ispecies],vperpmax[ispecies], nperp[ispecies])

	vpara2,vperp2=np.meshgrid(vpara,vperp)

	data=dist_bimax(vpara2, vperp2,dens[ispecies],mu[ispecies],beta_para[ispecies],beta_perp[ispecies],vdrift[ispecies])
	data=data.reshape(nperp[ispecies],npara[ispecies])

	data_new=np.zeros((nperp[ispecies],npara[ispecies]))

	for i in range(0,npara[ispecies]):
		for j in range(0,nperp[ispecies]):

			if(data[j,i]>limit):
				data_new[j,i]=data[j,i]*dens[ispecies]
			else:
				data_new[j,i]=0.0


	dat_fin = open(file_name, 'w')

	for i in range(0,npara[ispecies]):
		for j in range(0,nperp[ispecies]):
			dat_fin.write(str(vpara[i]))
			dat_fin.write(" ")
			dat_fin.write(str(vperp[j]))
			dat_fin.write(" ")
			dat_fin.write(str(data_new[j,i]))
			dat_fin.write("\n")
