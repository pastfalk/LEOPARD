import numpy
import pylab

#for python 2.7


#######################

Nspecies=2

#######################


npara=numpy.zeros(Nspecies,dtype='i4')
nperp=numpy.zeros(Nspecies,dtype='i4')

vparamin=numpy.zeros(Nspecies)
vparamax=numpy.zeros(Nspecies)

vperpmin=numpy.zeros(Nspecies)
vperpmax=numpy.zeros(Nspecies)

dens=numpy.zeros(Nspecies)
mu=numpy.zeros(Nspecies)
beta_para=numpy.zeros(Nspecies)
beta_perp=numpy.zeros(Nspecies)
vdrift=numpy.zeros(Nspecies)


########################

#species 1

npara[0]=255
nperp[0]=64

dens[0]=1.0
mu[0]=1.0
beta_para[0]=4.0
beta_perp[0]=2.0
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
beta_para[1]=1.0
beta_perp[1]=1.0
vdrift[1]=-0.0

vparamin[1]=-260.0
vparamax[1]=260.0

vperpmin[1]=0.0
vperpmax[1]=260.0


#############################

limit=10.0**(-300)


def dist_bimax((vpar,vper), n,m,beta_par,beta_per,drift):
	bimax=numpy.exp(-n*(vpar-drift)**2/beta_par/m -n*vper**2/beta_per/m)* n**1.5 /(m**1.5 *numpy.pi**1.5 *beta_per*numpy.sqrt(beta_par))
	return bimax.ravel()


for ispecies in range(0,Nspecies):

	file_name='distribution'+numpy.str(ispecies+1)+'.dat'

	vpara = numpy.linspace(vparamin[ispecies],vparamax[ispecies], npara[ispecies])
	vperp = numpy.linspace(vperpmin[ispecies],vperpmax[ispecies], nperp[ispecies])

	vpara2,vperp2=numpy.meshgrid(vpara,vperp)

	data=dist_bimax((vpara2, vperp2),dens[ispecies],mu[ispecies],beta_para[ispecies],beta_perp[ispecies],vdrift[ispecies])
	data=data.reshape(nperp[ispecies],npara[ispecies])

	data_new=numpy.zeros((nperp[ispecies],npara[ispecies]))

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
