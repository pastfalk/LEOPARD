import numpy
import pylab

#for python 2.7

file_name='distribution1.dat'

npara=255
nperp=128

limit=10.0**(-1776)

dens=1.0
beta_par=4.0
beta_per=2.0
vdrift=-0.0

vparamin=-12.0
vparamax=12.0

vperpmin=0.0
vperpmax=10.0


vpara = numpy.linspace(vparamin,vparamax, npara)
vperp = numpy.linspace(vperpmin,vperpmax, nperp)

vpara2,vperp2=numpy.meshgrid(vpara,vperp)


def dist_bimax((vpar,vper), beta_para,beta_perp,vdrift):

	bimax=numpy.exp(-(vpar-vdrift)**2/beta_para -vper**2/beta_perp)   /(numpy.pi**1.5 *beta_perp*numpy.sqrt(beta_para))

	return bimax.ravel()

data=dist_bimax((vpara2, vperp2),beta_par,beta_per,vdrift)
data=data.reshape(nperp,npara)

data_new=numpy.zeros((nperp,npara))

for i in range(0,npara):
	for j in range(0,nperp):

		if(data[j,i]>limit):
			data_new[j,i]=data[j,i]*dens
		else:
			data_new[j,i]=0.0
			

dat_fin = open(file_name, 'w')

for i in range(0,npara):
	for j in range(0,nperp):
		dat_fin.write(str(vpara[i]))
		dat_fin.write(" ")
		dat_fin.write(str(vperp[j]))
		dat_fin.write(" ")
		dat_fin.write(str(data_new[j,i]))
		dat_fin.write("\n")
