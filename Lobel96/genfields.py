""" GENFIELDS Generation fields routine
   This script implements the data generation for the Gradient Conjugated
   Method (Lobel et al., 1996). It computes the incident, total and
   scattered field as well as the Green function for both domains.

   Implemented by:
 
   Andre Costa Batista
   Universidade Federal de Minas Gerais

   REFERENCES

   Lobel, P., et al. "Conjugate gradient method for solving inverse
   scattering with experimental data." IEEE Antennas and Propagation
   Magazine 38.3 (1996): 48-51.
"""

# Importing general libraries
import numpy as np
import copy as cp
import pickle
from scipy import special as spl
from scipy.spatial import distance as dst

# Importing FDTD2D library
import fdtd2d as slv
import domain as dm
import boundary as bd
import waveform as wv
import source as sc
import probe as pb

# Model parameters
expname = 'basic'
I, J = 50, 50   # Number of cells in x,y-axis
N = I*J
dx, dy = 5e-3, 5e-3
epsrb, sigb = 1., 0.
sampled_frequencies = np.array([800e6])
wv_frequency = 800e6
dtheta = 12
M = round(360/dtheta)
Rs = 14e-2
ls_x, ls_y = 4*dx, 4*dy
magnitude = 1e1
time_window = 5e-8
lambda_b = 1/sampled_frequencies/np.sqrt(slv.mu_0*epsrb*slv.eps_0)
kb = 2*np.pi*sampled_frequencies*np.sqrt(slv.mu_0*slv.eps_0*epsrb)

Rx = pb.get_rx_array(Rs,dtheta,sampled_frequencies.size)
Tx = sc.Source(sc.get_tx_position(Rs,dtheta,0),
               wv.GaussianSignal(dx,dy,wv_frequency),magnitude)

mymodel = slv.FDTD2D(dm.Domain(dx,dy,I,J),epsrb,sigb,sampled_frequencies,ls_x,
                     ls_y,time_window,probes=Rx)

ei = np.zeros((I,J,M),dtype=complex)
es = np.zeros((M,M),dtype=complex)

for i in range(M):
   Tx = sc.Source(sc.get_tx_position(Rs,dtheta,i),
               wv.GaussianSignal(dx,dy,wv_frequency),magnitude)
   mymodel.run(Tx)
   ei[:,:,i] = np.squeeze(mymodel.get_intern_field())
   for j in range(M):
      es[j,i] = mymodel.probes[j].get_field_freq(0)

epsr = epsrb*np.ones((I,J))
sig = sigb*np.ones((I,J))
epsr[np.ix_(range(round(I/4)-round(.1*I),round(I/4)+round(.1*I)),
            range(round(J/4)-round(.1*J),round(J/4)+round(.1*J)))] = 5.0

et = np.zeros((I,J,M),dtype=complex)

for i in range(M):
   Tx = sc.Source(sc.get_tx_position(Rs,dtheta,i),
               wv.GaussianSignal(dx,dy,wv_frequency),magnitude)
   mymodel.run(Tx,epsr=epsr,sig=sig)
   et[:,:,i] = np.squeeze(mymodel.get_intern_field())
   for j in range(M):
      es[j,i] = mymodel.probes[j].get_field_freq(0)-es[j,i]

et = et.reshape((N,M))
ei = ei.reshape((N,M))
x, y = mymodel.get_intern_coordinates()
rx_xy = np.zeros((M,2))
for i in range(M):
   rx_xy[i,0], rx_xy[i,1] = Rx[i].position[0], Rx[i].position[1]

deltasn = dx*dy
an = np.sqrt(deltasn/np.pi)

xn, yn = np.meshgrid(x,y,indexing='ij')
xn, yn = xn.reshape(N), yn.reshape(N)
R = dst.cdist(rx_xy,np.stack((xn,yn),axis=1),'euclidean')

gs = np.zeros(R.shape,dtype=complex)
gs[R!=0] = 1j/2*np.pi*kb*an*spl.jv(1,kb*an)*spl.hankel2(0,kb*R[R!=0])
gs[R==0] = 1j/2*(np.pi*kb*an*spl.hankel2(1,kb*an)-2j)
gs = -gs

R = dst.cdist(np.stack((xn,yn),axis=1),np.stack((xn,yn),axis=1),'euclidean')
gd = np.zeros(R.shape,dtype=complex)
gd[R!=0] = 1j/2*np.pi*kb*an*spl.jv(1,kb*an)*spl.hankel2(0,kb*R[R!=0])
gd[R==0] = 1j/2*(np.pi*kb*an*spl.hankel2(1,kb*an)-2j)
gd = -gd

data = {
   'dx':dx,
   'dy':dy,
   'I':I,
   'J':J,
   'epsr':epsr,
   'sig':sig,
   'epsrb':epsrb,
   'sigb':sigb,
   'frequency':sampled_frequencies,
   'waveform_frequency':wv_frequency,
   'dtheta':dtheta,
   'Rs':Rs,
   'ls_x':ls_x,
   'ls_y':ls_y,
   'magnitude':magnitude,
   'time_window':time_window,
   'Rx':Rx,
   'kb':kb,
   'lambda_b':lambda_b,
   'ei':ei,
   'et':et,
   'es':es,
   'x':x,
   'y':y,
   'gs':gs,
   'gd':gd
}

with open(expname,'wb') as datafile:
   pickle.dump(data,datafile)