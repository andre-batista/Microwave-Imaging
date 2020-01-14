'''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   The Conjugated Gradient Method                        %
%                                                                         %
%       This scripts implements the Conjugated Gradient for inverse       %
%   scattering 2D TMz electromagnetic problems (Lobel et al., 1996).      %
%   Given the measurements of the scattering field in specific points of  %
%   a domain denoted by S, the incident field on a investigation domain D %
%   and the Green Function for both domains, the method recovers the      %
%   dielectric distribution within the region D.                          %
%                                                                         %
%   Inputs:                                                               %
%   - es: a M by L matrix with the M measured scattered fields for the L  %
%       sources [V/m]                                                     %
%   - ei: a N by L matrix with the N computed incident fields for the L   %
%       sources [V/m]                                                     %
%   - gd: a N by N matrix with Green function computed for all of the N   %
%       points of the mesh in respect to each of them                     %
%   - gs: a M by N matrix with Green function computed for all of the N   %
%       points of the mesh in respect to the M scattered field            %
%       measurements                                                      %
%                                                                         %
%   Data struct:                                                          %
%   - dx, dy: cell sizes [m]                                              %
%   - epsr, sig: correct information of the dielectric distribution of    %
%       the experiment (relative permittivity and conductivity [S/m])     %
%   - epsrb, sigb: relative permittivity and conductivity [S/m] of the    %
%       background                                                        %
%   - lambdab: wavelength of the background [m]                           %
%   - f: linear frequency of measurements [Hz]                            %
%                                                                         %
%   Output variables:                                                     %
%   - epsr: retrieved relative permittivity                               %
%   - sig: retrieved conductivity [S/m]                                   %
%                                                                         %
%   Implemented by:                                                       %
%                                                                         %
%   Andre Costa Batista                                                   %
%   Universidade Federal de Minas Gerais                                  %
%                                                                         %
%   References                                                            %
%                                                                         %
%   Lobel, P., et al. "Conjugate gradient method for solving inverse      %
%   scattering with experimental data." IEEE Antennas and Propagation     %
%   Magazine 38.3 (1996): 48-51.                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''

import numpy as np
import copy as cp
import pickle
import time
from scipy import sparse as sps
from numpy import linalg as lag
import matplotlib.pyplot as plt

def inner(v1,v2,d):
    return np.sum(d*v1*np.conj(v2))

def gsmethod(rho,v):

    gmean       = 0.618
    delta       = 1.0e-3
    [a, b]      = getinterval(rho,v)   
    xa          = b - gmean*(b-a)
    xb          = a + gmean*(b-a)
    fxa         = lag.norm(np.reshape(rho-xa*v,(-1,1)))**2
    fxb         = lag.norm(np.reshape(rho-xb*v,(-1,1)))**2
    
    while (b - a) > delta:
        if fxa > fxb:
            a = xa
            xa = xb
            fxa = fxb
            xb = a + gmean*(b - a)
            fxb = lag.norm(np.reshape(rho-xb*v,(-1,1)))**2
        else:
            b = xb
            xb = xa
            fxb = fxa
            xa = b - gmean*(b-a)
            fxa = lag.norm(np.reshape(rho-xa*v,(-1,1)))**2

    alpha = (b+a)/2
    return alpha
    
def getinterval(rho,v):
    
    step0           = 1.0e-03
    a               = 0
    Fa              = lag.norm(rho.reshape(-1))**2
    b               = step0
    Fb              = lag.norm(np.reshape(rho-step0*v,(-1,1)))**2
    stepsize        = step0
    acceleration    = 2
    
    while (Fa > Fb):
        stepsize = acceleration*stepsize
        b = b + stepsize
        Fb = lag.norm(np.reshape(rho-b*v,(-1,1)))**2
    return a,b

print('========== The Conjugated Gradient Method ==========')
expname = 'basic'

with open(expname,'rb') as datafile:
    data = pickle.load(datafile)

# Loading inputs
dx, dy = data['dx'], data['dy']
I, J = data['I'], data['J']
epsrb, sigb = data['epsrb'], data['sigb']
f = data['frequency']
kb, lambda_b = data['kb'], data['lambda_b']
ei, et, es = data['ei'], data['et'], data['es']
x, y = data['x'], data['y']
gs, gd = data['gs'], data['gd']
epsr, sig = data['epsr'], data['sig']

# General Parameters
maxit = 5            # Number of iterations
M, L = es.shape         # M measurements, L sources
N = ei.shape[0]         # N points within the mesh
dS = dx*dy              # Surface element [m^2]
eps0 = 8.85418782e-12   # Vaccum permittivity [F/m]
omega = 2*np.pi*f       # Angular frequency [rad/sec]

# How do you preffer the initial solution?
# 1 - Everything background
# 2 - Backpropagation method (Lobel et al., 1996)
# 3 - Exact solution
# 4 - Load last run
initopt = 2

if initopt is 1:
    C = sps.dia_matrix((N,N),dtype=complex)
    d = np.zeros((N,1),dtype=complex)
    g = np.ones((N,1),dtype=complex)

elif initopt is 2:
    gamma = lag.norm(np.reshape(gs.conj().T@es,(-1,1)))**2/lag.norm(np.reshape(gs@gs.conj().T@es,(-1,1)))**2
    w0 = gamma*gs.conj().T@es
    C = sps.dia_matrix(np.diag(1/L*np.sum(w0/ei,1)),dtype=complex)
    d = np.zeros((N,1),dtype=complex)
    g = np.ones((N,1),dtype=complex)

elif initopt is 3:
    C = sps.dia_matrix(np.diag(np.reshape((epsr-1j*sig/omega/eps0/epsrb)-(epsrb-1j*sigb/omega/eps0/epsrb),-1)),dtype=complex)
    d = np.zeros((N,1))
    g = np.ones((N,1))

else:
    pass
    # load ../../../../../../Documents/MATLAB/inverse-approximation/c.mat C g d   

# How do you preffer the choice of the alpha?
# 1 - (Lobel et al, 1996)
# 2 - Golden section method
alphaopt = 1

# Initializing variables
cnvg    = np.zeros((maxit+1,2))     # Convergence data
I       = sps.eye(N,dtype=complex)  # Identity matrix
LC      = lag.inv(I-gd@C)          # Initial inversion
rho     = np.asarray(es-gs@C@LC@ei)            # Initial residual

# Printing first solution
print('Iteration: 0 - Cost function: %.2e' %lag.norm(rho.reshape(-1))**2)

if initopt is not 2:
    cnvg[0,:] = np.array([lag.norm(rho.reshape(-1))**2,lag.norm(g)])
else:
    cnvg[0,:] = np.array([lag.norm(rho.reshape(-1))**2,.0])


totaltime = time.time()

# Iterations
for it in range(1,maxit):
    
    tic = time.time()
    
    # Computing the gradient
    gradJ = np.asarray(np.reshape(-2*np.conj(sps.spdiags(np.reshape(LC@ei,N*M),0,N*M,N*M) @ np.tile(LC,(L,1)))@gs.conj().T@rho,(N,-1)))
    gradJ = np.sum(gradJ[:,np.arange(0,L**2,L)+np.arange(0,L)],1)
    
    g_last = np.asarray(np.copy(g))
    g = np.asarray(-gradJ)
    g = g.reshape((-1,1))
    
    # Computing the optimum direction
    d = np.asarray(g) + inner(g,g-g_last,dS)/lag.norm(g_last)**2*d
    D = sps.spdiags(d.reshape(-1),0,N,N)

    # Computing v matrix
    v = np.asarray(gs@LC.T@D@LC@ei)
    
    # Computing step
    if alphaopt is 1:
        alpha = 0
        for l in range(L):
            alpha = alpha + inner(rho[:,l],v[:,l],dx)
        alpha = alpha/lag.norm(v.reshape(-1))**2
    else:
        alpha = gsmethod(rho,v)
      
    # Computing next contrast
    C = C + alpha*D
    
    # Computing the inverse matrix
    LC = lag.inv(I-gd*C)
    
    # Computing the residual
    rho = np.asarray(es-gs*C*LC*ei)
    
    # Computing the objective function
    J = lag.norm(rho.reshape(-1))**2
    t = time.time()-tic
    
    # Printing iteration
    print('Iteration: %d' %it
          + ' - Cost function: %.2e' %J
          + ' - norm(g): %.2e' %lag.norm(g)
          + ' - time: %.1f sec' %t)
    
    # Saving objetive function and gradient magnitude
    cnvg[it,:] = np.array([J,lag.norm(g)])


totaltime = time.time()-totaltime
print('Total time: %f' %totaltime + ' seconds')

# Recovering dielectric properties
tau     = np.reshape(np.diag(C),I,J)             # Constrast fuction
epsr    = np.real(tau) + epsrb           # Retrieved relative permittivity
sig     = -omega*eps0*epsrb*np.imag(tau) # Relative conductivity [S/m]

# % Plotting results
# figure
# load ./genfields/grid.mat


plt.imshow(epsr, extent = [x[0], x[-1], y[0], y[-1]])
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title(r'Relative Permittivity  - $f = $ %.1e [Hz]' %self.f[indx])
cbar = plt.colorbar()
cbar.set_label(r'$|\epsilon_r|$')
# plt.savefig(filename + '_f%d' %indx, format = fileformat)
plt.show()
plt.close()

# % Relative permittivity plot
# subplot(3,2,1)
# imagesc(y,x,epsr')
# set(gca,'YDir','normal')
# xlabel('x [m]')
# ylabel('y [m]')
# title('Relative permittivity')
# clb = colorbar;
# ylabel(clb, '\epsilon_r')

# % Conductivity plot
# subplot(3,2,2)
# imagesc(y,x,sig')
# set(gca,'YDir','normal')
# xlabel('x [m]')
# ylabel('y [m]')
# title('Conductivity')
# clb = colorbar;
# ylabel(clb, '\sigma [S/m]')

# % Gradient - Real
# subplot(3,2,3)
# imagesc(y,x,real(reshape(g,I,J))')
# set(gca,'YDir','normal')
# xlabel('x [m]')
# ylabel('y [m]')
# title('Gradient - Real')
# clb = colorbar;
# ylabel(clb, 'g')

# % Conductivity plot
# subplot(3,2,4)
# imagesc(y,x,imag(reshape(g,I,J))')
# set(gca,'YDir','normal')
# xlabel('x [m]')
# ylabel('y [m]')
# title('Gradient - Imaginary')
# clb = colorbar;
# ylabel(clb, 'g')

# % Convergence plot - Cost Function
# subplot(3,2,5)
# plot(0:maxit,cnvg(:,1),'linewidth',2)
# grid
# xlabel('Iterations')
# ylabel('J(C)')
# title('Cost Function')

# % Convergence plot - Gradient
# subplot(3,2,6)
# plot(0:maxit,cnvg(:,2),'linewidth',2)
# grid
# xlabel('Iterations')
# ylabel('|\nabla J(C)|')
# title('Gradient')

# savefig('lobel96fig.fig')

# % Saving solution
# save ../../../../../../Documents/MATLAB/inverse-approximation/lobel96.mat C cnvg totaltime -v7.3
