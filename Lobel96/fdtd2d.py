import copy as cp
import numpy as np
import pickle
import matplotlib.pyplot as plt
from scipy import interpolate

import domain as dm
import boundary as bd
import waveform as wv
import source as sc
import probe as pb

mu_0, eps_0 = 4*np.pi*1e-7, 8.854187817e-12
c = 1/np.sqrt(mu_0*eps_0)

class FDTD2D:

    domain = dm.Domain(.0, .0, 0, 0)
    boundary = [bd.BoundaryCondition(),bd.BoundaryCondition(),bd.BoundaryCondition(),bd.BoundaryCondition()]
    source = sc.Source((.0,.0),wv.Waveform(), .0)
    epsrb, sigb = float(), float()
    f = np.array([])
    probes = []
    time_window = float()
    ls_x, ls_y = float(), float()
    courant_factor = float()
    et, x, y = np.array([]), np.array([]), np.array([])
    epsr, sig = np.array([]), np.array([])
        
    def __init__(self, domain, epsrb, sigb, sampled_frequencies, ls_x, ls_y,
                 time_window, courant_factor = .9, probes=None, 
                 boundary = [bd.CPML('xn'), bd.CPML('xp'), 
                             bd.CPML('yn'), bd.CPML('yp')]):
        self.domain = cp.deepcopy(domain)
        self.boundary = cp.deepcopy(boundary)
        self.epsrb, self.sigb = epsrb, sigb
        self.f = sampled_frequencies
        self.time_window = time_window
        self.ls_x, self.ls_y = ls_x, ls_y
        self.courant_factor = courant_factor
        
        if probes is not None:
            self.probes = []
            for i in range(len(probes)):
                self.probes.append(cp.deepcopy(probes[i]))
                
    def run(self, source, epsr=None, sig=None, probes=None):
        
        courant_factor = self.courant_factor
        dx, dy = self.domain.dx, self.domain.dy
        
        if epsr is None and sig is None:
            epsr = self.epsrb*np.ones((self.domain.I,self.domain.J))
            sig = self.sigb*np.ones((self.domain.I,self.domain.J))
        elif epsr is None:
            epsr = self.epsrb*np.ones((self.domain.I,self.domain.J))
        elif sig is None:
            sig = self.sigb*np.ones((self.domain.I,self.domain.J))
        
        nxobj, nyobj = epsr.shape
        lx, ly = nxobj*dx, nyobj*dy
        dt = courant_factor/(c*np.sqrt((1/dx**2)+(1/dy**2)))
        
        if probes is not None:    
            for i in range(len(probes)):
                self.probes.append(cp.deepcopy(probes[i]))
                
        any_save_probe_signal = False
        for i in range(len(self.probes)):
            if self.probes[i].save_signal:
                any_save_probe_signal = True
                break
        
        impressed_J = []
        impressed_M = []
        impressed_J.append(cp.deepcopy(source))
        
        # Define source position
        pos_x = impressed_J[0].center_position[0]
        pos_y = impressed_J[0].center_position[1]
        impressed_J[0].min_x = pos_x-self.ls_x/2
        impressed_J[0].max_x = pos_x+self.ls_x/2
        impressed_J[0].min_y = pos_y-self.ls_y/2
        impressed_J[0].max_y = pos_y+self.ls_y/2
        
        pb_min_x, pb_max_x = -lx/2, lx/2
        pb_min_y, pb_max_y = -ly/2, ly/2
        for i in range(len(self.probes)):
            if pb_min_x > self.probes[i].position[0]:
                pb_min_x = self.probes[i].position[0]
            elif pb_max_x < self.probes[i].position[0]:
                pb_max_x = self.probes[i].position[0]
            if pb_min_y > self.probes[i].position[1]:
                pb_min_y = self.probes[i].position[1]
            elif pb_max_y < self.probes[i].position[1]:
                pb_max_y = self.probes[i].position[1]
        if pb_min_x > impressed_J[0].min_x:
            pb_min_x = impressed_J[0].min_x
        elif pb_max_x < impressed_J[0].max_x:
            pb_max_x < impressed_J[0].max_x
        if pb_min_y > impressed_J[0].min_y:
            pb_min_y = impressed_J[0].min_y
        elif pb_max_y < impressed_J[0].max_y:
            pb_max_y < impressed_J[0].max_y
        self.domain.min_x, self.domain.max_x = pb_min_x-5*dx, pb_max_x+5*dx
        self.domain.min_y, self.domain.max_y = pb_min_y-5*dy, pb_max_y+5*dy
        
        # Defining source parameters
        impressed_J[0].direction = 'zp'
        impressed_J[0].waveform_type = 'gaussian'
        
        # Determine the problem space boundaries including air buffers
        if self.boundary[0].bd_type is 'cpml':
            self.domain.min_x = (
                self.domain.min_x 
                - dx*self.boundary[0].air_buffer_number_of_cells)
        
        if self.boundary[1].bd_type is 'cpml':
            self.domain.max_x = (
                self.domain.max_x 
                + dx*self.boundary[1].air_buffer_number_of_cells)
        
        if self.boundary[2].bd_type is 'cpml':
            self.domain.min_y = (
                self.domain.min_y 
                - dy*self.boundary[2].air_buffer_number_of_cells)
        
        if self.boundary[3].bd_type is 'cpml':
            self.domain.max_y = (
                self.domain.max_y 
                + dy*self.boundary[3].air_buffer_number_of_cells)
        
        # Determine the problem space boundaries including cpml layers
        if (self.boundary[0].bd_type is 'cpml' 
            and self.boundary[0].n_cpml > 0):
            self.domain.min_x = (self.domain.min_x-dx*self.boundary[0].n_cpml)
            
        if (self.boundary[1].bd_type is 'cpml' 
            and self.boundary[1].n_cpml > 0):
            self.domain.max_x = (self.domain.max_x+dx*self.boundary[1].n_cpml)
            
        if (self.boundary[2].bd_type is 'cpml' 
            and self.boundary[2].n_cpml > 0):
            self.domain.min_y = (self.domain.min_y-dy*self.boundary[2].n_cpml)
            
        if (self.boundary[3].bd_type is 'cpml' 
            and self.boundary[3].n_cpml > 0):
            self.domain.max_y = (self.domain.max_y+dy*self.boundary[3].n_cpml)
        
        # Determining the problem space size
        self.domain.size_x = self.domain.max_x - self.domain.min_x
        self.domain.size_y = self.domain.max_y - self.domain.min_y

        # Number of cells in x and y directions
        nx = np.round(self.domain.size_x/dx).astype(int)
        ny = np.round(self.domain.size_y/dy).astype(int)
        
        # Adjust domain size by snapping to cells
        self.domain.size_x = nx * dx
        self.domain.size_y = ny * dy
        self.domain.max_x = self.domain.min_x + self.domain.size_x
        self.domain.max_y = self.domain.min_y + self.domain.size_y
        
        # Some frequently used auxiliary parameters
        nxp1, nxm1 = nx+1, nx-1
        nyp1, nym1 = ny+1, ny-1
        
        # Create arrays storing the center coordinates of the cells
        self.domain.cell_center_coordinates_x = np.zeros((nx, ny))
        self.domain.cell_center_coordinates_y = np.zeros((nx, ny))
        for ind in range(nx):
            self.domain.cell_center_coordinates_x[ind, :] = (
                (ind+1 - 0.5) * dx + self.domain.min_x)
        for ind in range(ny):
            self.domain.cell_center_coordinates_y[:, ind] = (
                (ind+1 - 0.5) * dy + self.domain.min_y)
            
        # x and y coordinates for epsr and Ez matrices
        xcoor = np.round(np.linspace(self.domain.min_x, 
                                     self.domain.max_x, nxp1), 4)
        ycoor = np.round(np.linspace(self.domain.min_y, 
                                     self.domain.max_y, nyp1), 4)
        
        # TMz components
        if impressed_J[0].direction[0] is 'z':
            eps_r_x = np.array([])
            eps_r_y = np.array([])
            eps_r_z = self.epsrb*np.ones((nxp1, nyp1))
            mu_r_x = np.ones((nxp1, ny))
            mu_r_y = np.ones((nx, nyp1))
            mu_r_z = np.array([])
            sigma_e_x = np.array([])
            sigma_e_y = np.array([])
            sigma_e_z = self.sigb*np.ones((nxp1, nyp1))
            sigma_m_x = np.zeros((nxp1, ny))
            sigma_m_y = np.zeros((nx, nyp1))
            sigma_m_z = np.array([])
        elif impressed_J[0].direction[0] is 'x':
            # Some day it'll be provided
            pass
        else: # y
             # Some day it'll be provided
            pass           
        
        # Adding user's domain
        eps_r_z[np.ix_(np.logical_and(xcoor > -lx/2, xcoor <= lx/2), 
                       np.logical_and(ycoor > -ly/2, ycoor <= ly/2))] = epsr
        sigma_e_z[np.ix_(np.logical_and(xcoor > -lx/2, xcoor <= lx/2), 
                         np.logical_and(ycoor > -ly/2, ycoor <= ly/2))] = sig

        # time array
        time = np.arange(dt, self.time_window+dt, dt)-.5*dt
        number_of_time_steps = time.size

        # Create and initialize field and current arrays
        Hx = np.zeros((nxp1, ny))
        Hy = np.zeros((nx, nyp1))
        Hz = np.array([])
        Ex = np.array([])
        Ey = np.array([])
        Ez = np.zeros((nxp1, nyp1))

        number_of_impressed_J = len(impressed_J)
        number_of_impressed_M = len(impressed_M)
        
        for i in range(number_of_impressed_J):
            impressed_J[i].signal.compute_waveform(time)
        
        # electric current sources
        for ind in range(number_of_impressed_J):
            i_s = np.floor((impressed_J[ind].min_x 
                            - self.domain.min_x)/dx).astype(int)
            j_s = np.floor((impressed_J[ind].min_y 
                            - self.domain.min_y)/dy).astype(int)
            i_e = np.floor((impressed_J[ind].max_x 
                            - self.domain.min_x)/dx).astype(int)
            j_e = np.floor((impressed_J[ind].max_y 
                            - self.domain.min_y)/dy).astype(int)
            
            impressed_J[ind].i_s = i_s
            impressed_J[ind].j_s = j_s
            impressed_J[ind].i_e = i_e
            impressed_J[ind].j_e = j_e

            if impressed_J[ind].direction[1] is 'n':
                j_magnitude_factor = -1*impressed_J[ind].magnitude
            else:
                j_magnitude_factor =  impressed_J[ind].magnitude

            # copy waveform of the waveform type to waveform of the source
            impressed_J[ind].signal.waveform = (
                j_magnitude_factor * impressed_J[ind].signal.waveform)

        # magnetic current sources
        for ind in range(number_of_impressed_M):
            i_s = np.round((impressed_M[ind].min_x 
                            - self.domain.min_x)/dx).astype(int)
            j_s = np.round((impressed_M[ind].min_y 
                            - self.domain.min_y)/dy).astype(int)
            i_e = np.round((impressed_M[ind].max_x 
                            - self.domain.min_x)/dx).astype(int)
            j_e = np.round((impressed_M[ind].max_y 
                            - self.domain.min_y)/dy).astype(int)
            
            impressed_M[ind].i_s = i_s
            impressed_M[ind].j_s = j_s
            impressed_M[ind].i_e = i_e
            impressed_M[ind].j_e = j_e

            if impressed_M[ind].direction[1] is 'n':
                m_magnitude_factor = -1*impressed_M[ind].magnitude
            else:
                m_magnitude_factor =  impressed_M[ind].magnitude

            # copy waveform of the waveform type to waveform of the source
            impressed_M[ind].signal.waveform = (
                m_magnitude_factor * impressed_M[ind].signal.waveform)
            
        # Coeffiecients updating Ez
        Ceze  =  (2*eps_r_z*eps_0-dt*sigma_e_z)/(2*eps_r_z*eps_0+dt*sigma_e_z)
        Cezhy =  (2*dt/dx)/(2*eps_r_z*eps_0 + dt*sigma_e_z)
        Cezhx = -(2*dt/dy)/(2*eps_r_z*eps_0 + dt*sigma_e_z)

        # Coeffiecients updating Hx
        Chxh  =  (2*mu_r_x*mu_0 - dt*sigma_m_x)/(2*mu_r_x*mu_0 + dt*sigma_m_x)
        Chxez = -(2*dt/dy)/(2*mu_r_x*mu_0 + dt*sigma_m_x)

        # Coeffiecients updating Hy
        Chyh  =  (2*mu_r_y*mu_0 - dt*sigma_m_y)/(2*mu_r_y*mu_0 + dt*sigma_m_y)
        Chyez =  (2*dt/dx)/(2*mu_r_y*mu_0 + dt*sigma_m_y)
        
        # Initialize coeffiecients for impressed current sources
        for ind in range(number_of_impressed_J):
            i_s = impressed_J[ind].i_s
            j_s = impressed_J[ind].j_s
            i_e = impressed_J[ind].i_e
            j_e = impressed_J[ind].j_e
            if impressed_J[ind].direction[0] is 'x':
                impressed_J[ind].Cexj = (
                    -(2*dt)
                    / (2*eps_r_x[np.ix_(range(i_s, i_e), range(j_s, j_e+1))]*eps_0 
                       + dt*sigma_e_x[np.ix_(range(i_s, i_e), range(j_s, j_e+1))])
                )
            elif impressed_J[ind].direction[0] is 'y':
                impressed_J[ind].Ceyj = (
                    -(2*dt)
                    / (2*eps_r_y[np.ix_(range(i_s, i_e+1), range(j_s, j_e))]*eps_0 
                       + dt*sigma_e_y[np.ix_(range(i_s, i_e+1), range(j_s, j_e))])
                )
            else:
                impressed_J[ind].Cezj = (
                    -(2*dt)
                    /(2*eps_r_z[np.ix_(range(i_s, i_e+1), range(j_s, j_e+1))]*eps_0 
                      + dt*sigma_e_z[np.ix_(range(i_s, i_e+1), range(j_s, j_e+1))])
                )

        for ind in range(number_of_impressed_M):
            i_s = impressed_M[ind].i_s
            j_s = impressed_M[ind].j_s
            i_e = impressed_M[ind].i_e
            j_e = impressed_M[ind].j_e
            if impressed_M[ind].direction[0] is 'x':
                impressed_M[ind].Chxm = (-(2*dt)
                                         /(2*mu_r_x[np.ix_(range(i_s, i_e+1), 
                                                           range(j_s, j_e))]
                                           * mu_0 + dt
                                           * sigma_m_x[np.ix_(range(i_s, i_e+1), 
                                                              range(j_s, j_e))]))
            elif impressed_M[ind].direction[0] is 'y':
                impressed_M[ind].Chym = (-(2*dt)
                                         / (2*mu_r_y[np.ix_(range(i_s, i_e), 
                                                            range(j_s, j_e+1))]
                                            * mu_0 + dt
                                            * sigma_m_y[np.ix_(range(i_s, i_e), 
                                                               range(j_s, j_e+1))]))
            else:
                impressed_M[ind].Chzm = (-(2*dt)
                                         / (2*mu_r_z[np.ix_(range(i_s, i_e), 
                                                            range(j_s, j_e))]
                                            * mu_0 + dt
                                            * sigma_m_z[np.ix_(range(i_s, i_e), 
                                                               range(j_s, j_e))]))
        
        # Initialize cpml for xn region
        if self.boundary[0].bd_type is 'cpml':
            self.boundary[0].set_parameters(dx, dy, dt, Cezhy, Chyez)
        
        # Initialize cpml for xp region
        if self.boundary[1].bd_type is 'cpml':
            self.boundary[1].set_parameters(dx, dy, dt, Cezhy, Chyez)
        
        # Initialize cpml for yn region
        if self.boundary[2].bd_type is 'cpml':
            self.boundary[2].set_parameters(dx, dy, dt, Cezhx, Chxez)
           
        # Initialize cpml for yp region 
        if self.boundary[3].bd_type is 'cpml':
            self.boundary[3].set_parameters(dx, dy, dt, Cezhx, Chxez)
        
        current_time = 0
        aux = (Ez.shape, self.f.size)
        self.et = np.zeros(aux[0]+(aux[1], ), dtype = complex)
        for i in range(len(self.probes)):
            self.probes[i].allocate_signal(time)
                
        # Starting the simulation
        for time_step in range(number_of_time_steps):

            current_time  = current_time + dt/2

            Hx [np.ix_(range(nxp1), range(ny))] = (
                Chxh[np.ix_(range(nxp1), range(ny))]
                * Hx[np.ix_(range(nxp1), range(ny))]
                + Chxez[np.ix_(range(nxp1), range(ny))]
                * (Ez[np.ix_(range(nxp1), range(1, nyp1))]
                - Ez[np.ix_(range(nxp1), range(ny))])
            )
            

            Hy[np.ix_(range(nx), range(nyp1))] = (
                Chyh[np.ix_(range(nx), range(nyp1))]
                * Hy[np.ix_(range(nx), range(nyp1))]
                + Chyez[np.ix_(range(nx), range(nyp1))]
                * (Ez[np.ix_(range(1, nxp1), range(nyp1))]
                - Ez[np.ix_(range(nx), range(nyp1))])
            )
            
            for ind in range(number_of_impressed_M):
                i_s = impressed_M[ind].i_s
                j_s = impressed_M[ind].j_s
                i_e = impressed_M[ind].i_e
                j_e = impressed_M[ind].j_e
                if impressed_M[ind].direction[0] is 'x':
                    Hx[np.ix_(np.arange(i_s, i_e+1), np.arange(j_s, j_e))] = (
                        Hx[i_s:i_e+1][:, j_s:j_e]
                        + impressed_M[ind].Chxm
                        * impressed_M[ind].signal.waveform[time_step]
                    )
                elif impressed_M[ind].direction[0] is 'y':
                    Hy[np.ix_(np.arange(i_s, i_e), np.arange(j_s, j_e+1))] = (
                        Hy[i_s:i_e][:, j_s:j_e+1]
                        + impressed_M[ind].Chym 
                        * impressed_M[ind].signal.waveform[time_step])
                else:
                    Hz[np.ix_(np.arange(i_s, i_e), np.arange(j_s, j_e))] = (
                        Hz[i_s:i_e][:, j_s:j_e]
                        + impressed_M[ind].Chzm
                        * impressed_M[ind].signal.waveform[time_step])
            
            # TMz
            # apply CPML to magnetic field components
            if self.boundary[0].bd_type is 'cpml':
                for i in range(self.boundary[0].n_cpml):
                    self.boundary[0].Psi_h[i, :] = (
                        self.boundary[0].cpml_b_m[i] 
                        * self.boundary[0].Psi_h[i, :] 
                        + self.boundary[0].cpml_a_m[i]*(Ez[i+1, :]-Ez[i, :]))
                Hy[:self.boundary[0].n_cpml, :] = (
                    Hy[:self.boundary[0].n_cpml, :] + self.boundary[0].CPsi_h 
                    * self.boundary[0].Psi_h)
            
            if self.boundary[1].bd_type is 'cpml':
                n_st = nx - self.boundary[1].n_cpml
                for i in range(self.boundary[1].n_cpml):
                    self.boundary[1].Psi_h[i, :] = (
                        self.boundary[1].cpml_b_m[i] 
                        * self.boundary[1].Psi_h[i, :] 
                        + self.boundary[1].cpml_a_m[i]
                        * (Ez[i+n_st+1, :]-Ez[i+n_st, :]))
                Hy[n_st:nx, :] = (Hy[n_st:nx, :] + self.boundary[1].CPsi_h 
                                  * self.boundary[1].Psi_h)
            
            if self.boundary[2].bd_type is 'cpml':                
                for i in range(self.boundary[2].n_cpml):
                    self.boundary[2].Psi_h[:, i] = (
                        self.boundary[2].cpml_b_m[i] 
                        * self.boundary[2].Psi_h[:, i] 
                        + self.boundary[2].cpml_a_m[i]*(Ez[:, i+1]-Ez[:, i]))
                Hx[:, :self.boundary[2].n_cpml] = (
                    Hx[:, :self.boundary[2].n_cpml] + self.boundary[2].CPsi_h 
                    * self.boundary[2].Psi_h)
            
            if self.boundary[3].bd_type is 'cpml':
                n_st = ny - self.boundary[3].n_cpml
                for i in range(self.boundary[3].n_cpml):
                    self.boundary[3].Psi_h[:, i] = (
                        self.boundary[3].cpml_b_m[i] 
                        * self.boundary[3].Psi_h[:, i] 
                        + self.boundary[3].cpml_a_m[i]
                        * (Ez[:, i+n_st+1] - Ez[:, i+n_st]))
                Hx[:, n_st:ny] = (Hx[:, n_st:ny] + self.boundary[3].CPsi_h 
                                  * self.boundary[3].Psi_h)
                
            current_time  = current_time + dt/2

            Ez[np.ix_(range(1, nx), range(1, ny))] = (
                Ceze[np.ix_(range(1, nx), range(1, ny))]
                * Ez[np.ix_(range(1, nx), range(1, ny))]
                + Cezhy[np.ix_(range(1, nx), range(1, ny))]
                * (Hy[np.ix_(range(1, nx), range(1, ny))]
                   - Hy[np.ix_(range(nxm1), range(1, ny))])
                + Cezhx[np.ix_(range(1, nx), range(1, ny))]
                * (Hx[np.ix_(range(1, nx), range(1, ny))]
                   - Hx[np.ix_(range(1, nx), range(nym1))])
            )
            
            for ind in range(number_of_impressed_J):
                i_s = impressed_J[ind].i_s
                j_s = impressed_J[ind].j_s
                i_e = impressed_J[ind].i_e
                j_e = impressed_J[ind].j_e
                if impressed_J[ind].direction[0] is 'x':
                    Ex[np.ix_(range(i_s, i_e), range(j_s, j_e+1))] = (
                        Ex[np.ix_(range(i_s, i_e), range(j_s, j_e+1))] 
                        + impressed_J[ind].Cexj 
                        * impressed_J[ind].signal.waveform[time_step]
                    )
                elif impressed_J[ind].direction[0] is 'y':
                    Ey[np.ix_(range(i_s, i_e+1), range(j_s, j_e))] = (
                        Ey[np.ix_(range(i_s, i_e+1), range(j_s, j_e))]
                        + impressed_J[ind].Ceyj * 
                        impressed_J[ind].signal.waveform[time_step]
                    )
                else:
                    Ez[np.ix_(range(i_s, i_e+1), range(j_s, j_e+1))] = (
                        Ez[np.ix_(range(i_s, i_e+1), range(j_s, j_e+1))]
                        + impressed_J[ind].Cezj
                        * impressed_J[ind].signal.waveform[time_step]
                    )
                    
            # update electric fields at the PML regions
            # TMz
            # apply CPML to electric field components
            if self.boundary[0].bd_type is 'cpml':
                for i in range(self.boundary[0].n_cpml):
                    self.boundary[0].Psi_e[i, :] = (
                        self.boundary[0].cpml_b_e[i] 
                        * self.boundary[0].Psi_e[i, :] 
                        + self.boundary[0].cpml_a_e[i]*(Hy[i+1, :] - Hy[i, :]))
                Ez[1:self.boundary[0].n_cpml+1, :] = (
                    Ez[1:self.boundary[0].n_cpml+1, :] 
                    + self.boundary[0].CPsi_e * self.boundary[0].Psi_e)

            if self.boundary[1].bd_type is 'cpml':
                n_st = nx - self.boundary[1].n_cpml
                for i in range(self.boundary[1].n_cpml):
                    self.boundary[1].Psi_e[i, :] = (
                        self.boundary[1].cpml_b_e[i] 
                        * self.boundary[1].Psi_e[i, :] 
                        + self.boundary[1].cpml_a_e[i]
                        * (Hy[i+n_st, :]-Hy[i+n_st-1, :]))
                Ez[n_st:nx, :] = (Ez[n_st:nx, :] + self.boundary[1].CPsi_e 
                                  * self.boundary[1].Psi_e)
                
            if self.boundary[2].bd_type is 'cpml':                
                for i in range(self.boundary[2].n_cpml):
                    self.boundary[2].Psi_e[:, i] = (
                        self.boundary[2].cpml_b_e[i] 
                        * self.boundary[2].Psi_e[:, i] 
                        + self.boundary[2].cpml_a_e[i]*(Hx[:, i+1] - Hx[:, i]))
                Ez[:, 1:self.boundary[2].n_cpml+1] = (
                    Ez[:, 1:self.boundary[2].n_cpml+1] 
                    + self.boundary[2].CPsi_e * self.boundary[2].Psi_e)
                
            if self.boundary[3].bd_type is 'cpml':                
                n_st = ny - self.boundary[3].n_cpml
                for i in range(self.boundary[3].n_cpml):
                    self.boundary[3].Psi_e[:, i] = (
                        self.boundary[3].cpml_b_e[i] 
                        * self.boundary[3].Psi_e[:, i] 
                        + self.boundary[3].cpml_a_e[i]
                        * (Hx[:, i+n_st]-Hx[:, i+n_st-1]))
                Ez[:, n_st:ny] = (Ez[:, n_st:ny] + self.boundary[3].CPsi_e 
                                  * self.boundary[3].Psi_e)
            
            if any_save_probe_signal:
                for i in range(len(self.probes)):
                    if self.probes[i].save_signal:
                        fp = interpolate.RectBivariateSpline(xcoor, ycoor, Ez)
                        self.probes[i].append(fp(self.probes[i].position[0],
                                                 self.probes[i].position[1]))
                
            for f_indx in range(self.f.size):
                self.et[:, :, f_indx] = (
                    self.et[:, :, f_indx] 
                    + Ez*np.exp(-1j*2*np.pi*self.f[f_indx]*(time_step+1)*dt)*dt)
                
  
        for f_indx in range(self.f.size):
            fint_real = interpolate.RectBivariateSpline(xcoor, ycoor, np.real(self.et[:, :, f_indx]))
            fint_imag = interpolate.RectBivariateSpline(xcoor, ycoor, np.imag(self.et[:, :, f_indx]))
            for irx in range(len(self.probes)):
                self.probes[irx].set_field_freq(fint_real(self.probes[i].position[0], 
                                                          self.probes[i].position[1]) 
                                                + 1j*fint_imag(self.probes[i].position[0], 
                                                               self.probes[i].position[1]),
                                                f_indx)
                
        # Setting the total electric field according to user's domain
        self.et = self.et[np.ix_(
            np.nonzero(np.logical_and(xcoor > 0, xcoor <= lx))[0], 
            np.nonzero(np.logical_and(ycoor > 0, ycoor <= ly))[0], 
            range(self.f.size))]
        self.x = xcoor[np.logical_and(xcoor > 0, xcoor <= lx)]
        self.y = ycoor[np.logical_and(ycoor > 0, ycoor <= ly)]
        
        self.epsr = np.copy(epsr)
        self.sig = np.copy(sig)
        
    def get_intern_field(self):
        return np.copy(self.et)
   
    def get_intern_coordinates(self):
        return np.copy(self.x), np.copy(self.y)
        
    def save_field_figure(self, filename, fileformat = '.eps'):
        for indx in range(self.f.size):
            plt.imshow(np.abs(self.et[:, :, indx]), extent = [self.x[0], self.x[-1], self.y[0], self.y[-1]])
            plt.xlabel('x [m]')
            plt.ylabel('y [m]')
            plt.title(r'Field Distribution - $f = $ %.1e [Hz]' %self.f[indx])
            cbar = plt.colorbar()
            cbar.set_label(r'$|E_z|$')
            plt.savefig(filename + '_f%d' %indx, format = fileformat)
            plt.close()
            
    def plot_probes(self,filename=None,fileformat=None,index=None):
        time = np.linspace(.0,self.time_window,self.probes[0].signal.size)
        if index is None and len(self.probes) > 1:
            fig = plt.figure()
            nrow = np.floor(np.sqrt(len(self.probes))).astype(int)
            ncol = np.ceil(len(self.probes)/nrow).astype(int)
            for ifig in range(len(self.probes)):
                ax = fig.add_subplot(nrow, ncol, ifig+1)
                ax.plot(time,self.probes[ifig].signal)
                ax.set_xlabel('Time [sec]')
                ax.set_ylabel('Electric Field Intensity [V/m]')
                ax.set_title('Probe at x = %.2e' %self.probes[ifig].position[0] 
                             + ' [m], y = %.2e' %self.probes[ifig].position[1] 
                             + ' [m]')
                plt.grid()
            plt.title('Probes\' signals')
        elif index is None:
            plt.plot(time,self.probes[0].signal)
            plt.title('Probe at x = %.2e' %self.probes[0].position[0] 
                      + ' [m], y = %.2e' %self.probes[0].position[1] + ' [m]')
            plt.xlabel('Time [sec]')
            plt.ylabel('Electric Field Intensity [V/m]')
            plt.grid()
        else:    
            plt.plot(time,self.probes[index].signal)
            plt.title('Probe at x = %.2e' %self.probes[index].position[0] 
                      + ' [m], y = %.2e' %self.probes[index].position[1] 
                      + ' [m]')
            plt.xlabel('Time [sec]')
            plt.ylabel('Electric Field Intensity [V/m]')
            plt.grid()            
            
        if filename is None:
            plt.show()
        else:
            if fileformat is None:
                plt.savefig(filename,'eps')
            else:
                plt.savefig(filename,fileformat)
        plt.close()            
            
    def save_data(self,filename):
        data = {'domain':self.domain,
                'boundary':self.boundary,
                'source':self.source,
                'epsrb':self.epsrb,
                'sigb':self.sigb,
                'f':self.f,
                'time_window':self.time_window,
                'dtheta':self.dtheta,
                'Rs':self.Rs,
                'ls_x':self.ls_x,
                'ls_y':self.ls_y,
                'courant_factor':self.courant_factor,
                'et':self.et,
                'x':self.x,
                'y':self.y,
                'probes':self.probes,
                'epsr':self.epsr,
                'sig':self.sig}
        
        with open(filename,'wb') as datafile:
            pickle.dump(data,datafile)

    def load_data(self,filename):
        with open(filename,'rb') as datafile:
            data = pickle.load(datafile)
        self.domain = data['domain']
        self.boundary = data['boundary']
        self.source = data['source']
        self.epsrb = data['epsrb']
        self.sigb = data['sigb']
        self.f = data['f']
        self.time_window = data['time_window']
        self.dtheta = data['dtheta']
        self.Rs = data['Rs']
        self.ls_x = data['ls_x']
        self.ls_y = data['ls_y']
        self.courant_factor = data['courant_factor']
        self.et = data['et']
        self.x = data['x']
        self.y = data['y']
        self.probes = data['probes']
        self.epsr = data['epsr']
        self.sig = data['sig']