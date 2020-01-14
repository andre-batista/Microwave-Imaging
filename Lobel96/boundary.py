import numpy as np

mu_0, eps_0 = 4*np.pi*1e-7, 8.854187817e-12
c = 1/np.sqrt(mu_0*eps_0)

class BoundaryCondition:
    
    bd_type = ''
    
    def __init__(self):
        pass
    
    def set_parameters_xn(self,dx,dy,dt):
        pass

class CPML(BoundaryCondition):
    
    bd_type = 'cpml'
    
    p_order = int()
    sigma_ratio = float()
    kappa_max = float()
    alpha_min = float()
    alpha_max = float
    side = ''
    air_buffer_number_of_cells = int()
    n_cpml = int()
    
    cpml_b_e, cpml_a_e = np.array([]), np.array([])
    cpml_b_m, cpml_a_m = np.array([]), np.array([])
    Psi_e, Psi_h = np.array([]), np.array([])
    CPsi_e, CPsi_h = np.array([]), np.array([])
    
    def __init__(self, side, cpml_order = 3, cpml_sigma_factor = 1.5, 
                 cpml_kappa_max = 7., cpml_alpha_min = .0, cpml_alpha_max = .05,
                 air_buffer_number_of_cells = 10,
                 cpml_number_of_cells = 8):
        
        self.side = side
        self.p_order = cpml_order
        self.sigma_ratio = cpml_sigma_factor
        self.kappa_max = cpml_kappa_max
        self.alpha_min, self.alpha_max = cpml_alpha_min, cpml_alpha_max
        self.air_buffer_number_of_cells = air_buffer_number_of_cells
        self.n_cpml = cpml_number_of_cells
        
    def set_parameters(self,dx,dy,dt,Ce,Ch):
        
        if self.side is 'xn':
        
            # define one-dimensional temporary cpml parameter arrays
            sigma_max = self.sigma_ratio  * (self.p_order+1)/(150*np.pi*dx)
            ncells = self.n_cpml
            rho_e = (np.arange(ncells,0,-1)-0.75)/ncells
            rho_m = (np.arange(ncells,0,-1)-0.25)/ncells
            sigma_pex_xn = sigma_max * rho_e**self.p_order
            sigma_pmx_xn = sigma_max * rho_m**self.p_order
            sigma_pmx_xn = (mu_0/eps_0) * sigma_pmx_xn
            kappa_ex_xn = 1 + (self.kappa_max - 1) * rho_e**self.p_order
            kappa_mx_xn = 1 + (self.kappa_max - 1) * rho_m**self.p_order
            alpha_ex_xn = self.alpha_min + (self.alpha_max - self.alpha_min) * (1-rho_e)
            alpha_mx_xn = self.alpha_min + (self.alpha_max - self.alpha_min) * (1-rho_m)
            alpha_mx_xn = (mu_0/eps_0) * alpha_mx_xn

            # define one-dimensional cpml parameter arrays
            self.cpml_b_e = np.exp((-dt/eps_0)*((sigma_pex_xn/kappa_ex_xn)+ alpha_ex_xn))
            self.cpml_a_e = (1/dx)*(self.cpml_b_e-1.0)* sigma_pex_xn/(kappa_ex_xn*(sigma_pex_xn+kappa_ex_xn*alpha_ex_xn))
            self.cpml_b_m = np.exp((-dt/mu_0)*((sigma_pmx_xn/kappa_mx_xn)+ alpha_mx_xn))
            self.cpml_a_m = (1/dx)*(self.cpml_b_m-1.0) * sigma_pmx_xn/(kappa_mx_xn*(sigma_pmx_xn+kappa_mx_xn*alpha_mx_xn))

            nyp1 = Ce.shape[1]

            # Create and initialize cpml convolution parameters
            self.Psi_e = np.zeros((ncells,nyp1))
            self.Psi_h = np.zeros((ncells,nyp1))

            # Create and initialize cpml convolution coefficients
            self.CPsi_e = Ce[1:ncells+1,:]*dx 
            self.CPsi_h = Ch[:ncells,:]*dx

            # Adjust FDTD coefficients in the CPML region
            for i in range(ncells):
                Ce[i+1,:] = Ce[i+1,:]/kappa_ex_xn[i]
                Ch[i,:] = Ch[i,:]/kappa_mx_xn[i]
            
        elif self.side is 'xp':

            # define one-dimensional temporary cpml parameter arrays
            sigma_max = self.sigma_ratio  * (self.p_order+1)/(150*np.pi*dx)
            ncells = self.n_cpml
            rho_e = (np.arange(1,ncells+1)-0.75)/ncells
            rho_m = (np.arange(1,ncells+1)-0.25)/ncells
            sigma_pex_xp = sigma_max * rho_e**self.p_order
            sigma_pmx_xp = sigma_max * rho_m**self.p_order
            sigma_pmx_xp = (mu_0/eps_0) * sigma_pmx_xp
            kappa_ex_xp = 1 + (self.kappa_max - 1) * rho_e**self.p_order
            kappa_mx_xp = 1 + (self.kappa_max - 1) * rho_m**self.p_order
            alpha_ex_xp = self.alpha_min + (self.alpha_max - self.alpha_min) * (1-rho_e)
            alpha_mx_xp = self.alpha_min + (self.alpha_max - self.alpha_min) * (1-rho_m)
            alpha_mx_xp = (mu_0/eps_0) * alpha_mx_xp

            # define one-dimensional cpml parameter arrays
            self.cpml_b_e = np.exp((-dt/eps_0)*((sigma_pex_xp/kappa_ex_xp)+alpha_ex_xp))
            self.cpml_a_e = (1/dx)*(self.cpml_b_e-1.0)*sigma_pex_xp/(kappa_ex_xp*(sigma_pex_xp+kappa_ex_xp*alpha_ex_xp))
            self.cpml_b_m = np.exp((-dt/mu_0)*((sigma_pmx_xp/kappa_mx_xp)+alpha_mx_xp))
            self.cpml_a_m = (1/dx)*(self.cpml_b_m-1.0) * sigma_pmx_xp/(kappa_mx_xp*(sigma_pmx_xp+kappa_mx_xp*alpha_mx_xp))

            nxp1, nyp1 = Ce.shape
            nx = nxp1-1 
            
            # Create and initialize cpml convolution parameters
            self.Psi_e = np.zeros((ncells,nyp1))
            self.Psi_h= np.zeros((ncells,nyp1))

            # Create and initialize cpml convolution coefficients
            self.CPsi_e = Ce[nxp1-ncells-1:nx,:]*dx
            self.CPsi_h = Ch[nxp1-ncells-1:nx,:]*dx

            # Adjust FDTD coefficients in the CPML region
            for i in range(ncells):
                Ce[nx-ncells+i,:] = Ce[nx-ncells+i,:]/kappa_ex_xp[i]
                Ch[nx-ncells+i,:] = Ch[nx-ncells+i,:]/kappa_mx_xp[i]
                
        elif self.side is 'yn':
            
            # define one-dimensional temporary cpml parameter arrays
            sigma_max = self.sigma_ratio  * (self.p_order+1)/(150*np.pi*dy)
            ncells = self.n_cpml
            rho_e = (np.arange(ncells,0,-1)-0.75)/ncells
            rho_m = (np.arange(ncells,0,-1)-0.25)/ncells
            sigma_pey_yn = sigma_max * rho_e**self.p_order
            sigma_pmy_yn = sigma_max * rho_m**self.p_order
            sigma_pmy_yn = (mu_0/eps_0) * sigma_pmy_yn
            kappa_ey_yn = 1 + (self.kappa_max - 1) * rho_e**self.p_order
            kappa_my_yn = 1 + (self.kappa_max - 1) * rho_m**self.p_order
            alpha_ey_yn = self.alpha_min + (self.alpha_max - self.alpha_min) * (1-rho_e)
            alpha_my_yn = self.alpha_min + (self.alpha_max - self.alpha_min) * (1-rho_m)
            alpha_my_yn = (mu_0/eps_0) * alpha_my_yn

            # define one-dimensional cpml parameter arrays
            self.cpml_b_e = np.exp((-dt/eps_0)*((sigma_pey_yn/kappa_ey_yn)+alpha_ey_yn))
            self.cpml_a_e = (1/dy)*(self.cpml_b_e-1.0)*sigma_pey_yn/(kappa_ey_yn*(sigma_pey_yn+kappa_ey_yn*alpha_ey_yn))
            self.cpml_b_m = np.exp((-dt/mu_0)*((sigma_pmy_yn/kappa_my_yn)+alpha_my_yn))
            self.cpml_a_m = (1/dy)*(self.cpml_b_m-1.0)*sigma_pmy_yn/(kappa_my_yn*(sigma_pmy_yn+kappa_my_yn*alpha_my_yn))

            nxp1 = Ce.shape[0]

            # Create and initialize cpml convolution parameters
            self.Psi_e = np.zeros((nxp1,ncells))
            self.Psi_h = np.zeros((nxp1,ncells))

            # Create and initialize cpml convolution coefficients
            self.CPsi_e = Ce[:,1:ncells+1]*dy
            self.CPsi_h = Ch[:,:ncells]*dy

            # Adjust FDTD coefficients in the CPML region
            for i in range(ncells):
                Ce[:,i+1] = Ce[:,i+1]/kappa_ey_yn[i]
                Ch[:,i] = Ch[:,i]/kappa_my_yn[i]
                
        else:

            # define one-dimensional temporary cpml parameter arrays
            sigma_max = self.sigma_ratio  * (self.p_order+1)/(150*np.pi*dy)
            ncells = self.n_cpml
            rho_e = (np.arange(1,ncells+1)-0.75)/ncells
            rho_m = (np.arange(1,ncells+1)-0.25)/ncells
            sigma_pey_yp = sigma_max * rho_e**self.p_order
            sigma_pmy_yp = sigma_max * rho_m**self.p_order
            sigma_pmy_yp = (mu_0/eps_0) * sigma_pmy_yp
            kappa_ey_yp = 1 + (self.kappa_max - 1) * rho_e**self.p_order
            kappa_my_yp = 1 + (self.kappa_max - 1) * rho_m**self.p_order
            alpha_ey_yp = self.alpha_min + (self.alpha_max - self.alpha_min) * (1-rho_e)
            alpha_my_yp = self.alpha_min + (self.alpha_max - self.alpha_min) * (1-rho_m)
            alpha_my_yp = (mu_0/eps_0) * alpha_my_yp

            # define one-dimensional cpml parameter arrays
            self.cpml_b_e = np.exp((-dt/eps_0)*((sigma_pey_yp/kappa_ey_yp)+alpha_ey_yp))
            self.cpml_a_e = (1/dy)*(self.cpml_b_e-1.0)*sigma_pey_yp/(kappa_ey_yp*(sigma_pey_yp+kappa_ey_yp*alpha_ey_yp))
            self.cpml_b_m = np.exp((-dt/mu_0)*((sigma_pmy_yp/kappa_my_yp)+alpha_my_yp))
            self.cpml_a_m = (1/dy)*(self.cpml_b_m-1.0)*sigma_pmy_yp/(kappa_my_yp*(sigma_pmy_yp+kappa_my_yp*alpha_my_yp))

            nxp1, nyp1 = Ce.shape
            ny = nyp1-1

            # Create and initialize cpml convolution parameters
            self.Psi_e = np.zeros((nxp1,ncells))
            self.Psi_h = np.zeros((nxp1,ncells))

            # Create and initialize cpml convolution coefficients
            self.CPsi_e = Ce[:,nyp1-ncells-1:ny]*dy
            self.CPsi_h = Ch[:,nyp1-ncells-1:ny]*dy

            # Adjust FDTD coefficients in the CPML region
            for i in range(ncells):
                Ce[:,ny-ncells+i] = Ce[:,ny-ncells+i]/kappa_ey_yp[i]
                Ch[:,ny-ncells+i] = Ch[:,ny-ncells+i]/kappa_my_yp[i]