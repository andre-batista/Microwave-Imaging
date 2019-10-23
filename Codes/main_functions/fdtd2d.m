function [et,x,y,rx_fields,rx_x,rx_y] = fdtd2d(DATA,PAR)

%FDTD2D Finite-Difference Time-Domain method for TMz-mode electromagnetic 
% simulation in two dimensions.
%   data = fdtd2d(data,par) simulates simulates the TMz-mode
%   electromagnetic problem using the Finite-Difference Time-Domain method.
%   The structs DATA and PAR specifies the domain and parameters of the
%   problem. The output is the same DATA struct with the electric field
%   distribution added. The sources are fixed in positions outside the 
%   given domain.
%
%   The DATA input struct must contain the following attributes:
%       dx:     Dimension of a unit cell in x direction [m].
%       dy:     Dimension of a unit cell in y direction [m].
%       epsrb:  Relative permittivity of the background media.
%       sigb:   Conductivity of the background media [S/m].
%       epsr:   Matrix containing the relative permittivity of each node of
%               your domain.
%       sig:    Matrix containing the conductivity of each node of your
%               domain [S/m].
%
%   The PAR input struct must contain the following attributes:
%       f:      Angular frequencies for Fourier Transform measurements [Hz].
%       fwv:    Angular frequency for waveform.
%       tw:     Time window [sec].
%       dtheta: Source angle discretization [degrees].
%       Rs:     Radius of source domain.
%       s_indx: Index of the source.
%       ls_x:   Length of source in x-axis
%       ls_y:   Length of source in y-axis
%
%   The DATA output struct constains the same attributes as the input,
%   plus:
%       dt:     Fixed time step [sec].
%       et:     Complex electric field distribution over the input domain for
%               each of the frequencies given. The dimension of this matrix is:
%               < number of cells in the x-direction, number of in the 
%               y-direction, number of given frequencies inputs >.
%       x,y:    Electric field domain coordinates
%
%   % Example:
%
%       % DATA struct
%       data.dx     = 0.1e-3;
%       data.dy     = 0.1e-3;
%       data.epsr   = ones(100,100);
%       data.sig    = zeros(size(data.epsr));
%       data.epsrb  = 1;
%       data.sigb   = 0;
%
%       % PAR struct
%       par.f       = [1e9;2e9;3e9];
%       par.tw      = 1e-8;
%       par.dtheta  = 20;
%       par.Rs      = 70 * 0.1e-3;
%       par.s_indx  = 4;
%       par.ls_x    = 4*dx;
%       par.ls_y    = 4*dy;
%
%       % Calling FDTD
%       [et,x,y,rx_fields,rx_x,rx_y] = fdtd2d(data,par);
%
%       % Plotting total electric field
%       figure
%       imagesc(y,x,abs(et(:,:,1))')
%       title('FDTD Method')
%       xlabel('x [m]')
%       ylabel('y [m]')
%
%   % Reference:
%
%   Elsherbeni, A. Z., & Demir, V. (2016). The finite-difference time-
%   domain method for electromagnetics with MATLAB simulations. The 
%   Institution of Engineering and Technology.

    % General parameters 
    courant_factor                  = 0.9;                                  % A factor that determines duration of a time step within the CFL limit.
    number_of_cells_per_wavelength  = 20;                                   % A factor determining the accuracy limit of FDTD results.
    dx                              = DATA.dx;                              % Dimensions of a unit cell in x direction (meters).
    dy                              = DATA.dy;                              % Dimensions of a unit cell in y direction (meters).
    [nxobj,nyobj]                   = size(DATA.epsr);                      % Size of relative permittivity array
    lx                              = nxobj*dx;                             % Size of the original domain in x-direction.
    ly                              = nyobj*dy;                             % Size of the original domain in y-direction.
    mu_0                            = 4*pi*1e-7;                            % Permeability of free space
    eps_0                           = 8.854187817e-12;                      % Permittivity of free space
    c                               = 1/sqrt(mu_0*eps_0);                   % Speed of light in free space
    dt                              = courant_factor/...                    % Time step [sec]
                                        (c*sqrt((1/dx^2)+(1/dy^2)));
    DATA.dt                         = dt;
    number_of_time_steps            = round(PAR.tw/dt);                     % Maximum number of time steps to run FDTD simulation.
    FLAG_RESOLUTION = false;

%     if 1/sqrt(mu_0*max(DATA.epsr(:))*eps_0)/max(PAR.f)/20 < min([dx,dy])
%         lambda_min = 1/sqrt(mu_0*max(DATA.epsr(:))*eps_0)/max(PAR.f);
%         dx = lambda_min/20;
%         dy = lambda_min/20;
%         old_dx = DATA.dx;
%         old_dy = DATA.dy;
%         [old_nxobj,old_nyobj] = size(DATA.epsr);
%         DATA = changemesh (DATA,dx,dy);
%         [nxobj,nyobj] = size(DATA.epsr);
%         dx = DATA.dx;
%         dy = DATA.dy;
%         lx = nxobj*dx;
%         ly = nyobj*dy;
%         dt = courant_factor/(c*sqrt((1/dx^2)+(1/dy^2)));
%         DATA.dt = dt;
%         FLAG_RESOLUTION = true;
%     end

    % Boundary condition parameters [CPML]
    boundary.cpml_order                     = 3;
    boundary.cpml_sigma_factor              = 1.5;
    boundary.cpml_kappa_max                 = 7;
    boundary.cpml_alpha_min                 = 0;
    boundary.cpml_alpha_max                 = 0.05;
    boundary.type_xn                        = 'cpml';
    boundary.type_xp                        = 'cpml';
    boundary.type_yn                        = 'cpml';
    boundary.type_yp                        = 'cpml';
    boundary.air_buffer_number_of_cells_xn  = 10;
    boundary.air_buffer_number_of_cells_xp  = 10;
    boundary.air_buffer_number_of_cells_yn  = 10;
    boundary.air_buffer_number_of_cells_yp  = 10;
    boundary.cpml_number_of_cells_xn        = 8;
    boundary.cpml_number_of_cells_xp        = 8;
    boundary.cpml_number_of_cells_yn        = 8;
    boundary.cpml_number_of_cells_yp        = 8;


    % Source initialization
    impressed_J = [];
    impressed_M = [];

    % Define source waveform types and parameters
    waveforms.gaussian(1).number_of_cells_per_wavelength ...
        = number_of_cells_per_wavelength;
    
    % Initializig domain sizes
    fdtd_domain.min_x = 0  -1.1*(PAR.Rs-lx/2);
    fdtd_domain.max_x = lx +1.1*(PAR.Rs-lx/2);
    fdtd_domain.min_y = 0  -1.1*(PAR.Rs-ly/2);
    fdtd_domain.max_y = ly +1.1*(PAR.Rs-ly/2);
    center_x = fdtd_domain.min_x+.5*(fdtd_domain.max_x-fdtd_domain.min_x);
    center_y = fdtd_domain.min_y+.5*(fdtd_domain.max_y-fdtd_domain.min_y);

    % Define source position
    theta = 0:PAR.dtheta:360;
    if theta(end) == 360
        theta(end) = [];
    end
    pos_x = center_x+PAR.Rs*cosd(theta(PAR.s_indx));
    pos_y = center_y+PAR.Rs*sind(theta(PAR.s_indx));
    impressed_J(1).min_x = pos_x-PAR.ls_x/2;
    impressed_J(1).max_x = pos_x+PAR.ls_x/2;
    impressed_J(1).min_y = pos_y-PAR.ls_y/2;
    impressed_J(1).max_y = pos_y+PAR.ls_y/2;

    % Defining source parameters
    impressed_J(1).direction = 'zp';
    impressed_J(1).magnitude = PAR.magnitude;
    impressed_J(1).waveform_type = 'gaussian';
    impressed_J(1).waveform_index = 1;
    
    % Defining Rx array data
    Nrx = round(360/PAR.dtheta);
    rx_fields = zeros(Nrx,numel(PAR.f));
    rx_x = center_x+PAR.Rs*cosd(theta);
    rx_y = center_y+PAR.Rs*sind(theta);

    % Determine the problem space boundaries including air buffers
    fdtd_domain.min_x = fdtd_domain.min_x ...
        - dx * boundary.air_buffer_number_of_cells_xn;
    fdtd_domain.min_y = fdtd_domain.min_y ...
        - dy * boundary.air_buffer_number_of_cells_yn;
    fdtd_domain.max_x = fdtd_domain.max_x ...
        + dx * boundary.air_buffer_number_of_cells_xp;
    fdtd_domain.max_y = fdtd_domain.max_y ...
        + dy * boundary.air_buffer_number_of_cells_yp;

    % Determine the problem space boundaries including cpml layers
    if strcmp(boundary.type_xn, 'cpml') && ...
            (boundary.cpml_number_of_cells_xn>0)
        fdtd_domain.min_x = fdtd_domain.min_x ...
            - dx * boundary.cpml_number_of_cells_xn;
    end
    if strcmp(boundary.type_xp, 'cpml') && ...
            (boundary.cpml_number_of_cells_xp>0)
        fdtd_domain.max_x = fdtd_domain.max_x ...
            + dx * boundary.cpml_number_of_cells_xp;
    end
    if strcmp(boundary.type_yn, 'cpml') && ...
            (boundary.cpml_number_of_cells_yn>0)
        fdtd_domain.min_y = fdtd_domain.min_y ...
            - dy * boundary.cpml_number_of_cells_yn;
    end
    if strcmp(boundary.type_yp, 'cpml') && ...
            (boundary.cpml_number_of_cells_yp>0)
        fdtd_domain.max_y = fdtd_domain.max_y ...
            + dy * boundary.cpml_number_of_cells_yp;
    end

    % Determining the problem space size
    fdtd_domain.size_x = fdtd_domain.max_x - fdtd_domain.min_x;
    fdtd_domain.size_y = fdtd_domain.max_y - fdtd_domain.min_y;

    % Number of cells in x and y directions
    nx = round(fdtd_domain.size_x/dx);
    ny = round(fdtd_domain.size_y/dy);

    % Adjust domain size by snapping to cells
    fdtd_domain.size_x = nx * dx;
    fdtd_domain.size_y = ny * dy;
    fdtd_domain.max_x = fdtd_domain.min_x + fdtd_domain.size_x;
    fdtd_domain.max_y = fdtd_domain.min_y + fdtd_domain.size_y;

    % Some frequently used auxiliary parameters
    nxp1 = nx+1;  nxm1 = nx-1; nxm2 = nx-2;
    nyp1 = ny+1;  nym1 = ny-1; nym2 = ny-2;

    % Create arrays storing the center coordinates of the cells
    fdtd_domain.cell_center_coordinates_x = zeros(nx,ny);
    fdtd_domain.cell_center_coordinates_y = zeros(nx,ny);
    for ind = 1:nx
        fdtd_domain.cell_center_coordinates_x(ind,:) = ...
            (ind - 0.5) * dx + fdtd_domain.min_x;
    end
    for ind = 1:ny
        fdtd_domain.cell_center_coordinates_y(:,ind) = ...
            (ind - 0.5) * dy + fdtd_domain.min_y;
    end

    % x and y coordinates for epsr and Ez matrices
    xcoor = round(linspace(fdtd_domain.min_x,fdtd_domain.max_x,nxp1),4);
    ycoor = round(linspace(fdtd_domain.min_y,fdtd_domain.max_y,nyp1),4);

    % TMz components
    eps_r_z     = DATA.epsrb*ones (nxp1, nyp1);
    mu_r_x      = ones (nxp1, ny);
    mu_r_y      = ones (nx  , nyp1);
    sigma_e_z   = DATA.sigb*ones(nxp1, nyp1);
    sigma_m_x   = zeros(nxp1, ny);
    sigma_m_y   = zeros(nx  , nyp1);
    
    % Adding user's domain
    eps_r_z     (xcoor>0 & xcoor<=lx,ycoor>0 & ycoor<=ly) = DATA.epsr;
    sigma_e_z   (xcoor>0 & xcoor<=lx,ycoor>0 & ycoor<=ly) = DATA.sig;

    % time array
    time = ((1:number_of_time_steps)-0.5)*dt;

    % Create and initialize field and current arrays
    Hx = zeros(nxp1,ny);
    Hy = zeros(nx,nyp1);
    Ez = zeros(nxp1,nyp1);

    number_of_impressed_J = size(impressed_J,2);
    number_of_impressed_M = size(impressed_M,2);

    % initialize waveforms
    % initialize sinusoidal waveforms
    if isfield(waveforms,'sinusoidal')
        for ind=1:size(waveforms.sinusoidal,2)
            waveforms.sinusoidal(ind).waveform = ...
                sin(2 * pi * waveforms.sinusoidal(ind).frequency * time);
        end
    end

    % initialize unit step waveforms
    if isfield(waveforms,'unit_step')
        for ind=1:size(waveforms.unit_step,2)
            start_index = waveforms.unit_step(ind).start_time_step;
            waveforms.unit_step(ind).waveform(1:number_of_time_steps) = 1;
            waveforms.unit_step(ind).waveform(1:start_index-1) = 0;
        end
    end

    % initialize Gaussian waveforms
    if isfield(waveforms,'gaussian')
        for ind=1:size(waveforms.gaussian,2)
            if waveforms.gaussian(ind).number_of_cells_per_wavelength == 0
                nc = number_of_cells_per_wavelength;
            else
                nc = waveforms.gaussian(ind).number_of_cells_per_wavelength;
            end
            waveforms.gaussian(ind).maximum_frequency = c/(nc*max([dx,dy]));
%             tau = (nc*max([dx,dy]))/(2*c);
            tau = sqrt(2.3)/pi/PAR.fwv;
%             tau = sqrt(2.3)/pi/max(PAR.f);
            waveforms.gaussian(ind).tau = tau;
            t_0 = 4.5 * waveforms.gaussian(ind).tau;
            waveforms.gaussian(ind).t_0 = t_0;
            waveforms.gaussian(ind).waveform = exp(-((time - t_0)/tau).^2);
        end
    end

    % initialize derivative of Gaussian waveforms
    if isfield(waveforms,'derivative_gaussian')
        for ind=1:size(waveforms.derivative_gaussian,2)
            if waveforms.derivative_gaussian(ind).number_of_cells_per_wavelength == 0
                nc = number_of_cells_per_wavelength;
            else
                nc = ...
                    waveforms.derivative_gaussian(ind).number_of_cells_per_wavelength;
            end
            waveforms.derivative_gaussian(ind).maximum_frequency = ...
                c/(nc*max([dx,dy]));
            tau = (nc*max([dx,dy]))/(2*c);
            waveforms.derivative_gaussian(ind).tau = tau;
            t_0 = 4.5 * waveforms.derivative_gaussian(ind).tau;
            waveforms.derivative_gaussian(ind).t_0 = t_0;
            waveforms.derivative_gaussian(ind).waveform = ...
                -(sqrt(2*exp(1))/tau)*(time - t_0).*exp(-((time - t_0)/tau).^2);
        end
    end

    % initialize cosine modulated Gaussian waveforms
    if isfield(waveforms,'cosine_modulated_gaussian')
        for ind=1:size(waveforms.cosine_modulated_gaussian,2)
            frequency = ...
                waveforms.cosine_modulated_gaussian(ind).modulation_frequency;
            tau = 0.966/waveforms.cosine_modulated_gaussian(ind).bandwidth;
            waveforms.cosine_modulated_gaussian(ind).tau = tau;
            t_0 = 4.5 * waveforms.cosine_modulated_gaussian(ind).tau;
            waveforms.cosine_modulated_gaussian(ind).t_0 = t_0;
            waveforms.cosine_modulated_gaussian(ind).waveform = ...
                cos(2*pi*frequency*(time - t_0)).*exp(-((time - t_0)/tau).^2);
        end
    end

    % electric current sources
    for ind = 1:number_of_impressed_J
        is = floor((impressed_J(ind).min_x - fdtd_domain.min_x)/dx)+1;
        js = floor((impressed_J(ind).min_y - fdtd_domain.min_y)/dy)+1;
        ie = floor((impressed_J(ind).max_x - fdtd_domain.min_x)/dx)+1;
        je = floor((impressed_J(ind).max_y - fdtd_domain.min_y)/dy)+1;
        impressed_J(ind).is = is;
        impressed_J(ind).js = js;
        impressed_J(ind).ie = ie;
        impressed_J(ind).je = je;

        if strcmp(impressed_J(ind).direction(2),'n')
            j_magnitude_factor = -1*impressed_J(ind).magnitude;
        else
            j_magnitude_factor =  impressed_J(ind).magnitude;
        end

        % copy waveform of the waveform type to waveform of the source
        wt_str = impressed_J(ind).waveform_type;
        wi_str = num2str(impressed_J(ind).waveform_index);
        eval_str = ['a_waveform = waveforms.' wt_str '(' wi_str ').waveform;'];
        eval(eval_str);
        impressed_J(ind).waveform = j_magnitude_factor * a_waveform;
    end

    % magnetic current sources
    for ind = 1:number_of_impressed_M
        is = round((impressed_M(ind).min_x - fdtd_domain.min_x)/dx)+1;
        js = round((impressed_M(ind).min_y - fdtd_domain.min_y)/dy)+1;
        ie = round((impressed_M(ind).max_x - fdtd_domain.min_x)/dx)+1;
        je = round((impressed_M(ind).max_y - fdtd_domain.min_y)/dy)+1;
        impressed_M(ind).is = is;
        impressed_M(ind).js = js;
        impressed_M(ind).ie = ie;
        impressed_M(ind).je = je;

        if strcmp(impressed_M(ind).direction(2),'n')
            m_magnitude_factor = -1*impressed_M(ind).magnitude;
        else
            m_magnitude_factor =  impressed_M(ind).magnitude;
        end

        % copy waveform of the waveform type to waveform of the source
        wt_str = impressed_M(ind).waveform_type;
        wi_str = num2str(impressed_M(ind).waveform_index);
        eval_str = ['a_waveform = waveforms.' wt_str '(' wi_str ').waveform;'];
        eval(eval_str);
        impressed_M(ind).waveform = m_magnitude_factor * a_waveform;
    end

    % Coeffiecients updating Ez
    Ceze  =  (2*eps_r_z*eps_0 - dt*sigma_e_z)./(2*eps_r_z*eps_0 + dt*sigma_e_z);
    Cezhy =  (2*dt/dx)./(2*eps_r_z*eps_0 + dt*sigma_e_z);
    Cezhx = -(2*dt/dy)./(2*eps_r_z*eps_0 + dt*sigma_e_z);

    % Coeffiecients updating Hx
    Chxh  =  (2*mu_r_x*mu_0 - dt*sigma_m_x)./(2*mu_r_x*mu_0 + dt*sigma_m_x);
    Chxez = -(2*dt/dy)./(2*mu_r_x*mu_0 + dt*sigma_m_x);

    % Coeffiecients updating Hy
    Chyh  =  (2*mu_r_y*mu_0 - dt*sigma_m_y)./(2*mu_r_y*mu_0 + dt*sigma_m_y);
    Chyez =  (2*dt/dx)./(2*mu_r_y*mu_0 + dt*sigma_m_y);

    % Initialize coeffiecients for impressed current sources
    for ind = 1:number_of_impressed_J
        is = impressed_J(ind).is;
        js = impressed_J(ind).js;
        ie = impressed_J(ind).ie;
        je = impressed_J(ind).je;
        switch (impressed_J(ind).direction(1))
            case 'x'
                impressed_J(ind).Cexj = -(2*dt)./ ...
                    (2*eps_r_x(is:ie-1,js:je)*eps_0 + dt*sigma_e_x(is:ie-1,js:je));
            case 'y'
                impressed_J(ind).Ceyj = -(2*dt)./ ...
                    (2*eps_r_y(is:ie,js:je-1)*eps_0 + dt*sigma_e_y(is:ie,js:je-1));
            case 'z'
                impressed_J(ind).Cezj = -(2*dt)./ ...
                    (2*eps_r_z(is:ie,js:je)*eps_0 + dt*sigma_e_z(is:ie,js:je));
        end
    end

    for ind = 1:number_of_impressed_M
        is = impressed_M(ind).is;
        js = impressed_M(ind).js;
        ie = impressed_M(ind).ie;
        je = impressed_M(ind).je;
        switch (impressed_M(ind).direction(1))
            case 'x'
                impressed_M(ind).Chxm = -(2*dt)./ ...
                    (2*mu_r_x(is:ie,js:je-1)*mu_0 + dt*sigma_m_x(is:ie,js:je-1));
            case 'y'
                impressed_M(ind).Chym = -(2*dt)./ ...
                    (2*mu_r_y(is:ie-1,js:je)*mu_0 + dt*sigma_m_y(is:ie-1,js:je));
            case 'z'
                impressed_M(ind).Chzm = -(2*dt)./ ...
                    (2*mu_r_z(is:ie-1,js:je-1)*mu_0 + dt*sigma_m_z(is:ie-1,js:je-1));
        end
    end

    is_cpml_xn = true;
    n_cpml_xn = abs(boundary.cpml_number_of_cells_xn);
    is_cpml_xp = true;
    n_cpml_xp = abs(boundary.cpml_number_of_cells_xp);
    is_cpml_yn = true;
    n_cpml_yn = abs(boundary.cpml_number_of_cells_yn);
    is_cpml_yp = true;
    n_cpml_yp = abs(boundary.cpml_number_of_cells_yp);
    is_any_side_cpml = true;

    % Call CPML initialization routine if any side is CPML
    if is_any_side_cpml
        % Initialize CPML boundary condition

        p_order = boundary.cpml_order; % order of the polynomial distribution
        sigma_ratio = boundary.cpml_sigma_factor;
        kappa_max = boundary.cpml_kappa_max;
        alpha_min = boundary.cpml_alpha_min;
        alpha_max = boundary.cpml_alpha_max;

        % Initialize cpml for xn region
        if is_cpml_xn

            % define one-dimensional temporary cpml parameter arrays
            sigma_max = sigma_ratio  * (p_order+1)/(150*pi*dx);
            ncells = n_cpml_xn;
            rho_e = ([ncells:-1:1]-0.75)/ncells;
            rho_m = ([ncells:-1:1]-0.25)/ncells;
            sigma_pex_xn = sigma_max * rho_e.^p_order;
            sigma_pmx_xn = sigma_max * rho_m.^p_order;
            sigma_pmx_xn = (mu_0/eps_0) * sigma_pmx_xn;
            kappa_ex_xn = 1 + (kappa_max - 1) * rho_e.^p_order;
            kappa_mx_xn = 1 + (kappa_max - 1) * rho_m.^p_order;
            alpha_ex_xn = alpha_min + (alpha_max - alpha_min) * (1-rho_e);
            alpha_mx_xn = alpha_min + (alpha_max - alpha_min) * (1-rho_m);
            alpha_mx_xn = (mu_0/eps_0) * alpha_mx_xn;

            % define one-dimensional cpml parameter arrays
            cpml_b_ex_xn = exp((-dt/eps_0) ...
                *((sigma_pex_xn./kappa_ex_xn)+ alpha_ex_xn));
            cpml_a_ex_xn = (1/dx)*(cpml_b_ex_xn-1.0).* sigma_pex_xn ...
                ./(kappa_ex_xn.*(sigma_pex_xn+kappa_ex_xn.*alpha_ex_xn));
            cpml_b_mx_xn = exp((-dt/mu_0) ...
                *((sigma_pmx_xn./kappa_mx_xn)+ alpha_mx_xn));
            cpml_a_mx_xn = (1/dx)*(cpml_b_mx_xn-1.0) .* sigma_pmx_xn ...
                ./(kappa_mx_xn.*(sigma_pmx_xn+kappa_mx_xn.*alpha_mx_xn));

            % Create and initialize cpml convolution parameters
            Psi_ezx_xn = zeros(ncells,nyp1);
            Psi_hyx_xn = zeros(ncells,nyp1);

            % Create and initialize cpml convolution coefficients
            CPsi_ezx_xn = Cezhy(2:ncells+1,:)*dx;
            CPsi_hyx_xn = Chyez(1:ncells,:)*dx;

            % Adjust FDTD coefficients in the CPML region
            for i = 1: ncells
                Cezhy(i+1,:) = Cezhy(i+1,:)/kappa_ex_xn(i);
                Chyez(i,:) = Chyez(i,:)/kappa_mx_xn(i);
            end

            % Delete temporary arrays. These arrays will not be used any more.
            clear sigma_pex_xn sigma_pmx_xn;
            clear kappa_ex_xn kappa_mx_xn;
            clear alpha_ex_xn alpha_mx_xn;
        end

        % Initialize cpml for xp region
        if is_cpml_xp

            % define one-dimensional temporary cpml parameter arrays
            sigma_max = sigma_ratio  * (p_order+1)/(150*pi*dx);
            ncells = n_cpml_xp;
            rho_e = ([1:ncells]-0.75)/ncells;
            rho_m = ([1:ncells]-0.25)/ncells;
            sigma_pex_xp = sigma_max * rho_e.^p_order;
            sigma_pmx_xp = sigma_max * rho_m.^p_order;
            sigma_pmx_xp = (mu_0/eps_0) * sigma_pmx_xp;
            kappa_ex_xp = 1 + (kappa_max - 1) * rho_e.^p_order;
            kappa_mx_xp = 1 + (kappa_max - 1) * rho_m.^p_order;
            alpha_ex_xp = alpha_min + (alpha_max - alpha_min) * (1-rho_e);
            alpha_mx_xp = alpha_min + (alpha_max - alpha_min) * (1-rho_m);
            alpha_mx_xp = (mu_0/eps_0) * alpha_mx_xp;

            % define one-dimensional cpml parameter arrays
            cpml_b_ex_xp = exp((-dt/eps_0) ...
                *((sigma_pex_xp./kappa_ex_xp)+ alpha_ex_xp));
            cpml_a_ex_xp = (1/dx)*(cpml_b_ex_xp-1.0).* sigma_pex_xp ...
                ./(kappa_ex_xp.*(sigma_pex_xp+kappa_ex_xp.*alpha_ex_xp));
            cpml_b_mx_xp = exp((-dt/mu_0) ...
                *((sigma_pmx_xp./kappa_mx_xp)+ alpha_mx_xp));
            cpml_a_mx_xp = (1/dx)*(cpml_b_mx_xp-1.0) .* sigma_pmx_xp ...
                ./(kappa_mx_xp.*(sigma_pmx_xp+kappa_mx_xp.*alpha_mx_xp));

            % Create and initialize cpml convolution parameters
            Psi_ezx_xp = zeros(ncells,nyp1);
            Psi_hyx_xp = zeros(ncells,nyp1);

            % Create and initialize cpml convolution coefficients
            CPsi_ezx_xp = Cezhy(nxp1-ncells:nx,:)*dx;
            CPsi_hyx_xp = Chyez(nxp1-ncells:nx,:)*dx;

            % Adjust FDTD coefficients in the CPML region
            for i = 1: ncells
                Cezhy(nx-ncells+i,:) = Cezhy(nx-ncells+i,:)/kappa_ex_xp(i);
                Chyez(nx-ncells+i,:) = Chyez(nx-ncells+i,:)/kappa_mx_xp(i);
            end

            % Delete temporary arrays. These arrays will not be used any more.
            clear sigma_pex_xp sigma_pmx_xp;
            clear kappa_ex_xp kappa_mx_xp;
            clear alpha_ex_xp alpha_mx_xp;
        end


        % Initialize cpml for yn region
        if is_cpml_yn

            % define one-dimensional temporary cpml parameter arrays
            sigma_max = sigma_ratio  * (p_order+1)/(150*pi*dy);
            ncells = n_cpml_yn;
            rho_e = ([ncells:-1:1]-0.75)/ncells;
            rho_m = ([ncells:-1:1]-0.25)/ncells;
            sigma_pey_yn = sigma_max * rho_e.^p_order;
            sigma_pmy_yn = sigma_max * rho_m.^p_order;
            sigma_pmy_yn = (mu_0/eps_0) * sigma_pmy_yn;
            kappa_ey_yn = 1 + (kappa_max - 1) * rho_e.^p_order;
            kappa_my_yn = 1 + (kappa_max - 1) * rho_m.^p_order;
            alpha_ey_yn = alpha_min + (alpha_max - alpha_min) * (1-rho_e);
            alpha_my_yn = alpha_min + (alpha_max - alpha_min) * (1-rho_m);
            alpha_my_yn = (mu_0/eps_0) * alpha_my_yn;

            % define one-dimensional cpml parameter arrays
            cpml_b_ey_yn = exp((-dt/eps_0) ...
                *((sigma_pey_yn./kappa_ey_yn)+ alpha_ey_yn));
            cpml_a_ey_yn = (1/dy)*(cpml_b_ey_yn-1.0).* sigma_pey_yn ...
                ./(kappa_ey_yn.*(sigma_pey_yn+kappa_ey_yn.*alpha_ey_yn));
            cpml_b_my_yn = exp((-dt/mu_0) ...
                *((sigma_pmy_yn./kappa_my_yn)+ alpha_my_yn));
            cpml_a_my_yn = (1/dy)*(cpml_b_my_yn-1.0) .* sigma_pmy_yn ...
                ./(kappa_my_yn.*(sigma_pmy_yn+kappa_my_yn.*alpha_my_yn));

            % Create and initialize cpml convolution parameters
            Psi_ezy_yn = zeros(nxp1,ncells);
            Psi_hxy_yn = zeros(nxp1,ncells);

            % Create and initialize cpml convolution coefficients
            CPsi_ezy_yn = Cezhx(:,2:ncells+1)*dy;
            CPsi_hxy_yn = Chxez(:,1:ncells)*dy;

            % Adjust FDTD coefficients in the CPML region
            for i = 1: ncells
                Cezhx(:,i+1) = Cezhx(:,i+1)/kappa_ey_yn(i);
                Chxez(:,i) = Chxez(:,i)/kappa_my_yn(i);
            end

            % Delete temporary arrays. These arrays will not be used any more.
            clear sigma_pey_yn sigma_pmy_yn;
            clear kappa_ey_yn kappa_my_yn;
            clear alpha_ey_yn alpha_my_yn;
        end

        % Initialize cpml for yp region
        if is_cpml_yp

            % define one-dimensional temporary cpml parameter arrays
            sigma_max = sigma_ratio  * (p_order+1)/(150*pi*dy);
            ncells = n_cpml_yp;
            rho_e = ([1:ncells]-0.75)/ncells;
            rho_m = ([1:ncells]-0.25)/ncells;
            sigma_pey_yp = sigma_max * rho_e.^p_order;
            sigma_pmy_yp = sigma_max * rho_m.^p_order;
            sigma_pmy_yp = (mu_0/eps_0) * sigma_pmy_yp;
            kappa_ey_yp = 1 + (kappa_max - 1) * rho_e.^p_order;
            kappa_my_yp = 1 + (kappa_max - 1) * rho_m.^p_order;
            alpha_ey_yp = alpha_min + (alpha_max - alpha_min) * (1-rho_e);
            alpha_my_yp = alpha_min + (alpha_max - alpha_min) * (1-rho_m);
            alpha_my_yp = (mu_0/eps_0) * alpha_my_yp;

            % define one-dimensional cpml parameter arrays
            cpml_b_ey_yp = exp((-dt/eps_0) ...
                *((sigma_pey_yp./kappa_ey_yp)+ alpha_ey_yp));
            cpml_a_ey_yp = (1/dy)*(cpml_b_ey_yp-1.0).* sigma_pey_yp ...
                ./(kappa_ey_yp.*(sigma_pey_yp+kappa_ey_yp.*alpha_ey_yp));
            cpml_b_my_yp = exp((-dt/mu_0) ...
                *((sigma_pmy_yp./kappa_my_yp)+ alpha_my_yp));
            cpml_a_my_yp = (1/dy)*(cpml_b_my_yp-1.0) .* sigma_pmy_yp ...
                ./(kappa_my_yp.*(sigma_pmy_yp+kappa_my_yp.*alpha_my_yp));

            % Create and initialize cpml convolution parameters
            Psi_ezy_yp = zeros(nxp1,ncells);
            Psi_hxy_yp = zeros(nxp1,ncells);

            % Create and initialize cpml convolution coefficients
            CPsi_ezy_yp = Cezhx(:,nyp1-ncells:ny)*dy;
            CPsi_hxy_yp = Chxez(:,nyp1-ncells:ny,:)*dy;

            % Adjust FDTD coefficients in the CPML region
            for i = 1: ncells
                Cezhx(:,ny-ncells+i) = Cezhx(:,ny-ncells+i)/kappa_ey_yp(i);
                Chxez(:,ny-ncells+i) = Chxez(:,ny-ncells+i)/kappa_my_yp(i);
            end

            % Delete temporary arrays. These arrays will not be used any more.
            clear sigma_pey_yp sigma_pmy_yp;
            clear kappa_ey_yp kappa_my_yp;
            clear alpha_ey_yp alpha_my_yp;
        end
    end

    
    start_time      = cputime;
    current_time    = 0;
    et              = zeros([size(Ez),length(PAR.f)]);
    
    % Constant for Fourier transform
    
    sample = [];
    
    % Starting the simulation
    for time_step = 1:number_of_time_steps

        current_time  = current_time + dt/2;

        Hx(1:nxp1,1:ny) = Chxh(1:nxp1,1:ny) .* Hx(1:nxp1,1:ny) ...
            + Chxez(1:nxp1,1:ny) .* (Ez(1:nxp1,2:nyp1)-Ez(1:nxp1,1:ny));

        Hy(1:nx,1:nyp1) = Chyh(1:nx,1:nyp1) .* Hy(1:nx,1:nyp1) ...
            + Chyez(1:nx,1:nyp1) .* (Ez(2:nxp1,1:nyp1)-Ez(1:nx,1:nyp1));

        for ind = 1:number_of_impressed_M
            is = impressed_M(ind).is;
            js = impressed_M(ind).js;
            ie = impressed_M(ind).ie;
            je = impressed_M(ind).je;
            switch (impressed_M(ind).direction(1))
                case 'x'
                    Hx(is:ie,js:je-1) = Hx(is:ie,js:je-1) ...
                        + impressed_M(ind).Chxm * ...
                        impressed_M(ind).waveform(time_step);
                case 'y'
                    Hy(is:ie-1,js:je) = Hy(is:ie-1,js:je) ...
                        + impressed_M(ind).Chym * ...
                        impressed_M(ind).waveform(time_step);
                case 'z'
                    Hz(is:ie-1,js:je-1) = Hz(is:ie-1,js:je-1) ...
                        + impressed_M(ind).Chzm * ...
                        impressed_M(ind).waveform(time_step);
            end
        end

        if is_any_side_cpml == false
            return;
        end

        % TMz
        % apply CPML to magnetic field components
        if is_cpml_xn
            for i = 1: n_cpml_xn
                Psi_hyx_xn(i,:,:) = cpml_b_mx_xn(i) * Psi_hyx_xn(i,:,:) ...
                    + cpml_a_mx_xn(i)*(Ez(i+1,:,:)-Ez(i,:,:));
            end
            Hy(1:n_cpml_xn,:,:) = Hy(1:n_cpml_xn,:,:) ...
                + CPsi_hyx_xn(:,:,:) .* Psi_hyx_xn(:,:,:);
        end

        if is_cpml_xp
            n_st = nx - n_cpml_xp;
            for i = 1:n_cpml_xp
                Psi_hyx_xp(i,:,:) = cpml_b_mx_xp(i) * Psi_hyx_xp(i,:,:) ...
                    + cpml_a_mx_xp(i)*(Ez(i+n_st+1,:,:)-Ez(i+n_st,:,:));
            end

            Hy(n_st+1:nx,:,:) = Hy(n_st+1:nx,:,:) ...
                + CPsi_hyx_xp(:,:,:) .* Psi_hyx_xp(:,:,:);
        end

        if is_cpml_yn
            for i = 1:n_cpml_yn
                Psi_hxy_yn(:,i,:) = cpml_b_my_yn(i) * Psi_hxy_yn(:,i,:) ...
                    + cpml_a_my_yn(i)*(Ez(:,i+1,:)-Ez(:,i,:));
            end
            Hx(:,1:n_cpml_yn,:) = Hx(:,1:n_cpml_yn,:) ...
                + CPsi_hxy_yn(:,:,:) .* Psi_hxy_yn(:,:,:);
        end

        if is_cpml_yp
            n_st = ny - n_cpml_yp;
            for i = 1:n_cpml_yp
                Psi_hxy_yp(:,i,:) = cpml_b_my_yp(i) * Psi_hxy_yp(:,i,:) ...
                    + cpml_a_my_yp(i)*(Ez(:,i+n_st+1,:)-Ez(:,i+n_st,:));
            end
            Hx(:,n_st+1:ny,:) = Hx(:,n_st+1:ny,:) ...
                + CPsi_hxy_yp(:,:,:) .* Psi_hxy_yp(:,:,:);
        end

        current_time  = current_time + dt/2;

        Ez(2:nx,2:ny) = ...
            Ceze(2:nx,2:ny).*Ez(2:nx,2:ny) ...
            + Cezhy(2:nx,2:ny) ...
            .* (Hy(2:nx,2:ny)-Hy(1:nxm1,2:ny)) ...
            + Cezhx(2:nx,2:ny) ...
            .* (Hx(2:nx,2:ny)-Hx(2:nx,1:nym1));

        for ind = 1:number_of_impressed_J
            is = impressed_J(ind).is;
            js = impressed_J(ind).js;
            ie = impressed_J(ind).ie;
            je = impressed_J(ind).je;
            switch (impressed_J(ind).direction(1))
                case 'x'
                    Ex(is:ie-1,js:je) = Ex(is:ie-1,js:je) ...
                        + impressed_J(ind).Cexj * ...
                        impressed_J(ind).waveform(time_step);
                case 'y'
                    Ey(is:ie,js:je-1) = Ey(is:ie,js:je-1) ...
                        + impressed_J(ind).Ceyj * ...
                        impressed_J(ind).waveform(time_step);
                case 'z'
                    Ez(is:ie,js:je) = Ez(is:ie,js:je) ...
                        + impressed_J(ind).Cezj * ...
                        impressed_J(ind).waveform(time_step);
            end
        end

        % update electric fields at the CPML regions
        if is_any_side_cpml == false
            return;
        end


        % update electric fields at the PML regions
        % TMz
        % apply CPML to electric field components
        if is_cpml_xn
            for i = 1:n_cpml_xn
                Psi_ezx_xn(i,:,:) = cpml_b_ex_xn(i) * Psi_ezx_xn(i,:,:) ...
                    + cpml_a_ex_xn(i)*(Hy(i+1,:,:)-Hy(i,:,:));
            end
            Ez(2:n_cpml_xn+1,:,:) = Ez(2:n_cpml_xn+1,:,:) ...
                + CPsi_ezx_xn .* Psi_ezx_xn;
        end

        if is_cpml_xp
            n_st = nx - n_cpml_xp;
            for i = 1:n_cpml_xp
                Psi_ezx_xp(i,:,:) = cpml_b_ex_xp(i) * Psi_ezx_xp(i,:,:) ...
                    + cpml_a_ex_xp(i)*(Hy(i+n_st,:,:)-Hy(i+n_st-1,:,:));
            end
            Ez(n_st+1:nx,:,:) = Ez(n_st+1:nx,:,:) ...
                + CPsi_ezx_xp .* Psi_ezx_xp;
        end

        if is_cpml_yn
            for i = 1:n_cpml_yn
                Psi_ezy_yn(:,i,:) = cpml_b_ey_yn(i) * Psi_ezy_yn(:,i,:) ...
                    + cpml_a_ey_yn(i)*(Hx(:,i+1,:)-Hx(:,i,:));
            end
            Ez(:,2:n_cpml_yn+1,:) = Ez(:,2:n_cpml_yn+1,:) ...
                + CPsi_ezy_yn .* Psi_ezy_yn;
        end

        if is_cpml_yp
            n_st = ny - n_cpml_yp;
            for i = 1:n_cpml_yp
                Psi_ezy_yp(:,i,:) = cpml_b_ey_yp(i) * Psi_ezy_yp(:,i,:) ...
                    + cpml_a_ey_yp(i)*(Hx(:,i+n_st,:)-Hx(:,i+n_st-1,:));
            end
            Ez(:,n_st+1:ny,:) = Ez(:,n_st+1:ny,:) ...
                + CPsi_ezy_yp .* Psi_ezy_yp;
        end


        for f_indx = 1:length(PAR.f)
            et(:,:,f_indx) =  et(:,:,f_indx) + ...
                Ez*exp(-1j*2*pi*PAR.f(f_indx)*...
                (time_step)*dt)*dt;
        end

    end
        
    for f_indx = 1:length(PAR.f)
        rx_fields(:,f_indx) = interp2(ycoor,xcoor,et(:,:,f_indx),rx_y,rx_x);
    end
    
    % Calculating the computational time
    end_time = cputime;
    total_time_in_minutes = (end_time - start_time)/60;

    % Setting the total electric field according to user's domain
    et = et(xcoor>0&xcoor<=lx,ycoor>0&ycoor<=ly,:);
    x = xcoor(xcoor>0&xcoor<=lx);
    y = ycoor(ycoor>0&ycoor<=ly);
    
    if FLAG_RESOLUTION == true
        
        newx = (1:old_nxobj)*old_dx;
        newy = (1:old_nyobj)*old_dy;
        [y2,x2] = meshgrid(newy,newx);
        
        [I,J] = size(et);
        x1 = (1:I)*dx;
        y1 = (1:J)*dy;
        
        for f_indx = 1:length(PAR.f)
            et(:,:,f_indx) = interp2(y1,x1,squeeze(et(:,:,f_indx)),y2,x2);
        end
        
        x = linspace(x(1),x(end),old_nxobj);
        y = linspace(y(1),y(end),old_nyobj);

    end

end

function [data] = changemesh (data,newdx,newdy)

    [I,J] = size(data.epsr);
    dx = data.dx;
    dy = data.dy;
    
    newdx = dx/(ceil(dx/newdx)+1);
    newdy = dy/(ceil(dy/newdy)+1);
    
    newx = dx:newdx:I*dx;
    newy = dy:newdy:J*dy;
    x = (1:I)*dx;
    y = (1:J)*dy;
    
    [yp,xp] = meshgrid(newy,newx);
    data.epsr = interp2(y,x,data.epsr,yp,xp);
    data.sig = interp2(y,x,data.sig,yp,xp);
    data.dx = newdx;
    data.dy = newdy;
end