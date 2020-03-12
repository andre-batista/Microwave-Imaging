function [model] = setmodel(varargin)
%SETMODEL Build inverse electromagnetic scattering model
%   model = setmodel(Nx,Ny,dx,dy,epsrb,sigb,frequencies,frequency_waveform,
%   Rs,dtheta,time_window,ls_x,ls_y,I,mi) creates a struct with essential 
%   parameters for solving an inverse electromagnetic scattering problem.
%   
%   In this case, the inputs mean:
%
%       - Nx, Ny: number of pixels in each axis.
%       - dx, dy: cell size in each axis [m]
%       - epsrb: background relative permittivity
%       - sigb: background conductivity [S/m]
%       - frequencies: measurement frequencies in which the equations will
%         be solved
%       - frequency_waveform: central frequency of the waveform used in the
%         time-domain foward solver.
%       - Rs: radius from the denter of the image to the circular array of 
%         sources [m]
%       - dtheta: angular space among sources in the circular array [deg]
%       - time_window: interval of simulation for time-domain foward solver
%         [sec]
%       - ls_x, ls_y: length of the square sources in each axis [m]
%       - I: impressed current for each source [A]
%       - mi [OPTIONAL]: an array containing all the indexes (linear way) 
%         of the image in which it is assumed background. Therefore, if in 
%         your experiment you have a priori knowledge of a set of pixels 
%         which are background, you may defined them here and the algorithm
%         will assume it when solving the equations. Default is [].
%
%   model = setmodel(model,nx,ny) creates a model with a new number of 
%   pixels. Given an initial model struct, it computes its mesh with new
%   number of elements in each axis (nx, ny). It also updates background
%   indexes, if present.
%
%   Implemented by Andre Costa Batista (Universidade Federal de Minas
%   Gerais, Brazil). Contact: andre-costa@ufmg.br

    if nargin == 14 || nargin == 15
    
        % Electromagnetic constants
        eps0 = 8.85418782e-12;
        mu0 = 4*pi*1e-7;
        
        % Mesh parameters
        model.I = varargin{1}; % Nx
        model.J = varargin{2}; % Ny
        model.dx = varargin{3}; % dx
        model.dy = varargin{4}; % dy;
        model.x = [];
        model.y = [];
        model.r = [];
        
        % Background parameters
        model.epsrb = varargin{5}; % epsrb
        model.sigb = varargin{6}; % sigb
        
        % Frequency parameters
        model.f = varargin{7}; % frequencies
        model.fwv = varargin{8}; % frequency_waveform
        model.lambda = 1./model.f/sqrt(mu0*model.epsrb*eps0);
        model.kb = 2*pi*model.f*sqrt(mu0*eps0*model.epsrb);
        
        % Source array parameters
        model.Rs = varargin{9}; % Rs
        model.dtheta = varargin{10}; % dtheta
        model.tw = varargin{11}; % time_window
        model.ls_x = varargin{12}; % ls_x
        model.ls_y = varargin{13}; % ls_y;
        model.magnitude = varargin{14}/(model.ls_x+model.dx)...
            /(model.ls_y+model.dy); % I
        
        if nargin == 15
            model.mi = varargin{15};
        else
            model.mi = [];
        end
        
    elseif nargin == 3
        
        % Inputs 
        model = varargin{1};
        nx = varargin{2};
        ny = varargin{3};
        
        % Recovering current
        current = model.magnitude*(model.ls_x+model.dx)*...
            (model.ls_y+model.dy);
        
        % Mesh computation
        newdx = model.dx*model.I/nx;
        newdy = model.dy*model.J/ny;
        
        % Updating background indexes
        if ~isempty(model.mi)
            background = zeros([model.I,model.J]);
            background(model.mi) = 1;
            x = model.dx:model.dx:model.I*model.dx;
            y = model.dy:model.dy:model.J*model.dy;
            newx = newdx:newdx:nx*newdx;
            newy = newdy:newdy:ny*newdy;
            [yp,xp] = meshgrid(newy,newx);
            background = interp2(y,x,background,yp,xp);
            background(background>0.5) = 1;
            model.mi = find(background==1);
        end
              
        % Updating mesh parameters
        model.dx = newdx;
        model.dy = newdy;
        model.x = [];
        model.y = [];
        model.I = nx;
        model.J = ny;
        
        % Updating source parameters
        model.ls_x = floor(model.ls_x/model.dx)*model.dx;
        model.ls_y = floor(model.ls_y/model.dy)*model.dy;
        model.magnitude = current/(model.ls_x+model.dx)...
            /(model.ls_y+model.dy); % I

        
%     elseif nargin == 3 || nargin == 5
%     
%         % Input variables
%         model = varargin{1};
%         nx = varargin{2};
%         ny = varargin{3};
%         epsilon_r = varargin{4};
%         sigma = varargin{5};
%         current = model.magnitude*(model.ls_x+model.dx)*(model.ls_y+model.dy);
%         
%         if ~isempty(model.mi)
%             background = zeros(size(epsilon_r));
%             background(model.mi) = 1;
%         end
%         
%         % Compute new grid
%         Lx = model.dx*model.I;
%         Ly = model.dy*model.J;
%         newdx = Lx/nx;
%         newdy = Ly/ny;
%         x = model.dx:model.dx:model.I*model.dx;
%         y = model.dy:model.dy:model.J*model.dy;
%         newx = newdx:newdx:nx*newdx;
%         newy = newdy:newdy:ny*newdy;
%         [yp,xp] = meshgrid(newy,newx);
%         epsilon_r = interp2(y,x,epsilon_r,yp,xp);
%         sigma = interp2(y,x,sigma,yp,xp);
%         model.dx = newdx;
%         model.dy = newdy;
%         model.x = [];
%         model.y = [];
%         model.I = nx;
%         model.J = ny;
%         
%         model.ls_x = floor(model.ls_x/model.dx)*model.dx;
%         model.ls_y = floor(model.ls_y/model.dy)*model.dy;
%         model.magnitude = current/(model.ls_x+model.dx)...
%             /(model.ls_y+model.dy); % I
%         
%         if nargin == 5 && ~isempty(model.mi)
%             background = interp2(y,x,background,yp,xp);
%             background(background>1) = 0;
%             model.mi = find(background==1);
% %             model.mi = find(epsilon_r(:)==model.epsrb & sigma(:)==model.sigb);
%         elseif nargin ~= 5
%             disp('ERROR SETMODEL: original maps are missing at the input!')
%         end 
        
    else
        disp('ERROR SETMODEL: Incorrect number of inputs!')
        return
    end
end