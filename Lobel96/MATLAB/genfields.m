%GENFIELDS Generation fields routine
%   This script implements the data generation for the Gradient Conjugated
%   Method (Lobel et al., 1996). It computes the incident, total and
%   scattered field as well as the Green function for both domains.
%
%   Implemented by:
% 
%   Andre Costa Batista
%   Universidade Federal de Minas Gerais
%
%   References
%
%   Lobel, P., et al. "Conjugate gradient method for solving inverse
%   scattering with experimental data." IEEE Antennas and Propagation
%   Magazine 38.3 (1996): 48-51.

clear, close all

% Domain parameters
I = 50;                                 % Number of cells in x-axis
J = 50;                                 % Number of cells in y-axis

% DATA struct
data.dx     = 5e-3;                     % Cell size in x-axis [m]
data.dy     = 5e-3;                     % Cell size in y-axis [m]
data.epsr   = ones(I,J);                % Relative permittivity data
data.sig    = zeros(size(data.epsr));   % Conductivity data [S/m]
data.epsrb  = 1;                        % Background relative permittivity
data.sigb   = 0;                        % Background conductivity [S/m]

% PAR struct
par.f       = 800e6;                    % Linear frequency of measurements [Hz]
par.nts     = 2^12;                     % Number of FDTD time steps
par.nsps    = 7;                        % Number of sources per side
data.lambda = 1/par.f/sqrt(4e-7*pi*...  % Wavelength [m]
    data.epsrb*8.854187817e-12);
data.kb     = 2*pi*par.f*sqrt(4e-7*...  % Wavenumber of background [1/m]
    pi*data.epsrb*8.854187817e-12);
data.f      = par.f;
data.I      = I;
data.J      = J;

% Incident field matrix
ei          = zeros(I,J,par.nsps*4);

% Grid variables
x = cell(4*par.nsps,1);
y = cell(4*par.nsps,1);

% Auxiliar parameter variable
aux_par = cell(4*par.nsps,1);
aux_par(:) = {par};

% Computing incident field
parfor s_indx = 1:par.nsps*4
    
    aux_par{s_indx}.s_indx = s_indx;
    [ei(:,:,s_indx),x{s_indx},y{s_indx}] = fdtd2d(data,aux_par{s_indx});

end

% Defining object
% data.epsr(round(I/2)-round(.1*I):round(I/2)+round(.1*I),...
%     round(J/2)-round(.1*J):round(J/2)+round(.1*J)) = 1.2;
data.epsr(round(I/4)-.1*I:round(I/4)+.1*I,...
   round(J/4)-.1*J:round(J/4)+.1*J) = 5.0;
% data.sig(round(I/4)-.1*I:round(I/4)+.1*I,...
%    round(J/4)-.1*J:round(J/4)+.1*J) = 0.02;

% Total electric field matrix
et          = zeros(I,J,par.nsps*4);

% Calling FDTD
parfor s_indx = 1:par.nsps*4
    
    aux_par{s_indx}.s_indx = s_indx;
    [et(:,:,s_indx),x{s_indx},y{s_indx}] = fdtd2d(data,aux_par{s_indx});

end

% Resetting grid variables
x = x{1};
y = y{1};

% Scattered field computation
es          = et-ei;
nmps        = 7;                        % Number of measurements per side

% Indexes computation
im          = round(linspace(1,I,nmps+2));
im([1,end]) = [];
im          = [im,I*ones(1,nmps),im,ones(1,nmps)]';
jm          = round(linspace(1,J,nmps+2));
jm([1,end]) = [];
jm          = [ones(1,nmps),jm,J*ones(1,nmps),jm]';
es          = reshape(permute(es,[3,1,2]),par.nsps*4,[]);
es          = permute(es(:,sub2ind([I,J],im,jm)),[2,1]);

% Measurements points [m]
r           = [x(im)',y(jm)'];
data.r = r;

% Model parameters
L       = par.nsps*4;                   % Number of sources
M       = nmps*4;                       % Number of measurements
N       = I*J;                          % Number of variables
deltasn = data.dx*data.dy;              % Surface element
an      = sqrt(deltasn/pi);             % Green function parameter
kb      = 2*pi*par.f*sqrt(pi*4e-7*...   % Wavenumber [1/m]
    8.854187817e-12*data.epsrb);

% Reshaping data
ei = permute(reshape(permute(ei,[3,1,2]),L,N),[2,1]);
et = permute(reshape(permute(et,[3,1,2]),L,N),[2,1]);

% Computing radius from S
R = zeros(M,N);
for m = 1:M
    n = 1;
    for j = 1:J
        for i = 1:I
            R(m,n) = sqrt((r(m,1)-x(i))^2+(r(m,2)-y(j))^2);
            n = n + 1;
        end
    end
end

% Computing Green function for S
gs          = zeros(size(R));
gs(R~=0)    = 1j/2*pi*kb*an*besselj(1,kb*an)*besselh(0,2,kb*R(R~=0));
gs(R==0)    = (1j/2)*(pi*kb*an*besselh(1,2,kb*an)-2j);
gs          = -gs;

% Computing radius from D
R = zeros(N,N);
m = 1;
for jm = 1:J
    for im = 1:I
        n = 1;
        for j = 1:J
            for i = 1:I
                R(m,n) = sqrt((x(im)-x(i))^2+(y(jm)-y(j))^2);
                n = n + 1;
            end
        end
        m = m + 1;
    end
end

% Computing Green function for D
gd          = zeros(size(R));
gd(R~=0)    = 1j/2*pi*kb*an*besselj(1,kb*an)*besselh(0,2,kb*R(R~=0));
gd(R==0)    = (1j/2)*(pi*kb*an*besselh(1,2,kb*an)-2j);
gd          = -gd;

% Others parameters
data.x = x;
data.y = y;

% Saving data
save data.mat   data par
save ei.mat     ei
save es.mat     es
save et.mat     et
save grid.mat   x y
save gs.mat     gs
save ../../../../../../../Documents/MATLAB/inverse-approximation/gd.mat gd -v7.3