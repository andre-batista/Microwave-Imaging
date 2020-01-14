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

disp(' ')
disp('========== The Conjugated Gradient Method ==========')
clear, close all

% Loading inputs
load ./genfields/es.mat
load ./genfields/ei.mat
load ../../../../../../Documents/MATLAB/inverse-approximation/gd.mat
load ./genfields/gs.mat
load ./genfields/data.mat

% General Parameters
maxit   = 1e3;                      % Number of iterations
[M,L]   = size(es);                 % M measurements, L sources
[N,~]   = size(ei);                 % N points within the mesh
dS      = data.dx*data.dy;          % Surface element [m^2]
eps0    = 8.85418782e-12;           % Vaccum permittivity [F/m]
omega   = 2*pi*data.f;              % Angular frequency [rad/sec]

% How do you preffer the initial solution?
% 1 - Everything background
% 2 - Backpropagation method (Lobel et al., 1996)
% 3 - Exact solution
% 4 - Load last run
initopt = 2;

switch initopt
    
    case 1
        C = sparse(diag(zeros([N,1])));
        d = zeros(N,1);
        g = ones(N,1);
    
    case 2
        gamma = norm(reshape(gs'*es,[],1))^2/...
            norm(reshape(gs*gs'*es,[],1))^2;
        w0 = gamma*gs'*es;
        C = sparse(diag(1/L*sum(w0./ei,2)));
        d = zeros(N,1);
        g = ones(N,1);
    
    case 3
        C = sparse(diag((...
            data.epsr(:)-1j*data.sig(:)/omega/eps0/data.epsrb)...
            -(data.epsrb-1j*data.sigb/omega/eps0/data.epsrb)));
        d = zeros(N,1);
        g = ones(N,1);
    
    otherwise
        load ../../../../../../Documents/MATLAB/inverse-approximation/f-ricardo.mat C g d   
end

% How do you preffer the choice of the alpha?
% 1 - (Lobel et al, 1996)
% 2 - Golden section method
alphaopt = 1;

% Auxiliar variables
dig = diag(diag(gd));
G_dig = gd-dig; % G - diagonal(G)
didi = dig*dig; % produto das diagonais de G
diGdi = dig*G_dig; % diagonal * (G - diagonal)
Gdidi = G_dig*dig;
Gdi2 = G_dig.^2;
GD = gd*C;
D2 = C*C;
diD = diag(C);

% Initializing variables
cnvg    = zeros(maxit+1,2);                     % Convergence data
I       = sparse(eye(N));                       % Identity matrix
LC      = inv(I-gd*C);                          % Initial inversion                    
rho     = es-gs*C*LC*ei;                        % Initial residual

% Printing first solution
disp(['Iteration: 0 - Cost function: ',...
    num2str(norm(rho(:))^2,'%.4e')])

if initopt ~= 2
    cnvg(1,:) = [norm(rho(:))^2,norm(g)];
else
    cnvg(1,:) = [norm(rho(:))^2,.0];
end

totaltime = cputime;

% Iterations
for it = 1:maxit
    
    tic
    
    % Computing the gradient    
    gradJ = reshape(-2*conj(sparse(1:N*L,1:N*L,reshape((LC*ei),[],1)) * ...
        repmat(LC,L,1))*gs'*rho,N,[]);
    gradJ = sum(gradJ(:,(1:L:L^2)+(0:L-1)),2);
    
    g_last = g;
    g = -gradJ;
    
    % Computing the optimum direction
    d = g + inner(g,g-g_last,dS)/norm(g_last)^2*d;
    D = sparse(diag(d));

    % Computing v matrix
	v = gs*LC.'*D*LC*ei;
    
    % Computing step
    switch alphaopt
        case 1
            alpha = 0;
            for l = 1:L
                alpha = alpha + inner(rho(:,l),v(:,l),data.dx);
            end
            alpha = alpha/norm(v(:))^2;
        otherwise
            alpha = gsmethod(rho,v);
    end
    
    % Computing next contrast
    C = C + alpha*D;
    
    % Computing the inverse matrix
    if mod(it,5) == 0
        LC = inv(I-gd*C);
    else
        GD = gd*C;
        D2 = C*C;
        diD = diag(C);
        LC      = I + GD ...
            + diag((Gdi2*diD).*diD) ...
            + C*diGdi*C ...
            + didi*D2 ...
            + Gdidi*D2;
    end
    
    % Computing the residual
    rho = es-gs*C*LC*ei;
    
    % Computing the objective function
    J = norm(rho(:))^2;
    t = toc;
    
    % Printing iteration
    disp(['Iteration: ',num2str(it),' - Cost function: ',num2str(J,'%.2e'),...
        ' - norm(g): ',num2str(norm(g),'%.2e'), ' - time: ',num2str(t,'%.1f'),' sec'])
    
    % Saving objetive function and gradient magnitude
    cnvg(it+1,:) = [J,norm(g)];
    
end

totaltime = cputime-totaltime;
disp(['Total time: ',num2str(totaltime),' seconds'])

% Recovering dielectric properties
[I,J]   = size(data.epsr);                  % Image size
tau     = reshape(diag(C),I,J);             % Constrast fuction
epsr    = real(tau) + data.epsrb;           % Retrieved relative permittivity
sig     = -omega*eps0*data.epsrb*imag(tau); % Relative conductivity [S/m]

% Plotting results
figure
load ./genfields/grid.mat

% Relative permittivity plot
subplot(3,2,1)
imagesc(y,x,epsr')
set(gca,'YDir','normal')
xlabel('x [m]')
ylabel('y [m]')
title('Relative permittivity')
clb = colorbar;
ylabel(clb, '\epsilon_r')

% Conductivity plot
subplot(3,2,2)
imagesc(y,x,sig')
set(gca,'YDir','normal')
xlabel('x [m]')
ylabel('y [m]')
title('Conductivity')
clb = colorbar;
ylabel(clb, '\sigma [S/m]')

% Gradient - Real
subplot(3,2,3)
imagesc(y,x,real(reshape(g,I,J))')
set(gca,'YDir','normal')
xlabel('x [m]')
ylabel('y [m]')
title('Gradient - Real')
clb = colorbar;
ylabel(clb, 'g')

% Conductivity plot
subplot(3,2,4)
imagesc(y,x,imag(reshape(g,I,J))')
set(gca,'YDir','normal')
xlabel('x [m]')
ylabel('y [m]')
title('Gradient - Imaginary')
clb = colorbar;
ylabel(clb, 'g')

% Convergence plot - Cost Function
subplot(3,2,5)
plot(0:maxit,cnvg(:,1),'linewidth',2)
grid
xlabel('Iterations')
ylabel('J(C)')
title('Cost Function')

% Convergence plot - Gradient
subplot(3,2,6)
plot(1:maxit,cnvg(2:end,2),'linewidth',2)
grid
xlabel('Iterations')
ylabel('|\nabla J(C)|')
title('Gradient')

savefig('f-ricardofig.fig')

% Saving solution
save ../../../../../../Documents/MATLAB/inverse-approximation/f-ricardo.mat C totaltime cnvg -v7.3