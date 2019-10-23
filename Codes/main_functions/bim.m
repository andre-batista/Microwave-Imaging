function [solution] = bim(varargin)
%BIM Born Iterative Method for microwave imaging.
%   solution = bim(ei,es,gs,model,parameters) solves the inverse
%   scattering problem reconstructing an image of dielectric properties
%   based on measurements of the incident field, scattered field, Green
%   function and a corresponding eletromagnetic model. The solution
%   sctructs contains the relative permittivity and conductivity maps.
%
%   solution = bim(...,initialguess) starts the algorithm with an initial
%   guess for contrast map. The variable must be a struct with fields
%   'epsilon_r' and 'sigma' and they must be a matrix with the same
%   dimensions of the model.
%
%   REFERENCES
%
%   Shah, Pratik, and Mahta Moghaddam. "A fast level set method for 
%   multimaterial recovery in microwave imaging." IEEE Transactions on 
%   Antennas and Propagation 66.6 (2018): 3017-3026.
%
%   Lobel, P., et al. "Conjugate gradient method for solving inverse 
%   scattering with experimental data." (1996).
%                                                                         
%   Implemented by Andre Costa Batista (Universidade Federal de Minas
%   Gerais, Brazil). Contact: andre-costa@ufmg.br           

    if nargin < 5 || nargin > 6
        disp('ERROR BIM: invalid number of inputs!')
        return
    end
    
    ei = varargin{1};
    es = varargin{2};
    gs = varargin{3};
    model = varargin{4};
    parameters = varargin{5};
    
    if nargin == 6
        initialguess = varargin{6};
    end
    
    MAX_IT = length(parameters.sfi)*parameters.misf;
   
    if parameters.piip ~= 0
        N = model.I*model.J;
        P = round(parameters.piip*N);
        indexes = randi(N,[P,1]);
        gd = getgreenfunction_d(model,indexes);
    else
        gd = [];
        indexes = [];
    end
    
    % Data size
    I = model.I; % Mesh size (x-axis)
    J = model.J; % Mesh size (y-axis)
    N = I*J; % Number of pixels        
    [~,L,F] = size(es); % Number of measurements, sources, freq.
    
    % Born approximation of first order
    et = ei;
    
    % Map initialization
    switch parameters.initialization
        
        case 1 % Shah & Moghaddam (2018) - normalized, multifrequency
            solution = initialsolution1 (es,ei,gs,model,parameters);
        case 2 % Shah & Moghaddam (2018) - multifrequency
            solution = initialsolution2 (es,ei,gs,model,parameters);
        case 3 % Lobel et al. (1996) - multifrequency
            solution = initialsolution3 (es,ei,gs,model,parameters);
        case 4 % Lower bound everywhere
            solution = initialsolution4 (model,parameters);
        case 5 % Upper bound everywhere
            solution = initialsolution5 (model,parameters);
        case 6 % Input initial guess
            solution.epsilon_r = initialguess.epsilon_r;
            solution.sigma = initialguess.sigma;
    end

    disp('========== The Born Iterative Method with QP ==========')

    auxpar = cell(L,1);
    zeta_r = zeros(MAX_IT,1);
    [data,par] = getfdtdstructs(solution.epsilon_r,solution.sigma,model);
    
    for l = 1:L
        auxpar{l} = par;
        auxpar{l}.s_indx = l;
    end
        
    tic
    i = 1;
    for j = 1:length(parameters.sfi)
        f = parameters.sfi{j};
        
        % Iterations
        for it = 1:parameters.misf
            
            et = zeros(I,J,F,L);
            fprintf(['Running iteration: ',num2str(i)])
            
            % Solving the foward problem
            parfor s_indx = 1:L
                [et(:,:,:,s_indx),~,~] = fdtd2d(data,auxpar{s_indx});
            end
            
            et = permute(reshape(permute(et,[3,4,1,2]),F,L,N),[3,2,1]);

            % Solving the inverse subproblem
            [solution,~] = recover(es,et,gs,model,parameters,...
                'piip',ei,gd,indexes,'frequencyindexes',f);
            
            % Computing residual
            [zeta_r(i),~] = computeresidual(solution,et,es,gs,model);
            fprintf([' - Residual: ',num2str(zeta_r(i),'%.2e'),'\n'])
            [data,~] = getfdtdstructs(solution.epsilon_r,solution.sigma,...
                model);
            i = i + 1;
        end
    end
    toc
    solution.et = et;
    solution.convergence = zeta_r;

end

function [solution] = initialsolution1 (es,et,gs,model,parameters)

    % Size variables
    [N,~,~] = size(et);
    [M,L,F] = size(es);
    
    % Auxiliary variables
    omega   = 2*pi*model.f;          % Angular frequency [rad/s]
    eps0    = 8.85418782e-12;
    minepsr = parameters.epsr_min;
    minsig  = parameters.sig_min;
    maxepsr = parameters.epsr_max;
    maxsig  = parameters.sig_max;
    
    % Coefficient matrix of linear system y=Ax
    A = zeros(2*M*L*F,2*N);
    
    % Setting A
    row = 1;
    n = 1:N;
    for f = 1:F
        for l = 1:L
            for m = 1:M
                
                % Real part of the equations
                A(row,n) = real(reshape(gs(m,n,f),1,[]).*...
                    reshape(et(n,l,f),1,[]));
                A(row,N+n) = imag(reshape(gs(m,n,f),1,[]).*...
                    reshape(et(n,l,f),1,[]))/omega(f)/model.epsrb/eps0;
                
                % Imaginary part of the equations
                A(row+1,n) = imag(reshape(gs(m,n,f),1,[]).*...
                    reshape(et(n,l,f),1,[]));
                A(row+1,N+n) = -real(reshape(gs(m,n,f),1,[]).*...
                    reshape(et(n,l,f),1,[]))/omega(f)/model.epsrb/eps0;
                
                row = row + 2;
            end
        end
    end

    % Setting Right-hand-side
    y = zeros(2*M*L*F,1);
    y(1:2:end) = real(es(:));
    y(2:2:end) = imag(es(:));

    chi = A'*y;
    chi_epsr = chi(1:N);
    chi_sig = chi(N+1:2*N);
    
    solution.epsilon_r = reshape(minepsr + (chi_epsr-min(chi_epsr(:)))...
        /(max(chi_epsr(:))-min(chi_epsr(:))) * (maxepsr-minepsr),...
        [model.I,model.J]);
    
    solution.sigma = reshape(minsig + (chi_sig-min(chi_sig(:)))...
        /(max(chi_sig(:))-min(chi_sig(:))) * (maxsig-minsig),...
        [model.I,model.J]);
    
    solution.epsilon_r(model.mi) = model.epsrb;
    solution.sigma(model.mi) = model.sigb;

end

function [solution] = initialsolution2 (es,et,gs,model,parameters)

    % Size variables
    [N,~,~] = size(et);
    [M,L,F] = size(es);
    
    % Auxiliary variables
    omega = 2*pi*model.f;
    eps0 = 8.85418782e-12;
     
    A = zeros(M*L,N,F);
    y = zeros(M*L,F);
    x = zeros(N,F);
    epsilon_r = zeros(N,F);
    sigma = zeros(N,F);
    for f = 1:F
        row = 1;
        for l = 1:L
            for m = 1:M
                A(row,:,f) = gs(m,:,f).'.*et(:,l,f);
                y(row,f) = es(m,l,f);
                row = row + 1;
            end
        end
        
        x(:,f) = A(:,:,f)'*y(:,f);
        epsilon_r(:,f) = real(x(:,f))+model.epsrb;
        sigma(:,f) = model.sigb-imag(x(:,f))*omega(f)*eps0;
    end
    
    epsilon_r = mean(epsilon_r,2);
    sigma = mean(sigma,2);
    solution.epsilon_r = reshape(epsilon_r,[model.I, model.J]);
    solution.sigma = reshape(sigma,[model.I,model.J]);
    
    solution.epsilon_r(solution.epsilon_r < parameters.epsr_min) ...
        = parameters.epsr_min;
    solution.sigma(solution.sigma < parameters.sig_min) ...
        = parameters.sig_min;
    solution.epsilon_r(solution.epsilon_r > parameters.epsr_max) ...
        = parameters.epsr_max;
    solution.sigma(solution.sigma > parameters.sig_max) ...
        = parameters.sig_max;
    
    solution.epsilon_r(model.mi) = model.epsrb;
    solution.sigma(model.mi) = model.sigb;

end

function [solution] = initialsolution3 (es,et,gs,model,parameters)

    % Size variables
    [N,~,~] = size(et);
    [~,L,F] = size(es);
    
    % Auxiliary variables
    omega = 2*pi*model.f;
    eps0 = 8.85418782e-12;
    
    epsilon_r = zeros(N,F);
    sigma = zeros(N,F);
    
    for fi = 1:F
    
        gs_aux = squeeze(gs(:,:,fi));
        es_aux = squeeze(es(:,:,fi));
        ei_aux = squeeze(et(:,:,fi));
    
        gamma = norm(reshape(gs_aux'*es_aux,[],1))^2 ...
            /norm(reshape(gs_aux*gs_aux'*es_aux,[],1))^2;
        w0 = gamma*gs_aux'*es_aux;
        C = sparse(diag(1/L*sum(w0./ei_aux,2)));
        tau = diag(C);
        epsilon_r(:,fi) = real(tau) + model.epsrb;
        sigma(:,fi) = -omega(fi)*eps0*model.epsrb*imag(tau); 
        
    end
    
    epsilon_r = mean(epsilon_r,2);
    sigma = mean(sigma,2);
    solution.epsilon_r = reshape(epsilon_r,[model.I, model.J]);
    solution.sigma = reshape(sigma,[model.I,model.J]);
    
    solution.epsilon_r(solution.epsilon_r < parameters.epsr_min) ...
        = parameters.epsr_min;
    solution.sigma(solution.sigma < parameters.sig_min) ...
        = parameters.sig_min;
    solution.epsilon_r(solution.epsilon_r > parameters.epsr_max) ...
        = parameters.epsr_max;
    solution.sigma(solution.sigma > parameters.sig_max) ...
        = parameters.sig_max;
    
    solution.epsilon_r(model.mi) = model.epsrb;
    solution.sigma(model.mi) = model.sigb;

end

function [solution] = initialsolution4 (model,parameters)

    solution.epsilon_r = parameters.epsr_min*ones(model.I,model.J);
    solution.sigma = parameters.sig_min*ones(model.I,model.J);
    solution.epsilon_r(model.mi) = model.epsrb;
    solution.sigma(model.mi) = model.sigb;

end

function [solution] = initialsolution5 (model,parameters)

    solution.epsilon_r = parameters.epsr_max*ones(model.I,model.J);
    solution.sigma = parameters.sig_max*ones(model.I,model.J);
    solution.epsilon_r(model.mi) = model.epsrb;
    solution.sigma(model.mi) = model.sigb;

end

function [data,par] = getfdtdstructs(epsilon_r,sigma,model)
%GETFDTDSTRUCTS Return the input structs for FDTD code.
%   [data,par] = getfdtdstructs(epsilon_r,sigma,model) returns the
%   necessary structs for running FDTD.

    % Mesh and dielectric properties
    data.epsr = epsilon_r;
    data.sig = sigma;
    data.dx = model.dx;
    data.dy = model.dy;
    data.epsrb = model.epsrb;
    data.sigb = model.sigb;
    
    % Source parameters
    par.f = model.f;
    par.fwv = model.fwv;
    par.tw = model.tw;
    par.dtheta = model.dtheta;
    par.Rs = model.Rs;
    par.s_indx = 1;
    par.ls_x = model.ls_x;
    par.ls_y = model.ls_y;
    par.magnitude = model.magnitude;

end