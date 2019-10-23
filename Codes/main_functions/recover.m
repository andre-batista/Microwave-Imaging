function [solution,fx] = recover(varargin)
%RECOVER Solves the inverse electromagnetic scattering subproblem.
%   [solution,fx] = recover(es,et,gs,model,parameters) solves the inverse
%   electromagnetic scattering subproblem with a given scattered field data
%   (es), a total eletric field map (et), the Green function map and
%   information about the electromagnetic model and parameters. The
%   functions returns a struct with the relative permittivity and
%   conductivity maps and the objective-function value.
%
%   [solution,fx] = recover(...,'piip',ei,gd,gd_indexes) solves the inverse
%   subproblem adding equations to the model computing the integral within
%   the imaging domain. You must give the incident field, Green function
%   computed for each integral and the indexes of the reference points
%   within the image domain.
%
%   [solution,fx] = recover(...,'frequencyindexes',f_indexes) solves the
%   inverse subproblem only for the frequencies specified by their indexes.
%
%   Implemented by Andre Costa Batista (Universidade Federal de Minas
%   Gerais, Brazil). Contact: andre-costa@ufmg.br
    
    % Checking the number of inputs
    if nargin < 5
        disp('ERROR RECOVER: invalid number of inputs!')
        return
    end
    
    % Reading standard inputs
    es = varargin{1};
    et = varargin{2};
    gs = varargin{3};
    model = varargin{4};
    parameters = varargin{5};
    
    % Allocating additional inputs
    ei = [];
    gd = [];
    gd_ind = [];
    P = 0;
    f_indx = 1:length(model.f);
    
    % Reading additional inputs
    i = 6;
    while i < nargin
        switch varargin{i}
            case 'piip'
                ei = varargin{i+1};
                gd = varargin{i+2};
                gd_ind = varargin{i+3};
                P = length(gd_ind);
                i = i + 4;
            case 'frequencyindexes'
                f_indx = varargin{i+1};
                i = i + 2;
        end
    end

    % Dimensions variables
    [N,~,~] = size(et);
    [M,L,~] = size(es);
    F = length(f_indx);
    Np = N-length(model.mi);
    dx = model.dx;
    dy = model.dy;
    
    % Auxiliary variables
    omega   = 2*pi*model.f;
    eps0    = 8.85418782e-12;
    minepsr = parameters.epsr_min;
    minsig  = parameters.sig_min;
    maxepsr = parameters.epsr_max;
    maxsig  = parameters.sig_max;
    epsb    = model.epsrb*eps0;
    REDUCED_MODEL_FLAG = parameters.coeff_epsr ~= 0 ...
        && parameters.coeff_sig ~= 0;
    
    % Number of variables and coefficient matrix
    if REDUCED_MODEL_FLAG
        nvar = Np+2*(M+P)*L*F;
        gurobimodel.A = zeros(2*(M+P)*L*F,Np);
    else
        nvar = 2*Np+2*(M+P)*L*F;
        gurobimodel.A = zeros(2*(M+P)*L*F,2*Np);
    end
    
    % Number of constraints, right-hand side and indexs
    nconst = 2*(M+P)*L*F;
    rhs = zeros(nconst,1);
    n = 1:N;
    n(model.mi) = [];
    
    % Build coefficient matrix for reduced model
    if REDUCED_MODEL_FLAG
        c = parameters.coeff_epsr;
        d = parameters.coeff_sig;
        a = 1/model.epsrb;
        row = 1;
        for f = f_indx
            for l = 1:L
                for m = 1:M

                    A = real(reshape(gs(m,n,f),1,[]).*...
                        reshape(et(n,l,f),1,[]));
                    B = imag(reshape(gs(m,n,f),1,[]).*...
                        reshape(et(n,l,f),1,[]));
                    b = c/omega(f)/epsb;
                    G = -1;
                    H = -d/omega(f)/epsb + model.sigb/omega(f)/epsb;
                    
                    % Real part of the equations
                    gurobimodel.A(row,:) = A.*a+B.*b;
                    rhs(row) = rhs(row) - sum(A.*G-B.*H);
                    
                    % Imaginary part of the equations
                    gurobimodel.A(row+1,:) = -A.*b + B.*a;
                    rhs(row+1) = rhs(row+1) - sum(A.*H+B.*G);

                    row = row + 2;
                end
            end
        end
    
        % Interior integration
        if P ~= 0
            for f = f_indx
                for l = 1:L
                    for m = 1:length(gd_ind)

                        A = real(reshape(gd(m,n,f),1,[]).*...
                            reshape(et(n,l,f),1,[]));
                        B = imag(reshape(gd(m,n,f),1,[]).*...
                            reshape(et(n,l,f),1,[]));
                        b = c/omega(f)/epsb;
                        G = -1;
                        H = -d/omega(f)/epsb + model.sigb/omega(f)/epsb;
                        
                        % Real part of the equations
                        gurobimodel.A(row,:) = A.*a+B.*b;
                        rhs(row) = rhs(row) - sum(A.*G-B.*H);
                        
                        % Imaginary part of the equations
                        gurobimodel.A(row+1,:) = -A.*b + B.*a;
                        rhs(row+1) = rhs(row+1) - sum(A.*H+B.*G);
                        
                        n2 = n2 + 1;
                        row = row + 2;
                    end
                end
            end
        end
        
    % Complete model
    else
        row = 1;
        for f = f_indx
            for l = 1:L
                for m = 1:M
                    
                    % Real part of the equations
                    gurobimodel.A(row,1:Np) = real(reshape(gs(m,n,f),1,[])...
                        .*reshape(et(n,l,f),1,[]));
                    gurobimodel.A(row,Np+1:end) = imag(reshape(gs(m,n,f),...
                        1,[]).*reshape(et(n,l,f),1,[]))/omega(f)/...
                        model.epsrb/eps0;
                    
                    % Imaginary part of the equations
                    gurobimodel.A(row+1,1:Np) = imag(reshape(gs(m,n,f),...
                        1,[]).*reshape(et(n,l,f),1,[]));
                    gurobimodel.A(row+1,Np+1:end) = -real(reshape(gs(m,n,f),....
                        1,[]).*reshape(et(n,l,f),1,[]))/omega(f)/...
                        model.epsrb/eps0;
                    
                    row = row + 2;
                end
            end
        end
        
        % Interior integration
        if P ~= 0
            for f = f_indx
                for l = 1:L
                    for m = 1:length(gd_ind)
%                         for n = 1:N
                            
                            % Real part of the equations
                            gurobimodel.A(row,1:Np) = real(reshape(gd(m,n,...
                                f),1,[]).*reshape(et(n,l,f),1,[]));
                            gurobimodel.A(row,Np+1:end) = imag(reshape(...
                                gd(m,n,f),1,[]).*reshape(et(n,l,f),1,[]))...
                                /omega(f)/model.epsrb/eps0;
                            
                            % Imaginary part of the equations
                            gurobimodel.A(row+1,1:Np) = imag(reshape(gd(m,n,f),1,[]).*...
                                reshape(et(n,l,f),1,[]));
                            gurobimodel.A(row+1,Np+1:end) = -real(reshape...
                                (gd(m,n,f),1,[]).*reshape(et(n,l,f),1,[]))...
                                /omega(f)/model.epsrb/eps0;
%                         end
                        row = row + 2;
                    end
                end
            end
        end
    end
    
    % Setting Right-hand-side   
    if P == 0
        rhs(1:2:end) = rhs(1:2:end) + reshape(real(es(:,:,f_indx)),[],1);
        rhs(2:2:end) = rhs(2:2:end) + reshape(imag(es(:,:,f_indx)),[],1);
    
    else
        u = 2*length(f_indx)*L*M;
        rhs(1:2:u) = rhs(1:2:u) + real(reshape(es(:,:,f_indx),[],1));
        rhs(2:2:u) = rhs(2:2:u) + imag(reshape(es(:,:,f_indx),[],1));
        rhs(u+1:2:end) = rhs(u+1:2:end) ...
            + reshape(real(et(gd_ind,:,f_indx)-ei(gd_ind,:,f_indx)),[],1);
        rhs(u+2:2:end) = rhs(u+2:2:end) ...
            + reshape(imag(et(gd_ind,:,f_indx)-ei(gd_ind,:,f_indx)),[],1);
    end

    % Regularizing model by adjusting coefficients
%     zeta = repmat(max(abs(gurobimodel.A),[],2),1,size(gurobimodel.A,2));
%     gurobimodel.A(abs(gurobimodel.A)<=1e-7*zeta) = 0;
%     gurobimodel.A = (gurobimodel.A).*1./zeta;
%     rhs = rhs./zeta(:,1);
    
    if REDUCED_MODEL_FLAG
        zeta = repmat(max(abs(gurobimodel.A),[],2),1,size(gurobimodel.A,2));
        gurobimodel.A(abs(gurobimodel.A)<=1e-7*zeta) = 0;
        gurobimodel.A = (gurobimodel.A).*1./zeta;
        rhs = rhs./zeta(:,1);
    else
        if min(abs(gurobimodel.A(gurobimodel.A(:)~=0))) <= 1e-13
            zeta = 1e-7/min(abs(gurobimodel.A(gurobimodel.A(:)~=0)));
        else
            zeta = 1;
        end
        gurobimodel.A = zeta*gurobimodel.A;
        rhs = zeta*rhs;
    end
                    
    % Model variables
    gurobimodel.Q     = sparse(nvar,nvar);          % Quadratic coefficient matrix
    gurobimodel.obj   = zeros(nvar,1);              % Linear coefficient matrix
    gurobimodel.A     = sparse([gurobimodel.A,...
        sparse(1:nconst,1:nconst,ones(nconst,1))]); % Coeffiecient matrix of y=Ax
    gurobimodel.rhs   = rhs;                        % Right-hand-side of Ax=y
    gurobimodel.sense = repmat('=',nconst,1);       % Type of equation
    gurobimodel.lb    = zeros(nvar,1);              % Lower bound for variables
    gurobimodel.ub    = zeros(size(gurobimodel.lb));% Upper bound for variables
    
    % Setting the basic quadratic function
    if REDUCED_MODEL_FLAG
        gurobimodel.Q(sub2ind(size(gurobimodel.Q),Np+1:nvar,Np+1:nvar)) = 1;
    else
        gurobimodel.Q(sub2ind(size(gurobimodel.Q),2*Np+1:nvar,2*Np+1:nvar)) = 1;
    end
    
    % Setting Tikhonov regularization: alpha*||x||
    if parameters.alpha ~= .0
        if REDUCED_MODEL_FLAG
            gurobimodel.Q(sub2ind([nvar,nvar],1:Np,1:Np)) ...
                = parameters.alpha;
        else
            gurobimodel.Q(sub2ind([nvar,nvar],1:2*Np,1:2*Np)) ...
                = parameters.alpha;
        end
    end
    
    % Setting Gradient-Integral regularization: integral of grad(x)
    if parameters.beta ~= .0
        
        I = model.I;
        J = model.J;
        
        if maxepsr-minepsr > 0 && maxsig-minsig > 0
            gamma = (maxepsr-minepsr) / (maxsig-minsig);
        else
            gamma = 1;
        end

        not_mi = 1:N;
        not_mi(model.mi) = [];
        
        for j = 2:J
            for i = 2:I
                
                beta = parameters.beta*zeta;
                
                % Foward finite difference
                ij = sub2ind([I,J],i,j);
                imj = sub2ind([I,J],i-1,j);
                ijm = sub2ind([I,J],i,j-1);
                
                if find(model.mi == ij,1) | find(model.mi == imj,1) ...
                        | find(model.mi == ijm,1)
                    continue
                end
                
                ij = find(not_mi == ij,1);
                imj = find(not_mi == imj,1);
                ijm = find(not_mi == ijm,1);

                gurobimodel.Q(ij,ij) = gurobimodel.Q(ij,ij)+1/gamma*beta*...
                    ((dx^2+dy^2)/dx/dy);
                gurobimodel.Q(ij,imj) = gurobimodel.Q(ij,imj)-1/gamma*...
                    beta*dy/dx;
                gurobimodel.Q(imj,ij) = gurobimodel.Q(imj,ij)-1/gamma*...
                    beta*dy/dx;
                gurobimodel.Q(ij,ijm) = gurobimodel.Q(ij,ijm)-1/gamma*...
                    beta*dx/dy;
                gurobimodel.Q(ijm,ij) = gurobimodel.Q(ijm,ij)-1/gamma*...
                    beta*dx/dy;
                gurobimodel.Q(imj,imj) = gurobimodel.Q(imj,imj)+1/gamma*...
                    beta*dy/dx;
                gurobimodel.Q(ijm,ijm) = gurobimodel.Q(ijm,ijm)+1/gamma*...
                    beta*dx/dy;
                
                if ~REDUCED_MODEL_FLAG
                    ij = ij + Np;
                    imj = imj + Np;
                    ijm = ijm + Np;
                    
                    gurobimodel.Q(ij,ij) = gurobimodel.Q(ij,ij) ...
                        + gamma*beta*((dx^2+dy^2)/dx/dy);
                    gurobimodel.Q(ij,imj) = gurobimodel.Q(ij,imj) ...
                        - gamma*beta*dy/dx;
                    gurobimodel.Q(imj,ij) = gurobimodel.Q(imj,ij) ...
                        - gamma*beta*dy/dx;
                    gurobimodel.Q(ij,ijm) = gurobimodel.Q(ij,ijm) ...
                        - gamma*beta*dx/dy;
                    gurobimodel.Q(ijm,ij) = gurobimodel.Q(ijm,ij) ...
                        - gamma*beta*dx/dy;
                    gurobimodel.Q(imj,imj) = gurobimodel.Q(imj,imj) ...
                        + gamma*beta*dy/dx;
                    gurobimodel.Q(ijm,ijm) = gurobimodel.Q(ijm,ijm) ...
                        + gamma*beta*dx/dy;
                end
            end
        end
    end
    
    [~,PSD_FLAG] = chol(gurobimodel.Q);
    if PSD_FLAG == 1 && ~REDUCED_MODEL_FLAG
        gurobimodel.Q = gurobimodel.Q*1e-4;
    end
    
    % Setting upper and lower bounds for relative permittivity
    if REDUCED_MODEL_FLAG
        if c*minepsr+d < 0
            minepsr = -d/c;
        end
        gurobimodel.lb(1:Np) = minepsr;
        gurobimodel.ub(1:Np) = maxepsr;
        gurobimodel.lb(Np+1:end) = -inf;
        gurobimodel.ub(Np+1:end) = inf;
        
    else
        gurobimodel.lb(1:Np) = minepsr/model.epsrb-1;
        gurobimodel.ub(1:Np) = maxepsr/model.epsrb-1;
        gurobimodel.lb(Np+1:2*Np) = minsig-model.sigb;
        gurobimodel.ub(Np+1:2*Np) = maxsig-model.sigb;
        gurobimodel.lb(2*Np+1:end) = -inf;
        gurobimodel.ub(2*Np+1:end) = inf;
    end
        
    gurobipar.OutputFlag = 0;               % Don't display data, Gurobi!
    rst = gurobi(gurobimodel,gurobipar);    % Run Gurobi!
    
    % Recovering dielectric information
    updated_cells = 1:N;
    updated_cells(model.mi) = [];
    if REDUCED_MODEL_FLAG 
        epsilon_r = model.epsrb*ones(N,1);
        epsilon_r(updated_cells) = rst.x(1:Np);
        epsilon_r = reshape(epsilon_r,[model.I,model.J]);
        sigma = c*epsilon_r + d;
        sigma(model.mi) = model.sigb;
        sigma(sigma<0) = 0;
        
    else
        epsilon_r = model.epsrb*ones(N,1);
        epsilon_r(updated_cells) = model.epsrb*(rst.x(1:Np)+1);
        epsilon_r = reshape(epsilon_r,[model.I,model.J]);
        sigma = model.sigb*ones(N,1);
        sigma(updated_cells) = model.sigb+rst.x(Np+1:2*Np);
        sigma = reshape(sigma,[model.I,model.J]);
    end
    solution.epsilon_r = epsilon_r;
    solution.sigma = sigma;
        
    % Objective function
    fx = rst.objval;

end