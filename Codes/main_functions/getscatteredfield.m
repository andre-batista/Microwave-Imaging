function [es,model] = getscatteredfield(contrastmap,model)
%GETSCATTEREDFIELD Simulates the scattered field for a given image.
%   [es,model] = getscatteredfield(contrastmap,model) computes the
%   scattered field for a given contrast map represented by a matrix of
%   relative permittivity and conductivity [S/m] and the model struct. The
%   contrast matrices must have the same dimension as indicated in the
%   model struct.
%
%   The output variable es stores the scattered field with the following
%   dimensions: M*L*F. Where M is the number of measurements per source (S)
%   and per frequency (F).
%
%   REFERENCES
%
%   Chew, Weng Cho. Waves and fields in inhomogeneous media. IEEE press,
%   1995.
%
%   Implemented by Andre Costa Batista (Universidade Federal de Minas
%   Gerais, Brazil). Contact: andre-costa@ufmg.br
    
    epsilon_r = contrastmap.epsilon_r;
    sigma = contrastmap.sigma;

    [I,J] = size(epsilon_r);
    if I ~= model.I && J ~= model.J
        disp(['ERROR GETSCATTERED FIELD: the size of the image is ',...
            'not compatible with the model!'])
        return
    end
    
    % Electromagnetic constant
    eps0 = 8.85418782e-12;
    
    % Dimensions
    L = round(360/model.dtheta);
    M = L;
    F = length(model.f);
    N = I*J;
    
    % Background constants
    epsilon_rb = model.epsrb;
    sigma_b = model.sigb;
    omega = 2*pi*model.f;

    % Contrast function
    tau = zeros(N,F);
    for f = 1:length(omega)
        tau(:,f) = (epsilon_r(:)-1j*sigma(:)/omega(f)/eps0)...
            ./(epsilon_rb-1j*sigma_b/omega(f)/eps0)-1;
    end
    
    % Structs to FDTD
    [data,sourcepar] = getfdtdstructs(epsilon_r,sigma,model);
    
    % Auxiliar variables
    et = zeros(I,J,F,L);
    x = cell(L,1);
    y = cell(L,1);
    xp = cell(L,1);
    yp = cell(L,1);
    par = cell(L,1);
    par(:) = {sourcepar};
    for l = 1:L
        par{l}.s_indx = l;
    end

    % Compute the total electric field at the measurement points
    parfor l = 1:L
        [et(:,:,:,l),xp{l},yp{l},~,x{l},y{l}] = fdtd2d(data,par{l});
    end
    et = permute(reshape(permute(et,[3,4,1,2]),F,L,N),[3,2,1]);

    % Save coordinates of image and measurement points
    model.r = [x{1}',y{1}'];
    model.x = xp{1};
    model.y = yp{1};
    
    % Compute the Green function
    gs = getgreenfunction_s(model);

    % Computes the Volume Integral Electric Field Equation (CHEW, 1995)
    es = zeros(M,L,F);
    for m = 1:M
        for l = 1:L
            for f = 1:F
                es(m,l,f) = sum(gs(m,:,f).'.*et(:,l,f).*tau(:,f));
            end
        end
    end
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