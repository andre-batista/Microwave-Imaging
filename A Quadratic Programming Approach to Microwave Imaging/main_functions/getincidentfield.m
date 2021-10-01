function [ei,model] = getincidentfield(model)
%GETINCIDENTFIELD Compute the incident electric field.
%   [ei,model] = getincidentfield(model) computes the incident field within
%   the image domain for a given model. The incident field is stored as a
%   matrix with the following dimensions: N*L*F, where N is the number of
%   pixels (linearly indexed), L is the number of sources and F is the
%   number of frequencies. In addition, the model struct is updated with
%   the coordinates of the pixels and measurement points.
%
%   Implemented by Andre Costa Batista (Universidade Federal de Minas
%   Gerais, Brazil). Contact: andre-costa@ufmg.br

    % Dimensions
    I = model.I;
    J = model.J;
    N = I*J;
    F = length(model.f);
    L = round(360/model.dtheta);

    % Allocation
    ei = zeros(I,J,F,L);
    x = cell(L,1);
    y = cell(L,1);
    xp = cell(L,1);
    yp = cell(L,1);
    
    % Background contrast
    epsilon_r = model.epsrb*ones(I,J);
    sigma = model.sigb*ones(I,J);
    [data,sourcepar] = getfdtdstructs(epsilon_r,sigma,model);

    % Auxiliar variables
    auxdata = data;
    par = cell(L,1);
    par(:) = {sourcepar};

    % Run simulations
    for l = 1:L
        par{l}.s_indx = l;
    end
    parfor l = 1:L
        [ei(:,:,:,l),xp{l},yp{l},~,x{l},y{l}] = fdtd2d(auxdata,par{l}); 
    end

    % Updating variables
    ei = permute(reshape(permute(ei,[3,4,1,2]),F,L,N),[3,2,1]);
    model.r = [x{1}',y{1}'];
    model.x = xp{1};
    model.y = yp{1};
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