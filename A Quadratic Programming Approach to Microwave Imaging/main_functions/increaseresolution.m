function [epsilon_r,sigma,model] = increaseresolution (map,model,nx,ny)
%INCREASERESOLUTION Increase the resolution of a map.
%   [epsilon_r,sigma,model] = increaseresolution(map,model,nx,ny)
%   increases the resolution of the contrast map (epsilon, sigma),
%   given its model struct, to a new number of elements in each axis
%   (nx, ny).
%
%   Implemented by Andre Costa Batista (Universidade Federal de Minas
%   Gerais, Brazil). Contact: andre-costa@ufmg.br

    epsilon_r = map.epsilon_r;
    sigma = map.sigma;
    x = model.x;
    y = model.y;
    newx = linspace(x(1),x(end),nx);
    newy = linspace(y(1),y(end),ny);
    [yp,xp] = meshgrid(newy,newx);
    epsilon_r = interp2(y,x,epsilon_r,yp,xp);
    sigma = interp2(y,x,sigma,yp,xp);
    model.dx = newx(2)-newx(1);
    model.dy = newy(2)-newy(1);
    model.x = newx;
    model.y = newy;
    model.I = nx;
    model.J = ny;
    
    if ~isempty(model.mi)
        model.mi = find(epsilon_r(:)==model.epsrb & sigma(:)==model.sigb);
    end
end