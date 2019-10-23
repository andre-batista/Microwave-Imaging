function [gs] = getgreenfunction_s(model)
%GETGREENFUNCTION Computes the Green Function.
%   gs = getgreenfunction_s(model) computes the 2D Green Function for the
%   measurement domain S. The input model must contain the coordinates of
%   the measurement points as well as each pixel location. The 3D matrix
%   output has the following dimension: M*N*F, where M is the number of
%   measurements, N is the number of pixels (linearly indexed) and F is the
%   number of frequencies considered.
%
%   Implemented by Andre Costa Batista (Universidade Federal de Minas
%   Gerais, Brazil). Contact: andre-costa@ufmg.br

    % Model constants
    r = model.r;
    x = model.x;
    y = model.y;
    dx = model.dx;
    dy = model.dy;
    f = model.f;
    kb = model.kb;
    deltasn = dx*dy;
    an = sqrt(deltasn/pi);

    % Dimensions
    F = length(f);
    
    % Compute radius
    [y,x] = meshgrid(y,x);
    R = pdist2(r,[x(:),y(:)]);

    % Compute the Green Function for each frequency
    gs = zeros([size(R),F]);
    for f = 1:F
        aux       = zeros(size(R));
        aux(R~=0) = 1j/2*pi*kb(f)*an*besselj(1,kb(f)*an)*besselh(0,2,kb(f)*R(R~=0));
        aux(R==0) = (1j/2)*(pi*kb(f)*an*besselh(1,2,kb(f)*an)-2j);
        gs(:,:,f) = -aux;
    end
end