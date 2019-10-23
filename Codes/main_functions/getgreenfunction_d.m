function [gd] = getgreenfunction_d(model,indexes)
%GETGREENFUNCTION_D Computes the Green function for image domain.
%   gd = getgreenfunction_d(model) computes the Green function for points
%   in the image domain. The indexes array represents the linear index of
%   the meausurement points in the image domain. They must be greater than
%   zero and at most equal to the number of pixels. The dimension of the
%   Green function matrix is P*N*F where P is the number of reference
%   points, N is the number of pixels an F is the number of frequencies
%   considered.
%
%   Implemented by Andre Costa Batista (Universidade Federal de Minas
%   Gerais, Brazil). Contact: andre-costa@ufmg.br

    % Mesh parameters
    dx = model.dx;
    dy = model.dy;
    kb = model.kb;
    x = model.x;
    y = model.y;
    F = length(model.f);

    % Green's constants
    deltasn = dx*dy;
    an = sqrt(deltasn/pi);

    % Compute radius
    [y,x] = meshgrid(y,x);
    R = pdist2([x(indexes),y(indexes)],[x(:),y(:)]);

    % Compute Green function
    gd = zeros([size(R),F]);
    for f = 1:F
        aux = zeros(size(R));
        aux(R~=0) = 1j/2*pi*kb(f)*an*besselj(1,kb(f)*an)*besselh(0,2,kb(f)*R(R~=0));
        aux(R==0) = (1j/2)*(pi*kb(f)*an*besselh(1,2,kb(f)*an)-2j);
        gd(:,:,f) = -aux;
    end
end