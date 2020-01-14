%GSMETHOD Golden Section Method
%   This function implements the Golden Section Method for unidimensional
%   optimization and it was designed for the Conjugated Gradient Method
%   applied to nonlinear inverse TMz electromagnetic scattering problem
%   proposed by Lobel et al., 1996.
%
%   Inputs:
%       - rho: residual M-by-L matrix
%       - v: auxiliar M-by-L matrix
%
%   Outputs:
%       - alpha: step size
%
%   Examples:
%
%   [alpha] = gsmethod(rand(100,20),rand(100,20));
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

function [alpha] = gsmethod (rho,v)

    gmean       = 0.618;
    delta       = 1.0e-3;
    [a, b]      = getinterval (rho,v);    
    xa          = b - gmean*(b - a);
    xb          = a + gmean*(b - a);
    fxa         = norm(reshape(rho-xa*v,[],1))^2;
    fxb         = norm(reshape(rho-xb*v,[],1))^2;
    
    while (b - a) > delta
        if fxa > fxb
            a = xa;
            xa = xb;
            fxa = fxb;
            xb = a + gmean*(b - a);
            fxb = norm(reshape(rho-xb*v,[],1))^2;
        else
            b = xb;
            xb = xa;
            fxb = fxa;
            xa = b - gmean*(b - a);
            fxa = norm(reshape(rho-xa*v,[],1))^2;
        end
    end

    alpha = (b + a)/2;
end

function [a,b] = getinterval (rho,v)
    
    step0           = 1.0e-03;
    a               = 0;
    Fa              = norm(rho(:))^2;
    b               = step0;
    Fb              = norm(reshape(rho-step0*v,[],1))^2;
    stepsize        = step0;
    acceleration    = 2;
    
    while (Fa > Fb)
        stepsize = acceleration*stepsize;
        b = b + stepsize;
        Fb = norm(reshape(rho-b*v,[],1))^2;
    end    
end