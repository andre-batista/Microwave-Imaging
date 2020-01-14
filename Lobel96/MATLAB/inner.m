%INNER Inner product
%   This function implements the inner product between two function defined
%   by Lobel et al., 1996.
%
%   Inputs:
%       - v1: fisrt function as a N-dimensional vector
%       - v2: second function as a N-dimensional vector
%       - d: element size of integration
%
%   Outputs:
%       - fx: inner product evaluation
%
%   Example:
%
%       [fx] = inner(rand(100,1),rand(100,1),rand);
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

function [fx] = inner(v1,v2,d)
    fx = sum(d.*v1.*conj(v2));