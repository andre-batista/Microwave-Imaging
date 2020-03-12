function [noise] = getnoise(SNR,data)
%GETNOISE Return noise array.
%   noise = getnoise(SNR,data) returns an array with noise to be added to a
%   data array controlled by a signal-to-noise ratio.

    alpha = 10^(-SNR/20);
    re = normrnd(0,1,size(data));
    im = normrnd(0,1,size(data));
    noise = abs(data).*alpha.*(re+1j*im)./abs(re+1j*im);