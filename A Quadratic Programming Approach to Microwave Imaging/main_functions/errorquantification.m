function [zeta_epsilon,zeta_sigma] = errorquantification(originalmap,...
    originalmodel,recoveredmap,recoveredmodel)

%ERRORQUANTIFICATION Compute the error between two maps.
%   [zeta_epsilon,zeta_sigma] = errorquantification(epsr_ori,sig_ori,
%   model_ori,epsr_rec,sig_rec,model_rec) computes the error between an
%   orignal map (epsr_ori, sig_ori, model_ori) and a recovered map
%   (epsr_rec, sig_rec, model_rec). The error on the relative permittivity
%   is defined as the average percentage deviation from the recovered map
%   to the original one. The error on the conductivity is defined in the
%   same way, however, it is normalized by the maximum conductivity in the
%   original map.
%
%   If the maps have different resolutions, the recovered one is
%   interpolate to same resolution as the original one.
%
%   Implemented by Andre Costa Batista (Universidade Federal de Minas
%   Gerais, Brazil). Contact: andre-costa@ufmg.br

    Io = originalmodel.I;
    Jo = originalmodel.J;
    Ir = recoveredmodel.I;
    Jr = recoveredmodel.J;
    
    epsr_ori = originalmap.epsilon_r;
    sig_ori = originalmap.sigma;
    epsr_rec = recoveredmap.epsilon_r;
    sig_rec = recoveredmap.sigma;
    
    mi = originalmodel.mi;

    if Io ~= Ir || Jo ~= Jr
        [epsr_rec,sig_rec,~] = increaseresolution (recoveredmap,...
            recoveredmodel,Io,Jo);                                             
    end

    if ~isempty(mi)
        epsr_rec(mi) = [];
        epsr_ori(mi) = [];
        sig_rec(mi) = [];
        sig_ori(mi) = [];
    end
    
    N = numel(epsr_ori);
    
    zeta_epsilon = 1/N*sqrt(sum(abs((epsr_rec(:)-epsr_ori(:))...
        ./epsr_rec(:)).^2))*100;
    zeta_sigma = 1/N*sqrt(sum(abs((sig_rec(:)-sig_ori(:))...
        ./max(sig_ori(:))).^2))*100;
end