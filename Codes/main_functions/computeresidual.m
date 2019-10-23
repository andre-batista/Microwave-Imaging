function [zeta_R,residual] = computeresidual(recoveredmap,et,es,gs,model)
%COMPUTERESIDUAL Compute the residual between data and approximation.
%   [zeta_R,residual] = computeresidual(epsilon_r,sigma,es,et,gs,model)
%   computes the error between the scattered field data and the
%   approximation done by BIM algorithm. The inputs epsilon_r and sigma
%   represents the recovered contrast map. The variable et is the intern 
%   total electric field estimated by BIM algorithm. The variables es and
%   gs are the scattered field data and Green function. The mean percentage
%   deviation of all equations is stored in zeta_R. The residual variable
%   stores each difference between data and approximation.
%
%   REFERENCES
%   
%   Chew, Weng Cho. Waves and fields in inhomogeneous media. IEEE press,
%   1995.
%
%   Pastorino, Matteo. Microwave imaging. Vol. 208. John Wiley & Sons, 2010.
%
%   Implemented by Andre Costa Batista (Universidade Federal de Minas
%   Gerais, Brazil). Contact: andre-costa@ufmg.br

    % Model constants and dimensions
    eps0 = 8.854187817e-12;  
    N = numel(recoveredmap.epsilon_r);
    [M,L,F] = size(es);
    omega = 2*pi*model.f;
    epsilon_rb = model.epsrb;
    sigma_b = model.sigb;
    epsilon_r = recoveredmap.epsilon_r;
    sigma = recoveredmap.sigma;

    % Contrast function
    tau = zeros(N,F);
    for f = 1:length(omega)
        tau(:,f) = (epsilon_r(:)-1j*sigma(:)/omega(f)/eps0)...
            ./(epsilon_rb-1j*sigma_b/omega(f)/eps0)-1;
    end
    
    % Volume Electric Field Integral Equation
    y = zeros(size(es));
    for m = 1:M
        for l = 1:L
            for f = 1:F
                y(m,l,f) = sum(gs(m,:,f).'.*et(:,l,f).*tau(:,f));
            end
        end
    end

    % Error calculation
    N = numel(y);
    residual = es-y;
    zeta_R = 1/N*sqrt(sum(abs((y(:)-es(:))./y(:)).^2))*100;
end