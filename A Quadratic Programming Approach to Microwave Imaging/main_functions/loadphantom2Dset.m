function [epsilon_r,sigma_eff,dx,dy,mi] = loadphantom2Dset(varargin)
%LOADPHANTOM2DSET Load a set of slices of breast phantoms.
%   [epsilon_r,sigma_eff,dx,dy,mi] = LOADPHANTOM2DSET(filepath,frequency,
%   epsilon_rb,sigma_b,nslices) return a set of relative 
%   permittivity maps (epsilon_r) and conductivity ones (sigma_eff) [S/m]
%   based on files provided by the UWCEM Numerical Breast Phantom
%   Repository. The string variable filepath must contain the path to where
%   the three necessary files are. Only one frequency of measurement is
%   allowed. The variables epsilon_rb and sigma_b represent the background
%   properties. The number of images collects is controlled by the variable
%   nslices. The output mi is a set of indexes in which the background is
%   present in the image.
%
%   [epsilon_r,sigma_eff,dx,dy,mi] = LOADPHANTOM2DSET(filepath,frequency,
%   epsilon_rb,sigma_b,nslices,new_res) returns images with a dimension
%   specified by the two-element integer array new_res.
%
%   REFERENCES
%
%   M. J. Burfeindt, T. J. Colgan, R. O. Mays, J. D. Shea, N. Behdad, B. D. 
%   Van Veen, and S. C. Hagness, "MRI-derived 3D-printed breast phantom for
%   microwave breast imaging validation," IEEE Antennas and Wireless 
%   Propagation Letters, vol. 11, pp. 1610-1613, 2012.
%
%   Lazebnik, M., McCartney, L., Popovic, D., Watkins, C.B., Lindstrom, 
%   M.J., Harter, J., Sewall, S., Magliocco, A., Booske, J.H., Okoniewski, 
%   M. and Hagness, S.C., 2007. A large-scale study of the ultrawideband 
%   microwave dielectric properties of normal breast tissue obtained from 
%   reduction surgeries. Physics in Medicine & Biology, 52(10), p.2637.
%
%   Implemented by Andre Costa Batista (Universidade Federal de Minas
%   Gerais, Brazil). Contact: andre-costa@ufmg.br
   
    if nargin < 5 || nargin > 6
        disp('ERROR LOADPHANTOM: invalid number of inputs!')
        return
    end
    
    filepath = varargin{1};
    frequency = varargin{2};
    epsilon_rb = varargin{3};
    sigma_b = varargin{4};
    nslices = varargin{5};
    
    if nargin == 6
        new_res = varargin{6};
    else
        new_res = [];
    end
    
    if frequency < 500e6 || frequency > 20e9
        disp('ERROR LOADPHANTOM: invalid frequency!')
        return
    end

    % Constants
    omega = 2*pi*frequency;
    dx = .5e-3; dy = .5e-3; dz = .5e-3;
    epsilon_0 = 8.854187817e-12;
    
    % Cole-Cole parameters - 500MHz < f < 20GHz
    epsilon_inf = [2.293, 2.908, 3.140, 4.031, 9.941, 7.821, 6.151, 1.];
    epsilon_delta = [.141, 1.2, 1.708, 3.654, 26.6, 41.48, 48.26, 66.31];
    tau = 1e-12 * [16.4, 16.88, 14.65, 14.12, 10.9, 10.66, 10.26, 7.585];
    alpha = [.251, .069, .061, .055, .003, .047, .049, .063];
    sigma_s = [.002, .02, .036, .083, .462, .713, .809, 1.37];

    % Open info file
    bi_file = fopen([filepath,'/breastInfo.txt'],'r');
    breastInfo = fscanf(bi_file,'%s');
    fclose(bi_file);

    % Open media index file
    m_file = fopen([filepath '/mtype.txt'],'r');
    mi = fscanf(m_file,'%f');
    fclose(m_file);
    
    % Open p-value file
    p_file = fopen([filepath '/pval.txt'],'r');
    p_i = fscanf(p_file,'%f');
    fclose(p_file);
    
    % Cole-Cole model
    epsilon_star = epsilon_inf ...
        + epsilon_delta./(1+(1j*omega*tau).^(1-alpha)) ...
        + sigma_s/1j/omega/epsilon_0;
    
    % Bounds for dielectric propertie
    bound_epsilon_r = real(epsilon_star);
    bound_sigma_eff = -imag(epsilon_star)*omega*epsilon_0;
    
    % Allocation
    epsilon_r = epsilon_rb * ones(size(mi));
    sigma_eff = sigma_b * ones(size(mi));
    tissues = [3.3, 3.2, 3.1, 2, 1.3, 1.2, 1.1];
    
    % Compute dielectric properties for each tissue
    for i = 1:length(tissues)
        epsilon_r(mi==tissues(i)) = p_i(mi==tissues(i))*bound_epsilon_r(i+1)...
            +(1-p_i(mi==tissues(i)))*bound_epsilon_r(i);
        sigma_eff(mi==tissues(i)) = p_i(mi==tissues(i))*bound_sigma_eff(i+1)...
            +(1-p_i(mi==tissues(i)))*bound_sigma_eff(i);     
    end
    
    % Muscle
    epsilon_r(mi==-4) = max(bound_epsilon_r);
    sigma_eff(mi==-4) = max(bound_sigma_eff);
    
    % Skin
    epsilon_r(mi==-2) = mean(bound_epsilon_r); % epsilon_rb;
    sigma_eff(mi==-2) = mean(bound_sigma_eff); % sigma_b;
    
    i = strfind(breastInfo,'s1=');
    j = strfind(breastInfo,'s2=');
    k = strfind(breastInfo,'s3=');
    l = strfind(breastInfo,'class');

    % Mesh size
    I = str2double(breastInfo(i+3:j-1));
    J = str2double(breastInfo(j+3:k-1));
    K = str2double(breastInfo(k+3:l-1));
    
    % Dielectric properties matrices
    epsilon_r = reshape(epsilon_r,[I,J,K]);
    sigma_eff = reshape(sigma_eff,[I,J,K]);
    
    % Image size
    Ly = J*dy;
    Lz = K*dz;
    
    % Interpolation for different resolution
    if ~isempty(new_res)
        if mod(new_res(1),1)==0 && mod(new_res(2),1)==0
            newj = new_res(1);
            newk = new_res(2);
            newdy = Ly/newj;
            newdz = Lz/newk;
        else
            newdy = new_res(1);
            newdz = new_res(2);
            newj = floor(Ly/newdy);
            newk = floor(Lz/newdz);
        end
        
        [y,x,z] = meshgrid((1:J)*dy,(1:I)*dx,(1:K)*dz);
        [newy,newx,newz] = meshgrid((1:newj)*newdy,(1:I)*dx,(1:newk)*newdz);
        
        epsilon_r = interp3(y,x,z,epsilon_r,newy,newx,newz);
        sigma_eff = interp3(y,x,z,sigma_eff,newy,newx,newz);
        dx = dx;
        dy = newdy;
        dz = newdz;
        
        epsilon_r(isnan(epsilon_r)) = epsilon_rb;
        sigma_eff(isnan(sigma_eff)) = sigma_b;
        
    end
    
    epsilon_r = epsilon_r(round(linspace(1,I-10,nslices+2)),:,:);
    sigma_eff = sigma_eff(round(linspace(1,I-10,nslices+2)),:,:);
    epsilon_r([1,end],:,:) = [];
    sigma_eff([1,end],:,:) = [];
    
    mi = cell(nslices,1);
    for i = 1:nslices
        mi{i} = find(squeeze(epsilon_r(i,:))==epsilon_rb ...
            & squeeze(sigma_eff(i,:))==sigma_b);
    end
    
    dx = dy;
    dy = dz;
end