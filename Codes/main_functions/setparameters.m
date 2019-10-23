function [par] = setparameters(varargin)
%SETPARAMETERS Set parameters for BIM algorithm
%   par = setparameters(epsr_bounds,sig_bounds,sfi,misf,...) builds a
%   struct which will be used as input for BIM algorithm. The variable
%   epsr_bounds may be a scalar or an array. If scalar, it will be assumed
%   that epsr_min=0 and epsr_max wil be this input. If array, the first
%   element is epsr_min and the second is epsr_max.
%
%   The input variable sfi must be a cell in which each element is an array
%   with the indexes of the available measurement frequencies. The
%   algorithm will address each set of frequencies at time, beginning with
%   the first array in the cell. For example, if sfi = {[1,3],[2,4]}, then
%   the algorithm will solve the equations beginning with the first and
%   third frequencies and it will address the second and the fourth latter.
%   
%   The input variable misf must be an integer and it means the number of
%   iterations for each set of frequencies. For example, if mist=5, then
%   the BIM algorithm will run 5 iterations for each set of frequencies.
%   Then the total number of iterations equals to misf times the number of
%   set of frequencies considered.
%   
%   The optional inputs must come in pairs. You must indicate the name in a
%   string with the desired value. The names and their meaning are:
%
%       - 'alpha': regularization coefficient for Tikhonov formulation.
%       Default is 0.
%       - 'beta': regularization coefficient for variational formulation.
%       Default is 0.
%       - 'piip': percentage of intern integral points. This means the
%       percentage of the points of the image domain that will be added in
%       the equations for integral. This is a kind of regularization and it
%       cost a lof of memory. Default is 0.
%       - 'initialization': an index for which method for defining the
%       initial guess of the contrast map. The options are:
%           1 - Born approximation defined according to Shah & Moghaddam
%               (2018) with normalization.
%           2 - Born approximation defined according to Shah & Moghaddam
%               (2018) with truncation.
%           3 - Back-propagation defined according to Lobel et al. (1996)
%               with truncation.
%           4 - Minimum contrast everywhere.
%           5 - Maximum contrast everywhere.
%       - 'coeff_epsr', 'coeff_sig': the conductivity may be represent as a
%       linear function in which the input is te relative permittivity,
%       such as sigma = A*epsilon_r + B. In this case, it can reduce the
%       complexity of inverse subproblem. Therefore, the parameter
%       'coeff_eprs' means the coefficient A in the previous relation, and
%       'coeff_sig', the coefficient B.
%   
%   REFERENCES
%
%   Shah, Pratik, and Mahta Moghaddam. "A fast level set method for 
%   multimaterial recovery in microwave imaging." IEEE Transactions on 
%   Antennas and Propagation 66.6 (2018): 3017-3026.
%
%   Lobel, P., et al. "Conjugate gradient method for solving inverse 
%   scattering with experimental data." (1996).
%                                                                         
%   Implemented by Andre Costa Batista (Universidade Federal de Minas
%   Gerais, Brazil). Contact: andre-costa@ufmg.br
    
    % Check the existence of the required inputs
    if nargin < 4
        disp('ERROR SETPARAMETERS: less inputs than required!')
        return 
    end
    
    % Set the bounds of the relative permittivity variables
    if length(varargin{1}) == 1
        par.epsr_min = 1.;
        par.epsr_max = varargin{1};
    else
        par.epsr_min = varargin{1}(1);
        par.epsr_max = varargin{1}(2);
    end
    
    % Set the bounds of the conductivity variables
    if isempty(varargin{2})
        par.sig_min = 0.;
        par.sig_max = 0;
    elseif length(varargin{2}) == 1
        par.sig_min = 0.;
        par.sig_max = varargin{2};
    else
        par.sig_min = varargin{2}(1);
        par.sig_max = varargin{2}(2);
    end
    
    % Set parameters corresponding to the set of frequencies
    par.sfi = varargin{3};
    par.misf = varargin{4};
    
    % Default optional parameters
    par.alpha = 0.;
    par.beta = 0.;
    par.piip = 0;
    par.initialization = 1;
    par.coeff_epsr = 0.;
    par.coeff_sig = 0.;
    
    i = 5;
    while i <= nargin
        
        switch varargin{i}
            case 'alpha'
                par.alpha = varargin{i+1};
            case 'beta'
                par.beta = varargin{i+1};
            case 'piip'
                par.piip = varargin{i+1};
            case 'initialization'
                par.initialization = varargin{i+1};
            case 'coeff_epsr'
                par.coeff_epsr = varargin{i+1};
            case 'coeff_sig'
                par.coeff_sig = varargin{i+1};
            otherwise
                disp(['WARNING SETPARAMETERS: invalid input ',varargin{i}])
        end
        
        i = i + 2;
    end

end

