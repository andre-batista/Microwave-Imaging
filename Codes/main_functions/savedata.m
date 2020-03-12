function [] = savedata(varargin)
%SAVEDATA Save experiment information
%   savedata(filepath,filename,solution,model,parameters,es,ei,gs,zeta_r,
%   zeta_e,zeta_s) saves these variables in a .mat file for experiments
%   with no testbench.
%
%   savedata(filepath,filename,testbench,solution,finemodel,coarsemodel,
%   parameters,es,ei,gs,zeta_r,zeta_e,zeta_s) saves these variables in a 
%   .mat file for experiments with testbench data.

    if nargin == 11
        filepath    = varargin{1};
        filename 	= varargin{2};
        solution    = varargin{3};
        model       = varargin{4};
        parameters  = varargin{5};
        es          = varargin{6};
        ei          = varargin{7};
        gs          = varargin{8};
        zeta_r      = varargin{9};
        zeta_e      = varargin{10};
        zeta_s      = varargin{11};
        
        save([filepath,filename,'.mat'],'solution','model','parameters',...
            'es','ei','gs','zeta_r','zeta_e','zeta_s')
        
    elseif nargin == 13
        filepath    = varargin{1};
        filename    = varargin{2};
        testbench   = varargin{3};
        solution    = varargin{4};
        finemodel   = varargin{5};
        coarsemodel = varargin{6};
        parameters  = varargin{7};
        es          = varargin{8};
        ei          = varargin{9};
        gs          = varargin{10};
        zeta_r      = varargin{11};
        zeta_e      = varargin{12};
        zeta_s      = varargin{13};
        
        save([filepath,filename,'.mat'],'testbench','solution','finemodel',...
            'coarsemodel','parameters','es','ei','gs','zeta_r','zeta_e',...
            'zeta_s')
    end
end
    
        