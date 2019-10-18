clc, clear, close all

% Adding directories
addpath('../../../main_functions')

% Checking parallel pool
if isempty(gcp('nocreate'))
    parpool;
end

% Model properties and parameters
fprintf             ('Setting model properties and parameters... ');
expname             = 'basic_exp01_twoellipses';% Experiment name
geometryname        = 'ellipses';           % Geometry name
Nx_high             = 60;                   % Number of pixels (x-axis), high resolution
Ny_high             = 60;                   % Number of pixels (y-axis), high resolution
Nx_low              = 48;                   % Number of pixels (x-axis), low resolution
Ny_low              = 48;                   % Number of pixels (y-axis), low resolution
dx                  = .6e-2;                % Pixel size (x-axis) [m]
dy                  = .6e-2;                % Pixel size (y-axis) [m]
epsilon_rb          = 1.;                   % Background relative permittivity
sigma_b             = 0;                    % Background conductivity [S/m]
epsilon_robj        = 2.;                   % Object relative permittivity
sigma_obj           = 0;                    % Object conductivity [S/m]
la                  = 7.5e-2/2;             % Semi-minor axis [m]
lb                  = 18e-2/2;              % Semi-major axis [m]
del                 = 6.75e-2;              % Displacement between two ellipses [m]
frequency           = [.8e9,1e9,2e9];       % Frequency of measurements [Hz]
frequency_waveform  = 2e9;                  % Central frequency of waveform [Hz]
time_window         = 5.2170e-08;           % Time window [sec]
dtheta              = 20;                   % Angular shift among sources [deg]
Rs                  = 25e-2;                % Source array radius [m]
ls_x                = 1.8e-2;               % Source size (x-axis) [m]
ls_y                = 1.8e-2;               % Source size (y-axis) [m]
current             = 5.76;                 % Impressed current [A]
sfi                 = {[1,2,3]};            % Set of frequency indexes
misf                = 10;                   % Number of iterations for sets
beta                = 1e-8;                 % Variational regularization
initialization      = 1;                    % Contrast initialization
fprintf             ('ok!\n')

% Set model struct
fprintf('Getting model struct... ')
finemodel = setmodel(Nx_high,Ny_high,dx,dy,epsilon_rb,sigma_b,frequency,...
    frequency_waveform,Rs,dtheta,time_window,ls_x,ls_y,current);
fprintf('ok!\n')

% Set parameters struct
fprintf('Getting parameters struct... ')
parameters = setparameters(epsilon_robj,sigma_obj,sfi,misf,'beta',beta,...
    'initialization',initialization);
fprintf('ok!\n')

% Set original contrast map
fprintf('Geting testbench map... ')
testbench = getmap(geometryname,finemodel,epsilon_robj,sigma_obj,la,lb,del);
fprintf('ok!\n')

% Get scattered field data
fprintf('Geting scattered field data... ')
[es,finemodel] = getscatteredfield(testbench,finemodel);
fprintf('ok!\n')

% View testbench map
fprintf('Plotting testbench... ')
viewresults(testbench,finemodel,[expname,'_testbench']);
fprintf('ok!\n')

% Set coarse model
fprintf('Setting a coarse model... ')
coarsemodel = setmodel(finemodel,Nx_low,Ny_low);
fprintf('ok!\n')

% Get incident field
fprintf('Getting incident field... ')
[ei,coarsemodel] = getincidentfield(coarsemodel);
fprintf('ok!\n')

% Get Green function data
fprintf('Getting Green function data... ')
gs = getgreenfunction_s(coarsemodel);
fprintf('ok!\n')

% Run BIM algorithm
solution = bim(ei,es,gs,coarsemodel,parameters);

% Plot solution
fprintf('Plotting solution... ')
viewresults(testbench,finemodel,solution,coarsemodel,expname);
fprintf('ok!\n')

% Error quantification
fprintf('Computing error... ')
[zeta_e,zeta_s] = errorquantification(testbench,finemodel,solution,...
    coarsemodel);
[zeta_r,~] = computeresidual(solution,solution.et,es,gs,coarsemodel);
fprintf('ok!\n')

% Save data
fprintf('Saving data... ')
savedata([],expname,testbench,solution,finemodel,coarsemodel,...
    parameters,es,ei,gs,zeta_r,zeta_e,zeta_s);
fprintf('ok!\n')