clc, clear, close all

% Adding directories
addpath('../../../main_functions/')

% Checking parallel pool
if isempty(gcp('nocreate'))
    parpool;
end

classnumber = '2';
samplenumber = '3';
cplanenumber = '5';
fprintf(['Model: Class ',classnumber,', Sample ',samplenumber,...
    ', Coronal plane ',cplanenumber,'.\n'])

% Model properties and parameters
fprintf             ('Setting model properties and parameters... ');
expname             = ['class',classnumber,...
                    'sample',samplenumber,...
                    'cplane',cplanenumber]; % Experiment name
filepath            = './maps';             % Path for breast files
breast_freq         = 1e9;                  % Sampling frequency for data collection
nslices             = 5;                    % Number of slices of the breast
Nx_high             = 150;                  % Number of pixels (x-axis), high resolution
Ny_high             = 150;                  % Number of pixels (y-axis), high resolution
Nx_low              = 100;                  % Number of pixels (x-axis), low resolution
Ny_low              = 100;                  % Number of pixels (y-axis), low resolution
epsilon_rb          = 10;                   % Background relative permittivity
sigma_b             = 0;                    % Background conductivity [S/m]
epsilon_rmax        = 58;                   % Maximum relative permittivity
epsilon_rmin        = 2.5;                  % Minimum relative permittivity
sigma_max           = 1.2;                  % Maximum conductivity [S/m]
sigma_min           = 0;                    % Minimum conductivity [S/m]
frequency           = linspace(3e8,...
                        1.15e9,15);         % Frequency of measurements [Hz]
frequency_waveform  = 1e9;                  % Central frequency of waveform [Hz]
time_window         = 3e-08;             % Time window [sec]
dtheta              = 20;                   % Angular shift among sources [deg]
Rs                  = 10e-2;                % Source array radius [m]
ls_x                = 1e-2;                 % Source size (x-axis) [m]
ls_y                = 1e-2;                 % Source size (y-axis) [m]
current             = 4.0;                  % Impressed current [A]
sfi                 = {1:3,4:6,7:9,...      % Set of frequency indexes
                        10:12,13:15};          
misf                = 5;                    % Number of iterations for sets
alpha               = 5e1;                  % Variational regularization
initialization      = 4;                    % Contrast initialization
coeff_epsr          = 0.0188;
coeff_sig           = -0.04711;
fprintf             ('ok!\n')

% Get breast model
fprintf('Get breast model... ')
load([filepath,'/class_',classnumber,'_',samplenumber,'_',cplanenumber,...
    '.mat'])
breast.epsilon_r = epsilon_r;
breast.sigma = sigma;
fprintf('ok!\n')

% Set model struct
fprintf('Getting model struct... ')
finemodel = setmodel(Nx_high,Ny_high,dx,dy,epsilon_rb,sigma_b,frequency,...
    frequency_waveform,Rs,dtheta,time_window,ls_x,ls_y,current,mi);
fprintf('ok!\n')

% Set parameters struct
fprintf('Getting parameters struct... ')
parameters = setparameters([epsilon_rmin,epsilon_rmax],[sigma_min,sigma_min],...
    sfi,misf,'alpha',alpha,'initialization',initialization,'coeff_epsr',...
    coeff_epsr,'coeff_sig',coeff_sig);
fprintf('ok!\n')

% Get scattered field data
fprintf('Geting scattered field data... ')
[es,finemodel] = getscatteredfield(breast,finemodel);
fprintf('ok!\n')

% View testbench map
fprintf('Plotting testbench... ')
viewresults(breast,finemodel,[expname,'_testbench']);
fprintf('ok!\n')

% Set coarse model
fprintf('Setting a coarse model... ')
coarsemodel = finemodel;
aux_x = finemodel.x;
aux_y = finemodel.y;
aux_xp = linspace(aux_x(1),aux_x(end),Nx_low);
aux_yp = linspace(aux_y(1),aux_y(end),Ny_low);
coarsemodel.dx = aux_xp(2)-aux_xp(1);
coarsemodel.dy = aux_yp(2)-aux_yp(1);
[aux_yp,aux_xp] = meshgrid(aux_yp,aux_xp);
aux_epsr = interp2(aux_y,aux_x,breast.epsilon_r,aux_yp,aux_xp);
aux_sig = interp2(aux_y,aux_x,breast.sigma,aux_yp,aux_xp);
coarsemodel.ls_x = .01;
coarsemodel.ls_y = .01;
coarsemodel.I = Nx_low;
coarsemodel.J = Ny_low;
coarsemodel.magnitude = current/(coarsemodel.ls_x+coarsemodel.dx)/...
    (coarsemodel.ls_y+coarsemodel.dy);
coarsemodel.mi = find(aux_epsr(:)==coarsemodel.epsrb ...
    & aux_sig(:)==coarsemodel.sigb);
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
viewresults(breast,finemodel,solution,coarsemodel,expname);
fprintf('ok!\n')

% Error quantification
fprintf('Computing error... ')
[zeta_e,zeta_s] = errorquantification(breast,finemodel,solution,...
    coarsemodel);
[zeta_r,~] = computeresidual(solution,solution.et,es,gs,coarsemodel);
fprintf('ok!\n')

% Save data
fprintf('Saving data... ')
savedata([],expname,breast,solution,finemodel,coarsemodel,...
    parameters,es,ei,gs,zeta_r,zeta_e,zeta_s);
fprintf('ok!\n')