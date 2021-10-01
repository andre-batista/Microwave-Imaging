expname = 'class_4_1_302_rec';

solution.sigma(:,:) = 0;
testbench.sigma(:,:) = 0;

% Plot solution
fprintf('Plotting solution... ')
% viewresults2(testbench,finemodel,solution,coarsemodel,expname);
% viewresults2(testbench,finemodel,expname);
viewresults2(solution,coarsemodel,expname);
fprintf('ok!\n')