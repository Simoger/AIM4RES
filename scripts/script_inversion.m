%% Anisotropic inversion

clear
close all

disp('Input data already existing ?')
disp('1. No: launch script_preparation_inversion.m to create synthetic data')
disp('2. Yes: import existing input data')
val = input('');

if val == 1
    script_preparation_inversion
elseif val == 2
    load synt_data_for_inversion
end

tic

gauss_newton_inversion_anis(param, XYZ)