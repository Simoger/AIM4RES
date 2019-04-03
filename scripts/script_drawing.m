%% Script that calls the suitable drawing script for the inverted results graphic representation

close all
clear

disp('Drawing inverted results:')
disp('1. Load anisotropic inversion of anisotropic model')
disp('2. Load anisotropic inversion of isotropic model')
disp('3. Load isotropic inversion of isotropic model')
prompt = '';

val = input(prompt);

if val == 1
    load inv_aniso_model_aniso
    script_drawing_aniso
elseif val == 2
    load inv_aniso_model_iso
    script_drawing_aniso
elseif val == 3
    load inv_iso_model_iso
    script_drawing_iso
else
    error('Please enter 1, 2 or 3')
end