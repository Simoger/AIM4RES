clear
close all

%% Model selection

% This script computes the synthetic data from model1 presented in the 
% Computer and Geosciences article: 
% model1 is a two layer model with ? = 0. Upper layer is 4m deep, with
% ?H = 100 ?.m and ?V = 400 ?.m. Bottom layer is a semi-infinite space
% with ?H = 10 ?.m and ?V = 40 ?.m.

%% Protocol design

load protocol.mat

[~,Ind] = sort(protocole_surface(:,4)-protocole_surface(:,1));
protocole_surface = protocole_surface(Ind(1:3:end),:);

p1t = [1:2:22 27:2:50]';
p2t = [2:2:23 28:2:50]';
c1t = (1:14)';
c2t = (2:15)';
c1 = kron(c1t,ones(length(p1t),1));
c2 = kron(c2t,ones(length(p1t),1));
p1 = kron(ones(length(c1t),1),p1t);
p2 = kron(ones(length(c1t),1),p2t);
protocole_surface_puits = [c1 c2 p1 p2];

XYZ.surface_electrode(:,1) = 0:50;
XYZ.surface_electrode(:,2) = 0;

XYZ.borehole1_electrode(:,2) = 0.5:.5:20;
xb1 = 25;
XYZ.borehole1_electrode(:,1) = xb1;

XYZ_C1_surface = [protocole_surface(:,1) zeros(length(protocole_surface),1)];
XYZ_C1_puits = [xb1*ones(length(protocole_puits),1) protocole_puits(:,1)];
XYZ_C1_surface_puits = [xb1*ones(length(protocole_surface_puits),1) protocole_surface_puits(:,1)];
XYZ.MEAS.C1 = [XYZ_C1_surface ; XYZ_C1_puits; XYZ_C1_surface_puits];

XYZ_C2_surface = [protocole_surface(:,2) zeros(length(protocole_surface),1)];
XYZ_C2_puits = [xb1*ones(length(protocole_puits),1) protocole_puits(:,2)];
XYZ_C2_surface_puits = [xb1*ones(length(protocole_surface_puits),1) protocole_surface_puits(:,2)];
XYZ.MEAS.C2 = [XYZ_C2_surface ; XYZ_C2_puits ; XYZ_C2_surface_puits];

XYZ_P1_surface = [protocole_surface(:,3) zeros(length(protocole_surface),1)];
XYZ_P1_puits = [xb1*ones(length(protocole_puits),1) protocole_puits(:,3)];
XYZ_P1_surface_puits = [protocole_surface_puits(:,3) zeros(length(protocole_surface_puits),1)];
XYZ.MEAS.P1 = [XYZ_P1_surface ; XYZ_P1_puits ; XYZ_P1_surface_puits];

XYZ_P2_surface = [protocole_surface(:,4) zeros(length(protocole_surface),1)];
XYZ_P2_puits = [xb1*ones(length(protocole_puits),1) protocole_puits(:,4)];
XYZ_P2_surface_puits = [protocole_surface_puits(:,4) zeros(length(protocole_surface_puits),1)];
XYZ.MEAS.P2 = [XYZ_P2_surface ; XYZ_P2_puits ; XYZ_P2_surface_puits];

geo_factor_temp = 10*double(XYZ.MEAS.C1(:,2)==0) + double(XYZ.MEAS.C2(:,2)==0); % XYZ.MEAS.C1(:,2)==0 true if C1 at surface
geo_factor_temp(geo_factor_temp == 10) = 4; % 10 : C1 surface
geo_factor_temp(geo_factor_temp == 01) = 3; % 01 : C2 surface
geo_factor_temp(geo_factor_temp == 11) = 2; % 11 : C1 & C2 borehole
geo_factor_temp(geo_factor_temp == 00) = 1; % 00 : C1 & C2 surface

param.flag.geo_factor = geo_factor_temp;

Elec_pos_mat = [XYZ.MEAS.C1 XYZ.MEAS.C2 XYZ.MEAS.P1 XYZ.MEAS.P2];

%% Grid preparation and model building

param.cell_size = [4 4];

hx = 1/param.cell_size(1) * ones(param.cell_size(1)*ceil(max(max(Elec_pos_mat(:,[1 3 5 7])))), 1); 
hz = 1/param.cell_size(2) * ones(param.cell_size(2)*ceil(max(max(Elec_pos_mat(:,[2 4 6 8])))), 1);
    if isempty(hz)
        hz = ones(1,15);
    end
    
param.nb_pad_bloc = 17;
param.fact = 1.5;
param.nb_raff = [1 1];
param.nb_surr = 4*param.cell_size(1);

[param.grille,param.h_x,param.h_z] = grille2d_elect(hx,hz,param.fact,param.nb_pad_bloc,param.nb_raff,param.nb_surr,0);

param.nb_ligne = length(param.h_z)+1;
param.nb_col = length(param.h_x)+1;
k1 = 8; % wavenumber 
nb_dondes = 13;

[param.k,param.wk,param.nbk1]=find_k(nb_dondes,k1);

if ~isfield(XYZ,'surface_electrode')
    XYZ.surface_electrode(:,1) = NaN;
    XYZ.surface_electrode(:,2) = NaN;
end

if ~isfield(XYZ,'borehole1_electrode')
    XYZ.borehole1_electrode(:,1) = NaN;
    XYZ.borehole1_electrode(:,2) = NaN;
end

if ~isfield(XYZ,'borehole2_electrode')
    XYZ.borehole2_electrode(:,1) = NaN;
    XYZ.borehole2_electrode(:,2) = NaN;
end

[param.pos,param.ind,param.L,param.sign]= tri_electr_pos(XYZ.MEAS.C1 , XYZ.MEAS.C2,...
    XYZ.borehole1_electrode , XYZ.borehole2_electrode , XYZ.surface_electrode , XYZ.MEAS.P1 , XYZ.MEAS.P2 , param.grille); 

param.rho.xx = 10*ones(param.nb_ligne,param.nb_col);
param.rho.xx(1:4*param.cell_size(2),:) = 100;
param.rho.xz = zeros(param.nb_ligne,param.nb_col);
param.rho.zz = 40*ones(param.nb_ligne,param.nb_col);
param.rho.zz(1:4*param.cell_size(2),:) = 400;
param.rho.yy = param.rho.xx;


param.flag.inv.p = 2; % no xz component
param.data_nb = length(XYZ.MEAS.C1(:,1));

%% Plots

txtnum1 = num2str(param.rho.xx(1));
    txtnum11 = num2str(param.rho.xx(end));
txtnum2 = num2str(param.rho.zz(1));
    txtnum22 = num2str(param.rho.zz(end));
txtnum3 = num2str(param.rho.xz(1));
    txtnum33 = num2str(param.rho.xz(end));

figure

subplot(2,2,1)
    imagesc(param.rho.xx)
    colormap jet
    colorbar
    title('\rho_x_x')
    xlabel('x')
    ylabel('z')
    set(gca,'FontSize',20)
    caxis([0 400])
    text(100,10,strcat('\rho_x_x=',txtnum1,'\Omega.m'),'fontsize',20,'color','w')
    text(100,50,strcat('\rho_x_x=',txtnum11,'\Omega.m'),'fontsize',20,'color','w')
    axis square
subplot(2,2,2)
    imagesc(param.rho.zz)
    colormap jet
    colorbar
    title('\rho_z_z')
    xlabel('x')
    ylabel('z')
    set(gca,'FontSize',20)
    caxis([0 400])
    text(100,10,strcat('\rho_z_z=',txtnum2,'\Omega.m'),'fontsize',20,'color','w')
    text(100,50,strcat('\rho_z_z=',txtnum22,'\Omega.m'),'fontsize',20,'color','w')
    axis square
subplot(2,2,3.5)
    imagesc(param.rho.xz)
    colormap jet
    colorbar
    title('\rho_x_z')
    xlabel('x')
    ylabel('z')
    set(gca,'FontSize',20)
    text(50,50,strcat('\rho_x_z=',txtnum3,'\Omega.m'),'fontsize',20,'color','w')
    caxis([0 400])
    axis square

%% Forward modeling: apparent resistivity data to inverse

tic
[u,~] = calcul_u_S_anis(param,XYZ,param.rho,0);
toc

K = geometric_factor(XYZ,param.flag.geo_factor);
param.K = K(:);
param.MEAS.Res = u(:);

rho_calc = param.K.*param.MEAS.Res';


%% Inversion parameters
% parameters needed for the inversion. They should be edited by the user in
% order to improve the inversion results according to the considered 
% inverse problem

% Initial anisotropy
param.anis_init = 4;

% Constraint vector 
constr_col = 0:1;
b_depth = 15; % borehole depth
I = kron( [ find(sum(unique(param.grille(:,1)) == xb1-constr_col/(param.nb_raff(1)*param.cell_size(1)),2)) ; ...
            find(sum(unique(param.grille(:,1)) == xb1+constr_col/(param.nb_raff(1)*param.cell_size(1)),2)) ] , ...
            ones(length(find(unique(param.grille(:,2)) <= b_depth)),1) ...
        ); %
J = kron(ones(2*length(constr_col),1), find(unique(param.grille(:,2)) <= b_depth));
param.const_ind = sub2ind([param.nb_col param.nb_ligne],I,J);
const_1 = param.rho.xx(1:10*param.nb_raff(2)+1 , param.nb_surr+param.nb_pad_bloc-1 : param.nb_surr+param.nb_pad_bloc+1);
const_1 = const_1(:);
const_2 = param.rho.xx(1:10*param.nb_raff(2)+1 , param.nb_surr+param.nb_pad_bloc+8*param.nb_raff(2)-1 : param.nb_surr+param.nb_pad_bloc+8*param.nb_raff(2)+1);
const_2 = const_2(:);
param.const_vect = [const_1 ; const_2];
param.const_TrueFalse = false; % true: constraint added to the inversion; false: without constraint

param.angles = 0; % theta = 0

% constraint variables (if not assigned, default values will be used)
param.const.rho_init = []; % Initial values of rho in all cells
param.const.angles_init = param.angles;
param.const.rho_ref = []; % Reference model used
param.const.rho_min = 1e-3; % Min attributable value of rho in in the grid, scalar
param.const.rho_max = []; % Max attributable value of rho in in the grid, scalar
param.const.WGx = []; % weights on x
param.const.WGz = []; % weights on z

% inverse modeling features
param.inv.invparam = 'log resistivity'; % 'log resistivity' or 'resistivity' or 'resistance'
param.inv.fct_reg  = 'flatness'; % ?   % 'flatness' or 'smoothness'
if strcmp(param.inv.fct_reg,'flatness')
    param.flag.reg_fct = 1;
elseif strcmp(param.inv.fct_reg,'smoothness')
    param.flag.reg_fct = 2;
end
param.inv.appl_fct_reg = 'model'; % 'model' or 'model perturbation'

% data weight
param.inv.dataweight.type = 'const. weight'; % 'const. weight' or 'W-estimator' or 'L1'
param.inv.dataweight.b = 0.01; % default = 0.01
param.inv.dataweight.a = 0.05; % 5% d'erreur sur les données, default = 0.03

param.inv.alx = 1; % default = 1
param.inv.alz = 1; % default = 1
param.inv.als = 0.01; % default = 0.01
param.flag.sen = 1; % default = 1

% regularization coefficient options
param.inv.BETA = 10; % default = no default value. Depends highly on the considered inverse problem

% convergence criteria
param.inv.rms_model = 1; % default = 1
param.inv.tol = 1e-3; % default = 1e-3
param.inv.maxit = 20; % nbr d'iteration, default = 20

param.inv.weight = false;
param.inv.weightMinit = 1;
param.inv.weightMaxit = param.inv.maxit;
param.inv.weightFun = 'distance';

% Line search option
param.inv.alp = 1e-4; % Armijo coefficient, default = 1e-4

%% Save synthetic data

data_save_path = './data/synt_data_for_inversion.mat';
save(data_save_path,'param','XYZ')

