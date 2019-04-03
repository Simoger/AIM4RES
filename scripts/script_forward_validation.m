clear
close all

%% Model selection

% This script computes the two synthetic models presented in the Computer
% and Geosciences article: 
%   - model1 is a two layer model with ? = 0. Upper layer is 4m deep, with
%     ?H = 100 ?.m and ?V = 400 ?.m. Bottom layer is a semi-infinite space
%     with ?H = 10 ?.m and ?V = 40 ?.m.
%   - model2 is a homogeneous semi-infinite space model with theta = 30°. 
%     ?1 = 100 ?.m and ?3 = 400 ?.m.

disp('Model to consider:')
disp('1. 2 layer anisotropic model')
disp('2. homogeneous anisotropic model with theta = 30°')
prompt = '';
model = input(prompt);

%% Protocol design

% Model 1 protocol
if model == 1

    % protocol creation
    a = 1:50;
    C1 = 75 - 3/2 * a;
    P1 = 75 - a/2; 
    P2 = 75 + a/2;
    C2 = 75 + 3/2*a;

    % electrode positions
    XYZ.surface_electrode(:,1) = 0:.5:150;
    XYZ.surface_electrode(:,2) = 0;

    % electrode position structure
    XYZ.MEAS.C1 = [C1' zeros(length(C1),1)];
    XYZ.MEAS.C2 = [C2' zeros(length(C2),1)];
    XYZ.MEAS.P1 = [P1' zeros(length(P1),1)];
    XYZ.MEAS.P2 = [P2' zeros(length(P2),1)];

    % geometric factor calculation
    geo_factor_temp = 10*double(XYZ.MEAS.C1(:,2)==0) + double(XYZ.MEAS.C2(:,2)==0); % XYZ.MEAS.C1(:,2)==0 true if C1 at surface
    geo_factor_temp(geo_factor_temp == 10) = 4; % 10 : C1 surface
    geo_factor_temp(geo_factor_temp == 01) = 3; % 01 : C2 surface
    geo_factor_temp(geo_factor_temp == 11) = 2; % 11 : C1 & C2 borehole
    geo_factor_temp(geo_factor_temp == 00) = 1; % 00 : C1 & C2 surface
    param.flag.geo_factor = geo_factor_temp;

    Elec_pos_mat = [XYZ.MEAS.C1 XYZ.MEAS.C2 XYZ.MEAS.P1 XYZ.MEAS.P2];

    % Get quadrilaterals area
    quadrils_points = Elec_pos_mat; % All my quadrilateres
    quadrils_points = mat2cell(quadrils_points,ones(length(quadrils_points),1));
    normRatio = 60/8; % value chosen such as the maximum inline quadrilaterals value is equivalent to the max quadrilaterals area
    XYZ.MEAS.areas = cellfun(@(x) getAperture(x,normRatio), quadrils_points); % areas of all the quadrilaterals

% Model 2 protocol
elseif model == 2

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
    
end

%% Grid preparation and model building

hx = ones(ceil(max(max(Elec_pos_mat(:,[1 3 5 7])))), 1); 
hz = ones(ceil(max(max(Elec_pos_mat(:,[2 4 6 8])))), 1);
    if isempty(hz)
        hz = ones(1,15);
    end
    
param.nb_pad_bloc = 17;
param.nb_raff = [4 4];
param.nb_surr = 4;
param.fact = 1.5;

[param.grille,param.h_x,param.h_z] = grille2d_elect(hx,hz,param.fact,param.nb_pad_bloc,param.nb_raff,param.nb_surr,0);

param.nb_ligne = length(param.h_z)+1; % model row number
param.nb_col = length(param.h_x)+1; % model column number
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

% Model 1
if model == 1
    param.rho.xx = 10*ones(param.nb_ligne,param.nb_col);
    param.rho.xx(1:4*param.nb_raff(2),:) = 100;
    param.rho.xz = zeros(param.nb_ligne,param.nb_col);
    param.rho.zz = 40*ones(param.nb_ligne,param.nb_col);
    param.rho.zz(1:4*param.nb_raff(2),:) = 400;
    param.rho.yy = param.rho.xx;
    
    param.flag.inv.p = 2;

% Model 2   
elseif model == 2
    
    angle_deg = 30;
    
    rhoxx = 100*cosd(angle_deg)^2+400*sind(angle_deg)^2;
    rhozz = 100*sind(angle_deg)^2+400*cosd(angle_deg)^2;
    rhoxz = (400-100)*cosd(angle_deg)*sind(angle_deg);

    param.rho.xx = rhoxx*ones(param.nb_ligne,param.nb_col);
    param.rho.zz = rhozz*ones(param.nb_ligne,param.nb_col);
    param.rho.xz = rhoxz*ones(param.nb_ligne,param.nb_col);
    param.rho.yy = 100*ones(param.nb_ligne,param.nb_col);
    
    param.flag.inv.p = 3;
end

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
    if strcmp(model,'model1')
        text(100,10,strcat('\rho_x_x=',txtnum1,'\Omega.m'),'fontsize',20,'color','w')
        text(100,50,strcat('\rho_x_x=',txtnum11,'\Omega.m'),'fontsize',20,'color','w')
    elseif strcmp(model,'model2')
        text(50,50,strcat('\rho_x_x=',txtnum1,'\Omega.m'),'fontsize',20,'color','w')
    end
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
    if strcmp(model,'model1')
        text(100,10,strcat('\rho_z_z=',txtnum2,'\Omega.m'),'fontsize',20,'color','w')
        text(100,50,strcat('\rho_z_z=',txtnum22,'\Omega.m'),'fontsize',20,'color','w')
    elseif strcmp(model,'model2')
        text(50,50,strcat('\rho_z_z=',txtnum2,'\Omega.m'),'fontsize',20,'color','w')
    end
    axis square
subplot(2,2,3.5)
    imagesc(param.rho.xz)
    colormap jet
    colorbar
    title('\rho_x_z')
    xlabel('x')
    ylabel('z')
    set(gca,'FontSize',20)
    if strcmp(model,'model1')
        text(100,50,strcat('\rho_x_z=',txtnum33,'\Omega.m'),'fontsize',20,'color','w')
    elseif strcmp(model,'model2')
        text(50,50,strcat('\rho_x_z=',txtnum3,'\Omega.m'),'fontsize',20,'color','w')
    end
    caxis([0 400])
    axis square

%% Forward modeling

tic
[u,~] = calcul_u_S_anis(param,XYZ,param.rho,0);
toc

K = geometric_factor(XYZ,param.flag.geo_factor);
param.K = K(:);
rho_calc = K.*u';
u_num = u;

if model == 1
    script_analytical_anis_2layers
elseif model == 2
    script_analytical_anis_semiInfiniteSpace
end

