function [u,S] = calcul_u_S_anis(param,XYZ,rho,sens_calc)

%%========================================================================%
%                                                                         %
% Calculation of potential u and sensitivity S in anisotropic media       %
%                                                                         %
%%========================================================================%
%                                                                         %
% In:                                                                     %
% -----------                                                             %
%   sens_calc: sensitivity calculation                                    %
%     sens_calc == 0: don't calculate sensitivity matrix                  %
%     sens_calc == 1: calculate sensistivity matrix                       %
%   param: parameters structure                                           %
%   XYZ: electrodes position structure                                    %
%                                                                         %
% Out:                                                                    %
% -----------                                                             %
%   u: calculated potential                                               %
%       I = 1 == > u = resistance (= u/I = u/1)                           %
%   S: sensitivity matrix (if sens_calc == 1)                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%========================================================================%
% Copyright (C) 2019 Simon GERNEZ and Abderrezak BOUCHEDDA                %
%%========================================================================%
%                                                                         %
% Contacts:                                                               %
%                                                                         %
% Simon GERNEZ                                                            %
%     simon.gernez@ete.inrs.ca                                            %
%     Institut National de la Recherche Scientifique                      %
%     Centre Eau-Terre-Environnement                                      %
%     http://www.ete.inrs.ca/                                             %
%                                                                         %
% Abderrezak BOUCHEDDA                                                    %
%     Abderrezak.Bouchedda@ete.inrs.ca                                    %
%     Institut National de la Recherche Scientifique                      %
%     Centre Eau-Terre-Environnement                                      %
%     http://www.ete.inrs.ca/                                             %
%                                                                         %
% This program is free software; you can redistribute it and/or modify    %
% it under the terms of the GNU General Public License as published by    %
% the Free Software Foundation; either version 2 of the License, or       %
% (at your option) any later version.                                     %
%                                                                         %
% This program is distributed in the hope that it will be useful,         %
% but WITHOUT ANY WARRANTY; without even the implied warranty of          %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           %
% GNU General Public License for more details.                            %
%                                                                         %
% You should have received a copy of the GNU General Public License       %
% along with this program; if not, write to the Free Software             %
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA%

%%========================================================================%

%% Initialisation

    S = [];
    eps = 1e-5;
    position = param.pos;
    ind = param.ind;
    I = 1; % current intensity
    nb_ligne = max(size(param.h_z))+1;
    nb_col = max(size(param.h_x))+1;
    nb_t = nb_ligne*nb_col;
    nb_electrode1 = size(position.C,1);
    nb_electrode2 = size(position.P);
        nb_electrode2 = nb_electrode2(1,1);
    nb_electrode = nb_electrode1+nb_electrode2; % number of electrodes
    nb_k = length(param.k); % wave number
    nb_meas = param.data_nb; % measures number
    position1 = [position.C ; position.P];
    u_P1 = zeros(2*nb_meas,1);
    u_P2 = zeros(2*nb_meas,1);

    if sens_calc ==  0 % no sensitivity
        param.flag.sen = 0; % force to no sensitivity
    end

    if param.flag.sen ==  1
        u1 = zeros(nb_col*nb_ligne,nb_electrode*nb_k);
    end

    qq = sparse(nb_col*nb_ligne,nb_electrode);
    ind_source = [];
    u = zeros(nb_meas,nb_k);

%% Matrix coeff calculation

    [A1,Ak] = matrix_coeff_anis(param.h_x,param.h_z,param.grille,rho,param.k);

    % Ak shouldn't be complex
    try 
        isreal(Ak) == false;
    catch
        error('complex Ak: Ak should be real')
    end

    for ii = 1:nb_electrode
        iind = find(abs(position1(ii,1)-param.grille(:,1)) < eps & abs(position1(ii,2)-param.grille(:,2)) < eps);
        qq(iind,ii) = I/2;
        ind_source = [ind_source;iind];
    end

    %% Forward problem resolution

    if nb_t > 10^6 % if huge parameters number to estimate

        for j = 1:nb_k % loop for k

            A = A1 + sparse(1:nb_t,1:nb_t,Ak(:,j)); 

            % Symmetric reverse Cuthill-McKee permutation
            p = symrcm(A);  
            A = A(p,p);
            [~,up] = sort(p);
            %--------
            setup.droptol = 10^-5;
            % LU precondit. 
            [L1,U1] = ilu(A,setup);
            %--------
            l1 = 0;
            l2 = 0;

            for i = 1:nb_electrode1 % loop for current electrode position 

                q = qq(:,i);
                q = q(p);

                % resoudre Au = q
                [uu,~] = bicgstab(A,q,10^-9,20,L1,U1);            

                uu = uu(up); % rearrangement after Cuthill-McKee permutation

                l2 = l2+param.L(i);

                u1(:,i+(j-1)*nb_electrode) = uu;

                u_P1((1+l1:l2)) = uu(ind.P1(1+l1:l2));
                u_P2((1+l1:l2)) = uu(ind.P2(1+l1:l2));

                l1 = l2;

            end

            u_P1C1(ind.C1C2(param.sign == 1))  = u_P1(param.sign == 1);
            u_P1C2(ind.C1C2(param.sign == -1)) = u_P1(param.sign == -1);
            u_P2C1(ind.C1C2(param.sign == 1))  = u_P2(param.sign == 1);
            u_P2C2(ind.C1C2(param.sign == -1)) = u_P2(param.sign == -1);

            u(:,j) = (u_P1C1-u_P1C2)-(u_P2C1-u_P2C2);

        end

        % integration into the spatial domain
        u = integration_U(u,param.wk,param.k,param.nbk1);

    else

        for j = 1:nb_k % loop for k

            A = A1 + sparse(1:nb_t,1:nb_t,Ak(:,j));
            
            if param.flag.inv.p == 3
                [L,U] = lu(A);
            elseif param.flag.inv.p == 2
                % symmetric coefficient matrix if only 2 components, chol
                % more efficient than lu in that case
                R = chol(A);
                RT = R';
            end
            
            l1 = 0;
            l2 = 0;

            q = qq;

            for i = 1:nb_electrode1 % loop for current electrode position 
                if param.flag.inv.p == 3
                    uu = (U\(L\(q(:,i))));
                elseif param.flag.inv.p == 2
                    uu = (R\(RT\(q(:,i))));
                end

                l2 = l2+param.L(i);

                if param.flag.sen ==  1
                    u1(:,i+(j-1)*nb_electrode) = uu;
                end

                u_P1((1+l1:l2)) = uu(ind.P1(1+l1:l2));
                u_P2((1+l1:l2)) = uu(ind.P2(1+l1:l2));

                l1 = l2;
            end

            u_P1C1(ind.C1C2(param.sign ==  1)) = u_P1(param.sign ==  1);
            u_P1C2(ind.C1C2(param.sign == -1)) = u_P1(param.sign == -1);
            u_P2C1(ind.C1C2(param.sign ==  1)) = u_P2(param.sign ==  1);
            u_P2C2(ind.C1C2(param.sign == -1)) = u_P2(param.sign == -1);

            if isnan(XYZ.MEAS.C2(1,1)) && isnan(XYZ.MEAS.P2(1,1)) % pole-pole
                u_P2C1 = 0; 
                u_P2C2 = 0; 
                u_P1C2 = 0;
            elseif isnan(XYZ.MEAS.C2(1,1)) && ~isnan(XYZ.MEAS.P2(1,1)) % pole-dipole
                u_P2C2 = 0; 
                u_P1C2 = 0;
            end

            u_k = (u_P1C1-u_P1C2) - (u_P2C1-u_P2C2);
            u(:,j) = u_k(:);

        end

        % integration into the spatial domain
        u = integration_U(u,param.wk,param.k,param.nbk1);

    end

%% Sensitivity matrix calculation

    % we start with the calculation of the reciprocal potential pos. misses

    if sens_calc == 1 % sensitivity calculation
        
        if param.flag.inv.p == 2
            if nb_electrode2 ~= 0 % calculate reciprocal potential if P1 or P2 are not used as C1 or C2
                for j = 1:nb_k
                    A = A1 + sparse(1:nb_t,1:nb_t,Ak(:,j));              
                    [R,~,p] = chol(A,'vector'); % Cholesky decomposition with ordering
                    [~,up] = sort(p); % ordering indices
                    for i = 1:nb_electrode2    
                        q = qq(:,i+nb_electrode1);
                        % Au = q
                        uu = (R\(R'\(q(p))));
                        uu = uu(up);
                        u1(:,i+nb_electrode1+(j-1)*nb_electrode) =  uu;
                    end
                end
            end
            
        elseif param.flag.inv.p == 3
            if nb_electrode2 ~= 0 % calculate reciprocal potential if P1 or P2 are not used as C1 or C2
                for j = 1:nb_k
                    A = A1 + sparse(1:nb_t,1:nb_t,Ak(:,j));              
                    [L,U,p] = lu(A,'vector');
                    [~,up] = sort(p); % ordering indices
                    for i = 1:nb_electrode2    
                        q = qq(:,i+nb_electrode1);
                        % Au = q
                        uu = (U\(L\(q(p))));
                        uu = uu(up);
                        u1(:,i+nb_electrode1+(j-1)*nb_electrode) =  uu;
                    end
                end
            end
            
        end

        n = size(ind.meas,1); 

        rho.xx = rho.xx';
        rho.xz = rho.xz';
        rho.zz = rho.zz';
        rho.yy = rho.yy';

        rho.xx = rho.xx(:)'; % Mandatory step so that cells are considered 
        rho.xz = rho.xz(:)'; % as an arranged vector as for inversion 
        rho.zz = rho.zz(:)';
        rho.yy = rho.yy(:)';

        S.xx = zeros(nb_meas,nb_ligne*nb_col);
        S.zz = zeros(nb_meas,nb_ligne*nb_col);	
        S.xz = zeros(nb_meas,nb_ligne*nb_col);
        % % Uncomment if lambda sensitivity is needed (if you don't know wether it is needed, assume it is not)
        %  Slambda = zeros(nb_meas,nb_ligne*nb_col);

        rho.rho_1 = rho.rho_1';
        rho.rho_2 = rho.rho_2';
        rho.angle = rho.angle';
        rho.rho_1 = rho.rho_1(:)';
        rho.rho_2 = rho.rho_2(:)';
        rho.angle = rho.angle(:)';

        if strcmp(param.inv.invparam,'log resistivity') 

            % dipole-dipole
            if ~isnan(XYZ.MEAS.C2(1,1)) && ~isnan(XYZ.MEAS.P2(1,1))

                parfor i = 1:n

                    u_1  = (u1(:,ind.meas(i,1):nb_electrode:end)-u1(:,ind.meas(i,2):nb_electrode:end));% u(c1)-u(c2)
                    u_1R = (u1(:,ind.meas(i,3):nb_electrode:end)-u1(:,ind.meas(i,4):nb_electrode:end));% u(p1)-u(p2)

                    % Gradient calculation
                    g = gradient_product_anis(u_1R, u_1, param.h_x, param.h_z, param.k);

                    % xx
                    g.xx = integration_U(g.xx,param.wk,param.k,param.nbk1); % Inverse Fourier transform
                    g.xx = reshape(g.xx,nb_ligne,nb_col);
                    g.xx = reshape(g.xx',(nb_ligne)*(nb_col),1)'; % horizontally arranged

                    % xz
                    g.xz = integration_U(g.xz,param.wk,param.k,param.nbk1); % Inverse Fourier transform
                    g.xz = reshape(g.xz,nb_ligne,nb_col);
                    g.xz = reshape(g.xz',(nb_ligne)*(nb_col),1)'; % horizontally arranged

                    % zz
                    g.zz = integration_U(g.zz,param.wk,param.k,param.nbk1); % Inverse Fourier transform
                    g.zz = reshape(g.zz,nb_ligne,nb_col);
                    g.zz = reshape(g.zz',(nb_ligne)*(nb_col),1)'; % horizontally arranged

                    % yy
                    g.yy = integration_U(g.yy,param.wk,param.k,param.nbk1); % tr. Fourier inv.
                    g.yy = reshape(g.yy,nb_ligne,nb_col);
                    g.yy = reshape(g.yy',(nb_ligne)*(nb_col),1)'; % horizontally arranged

                    % Calculation of the sensitivity (ln(resistivity))
                    Sxx(i,:) = (1./(rho.rho_1*u(i))).*(g.xx.*(cosd(rho.angle)).^2 + ...
                                                       g.zz.*(sind(rho.angle)).^2 - ...
                                                       g.xz.*cosd(rho.angle).*sind(rho.angle) + ...
                                                       g.yy);

                    Szz(i,:) = (1./(rho.rho_2*u(i))).*(g.xx.*(sind(rho.angle)).^2 + ...
                                                       g.zz.*(cosd(rho.angle)).^2 + ...
                                                       g.xz.*cosd(rho.angle).*sind(rho.angle));

                    Sxz(i,:) = -(rho.angle./u(i)).*(g.xx.*(1./rho.rho_2-1./rho.rho_1).*sind(2*rho.angle) + ...
                                                    g.zz.*(1./rho.rho_1-1./rho.rho_2).*sind(2*rho.angle) + ...
                                                    g.xz.*(1./rho.rho_2-1./rho.rho_1).*cosd(2*rho.angle));

                    % Uncomment if lambda sensitivity is needed (if you don't know wether it is needed, assume it is not)
                        % s_m = sqrt(rho.rho_1.*rho.rho_2);
                        % lambd = sqrt(rho.rho_2./rho.rho_1);
                        % Slambda(i,:) = - ( d_phi.xx.*s_m.*(cosd(rho.angle).^2 - 1./lambd.^2.*sind(rho.angle).^2) ...
                        %                  + d_phi.zz.*s_m.*(sind(rho.angle).^2 - 1./lambd.^2.*cosd(rho.angle).^2) ...
                        %                  - d_phi.xz.*0.5.*s_m.*(1 + 1./lambd.^2).*sind(2*rho.angle));

                end

                S.xx = Sxx;
                S.xz = Sxz;
                S.zz = Szz;
                % Uncomment if lambda sensitivity is needed
                    % S.lambda = Slambda;

            % pole-pole
            elseif isnan(XYZ.MEAS.C2(1,1)) && isnan(XYZ.MEAS.P2(1,1))

                parfor i = 1:n    

                    u_1  = u1(:,ind.meas(i,1):nb_electrode:end); % u(c1)
                    u_1R = u1(:,ind.meas(i,2):nb_electrode:end); % u(p1)

                    % Gradient calculation
                    g = gradient_product_anis(u_1R, u_1, param.h_x, param.h_z, param.k, param.flag.inv.p);

                    % xx
                    g.xx = integration_U(g.xx,param.wk,param.k,param.nbk1); % Inverse Fourier transform
                    g.xx = reshape(g.xx,nb_ligne,nb_col);
                    g.xx = reshape(g.xx',(nb_ligne)*(nb_col),1)'; % horizontally arranged

                    % xz
                    g.xz = integration_U(g.xz,param.wk,param.k,param.nbk1); % Inverse Fourier transform
                    g.xz = reshape(g.xz,nb_ligne,nb_col);
                    g.xz = reshape(g.xz',(nb_ligne)*(nb_col),1)'; % horizontally arranged

                    % zz
                    g.zz = integration_U(g.zz,param.wk,param.k,param.nbk1); % Inverse Fourier transform
                    g.zz = reshape(g.zz,nb_ligne,nb_col);
                    g.zz = reshape(g.zz',(nb_ligne)*(nb_col),1)'; % horizontally arranged

                    % yy
                    g.yy = integration_U(g.yy,param.wk,param.k,param.nbk1); % tr. Fourier inv.
                    g.yy = reshape(g.yy,nb_ligne,nb_col);
                    g.yy = reshape(g.yy',(nb_ligne)*(nb_col),1)'; % horizontally arranged

                    % Calculation of the sensitivity (ln(resistivity))
                    Sxx(i,:) = (1./(rho.rho_1*u(i))).*(g.xx.*(cosd(rho.angle)).^2 + ...
                                                       g.zz.*(sind(rho.angle)).^2 - ...
                                                       g.xz.*cosd(rho.angle).*sind(rho.angle) + ...
                                                       g.yy);

                    Szz(i,:) = (1./(rho.rho_2*u(i))).*(g.xx.*(sind(rho.angle)).^2 + ...
                                                       g.zz.*(cosd(rho.angle)).^2 + ...
                                                       g.xz.*cosd(rho.angle).*sind(rho.angle));

                    Sxz(i,:) = -(rho.angle./u(i)).*(g.xx.*(1./rho.rho_2-1./rho.rho_1).*sind(2*rho.angle) + ...
                                                    g.zz.*(1./rho.rho_1-1./rho.rho_2).*sind(2*rho.angle) + ...
                                                    g.xz.*(1./rho.rho_2-1./rho.rho_1).*cosd(2*rho.angle));

                    % Uncomment if lambda sensitivity is needed (if you don't know wether it's needed, assume it is not)
                        % s_m = sqrt(rho.rho_1.*rho.rho_2);
                        % lambd = sqrt(rho.rho_2./rho.rho_1);
                        % Slambda(i,:) = - ( d_phi.xx.*s_m.*(cosd(rho.angle).^2 - 1./lambd.^2.*sind(rho.angle).^2) ...
                        %                  + d_phi.zz.*s_m.*(sind(rho.angle).^2 - 1./lambd.^2.*cosd(rho.angle).^2) ...
                        %                  - d_phi.xz.*0.5.*s_m.*(1 + 1./lambd.^2).*sind(2*rho.angle));

                end

                S.xx = Sxx;
                S.xz = Sxz;
                S.zz = Szz;
                % Uncomment if lambda sensitivity is needed
                    % S.lambda = Slambda;

            % pole-dipole
            elseif isnan(XYZ.MEAS.C2(1,1)) && ~isnan(XYZ.MEAS.P2(1,1)) 

                parfor i = 1:n   

                    u_1  = u1(:,ind.meas(i,1):nb_electrode:end); % u(c1)
                    u_1R = (u1(:,ind.meas(i,2):nb_electrode:end)-u1(:,ind.meas(i,3):nb_electrode:end)); % u(p1)-u(p2)

                    % Gradient calculation
                    g = gradient_product_anis(u_1R, u_1, param.h_x, param.h_z, param.k, param.flag.inv.p);

                    % xx
                    g.xx = integration_U(g.xx,param.wk,param.k,param.nbk1); % Inverse Fourier transform
                    g.xx = reshape(g.xx,nb_ligne,nb_col);
                    g.xx = reshape(g.xx',(nb_ligne)*(nb_col),1)'; % horizontally arranged

                    % xz
                    g.xz = integration_U(g.xz,param.wk,param.k,param.nbk1); % Inverse Fourier transform
                    g.xz = reshape(g.xz,nb_ligne,nb_col);
                    g.xz = reshape(g.xz',(nb_ligne)*(nb_col),1)'; % horizontally arranged

                    % zz
                    g.zz = integration_U(g.zz,param.wk,param.k,param.nbk1); % Inverse Fourier transform
                    g.zz = reshape(g.zz,nb_ligne,nb_col);
                    g.zz = reshape(g.zz',(nb_ligne)*(nb_col),1)'; % horizontally arranged

                    % yy
                    g.yy = integration_U(g.yy,param.wk,param.k,param.nbk1); % tr. Fourier inv.
                    g.yy = reshape(g.yy,nb_ligne,nb_col);
                    g.yy = reshape(g.yy',(nb_ligne)*(nb_col),1)'; % horizontally arranged

                    % Calculation of the sensitivity (ln(resistivity))
                    Sxx(i,:) = (1./(rho.rho_1*u(i))).*(g.xx.*(cosd(rho.angle)).^2 + ...
                                                       g.zz.*(sind(rho.angle)).^2 - ...
                                                       g.xz.*cosd(rho.angle).*sind(rho.angle) + ...
                                                       g.yy);

                    Szz(i,:) = (1./(rho.rho_2*u(i))).*(g.xx.*(sind(rho.angle)).^2 + ...
                                                       g.zz.*(cosd(rho.angle)).^2 + ...
                                                       g.xz.*cosd(rho.angle).*sind(rho.angle));

                    Sxz(i,:) = -(rho.angle./u(i)).*(g.xx.*(1./rho.rho_2-1./rho.rho_1).*sind(2*rho.angle) + ...
                                                    g.zz.*(1./rho.rho_1-1./rho.rho_2).*sind(2*rho.angle) + ...
                                                    g.xz.*(1./rho.rho_2-1./rho.rho_1).*cosd(2*rho.angle));

                    % Uncomment if lambda sensitivity is needed (if you don't know wether it's needed, assume it is not)
                        % s_m = sqrt(rho.rho_1.*rho.rho_2);
                        % lambd = sqrt(rho.rho_2./rho.rho_1);
                        % Slambda(i,:) = - ( d_phi.xx.*s_m.*(cosd(rho.angle).^2 - 1./lambd.^2.*sind(rho.angle).^2) ...
                        %                  + d_phi.zz.*s_m.*(sind(rho.angle).^2 - 1./lambd.^2.*cosd(rho.angle).^2) ...
                        %                  - d_phi.xz.*0.5.*s_m.*(1 + 1./lambd.^2).*sind(2*rho.angle));

                end

                S.xx = Sxx;
                S.xz = Sxz;
                S.zz = Szz;
                % Uncomment if lambda sensitivity is needed
                    % S.lambda = Slambda;
            end
        end
    end
    
end
