function [A,Ak] = matrix_coeff_anis(h_x,h_z,grille,rho,k)

%%========================================================================%
%                                                                         %
% capacitance matrix calculation                                          %
%                                                                         %
%%========================================================================%
%                                                                         %
% In:                                                                     %
% -------                                                                 %
%                                                                         %
%   h_x : x-steps vector                                                  %
%   h_z : z-steps vector                                                  %
%   rho : conductivity matrix (structure) defined as follows              %
%           |¯               ¯|                                           %
%   rho =   | rho.xx   rho.xz |                                           %
%           |                 |                                           %
%           | rho.xz   rho.zz |                                           %
%           |_               _|                                           %
%   rho.yy = max(eig(rho)): mandatory for the 2.5D calculus               %
%   k : wave number vector                                                %
%   grille :                                                              %
%      grille(:,1) = x_position                                           %
%      grille(:,2) = z_position                                           %
%                                                                         %
% Out:                                                                    %
% -----------                                                             %
%   A and Ak such as G = A + diag(Ak(:,nb_k))                             %
%   G u = q                                                               %
%                                                                         %
% test :                                                                  %
%  --------                                                               %
%   hz=ones(10,1);hx=ones(10,1);                                          %
%   [X,Y] = meshgrid(0:9,0:9);                                            %
%   x=X(:);                                                               %
%   y=Y(:);grille=[];                                                     %
%   grille(:,1)=x;                                                        %
%   grille(:,2)=y;                                                        %
%   sigma=ones(10,10);                                                    %
%   [grille,h_x,h_z] = grille2d_elect(hx,hz,1.3,10,[1;1],0,1);            %
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
%     bouchedda@geo.polymtl.ca                                            %
%     Ecole polytechnique de Montreal                                     %
%     Departement de géophysique appliquée                                %
%     http://geo.polymtl.ca                                               %
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

%% Calculation

    assert(isequal(size(rho.xx) , size(rho.zz)),'The sigma tensor is inconsistent');
    assert(isequal(size(rho.xx) , size(rho.xz)),'The sigma tensor is inconsistent');
    % assert(isequal(size(rho.xx) , size(rho.yy)),'The sigma tensor is inconsistent');

    % initialisation
    %--------------------------
    [n,m]=size(rho.xx);
    nb_k=length(k);

    C_L  = zeros(n*m,1);%m-n
    C_R  = zeros(n*m,1);%m-n
    C_B  = zeros(n*m,1);%m-1
    C_T  = zeros(n*m,1);%m-n
    C_TL  = zeros(n*m,1);%m-n
    C_TR  = zeros(n*m,1);%m-n
    C_BL  = zeros(n*m,1);%m-1
    C_BR  = zeros(n*m,1);%m-n


    % C_P1 = zeros(n*xm,1);
    % A    = zeros(n*m,1);
    Ak   = zeros(n*m,nb_k);
    h_x  = h_x(:);
    h_z  = h_z(:);

    rho.xx = rho.xx(:);
    rho.xz = rho.xz(:);
    rho.zz = rho.zz(:);

    for i = 1:length(rho.xx)
        sigma_temp = pinv([rho.xx(i) rho.xz(i); rho.xz(i) rho.zz(i)]);
        sigma.xx(i) = sigma_temp(1);
        sigma.xz(i) = sigma_temp(2);
        sigma.zz(i) = sigma_temp(4);
        sigma.yy(i) = max(eig(sigma_temp));
    end

    sigma.xx = sigma.xx(:);
    sigma.xz = sigma.xz(:);
    sigma.zz = sigma.zz(:);

    rho.yy = 1./sigma.yy;
    rho.yy = rho.yy(:);

    sigma.yy = sigma.yy(:);

    k=k(:)';

    ind_c = n*floor(m/2)+1; % consider the middle surface point
    X = (grille(:,1)-grille(ind_c,1)); 
    Z = (grille(:,2)-grille(ind_c,2)); 

    % index nodes for boundary condition (bc)
    %------------------------------------
    ind_Rtc = sub2ind([n m],1,m);                  % index right top corner
    ind_Ltc = sub2ind([n m],1,1);                  % index left top corner
    ind_Rbc = sub2ind([n m],n,m);                  % index right bottom corner
    ind_Lbc = sub2ind([n m],n,1);                  % index left bottom corner
    ind_T   = sub2ind([n m],ones(1,m-2),2:m-1);    % index top edge
    ind_B   = sub2ind([n m],n*ones(1,m-2),2:m-1);  % index bottom edge
    ind_R   = sub2ind([n m],2:n-1,m*ones(1,n-2));  % index right edge
    ind_L   = sub2ind([n m],2:n-1,ones(1,n-2));    % index left edge
    ind     = 1:n*m;                               % index for all nodes
    ind_wbc = ind;
    ind_wbc ([ind_Ltc  ind_Lbc ind_L ind_R ind_Rtc ind_Rbc ind_B ind_T]) = []; % index without bc

    dz = repmat([h_z;1],m,1);
    dx = repmat([h_x' 1],n,1);
    dx=dx(:);

    % coeff of non-boundary cd nodes
    %==========================================================================

    C_L(ind_wbc) = ( - dz(ind_wbc).*sigma.xx(ind_wbc-n) - dz(ind_wbc-1).*sigma.xx(ind_wbc-n-1)) ./ (2*dx(ind_wbc-n));
    C_R(ind_wbc) = ( - dz(ind_wbc).*sigma.xx(ind_wbc)   - dz(ind_wbc-1).*sigma.xx(ind_wbc-1))   ./ (2*dx(ind_wbc));
    C_T(ind_wbc) = ( - dx(ind_wbc).*sigma.zz(ind_wbc-1) - dx(ind_wbc-n).*sigma.zz(ind_wbc-n-1)) ./ (2*dz(ind_wbc-1));
    C_B(ind_wbc) = ( - dx(ind_wbc).*sigma.zz(ind_wbc)   - dx(ind_wbc-n).*sigma.zz(ind_wbc-n))   ./ (2*dz(ind_wbc));

    C_P1 = -(C_L + C_R + C_T + C_B);

    ddx = 0.25./(dx(ind_wbc)+dx(ind_wbc-n));
    ddz = 0.25./(dz(ind_wbc)+dz(ind_wbc-1));

    C_L(ind_wbc) = C_L(ind_wbc) + ddx.*(...
        dx(ind_wbc).*sigma.xz(ind_wbc)+dx(ind_wbc-n).*sigma.xz(ind_wbc-n)-...
        dx(ind_wbc).*sigma.xz(ind_wbc-1)-dx(ind_wbc-n).*sigma.xz(ind_wbc-n-1));

    C_R(ind_wbc) = C_R(ind_wbc) + ddx.*(...
        dx(ind_wbc).*sigma.xz(ind_wbc-1)+dx(ind_wbc-n).*sigma.xz(ind_wbc-n-1)-...
        dx(ind_wbc).*sigma.xz(ind_wbc)-dx(ind_wbc-n).*sigma.xz(ind_wbc-n));


    C_T(ind_wbc) = C_T(ind_wbc) + ddz.*(...
        dz(ind_wbc).*sigma.xz(ind_wbc)+dz(ind_wbc-1).*sigma.xz(ind_wbc-1)-...
        dz(ind_wbc).*sigma.xz(ind_wbc-n)-dz(ind_wbc-1).*sigma.xz(ind_wbc-n-1));

    C_B(ind_wbc) = C_B(ind_wbc) + ddz.*(...
        -dz(ind_wbc).*sigma.xz(ind_wbc)-dz(ind_wbc-1).*sigma.xz(ind_wbc-1)+...
        dz(ind_wbc).*sigma.xz(ind_wbc-n)+dz(ind_wbc-1).*sigma.xz(ind_wbc-n-1));

    C_TL(ind_wbc) = -ddx.*(dx(ind_wbc).*sigma.xz(ind_wbc-1)+dx(ind_wbc-n).*sigma.xz(ind_wbc-n-1))...
                    -ddz.*(dz(ind_wbc).*sigma.xz(ind_wbc-n)+dz(ind_wbc-1).*sigma.xz(ind_wbc-n-1));

    C_BL(ind_wbc) = ddx.*(dx(ind_wbc).*sigma.xz(ind_wbc)+dx(ind_wbc-n).*sigma.xz(ind_wbc-n))...
                    +ddz.*(dz(ind_wbc).*sigma.xz(ind_wbc-n)+dz(ind_wbc-1).*sigma.xz(ind_wbc-n-1));

    C_TR(ind_wbc) = ddx.*(dx(ind_wbc).*sigma.xz(ind_wbc-1)+dx(ind_wbc-n).*sigma.xz(ind_wbc-n-1))...
                    +ddz.*(dz(ind_wbc-1).*sigma.xz(ind_wbc-1)+dz(ind_wbc).*sigma.xz(ind_wbc));

    C_BR(ind_wbc) = -ddx.*(dx(ind_wbc).*sigma.xz(ind_wbc)+dx(ind_wbc-n).*sigma.xz(ind_wbc-n))...
                    -ddz.*(dz(ind_wbc).*sigma.xz(ind_wbc)+dz(ind_wbc-1).*sigma.xz(ind_wbc-1));


    % Use of sigma.yy in this part
    A = (sigma.yy(ind_wbc-n-1).*dx(ind_wbc-n).*dz(ind_wbc-1)/4 + sigma.yy(ind_wbc-1).*dx(ind_wbc).*dz(ind_wbc-1)/4 +...
                  sigma.yy(ind_wbc).*dx(ind_wbc).*dz(ind_wbc)/4 + sigma.yy(ind_wbc-n).*dx(ind_wbc-n).*dz(ind_wbc)/4);
    la = length(A);
    k1 = repmat(k,la,1); k1=k1.^2;
    A = repmat(A,1,nb_k);
    Ak(ind_wbc,:) = k1.*A;

    % coeff of boundary cd nodes
    %==========================================================================

    % nodes located on the ground surface
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    C_L(ind_T) = ( - dz(ind_T).*sigma.xx(ind_T-n))./(2*dx(ind_T-n));

    % C_T(ind_T) = 0;

    C_R(ind_T) = ( - dz(ind_T).*sigma.xx(ind_T))./(2*dx(ind_T)); 

    C_B(ind_T) = (- dx(ind_T).*sigma.zz(ind_T) - dx(ind_T-n).*sigma.zz(ind_T-n))./(2*dz(ind_T));

    C_P1(ind_T) = - (C_L(ind_T) + C_R(ind_T) + C_B(ind_T));


    ddx = 0.25./(dx(ind_T)+dx(ind_T-n));


    C_L(ind_T) = C_L(ind_T) - 0.125*sigma.xz(ind_T-n) + ddx.*(...
        dx(ind_T).*sigma.xz(ind_T)+dx(ind_T-n).*sigma.xz(ind_T-n));

    C_R(ind_T) = C_R(ind_T) + 0.125*sigma.xz(ind_T) + ddx.*(...
        -dx(ind_T).*sigma.xz(ind_T)-dx(ind_T-n).*sigma.xz(ind_T-n));

    C_B(ind_T) = C_B(ind_T) +0.125*(-sigma.xz(ind_T) + sigma.xz(ind_T-n));

    C_BL(ind_T) = 0.125*sigma.xz(ind_T-n)+ ddx.*(dx(ind_T).*sigma.xz(ind_T)+dx(ind_T-n).*sigma.xz(ind_T-n));

    C_BR(ind_T) = -0.125*sigma.xz(ind_T)-ddx.*(dx(ind_T).*sigma.xz(ind_T)+dx(ind_T-n).*sigma.xz(ind_T-n));

    C_P1(ind_T)= C_P1(ind_T) + 0.125*(sigma.xz(ind_T) - sigma.xz(ind_T-n));

    A = (sigma.yy(ind_T).*dx(ind_T).*dz(ind_T)/4 + sigma.yy(ind_T-n).*dx(ind_T-n).*dz(ind_T)/4);
    la = length(A);
    k1 = repmat(k,la,1); k1=k1.^2;
    A = repmat(A,1,nb_k);
    Ak(ind_T,:) = k1.*A; 

    % nodes located on the top left corner
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %C_L (ind_Ltc) = 0;
    C_R (ind_Ltc) = (- dz(ind_Ltc).*sigma.xx(ind_Ltc))./(2*dx(ind_Ltc));
    C_B (ind_Ltc) = ( - dx(ind_Ltc).*sigma.zz(ind_Ltc))./(2*dz(ind_Ltc));

    C_P1(ind_Ltc) = -(C_R(ind_Ltc) + C_B(ind_Ltc)) + 0.5*sigma.xz(ind_Ltc);

    C_BR(ind_Ltc) = -0.5*sigma.xz(ind_Ltc);

    A  = (sigma.yy(ind_Ltc).*dx(ind_Ltc).*dz(ind_Ltc)/4);
    la = length(A);
    k1 = repmat(k,la,1); k1=k1.^2;
    A = repmat(A,1,nb_k);

    XLtc = X(ind_Ltc,1);
    ZLtc = Z(ind_Ltc,1);
    pos = [XLtc ZLtc];

    sigma_Ltc.xx = sigma.xx(ind_Ltc);
    sigma_Ltc.xz = sigma.xz(ind_Ltc);
    sigma_Ltc.zz = sigma.zz(ind_Ltc);
    sigma_Ltc.yy = sigma.yy(ind_Ltc);
    rho_Ltc.xx = rho.xx(ind_Ltc);
    rho_Ltc.xz = rho.xz(ind_Ltc);
    rho_Ltc.zz = rho.zz(ind_Ltc);
    rho_Ltc.yy = rho.yy(ind_Ltc);
    rho_Ltc_M = [rho_Ltc.xx rho_Ltc.xz
                   rho_Ltc.xz rho_Ltc.zz];
    a = 1./rho_Ltc.yy.*pos*rho_Ltc_M*pos';
    assert(a>0,'Anisotropic tensor is not positive defined !')

    normale = [-1;0];

    alpha_L = calcul_alpha_anis(a,k,rho_Ltc,sigma_Ltc,pos,normale)*dz(ind_Ltc)/2; % boundary term
    Ak(ind_Ltc,:) = k1.*A + alpha_L;


    % nodes located on the top right corner
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    C_L (ind_Rtc) = ( - dz(ind_Rtc).*sigma.xx(ind_Rtc-n))./(2*dx(ind_Rtc-n));
    %C_R (ind_Rtc) = 0;
    %C_T (ind_Rtc) = 0;
    C_B (ind_Rtc) = ( - dx(ind_Rtc-n).*sigma.zz(ind_Rtc-n))./(2*dz(ind_Rtc));

    C_P1 (ind_Rtc) = -(C_L (ind_Rtc) + C_B (ind_Rtc))-0.5*sigma.xz(ind_Rtc-n);

    C_BL(ind_Rtc) = -0.5*sigma.xz(ind_Rtc-n);

    A = sigma.yy(ind_Rtc-n).*dx(ind_Rtc-n).*dz(ind_Rtc)/4;
    la = length(A);
    k1 = repmat(k,la,1); k1=k1.^2;          
    A = repmat(A,1,nb_k);

    XRtc = X(ind_Rtc,1); % /!\ nodes position
    ZRtc = Z(ind_Rtc,1); % /!\ nodes position
    pos = [XRtc ZRtc];

    sigma_Rtc.xx = sigma.xx(ind_Rtc-n); % /!\ cells position =/= nodes position
    sigma_Rtc.xz = sigma.xz(ind_Rtc-n); % /!\ cells position =/= nodes position
    sigma_Rtc.zz = sigma.zz(ind_Rtc-n); % /!\ cells position =/= nodes position
    sigma_Rtc.yy = sigma.yy(ind_Rtc-n); % /!\ cells position =/= nodes position
    rho_Rtc.xx = rho.xx(ind_Rtc-n); % /!\ cells position =/= nodes position
    rho_Rtc.xz = rho.xz(ind_Rtc-n); % /!\ cells position =/= nodes position
    rho_Rtc.zz = rho.zz(ind_Rtc-n); % /!\ cells position =/= nodes position
    rho_Rtc.yy = rho.yy(ind_Rtc-n); % /!\ cells position =/= nodes position
    rho_Rtc_M = [rho_Rtc.xx rho_Rtc.xz
                   rho_Rtc.xz rho_Rtc.zz];
    a = 1./rho_Rtc.yy.*pos*rho_Rtc_M*pos';
    assert(a>0,'Anisotropic tensor is not positive defined !')

    normale = [1;0];

    alpha_R = calcul_alpha_anis(a,k,rho_Rtc,sigma_Rtc,pos,normale)*dz(ind_Rtc)/2; 
    Ak(ind_Rtc,:) = k1.*A + alpha_R;


    % % nodes located on the bottom left corner 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %C_L (ind_Lbc) = 0;
    C_R (ind_Lbc) = ( - dz(ind_Lbc-1).*sigma.xx(ind_Lbc-1))./(2*dx(ind_Lbc));
    C_T (ind_Lbc) = ( - dx(ind_Lbc).*sigma.zz(ind_Lbc-1))./(2*dz(ind_Lbc-1));
    %C_B (ind_Lbc) = 0;

    C_P1(ind_Lbc) = -(C_R(ind_Lbc) + C_T(ind_Lbc))-0.5*sigma.xz(ind_Lbc-1);

    C_TR(ind_Lbc) = 0.5*sigma.xz(ind_Lbc-1);

    A = (sigma.yy(ind_Lbc-1).*dx(ind_Lbc).*dz(ind_Lbc-1)/4);
    la = length(A);
    k1 = repmat(k,la,1); k1=k1.^2;
    A = repmat(A,1,nb_k);

    XLbc = X(ind_Lbc,1);
    ZLbc = Z(ind_Lbc,1);
    pos = [XLbc ZLbc];

    sigma_Lbc.xx = sigma.xx(ind_Lbc-1);
    sigma_Lbc.xz = sigma.xz(ind_Lbc-1);
    sigma_Lbc.zz = sigma.zz(ind_Lbc-1);
    sigma_Lbc.yy = sigma.yy(ind_Lbc-1);
    rho_Lbc.xx = rho.xx(ind_Lbc-1);
    rho_Lbc.xz = rho.xz(ind_Lbc-1);
    rho_Lbc.zz = rho.zz(ind_Lbc-1);
    rho_Lbc.yy = rho.yy(ind_Lbc-1);
    rho_Lbc_M = [rho_Lbc.xx rho_Lbc.xz
                   rho_Lbc.xz rho_Lbc.zz];
    a = 1/rho_Lbc.yy*pos*rho_Lbc_M*pos';
    assert(a>0,'Anisotropic tensor is not positive defined !')

    % Left part
    normale = [-1;0];
    alpha_L1 = calcul_alpha_anis(a,k,rho_Lbc,sigma_Lbc,pos,normale)*dz(ind_Lbc-1)/2;

    % Bottom part
    normale = [0;1];
    alpha_L2 = calcul_alpha_anis(a,k,rho_Lbc,sigma_Lbc,pos,normale)*dx(ind_Lbc)/2;

    Ak(ind_Lbc,:) = k1.*A + alpha_L1 + alpha_L2;

    % % nodes located on the bottom right corner 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    C_L (ind_Rbc) = ( - dz(ind_Rbc-1).*sigma.xx(ind_Rbc-n-1))./(2*dx(ind_Rbc-n));
    % C_R (ind_Rbc) = 0;
    C_T (ind_Rbc) = ( - dx(ind_Rbc-n).*sigma.zz(ind_Rbc-n-1))./(2*dz(ind_Rbc-1));
    % C_B (ind_Rbc) = 0;


    C_P1(ind_Rbc) =  -(C_L(ind_Rbc) + C_T(ind_Rbc))+0.5*sigma.xz(ind_Rbc-n-1);

    C_TL(ind_Rbc) = 0.5*sigma.xz(ind_Rbc-n-1);

    A = (sigma.yy(ind_Rbc-n-1).*dx(ind_Rbc-n).*dz(ind_Rbc-1)/4);
    la = length(A);
    k1 = repmat(k,la,1); k1=k1.^2;          
    A = repmat(A,1,nb_k);

    XRbc = X(ind_Rbc,1);
    ZRbc = Z(ind_Rbc,1);
    pos = [XRbc ZRbc];

    sigma_Rbc.xx = sigma.xx(ind_Rbc-n-1);
    sigma_Rbc.xz = sigma.xz(ind_Rbc-n-1);
    sigma_Rbc.zz = sigma.zz(ind_Rbc-n-1);
    sigma_Rbc.yy = sigma.yy(ind_Rbc-n-1);
    rho_Rbc.xx = rho.xx(ind_Rbc-n-1);
    rho_Rbc.xz = rho.xz(ind_Rbc-n-1);
    rho_Rbc.zz = rho.zz(ind_Rbc-n-1);
    rho_Rbc.yy = rho.yy(ind_Rbc-n-1);
    rho_Rbc_M = [rho_Rbc.xx rho_Rbc.xz
                 rho_Rbc.xz rho_Rbc.zz];
    a = (1/rho_Rbc.yy).*pos*rho_Rbc_M*pos';
    assert(a>0,'Anisotropic tensor is not positive defined !')

    % Right part
    normale = [1;0];
    alpha_L1 = calcul_alpha_anis(a,k,rho_Rbc,sigma_Rbc,pos,normale)*dz(ind_Rbc-1)/2;

    % Bottom part
    normale = [0;1];
    alpha_L2 = calcul_alpha_anis(a,k,rho_Rbc,sigma_Rbc,pos,normale)*dx(ind_Rbc-n)/2; 

    Ak(ind_Rbc,:) = k1.*A + alpha_L1 + alpha_L2;

    % nodes located on the bottom edge 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    C_L (ind_B) = ( - dz(ind_B-1).*sigma.xx(ind_B-n-1))./(2*dx(ind_B-n));
    C_R (ind_B) = ( - dz(ind_B-1).*sigma.xx(ind_B-1))./(2*dx(ind_B));
    C_T (ind_B) = ( - dx(ind_B).*sigma.zz(ind_B-1) - ...
        dx(ind_B-n).*sigma.zz(ind_B-n-1))./(2*dz(ind_B-1));
    % C_B (ind_B) = 0;

    C_P1 (ind_B) = -(C_L (ind_B) + C_R (ind_B) + C_T (ind_B));

    ddx = 0.25./(dx(ind_B)+dx(ind_B-n));
    % ddz = 0.25./(dz(ind_B-1)+dz(ind_B-1));

    C_L(ind_B) = C_L(ind_B) + 0.25*sigma.xz(ind_B-n-1)+ ddx.*(...
        -dx(ind_B).*sigma.xz(ind_B-1)-dx(ind_B-n).*sigma.xz(ind_B-n-1));

    C_R(ind_B) = C_R(ind_B) - 0.25*sigma.xz(ind_B-1)+ ddx.*(...
        dx(ind_B).*sigma.xz(ind_B-1)+dx(ind_B-n).*sigma.xz(ind_B-n-1));

    C_T(ind_B) = C_T(ind_B) + 0.25*(sigma.xz(ind_B-1)-sigma.xz(ind_B-n-1));

    C_P1 (ind_B)=C_P1 (ind_B) - 0.25*(sigma.xz(ind_B-1)-sigma.xz(ind_B-n-1));

    C_TL(ind_B) = -0.25*sigma.xz(ind_B-n-1)-ddx.*(dx(ind_B).*sigma.xz(ind_B-1)...
                                                  +dx(ind_B-n).*sigma.xz(ind_B-n-1));

    C_TR(ind_B) = 0.25*sigma.xz(ind_B-n-1)+ ddx.*(dx(ind_B).*sigma.xz(ind_B-1)...
                                                 +dx(ind_B-n).*sigma.xz(ind_B-n-1));

    A = (sigma.yy(ind_B-n-1).*dx(ind_B-n).*dz(ind_B-1)/4 + sigma.yy(ind_B-1).*dx(ind_B).*dz(ind_B-1)/4);
    la = length(A);
    k1 = repmat(k,la,1); k1=k1.^2;          
    A = repmat(A,1,nb_k);

    XB = X(ind_B,1);
    ZB = Z(ind_B,1);
    pos = [XB ZB];

    normale = [0;1];

    % Left part
    sigma_B.xx = sigma.xx(ind_B-n-1); 
    sigma_B.xz = sigma.xz(ind_B-n-1);
    sigma_B.zz = sigma.zz(ind_B-n-1);
    sigma_B.yy = sigma.yy(ind_B-n-1);
    rho_B.xx = rho.xx(ind_B-n-1); 
    rho_B.xz = rho.xz(ind_B-n-1);
    rho_B.zz = rho.zz(ind_B-n-1);
    rho_B.yy = rho.yy(ind_B-n-1);
    % a formulation is expanded because pos is a matrix and not a 2-uplet
    % anymore
    a = 1./rho_B.yy.*(XB.^2.*rho_B.xx + 2*XB.*ZB.*rho_B.xz + ZB.^2.*rho_B.zz);
    assert(isempty(find(a<=0)),'Anisotropic tensor is not positive defined !')
    alpha_B1 = calcul_alpha_anis(a,k,rho_B,sigma_B,pos,normale);
    dxL = repmat(dx(ind_B-n)/2,1,nb_k);
    alpha_B1 = alpha_B1.*dxL; 


    % Right part
    sigma_B.xx = sigma.xx(ind_B-1); 
    sigma_B.xz = sigma.xz(ind_B-1);
    sigma_B.zz = sigma.zz(ind_B-1);
    sigma_B.yy = sigma.yy(ind_B-1);
    rho_B.xx = rho.xx(ind_B-1); 
    rho_B.xz = rho.xz(ind_B-1);
    rho_B.zz = rho.zz(ind_B-1);
    rho_B.yy = rho.yy(ind_B-1);
    a = 1./rho_B.yy.*(XB.^2.*rho_B.xx + 2*XB.*ZB.*rho_B.xz + ZB.^2.*rho_B.zz);
    assert(isempty(find(a<=0)),'Anisotropic tensor is not positive defined !')
    alpha_B2 = calcul_alpha_anis(a,k,rho_B,sigma_B,pos,normale); 
    dxR = repmat(dx(ind_B)/2,1,nb_k);
    alpha_B2 = alpha_B2.*dxR; 

    Ak(ind_B,:) = k1.*A + alpha_B1 + alpha_B2;

    %  nodes located on the left edge
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %C_L (ind_L) = 0;
    C_R (ind_L) = ( - dz(ind_L).*sigma.xx(ind_L) - dz(ind_L-1).*sigma.xx(ind_L-1))./(2*dx(ind_L));
    C_T (ind_L) = ( - dx(ind_L).*sigma.zz(ind_L-1))./(2*dz(ind_L-1));
    C_B (ind_L) = ( - dx(ind_L).*sigma.zz(ind_L))./(2*dz(ind_L));

    C_P1(ind_L) = -(C_R(ind_L) + C_T(ind_L) + C_B(ind_L));

    ddz = 0.25./(dz(ind_L)+dz(ind_L-1));

    C_R(ind_L) = C_R(ind_L) + 0.125*(sigma.xz(ind_L-1)-sigma.xz(ind_L));

    C_P1(ind_L) = C_P1(ind_L) - 0.125*(sigma.xz(ind_L-1)-sigma.xz(ind_L));

    C_T(ind_L) = C_T(ind_L) - 0.125*sigma.xz(ind_L-1)+ ddz.*(...
        dz(ind_L).*sigma.xz(ind_L)+dz(ind_L-1).*sigma.xz(ind_L-1));

    C_B(ind_L) = C_B(ind_L) + 0.125*sigma.xz(ind_L)+ ddz.*(...
        -dz(ind_L).*sigma.xz(ind_L)-dz(ind_L-1).*sigma.xz(ind_L-1));

    C_TR(ind_L) =  0.125*sigma.xz(ind_L-1) + ddz.*(dz(ind_L-1).*sigma.xz(ind_L-1)+...
                                                              dz(ind_L).*sigma.xz(ind_L));

    C_BR(ind_L) = -0.125*sigma.xz(ind_L)-ddz.*(dz(ind_L).*sigma.xz(ind_L)+...
                                              dz(ind_L-1).*sigma.xz(ind_L-1));


    A = (sigma.yy(ind_L-1).*dx(ind_L).*dz(ind_L-1)/4 + sigma.yy(ind_L).*dx(ind_L).*dz(ind_L)/4);
    la = length(A);
    k1 = repmat(k,la,1); k1=k1.^2;          
    A = repmat(A,1,nb_k);

    XL = X(ind_L,1);
    ZL = Z(ind_L,1);
    pos = [XL ZL];
    normale = [-1;0];

    % Upper part
    sigma_U.xx = sigma.xx(ind_L-1); 
    sigma_U.xz = sigma.xz(ind_L-1);
    sigma_U.zz = sigma.zz(ind_L-1); 
    sigma_U.yy = sigma.yy(ind_L-1);
    rho_U.xx = rho.xx(ind_L-1); 
    rho_U.xz = rho.xz(ind_L-1);
    rho_U.zz = rho.zz(ind_L-1); 
    rho_U.yy = rho.yy(ind_L-1);
    % a formulation is expanded
    a = (1./rho_U.yy).*(XL.^2.*rho_U.xx + 2*XL.*ZL.*rho_U.xz + ZL.^2.*rho_U.zz);
    assert(isempty(find(a<=0)),'Anisotropic tensor is not positive defined !')
    alpha_L1 = calcul_alpha_anis(a,k,rho_U,sigma_U,pos,normale);
    dzU = repmat(dz(ind_L-1)/2,1,nb_k);
    alpha_L1 = alpha_L1.*dzU;

    % Lower part
    sigma_Lo.xx = sigma.xx(ind_L); 
    sigma_Lo.xz = sigma.xz(ind_L); 
    sigma_Lo.zz = sigma.zz(ind_L); 
    sigma_Lo.yy = sigma.yy(ind_L);
    rho_Lo.xx = rho.xx(ind_L); 
    rho_Lo.xz = rho.xz(ind_L); 
    rho_Lo.zz = rho.zz(ind_L); 
    rho_Lo.yy = rho.yy(ind_L);
    a = 1./rho_Lo.yy.*(XL.^2.*rho_Lo.xx + 2*XL.*ZL.*rho_Lo.xz + ZL.^2.*rho_Lo.zz);
    assert(isempty(find(a<=0)),'Anisotropic tensor is not positive defined !')
    alpha_L2 = calcul_alpha_anis(a,k,rho_Lo,sigma_Lo,pos,normale); 
    dzLo = repmat(dz(ind_L)/2,1,nb_k);
    alpha_L2 = alpha_L2.*dzLo;

    Ak(ind_L,:) = k1.*A + alpha_L1 + alpha_L2;


    %  nodes located on the right edge
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    C_L (ind_R) = (dx(ind_R-n).*(-sigma.xz(ind_R-n-1)) - ...
        dz(ind_R).*sigma.xx(ind_R-n) - dz(ind_R-1).*sigma.xx(ind_R-n-1))./(2*dx(ind_R-n));
    %C_R (ind_R) = 0;
    C_T (ind_R) = ( - dx(ind_R-n).*sigma.zz(ind_R-n-1))./(2*dz(ind_R-1));
    C_B (ind_R) = (- dx(ind_R-n).*sigma.zz(ind_R-n))./(2*dz(ind_R));

    C_P1(ind_R) = -(C_L(ind_R) + C_T (ind_R) + C_B (ind_R));


    ddz=0.25./(dz(ind_R)+dz(ind_R-1));

    C_L(ind_R) = C_L(ind_R) + 0.125*(sigma.xz(ind_R-n)-sigma.xz(ind_R-n-1));

    C_P1(ind_R) = C_P1(ind_R) - 0.125*(sigma.xz(ind_R-n)-sigma.xz(ind_R-n-1));

    C_T(ind_R) = C_T(ind_R) + 0.125*sigma.xz(ind_R-n-1)+ ddz.*(...
        -dz(ind_R).*sigma.xz(ind_R-n)-dz(ind_R-1).*sigma.xz(ind_R-n-1));

    C_B(ind_R) = C_B(ind_R) - 0.125*sigma.xz(ind_R-n) + ddz.*(dz(ind_R).*sigma.xz(ind_R-n)+...
                                                             dz(ind_R-1).*sigma.xz(ind_R-n-1));

    C_TL(ind_R) = -0.125*sigma.xz(ind_R-n-1)-ddz.*(dz(ind_R).*sigma.xz(ind_R-n)+...
                                                  dz(ind_R-1).*sigma.xz(ind_R-n-1));

    C_BL(ind_R) = 0.125*sigma.xz(ind_R-n)+ddz.*(dz(ind_R).*sigma.xz(ind_R-n)+...
                                               dz(ind_R-1).*sigma.xz(ind_R-n-1));


    A = (sigma.yy(ind_R-n-1).*dx(ind_R-n).*dz(ind_R-1)/4 + sigma.yy(ind_R-n).*dx(ind_R-n).*dz(ind_R)/4);
    la = length(A);
    k1 = repmat(k,la,1); k1=k1.^2;          
    A = repmat(A,1,nb_k);

    XR = X(ind_R,1);
    ZR = Z(ind_R,1);
    pos = [XR ZR];
    normale = [1;0];

    % Upper part
    sigma_U.xx = sigma.xx(ind_R-n-1); 
    sigma_U.xz = sigma.xz(ind_R-n-1);
    sigma_U.zz = sigma.zz(ind_R-n-1);
    sigma_U.yy = sigma.yy(ind_R-n-1);
    rho_U.xx = rho.xx(ind_R-n-1); 
    rho_U.xz = rho.xz(ind_R-n-1);
    rho_U.zz = rho.zz(ind_R-n-1);
    rho_U.yy = rho.yy(ind_R-n-1);
    a = 1./rho_U.yy.*(XR.^2.*rho_U.xx + 2*XR.*ZR.*rho_U.xz + ZR.^2.*rho_U.zz);
    assert(isempty(find(a<=0)),'Anisotropic tensor is not positive defined !')
    alpha_L1 = calcul_alpha_anis(a,k,rho_U,sigma_U,pos,normale);
    dzU = repmat(dz(ind_R-1)/2,1,nb_k);
    alpha_L1 = alpha_L1.*dzU;

    % Lower part
    sigma_Lo.xx = sigma.xx(ind_R-n);
    sigma_Lo.xz = sigma.xz(ind_R-n);
    sigma_Lo.zz = sigma.zz(ind_R-n);
    sigma_Lo.yy = sigma.yy(ind_R-n);
    rho_Lo.xx = rho.xx(ind_R-n);
    rho_Lo.xz = rho.xz(ind_R-n);
    rho_Lo.zz = rho.zz(ind_R-n);
    rho_Lo.yy = rho.yy(ind_R-n);
    a = 1./rho_Lo.yy.*(XR.^2.*rho_Lo.xx + 2*XR.*ZR.*rho_Lo.xz + ZR.^2.*rho_Lo.zz);
    assert(isempty(find(a<=0)),'Anisotropic tensor is not positive defined !')
    alpha_L2 = calcul_alpha_anis(a,k,rho_Lo,sigma_Lo,pos,normale); 
    dzLo = repmat(dz(ind_R)/2,1,nb_k);
    alpha_L2 = alpha_L2.*dzLo;

    Ak(ind_R,:) = k1.*A + alpha_L1 + alpha_L2;

    % Bulding up the entire C matrix
    %==========================================================================
    C_L = C_L(n+1:end);
    C_R = C_R(1:end-n);
    C_T = C_T(2:end);
    C_B = C_B(1:end-1);

    C_TL = C_TL(n+2:end);
    C_BL = C_BL(n:end);
    C_BR = C_BR(1:end-n-1);
    C_TR = C_TR(1:end-n+1);

    ind = 1:n*m;
    A = sparse(ind,ind,C_P1,n*m,n*m);
    ind1 = n+1:n*m; 
    ind2 = 1:n*m-n;
    A = A + sparse(ind1,ind2,C_L,n*m,n*m);
    A = A + sparse(ind2,ind1,C_R,n*m,n*m);
    ind1 = 2:n*m;
    ind2 = 1:n*m-1;
    A = A + sparse(ind1,ind2,C_T,n*m,n*m);
    A = A + sparse(ind2,ind1,C_B,n*m,n*m);


    ind1=n+2:n*m; 
    ind2 = 1:n*m-n-1;
    A = A + sparse(ind1,ind2,C_TL,n*m,n*m);

    ind1=n:n*m; 
    ind2 = 1:n*m-n+1;
    A = A + sparse(ind1,ind2,C_BL,n*m,n*m);

    ind1=n:n*m; 
    ind2 = 1:n*m-n+1;
    A = A + sparse(ind2,ind1,C_TR,n*m,n*m);

    ind1=n+2:n*m; 
    ind2 = 1:n*m-n-1;
    A = A + sparse(ind2,ind1,C_BR,n*m,n*m);

end