function gpa = gradient_product_anis(G,u,h_x,h_z,k)

%%========================================================================%
%                                                                         %
% Calculation of the finite-differences gradient needed for sensitivity   %
% in anisotropic media using the adjoint equation approach                %
% for each discrete frequency                                             %
%                                                                         %
%%========================================================================%
%                                                                         %
% In:                                                                     %
% -----------                                                             %
%   G : adjoint Green's function.                                         %
%       For the electrical resistivity problem, G is the reciprocal       %
%       potentials                                                        %
%   u : potential data                                                    %
%   h_x : x step vector                                                   %
%   h_z : z step vector                                                   %
%   k : wave number vector                                                %
%                                                                         %
% Out:                                                                    %
% -----------                                                             %
%   g : anisotropic gradient product                                      %
%                                                                         %
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

    n = length(h_z)+1;
    m = length(h_x)+1;
    la = n*m;
    nb_k = length(k);
    gpa.xx = zeros(n*m,nb_k); %m-n ; phi : McGillivray's notation
    gpa.xz = zeros(n*m,nb_k);
    gpa.zz = zeros(n*m,nb_k);
    gpa.yy = zeros(n*m,nb_k);
    %A   = zeros(n*m,nb_k);
    %Sk = zeros(n*m,nb_k);
    h_x = h_x(:);
    h_z = h_z(:);
    k = k(:)';
    nb_k = length(k);

    % Index nodes for boundary condition (bc)
    ind_Rtc = sub2ind([n m],1,m);                  % index right top corner
    ind_Ltc = sub2ind([n m],1,1);                  % index left top corner
    ind_Rbc = sub2ind([n m],n,m);                  % index right bottom corner
    ind_Lbc = sub2ind([n m],n,1);                  % index left bottom corner
    %ind_T  = sub2ind([n m],ones(1,m-2),[2:m-1]);  % index top edge
    ind_B   = sub2ind([n m],n*ones(1,m-2),2:m-1);  % index bottom edge
    ind_R   = sub2ind([n m],2:n-1,m*ones(1,n-2));  % index right edge
    ind_L   = sub2ind([n m],2:n-1,ones(1,n-2));    % index left edge
    ind = 1:n*m; % index for all nodes
    ind_wbc = ind;
    ind_wbc ([ind_Ltc  ind_Lbc ind_L ind_R ind_Rtc ind_Rbc ind_B]) = []; % index without bc

    dz = repmat([h_z;h_z(end)],m,1);
    dx = repmat([h_x' h_x(end)],n,1);
    dx = dx(:);

    dx = dx(ind_wbc,ones(1,nb_k));
    dz = dz(ind_wbc,ones(1,nb_k));

    U = zeros((n+2)*(m+1),nb_k);
    g = U;

    for i=1:nb_k
        u1 = reshape(u(:,i),n,m);
        u1 = [u1(2,:);u1;u1(end,:)];
        u1 = [u1 u1(:,end)];
        U(:,i) = u1(:);

        G1 = reshape(G(:,i),n,m);
        G1 = [G1(2,:) ; G1 ; G1(end,:)];
        G1 = [G1 G1(:,end)];
        g(:,i) = G1(:); 
    end

    G = g;
    u = U;
    clear G1 u1 U g

    ind_wbc1 = ind_wbc;
    n = n+2;
    m = m+1;
    ind_wbc = 1:n*m; % index for all nodes
    ind_wbc = reshape(ind_wbc,n,m);
    ind_wbc = ind_wbc(2:end-2,2:end-2);
        ind_wbc = ind_wbc(:);

%% Calculation
        
    gpa.xx(ind_wbc1,:) =   ( ...
                             (u(ind_wbc+n,:)-u(ind_wbc,:)) .* (G(ind_wbc+n,:)-G(ind_wbc,:)) ...
                           + (u(ind_wbc+n+1,:)-u(ind_wbc+1,:)) .* (G(ind_wbc+n+1,:)-G(ind_wbc+1,:)) ...
                           + ((u(ind_wbc+n,:)+u(ind_wbc+n+1,:)-u(ind_wbc-n,:)-u(ind_wbc-n+1,:))/4).*((G(ind_wbc+n,:)+G(ind_wbc+n+1,:)-G(ind_wbc-n,:)-G(ind_wbc-n+1,:))/4) ...
                           + ((u(ind_wbc+2*n,:)+u(ind_wbc+2*n+1,:)-u(ind_wbc,:)-u(ind_wbc+1,:))/4).*((G(ind_wbc+2*n,:)+G(ind_wbc+2*n+1,:)-G(ind_wbc,:)-G(ind_wbc+1,:))/4) ...
                             ).*dz./dx;

    gpa.zz(ind_wbc1,:) =   ( ...
                             (u(ind_wbc+1,:)-u(ind_wbc,:)).*(G(ind_wbc+1,:)-G(ind_wbc,:)) ...
                           + (u(ind_wbc+n+1,:)-u(ind_wbc+n,:)).*(G(ind_wbc+n+1,:)-G(ind_wbc+n,:)) ...
                           + ((u(ind_wbc+1,:)+u(ind_wbc+n+1,:)-u(ind_wbc-1,:)-u(ind_wbc+n-1,:))/4).*((G(ind_wbc+1,:)+G(ind_wbc+n+1,:)-G(ind_wbc-1,:)-G(ind_wbc+n-1,:))/4) ...
                           + ((u(ind_wbc+2,:)+u(ind_wbc+n+2,:)-u(ind_wbc,:)-u(ind_wbc+n,:))/4).*((G(ind_wbc+2,:)+G(ind_wbc+n+2,:)-G(ind_wbc,:)-G(ind_wbc+n,:))/4) ...
                             ).*dx./dz;    

    gpa.xz(ind_wbc1,:) =   ( ...
                           + ((u(ind_wbc+n,:)-u(ind_wbc,:))).* ...
                                    ((G(ind_wbc+1,:)+G(ind_wbc+n+1,:)-G(ind_wbc-1,:)-G(ind_wbc+n-1,:))/4) ... 
                           + ((u(ind_wbc+n+1,:)-u(ind_wbc+1,:))).* ...
                                    ((G(ind_wbc+2,:)+G(ind_wbc+n+2,:)-G(ind_wbc,:)-G(ind_wbc+n,:))/4) ...
                           + (u(ind_wbc+1,:)-u(ind_wbc,:)).* ...
                                    ((G(ind_wbc+n,:)+G(ind_wbc+n+1,:)-G(ind_wbc-n,:)-G(ind_wbc-n+1,:))/4) ...                           
                           + (u(ind_wbc+n+1,:)-u(ind_wbc+1,:)).* ...
                                    ((G(ind_wbc+n,:)+G(ind_wbc+n+1,:)-G(ind_wbc-n,:)-G(ind_wbc-n+1,:))/4) ...
                           + ((G(ind_wbc+n,:)-G(ind_wbc,:))).* ...
                                    ((u(ind_wbc+1,:)+u(ind_wbc+n+1,:)-u(ind_wbc-1,:)-u(ind_wbc+n-1,:))/4) ... 
                           + ((G(ind_wbc+n+1,:)-G(ind_wbc+1,:))).* ...
                                    ((u(ind_wbc+2,:)+u(ind_wbc+n+2,:)-u(ind_wbc,:)-u(ind_wbc+n,:))/4) ...
                           + (G(ind_wbc+1,:)-G(ind_wbc,:)).* ...
                                    ((u(ind_wbc+n,:)+u(ind_wbc+n+1,:)-u(ind_wbc-n,:)-u(ind_wbc-n+1,:))/4) ...                           
                           + (G(ind_wbc+n+1,:)-G(ind_wbc+1,:)).* ...
                                    ((u(ind_wbc+n,:)+u(ind_wbc+n+1,:)-u(ind_wbc-n,:)-u(ind_wbc-n+1,:))/4) ...
                             );

    k1 = (k.^2)/4;
    gpa.yy(ind_wbc1,:) = 2*( G(ind_wbc,:).*u(ind_wbc,:) + G(ind_wbc+n,:).*u(ind_wbc+n,:) +...
                             G(ind_wbc+1,:).*u(ind_wbc+1,:) + G(ind_wbc+n+1,:).*u(ind_wbc+n+1,:) ...
                           ).*dx.*dz;
    gpa.yy = k1(ones(la,1),:).*gpa.yy;

    return

end