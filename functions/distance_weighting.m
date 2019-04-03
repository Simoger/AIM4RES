function Q=distance_weighting(param,beta)

%%========================================================================%
%                                                                         %
% Calculation of distance weigthing values                                %
%                                                                         %
%%========================================================================%
%                                                                         %
% In:                                                                     %
% -----------                                                             %
%   param: parameters structure                                           %
%   beta: weigthing coefficient                                           %
%                                                                         %
% Out:                                                                    %
% -----------                                                             %
%   Q = weigthing values                                                  %
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

    % Li and Oldenburg, 2000,Joint inversion of surface and three-component
    % borehole magnetic data, Geophysics.

    X = param.grille(:,1);
    X = reshape(X,param.nb_ligne,param.nb_col);
    X = X';
    X = X(:);

    Z = param.grille(:,2);
    Z = reshape(Z,param.nb_ligne,param.nb_col);
    Z = Z';
    Z = Z(:);
    
    h_x = [param.h_x ; param.h_x(end)];
    h_z = [param.h_z ; param.h_z(end)];

    h_xr = repmat(h_x(:)',length(h_z),1);
    h_zr = repmat(h_z(:),1,length(h_x));
    dV = (h_xr.*h_zr)';
        dV = dV(:);

    R0=min(dV)^(1/3);
    Q=0;
    pos = [param.pos.C ; param.pos.P];
%         pos = pos(pos(:,2)==0 , :); 
%         pos(pos(:,2)==0 , :) = []; 
    n=length(pos);

    for i=1:n
        R = sqrt((X-pos(i,1)).^2 + (Z-pos(i,2)).^2);

        R=(R+R0).^3;
        R=(dV./R).^2;
        Q=Q+R;
    end

    Q=Q.^(beta/4);
%     Q(Q<1e-3) = 1e-3;
    
    Q = reshape(Q,param.nb_col,param.nb_ligne)';
    Q(:,[1:param.nb_pad_bloc+param.nb_surr-3 end-(param.nb_pad_bloc+param.nb_surr)+3:end]) = 1;
    Q(end-(param.nb_pad_bloc+param.nb_surr)+3:end,:) = 1;
end