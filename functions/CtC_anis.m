function[CTC,Cx,Cz] = CtC_anis(wt,param,nb_composantes)

%%========================================================================%
%                                                                         %
% CtC calculation - the model regularization matrix                       %
%                                                                         %
%%========================================================================%
%                                                                         %
% In:                                                                     %
% -----------                                                             %
%   wt: weigths applied on the regularization matrix                      %
%   param.param.inv.alx: smoothing parameter in x direction               %
%   param.param.inv.alz: smoothing parameter in z direction               %
%   param.h_x: x step                                                     %
%   param.h_z: z step                                                     %
%                                                                         %
% Out:                                                                    %
% -----------                                                             %
%   CTC: model weighting matrix                                           %
%   Cx: model weighting matrix in x direction                             %
%   Cx: model weighting matrix in z direction                             %
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

    param.h_x = [param.h_x ; param.h_x(end)];
    param.h_z = [param.h_z ; param.h_z(end)];

    %smoothing parammeters
    alx = param.inv.alx;
    alz = param.inv.alz;

    nx = param.nb_col;
    nz = param.nb_ligne;

    %Smallness parameter
    als = param.inv.als;
    wts = wt(:);

    % derivative matrix
    [Cx,Cz] = deriv2D(param.h_x,param.h_z,param.flag.reg_fct);

    if ~isempty(param.const.WGx)
        Cx = spdiags(param.const.WGx,0,nx*nz,nx*nz)*Cx;
    end

    if ~isempty(param.const.WGz)
        Cz = spdiags(param.const.WGz,0,nx*nz,nx*nz)*Cz;
    end

    if isfield(param.const,'Wti')
        Ci = spdiags(param.const.Wti,0,nx*nz,nx*nz);
    else
        Ci = speye(nx*nz,nx*nz);
    end

%% Create a weighted smallness term

    % Weights certain points more than others

    Wts_1 = spdiags(wts(1:nx*nz), 0, nx*nz, nx*nz); 
    Wts_2 = spdiags(wts(nx*nz+1:2*nx*nz), 0, nx*nz, nx*nz); 

    if nb_composantes == 3
        Wts_3 = spdiags(wts(2*nx*nz+1:end), 0, nx*nz, nx*nz);
    end

    % Assemble the Anisotropic derivative operator
    Gs = [alx*Cx ; alz*Cz];

    % Assemble the 2d weighting matrix
    CTC = (Gs'*Gs + als * Ci);
    % if flag == 0
    %     CTC = (Gs'*Gs + als * Ci);
    % else
    %     CTC = (Gs'*Gs + als * speye(nx*nz,nx*nz));
    % end

    if nb_composantes == 2
        CTC = sparse( blkdiag(Wts_1',Wts_2') * blkdiag(CTC,CTC) * blkdiag(Wts_1,Wts_2) );
    elseif nb_composantes == 3
        CTC = sparse( blkdiag(Wts_1',Wts_2',Wts_3') * blkdiag(CTC,CTC,CTC) * blkdiag(Wts_1,Wts_2,Wts_3) );
    end
end