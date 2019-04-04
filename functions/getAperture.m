function aperture = getAperture(quadril_points, linearArea_coef)

%%========================================================================%
%                                                                         %
% Gauss Newton inversion                                                  %
%                                                                         %
%%========================================================================%
%                                                                         %
% In :                                                                    %
% -----------                                                             %
%   quadril_points: 4 points forming the quadrilateral of which we want   %
%                   to determine the surface.                             %
%                   (row vector) -> 4 points, 8 coordinates               %
%                   [(xC1,zC1), (xC2,zC2), (xP1,zP1), (xP2,zP2)]          %
%   linearArea_coef: coefficient applied on the inline quadrilateral      % 
%                    to make its length comparable to a surface           %
%                                                                         %
% Out :                                                                   %
% -----------                                                             %
%   aperture: area of the quadrilateral                                   %
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

%% Calculation

    quadril_points = reshape(quadril_points,2,4)'; % 4x2 matrix : c1 : x coord ; c2 : z coord

    % Lengths
    
%     if length(unique(quadril_points(:,1))) == 1 || length(unique(quadril_points(:,2))) == 1 % Test colinearity vectors better
    % Inline electrodes: no proper area calculable
    if length(unique(quadril_points(:,1))) == 1
        Pmax = max(quadril_points(:,2));
        Pmin = min(quadril_points(:,2));
        aperture = (Pmax - Pmin)*linearArea_coef;
        return
    elseif length(unique(quadril_points(:,2))) == 1
        Pmax = max(quadril_points(:,1));
        Pmin = min(quadril_points(:,1));
        aperture = (Pmax - Pmin)*linearArea_coef;
        return
    end
    
    % Areas

    indP = convhull(quadril_points);    

    quadril_points = quadril_points(indP(1:end-1),:);

    x = quadril_points(:,1);
    y = quadril_points(:,2);

    % area : http://www.mathopenref.com/coordpolygonarea.html

    x2 = circshift(x,-1);
    y2 = circshift(y,-1);

    aperture = abs(sum(x.*y2 - y.*x2)/2);   
    
    
end