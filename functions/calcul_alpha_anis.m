function alpha = calcul_alpha_anis(a,k,rho,sigma,pos,normale)

%%========================================================================%
%                                                                         %
% calculation of nu for mixed boundary condition                          %
%                                                                         %
%%========================================================================%
%                                                                         %
% In:                                                                     %
% -------                                                                 %
%   r: distance between current source frontiers points                   %
%   k: wave number vector                                                 %
%   rho: resistivity values                                               %
%   sigma: conductivity values                                            %
%   pos: positions vector                                                 %
%   normale: normale direction vector                                     %
%                                                                         %
% Out:                                                                    %
% -----------                                                             %
%   nu: mixed boundary condition coefficient                              %
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

    eps = 10^-10;
    nk = length(k);
    nr = length(a);
    alpha = zeros(nr,nk);
    normale = normale(:);

    x = pos(:,1);
    z = pos(:,2);

    
%     scal = 2*[x z]; % [sigma . grad(a)] dot product
    scal = [sigma.xx.*(2.*x.*rho.xx + 2.*z.*rho.xz) + sigma.xz.*(2.*x.*rho.xz + 2.*z.*rho.zz), ...
            sigma.xz.*(2.*x.*rho.xx + 2.*z.*rho.xz) + sigma.zz.*(2.*x.*rho.xz + 2.*z.*rho.zz)];
    scal_norm = scal*normale;
    scal_norm = (1./rho.yy).*scal_norm;

    for j=1:nk
    %     K0 = besselk(0,(k(j)*sqrt(a)),0);
    %     K1 = besselk(1,(k(j)*sqrt(a)),0);
        K0 = besselk(0,(k(j)*sqrt(a)));
        K1 = besselk(1,(k(j)*sqrt(a)));
        alpha(:,j) = k(j)./(2*sqrt(a)).*(K1./(K0+eps)).*scal_norm;
    end

end
    
    
    
    
    