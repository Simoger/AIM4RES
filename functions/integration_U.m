function u = integration_U(U,wk,k,nbk1)

%%========================================================================%
%                                                                         %
% Integration of potential u, inverse Fourier transform                   %
%                                                                         %
% ======================================================================= %
%                                                                         %
% In:                                                                     %
% -----------                                                             %
%   U : potential in wave number domain (matrix or vector)                %
%                                                                         %
%     for one current electrode                                           %
%     U = [U1_k1 U1_k2 ..... U1_kn]                                       %
%         [U2_k1 U2_k2 ..... U2_kn]                                       %
%         [.      .          .    ]                                       %
%         [.      .          .    ]                                       %
%         [Um_k1 U2_k2 ..... Um_kn]                                       %
%                                                                         %
%     for nb current electrodes                                           %
%     U = [U11_k1 U12_k1 ...U1nb_k1 U11_k2 U12_k2 ..... U1nb_kn]          %
%         [U21_k1 U22_k1 ...U2nb_k1 U21_k2 U22_k2 ..... U2nb_kn]          %
%         [.      .          .                                 ]          %
%         [.      .          .                                 ]          %
%         [Um1_k1 Um2_k1 ...Umnb_k1 Um1_k2 Um2_k2 ..... Umnb_kn]          %
%                                                                         %
%   k : wave number                                                       %
%       k = [k1 k2 .... kn]                                               %
%                                                                         %
% Out:                                                                    %
% -----------                                                             %
%   u : integrated potential into the spatial domain                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2007 Abderrezak BOUCHEDDA                                 %
%=====oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo====%
%          contact: ---------->\\\////                                    %
%                               |_ _|                                     %
%                               (@ @)                                     %
%               **********oooO***(_)***Oooo**********                     %
%               * -----> Abderrezak BOUCHEDDA<----- *                     %
%               *  Abderrezak.Bouchedda@ete.inrs.ca *                     %
%               *  INRS-ETE                         *                     %
%               *  http://www.ete.inrs.ca/ete       *                     %
%               *************************************                     %
%                             |_______|                                   %
%                              |__|__|                                    %
%                               () ()                                     %
%                              ooO Ooo                                    %
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
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%========================================================================%

%% Calculation

    [h,l]=size(U);
    k1=sqrt(k(1:nbk1));
    k1=k1(:)'; 

    % potential integration using gauss quadrature
    u = (4/pi) * ( wk(1:nbk1) * (U(:,1:nbk1).*k1(ones(h,1),:))' ) ...
        + (2/pi) * ( wk(nbk1+1:end) * U(:,1+nbk1:end)' );
    
end
    
    
