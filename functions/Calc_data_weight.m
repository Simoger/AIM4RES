function [D,dtw_pre] = Calc_data_weight(param,e,dtw_pre)

%%========================================================================%
%                                                                         %
% Data weighting matrix D creation                                        %
%                                                                         %
%%========================================================================%
%                                                                         %
% In:                                                                     %
% -----------                                                             %
%   param.MEAS.Res : measured resistance                                  %
%   param.inv.dataweight.a : minimum resistance error (positive constant) % 
%   param.inv.dataweight.b : percentage error on measured resistance      %
%                            (value between 0 & 1)                        %
%                                                                         %
%   data noise model is defined by :                                      %                      
%     |dR| = a*|Res|+b;                                                   %
%   dR have zero mean and are normally distributed                        %
%   a and b can be obtained by reciprocity error measurements. for more   %
%   details see :                                                         %
%   Binley, A.M., A.Ramirez and W.Daily, 1995, Regularised Image Recons-  %
%   truction of Noisy Electrical resistance Tomography Data, In: Process  %
%   Tomography - 1995, by Beck, M.S. et al. (Eds.), Proc. Fourth Workshop %
%   of the European Concerted Action on Process Tomography, Bergen,       %
%   Norway, 401-410.                                                      %
%                                                                         %
%   e: data misfit ( e= D*(rho-rho_obs) ) at iteration j                  %
%                                                                         %
%   dtw_pre: previous data weighting vector (at iteration j-1)            %
%                                                                         %
%   param.inv.dataweight.type: 2 robust esimators (reweighted weight)     %
%                              'W-estimator' or 'L1'                      %
%                                                                         %
%   for 'L1' type see :                                                   %
%   ----------------------------                                          %
%   Claerbout and Muir, 1973 J.F. Claerbout and F. Muir, Robust modeling  % 
%   with erratic data, Geophysics 38 (1973), pp. 826–844.                 %
%   
%   for 'W-estimator' type see :                                          %
%   -------------------------------------                                 %
%   Daily, W., A. Ramirez, A. Binley and D. LaBrecque, Electrical         %
%   Resistance Tomography - Theory and Practice, Near Surface Geophysics, %
%   Investigationsin Geophysics No. 13, Society of Exploration            %
%   Geophysicists, 525-550.                                               %
%                                                                         %
% Out:                                                                    %
% -----------                                                             %  
%   D: data weighting matrix (D is diagonal matrix(uncorelated error))    %
%                                                                         %
%=========================================================================%
% Copyright (C) 2008 Abderrezak BOUCHEDDA                                 %
%=====oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo====%
%          contact: ---------->\\\////                                    %
%                               |_ _|                                     %
%                               (@ @)                                     %
%               **********oooO***(_)***Oooo**********                     %
%               * -----> Abderrezak BOUCHEDDA<----- *                     %
%               *      bouchedda@geo.polymtl.ca     *                     %
%               *  Ecole polytechnique de Montreal  *                     %
%               *  depart. de géophysique appliquée *                     %
%               *  http://geo.polymtl.ca            *                     %
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
%%========================================================================%

%% Calculation

    if nargin < 2
        e = []; 
        dtw_pre = [];
    elseif nargin < 3
        dtw_pre = [];
    end

    if isempty (dtw_pre) || strcmp(param.inv.dataweight.type,'const. weight')

         %get the absolute error for each mesured resistance datum, plus minimum error b
         dtw = (param.inv.dataweight.a * abs(param.MEAS.Res) + param.inv.dataweight.b);

          if strcmp(param.inv.invparam,'log resistivity') 

              dtw = log(1+(dtw./abs(param.MEAS.Res)));
              dtw = 1./dtw;

          elseif strcmp(param.inv.invparam,'resistivity') 

              dtw = abs(param.K).*dtw;
              dtw = 1./dtw;
          else    
              dtw = 1./dtw;
          end


          %Normalize things to 1
          %dtw = dtw./max(dtw);
    end



    if strcmp(param.inv.dataweight.type,'L1') && ~isempty(e)
       % L1 robust estimator
       wd = sum(e.^2)./(abs(e)*sum(abs(e)));
       wd(wd > 1) = 1; wd (wd < .05)=.05; 
       dtw = wd.*dtw_pre;
       dtw_pre = dtw;

    elseif strcmp(param.inv.dataweight.type,'W-estimator') && ~isempty(dtw_pre)
        % W-estimator 
        dtw = dtw_pre; 
        dtw_int = dtw_pre./(sqrt(abs(e.*dtw_pre)));
        dtw_int1 = dtw_int*(sum(abs(e.*dtw_pre))/sum(abs(e.*dtw_int)));
        dtw(dtw_int1 < dtw_pre) = dtw_int1(dtw_int1 < dtw_pre);

        dtw_pre = dtw;

    elseif isempty(e) 
         dtw_pre = dtw;
    else
        errordlg('Calc_data_weight function : data weight type is not defined or not suitable',...
            'data weight Error');
    end


    % data weighting matrix
    D = spdiags(dtw, 0, length(dtw), length(dtw));
    
end
