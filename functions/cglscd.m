function [x, error, iter] = cglscd(J, x, b, BETA, CTC, dxc, D, P, max_it, tol, appl_fct_reg)

%%========================================================================%
%                                                                         %
% This function calculate resistivity perturabation by cg                 % 
% it resolves :                                                           %
%  for model reg. (J'*D'*D*J+ BETA CTC) x= J'*D'*D (d-dobs)-BETA(CTC*dxc) %
%  for model pert. reg. (J'*D'*D*J+ BETA CTC) x= J'*D'*D (d-dobs)         %
%                                                                         %
%%========================================================================%
%                                                                         %
% In:                                                                     %
% -----------                                                             %
%   J : Jacobian matrix                                                   %
%   x : initial guess vector                                              % 
%   b : d-dobs                                                            %
%   BETA : regularization coefficient                                     %
%   CTC : C'*C where C is regularization matrix                           %
%   dxc : rho-rho_ref                                                     %
%   D : data weigthing matrix                                             %
%   P : constraints vector to fix resistivity ( 0 fix, 1 unfix)           %
%   max_it : maximum number of iteration                                  %
%   tol :  error tolerance                                                %
%   appl_fct_reg : applied regularization ( model or model perturbation)  %
%                                                                         %
% Out:                                                                    %
% -----------                                                             %
%   x        REAL solution vector                                         %
%   error    REAL error norm                                              %
%   iter     INTEGER number of iterations performed                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    iter = 0;

    if isempty(P)
        P = ones(size(x));
    end
 
    x = P.*x;
    % M : jacoby preconditioner 
    M = (sum((D*J).^2)'+BETA*diag(CTC));
    M(M==0) = eps;
    M = 1./M; 
    %
    zz = D*(b - J*x); % gradient initialisation 

    if strcmp(appl_fct_reg,'model') 

        b1 = (((D.^2)*b)'*J)' - BETA*CTC*(dxc); 
        bnrm2 = norm( b1 );

        if  ( bnrm2 == 0.0 )
            bnrm2 = 1.0; 
        end

        r= P.*(((zz' * D )*J)'- BETA*CTC * (x+dxc));

        error = norm( r ) / bnrm2; % error initialisation

        if ( error < tol ) 
            return
        end

        for iter = 1:max_it                       

            z  = M.* r; % modified gradient 
            rho = (r'*z);

            if ( iter > 1 )                   
                beta = rho / rho_1; % beta_k coefficient
                p = z + beta*p; % direction
            else
                p = z; % direction initialisation
            end

            q = D*(J*p);
            alpha = rho / (q'*q + (p'*BETA*CTC*p)); % step
            x = x + alpha * p; % descent                  
            zz = zz - alpha*q;
            r = P.*(((zz'*D)*J)' - BETA*CTC*(x+dxc));   % gradient                  
            error = norm( r ) / bnrm2; 

            if ( error <= tol ) 
                break 
            end % stop criterion

            rho_1 = rho;

        end


    elseif strcmp(appl_fct_reg,'model perturbation') 

        r= P.*((zz' * D )*J)'- BETA*CTC * x;
        b1=(((D.^2)*b)'*J)'; 
        bnrm2 = norm( b1 );
        error = norm( r ) / bnrm2; % error initialisation

        if ( error < tol ) 
            return 
        end

        for iter = 1:max_it                       

        z  = M.* r; % modified gradient
        rho = (r'*z);

        if ( iter > 1 )                      
            beta = rho / rho_1; % beta_k coefficient
            p = z + beta*p; % direction
        else
            p = z; % direction initialisation
        end

        q = D*J*p;
        alpha = rho / (q'*q + BETA *(p'*CTC*p)); % step
        x = x + alpha * p; % descent                  
        zz = zz - alpha*q;
        r = P.*((zz'*D)*J)'- BETA * CTC * x;  % gradient                  
        error = norm( r ) / bnrm2;            
        if ( error <= tol )
            break
        end % stop criterion

        rho_1 = rho;

        end

    end
end