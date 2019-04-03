function [k,wk,N] = find_k(N,k2)

%%========================================================================%
%                                                                         %
% 2.5D frequency process                                                  %
%                                                                         %
%%========================================================================%
%                                                                         %
% In:                                                                     %
% -----------                                                             %
%   N : number of frequencies considered                                  %
%   k2 : reference frequency                                              %
%                                                                         %
% Out:                                                                    %
% -----------                                                             %
%   k : number of frequencies considered                                  %
%   wk : reference frequency                                              %
%   N : number of frequencies considered                                  %
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

%% Code

    kk2=k2;
    k2=sqrt(k2);
    k1=0;

    [xk1,wk1]=Gausslp(N); 
    [xk2,wk2]=Gausslgp(N-3);

    k=[(((k2-k1)*xk1+k1+k2)/2).^2 xk2+kk2]';
    wk=[(k2-k1)*wk1/2 wk2];

end

    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function p=Lgndrp(N) % Legendre polynomial
    if N<=0
        p=1; % n*Ln(t)=(2n-1)t Ln-1(t)-(n-1)Ln-2(t)
    elseif N==1
        p=[1 0];
    else
        p=((2*N-1)*[Lgndrp(N-1) 0]-(N-1)*[0 0 Lgndrp(N-2)])/N;
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function p=Laguerp(N) %Legendre polynomial
    M=N:-1:0;
    p=factorial(N)*(((-1).^M)./factorial(M)).*((factorial(N))./(factorial(N-M).*factorial(M)));
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,w]=Gausslp(N)
    if N<0 
        error('Gauss-Legendre polynomial of negative order?');
    end 
    t=roots(Lgndrp(N))'; % make it a row vector
    A(1,:)= ones(1,N);  b(1)=2;
    for n=2:N
       A(n,:)=A(n-1,:).*t; 
       if mod(n,2)==0
           b(n)=0;
       else
           b(n)=2/n;
       end
    end
    w= b/A';
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,w]=Gausslgp(N)
    if N<0
        error('Gauss-Laguerre polynomial of negative order?'); 
    end 
    t=roots(Laguerp(N))'; % make it a row vector
    A(1,:)= ones(1,N);  b(1)=1;
    for n=2:N
       A(n,:)=A(n-1,:).*t; 
       b(n)=(n-1)*b(n-1); 
    end
    w= b/A';
end