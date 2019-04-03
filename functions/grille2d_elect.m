function [grille,h_x,h_z] =grille2d_elect(hx,hz,fact,nb_pad_bloc,nb_raff,nb_surr,grid_plot)

%%========================================================================%
%                                                                         %
%  % In:                                                                  %
% -----------                                                             %
%    hx : bloc thick in x direction of interesting region                 %
%    hz : bloc thick in z direction of interesting region                 %
%    fact : pad factor                                                    %
%    nb_pad_bloc : nb of padded blocs                                     %
%    nb_raff : refine the blocs in the intersting region by nb_aff(1)     %
%              nb_aff(2) in z direction                                   %
%    nb_surr : number of surrounding blocs                                %
%                                                                         %
%    such as:                                                             %
%                                                                         %
%     nb_pad_bloc  nb_surr  hx1  hx2  hx1 nb_surr  nb_pad_bloc            %
%    _____________________x0________________________________________      %
%    |     |      |surr b |    |    |    |surr b |       |       |        %
%    |_____|______|______ |____|____|____|_______|_______|_______|hz1     %
%    |     |      |surr b |    |    |    |surr b |       |       |        %
%    |_____|______|_______|____|____|____|_______|_______|_______|hz2     %
%    |     |      |surr b |    |    |    |surr b |       |       |        %
%    |_____|______|_______|____|____|____|_______|_______|_______|hz3     %
%    |     |      |surr   |surr|surr|surr| surr  |       |       |        %
%    |padded blocs|  b    | b  |  b |  b |   b   |padded blocs   |        %
%    |_____|______|_______|____|____|____|_______|_______|_______|        %
%    |     |      |       |padded blocs  |       |       |       |        %
%    |_____|______|_______|____|____|____|_______|_______|_______|        %
%    |     |      |       |    |    |    |       |       |       |        %
%    |_____|______|_______|____|____|____|_______|_______|_______|        %
%                                                                         %
%  % Out:                                                                 %
% -----------                                                             %
%    grille : grid ( ordered column by column                             %
%    h_x : mesh bloc thick   in x direction                               %
%    h_z : mesh bloc thick  in z direction                                %
%                                                                         %
%  % test:                                                                %
% -----------                                                             %
%    hz=ones(10,1);                                                       %
%    hx=ones(10,1);                                                       %
%    [X,Y] = meshgrid(0:9,0:9);                                           %
%    x=X(:);                                                              %
%    y=Y(:);grille=[];                                                    %
%    grille(:,1)=x;                                                       %
%    grille(:,2)=y;                                                       %
%    sigma=ones(10,10);                                                   %
%    [grille,h_x,h_z] = grille2d_elect(hx,hz,1.3,10,[1;1],0,1);           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2007 Abderrezak BOUCHEDDA
%=====oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo====%
%   contact:                   \\\////                                    %
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%========================================================================%

%% Calculation

    ind_x0 = nb_raff(1)*(nb_surr) + nb_pad_bloc + 1; % first electrode position in grid
    hx=hx(:);
    hz=hz(:);
    pwr = 1:nb_pad_bloc;
    factor = fact.^pwr;
    factor=factor(:);

    if isempty(nb_raff)
       nb_raff=[1;1];
    end

    hx1 = hx(1);
    hx2 = hx(end);
    hz1 = hz(end);
    ind_x0 = ind_x0-3;
    hx=[hx1*ones(nb_surr,1);
        hx;
        hx2*ones(nb_surr,1)];
    hz=[hz;hz1*ones(nb_surr,1)];

    hx=repmat(hx',nb_raff(1),1); 
        hx=hx(:)/nb_raff(1);
    hz=repmat(hz',nb_raff(2),1); 
    hz=hz(:)/nb_raff(2);

    h_x = [hx1*factor(end:-1:1);hx;hx2*factor];
    h_z = [hz;hz1*factor]; 
    hx = [0;h_x];
    hz = [0;h_z];

    hx=cumsum(hx);
    hz=cumsum(hz);

    [X,Z]=meshgrid(hx,hz);
    x0=X(1,ind_x0);
    grille(:,1)=X(:)-x0; % x=0 position of first electrode (on the left)
    grille(:,2)=Z(:);

    if grid_plot==1
        plot(grille(:,1),-grille(:,2),'s','MarkerFaceColor','g','MarkerSize',2);
        gridxy(hx-x0,-hz);
        axis([min(grille(:,1)) max(grille(:,1)) min(-grille(:,2)) max(-grille(:,2))])
    end
end
