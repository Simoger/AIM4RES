function [pos,ind,L,sign]= tri_electr_pos(meas_posC1 , meas_posC2 , elect_bh1 , elect_bh2 , elect_sur , meas_posP1 , meas_posP2 , grille)

%%========================================================================%
%                                                                         %
% Electrode sorting function                                              % 
%                                                                         %
%%========================================================================%
%                                                                         %
% Sorts electrodes in order to model with a pole-pole, which allows for   %
% potential recombination in order to obtain the potential of a dipolar   %
% device                                                                  %
%                                                                         %
% In:                                                                     %
% -----------                                                             %
%   meas_posC1: C1 measurement position                                   %
%   meas_posC2: C2 measurement position                                   %
%   elect_bh1 : borehole 1 electrodes position                            %
%   elect_bh2 : borehole 2 electrodes position                            %
%   elect_sur : surface electrodes position                               %
%   grille    : grid nodes position                                       %
%                                                                         %
% Out:                                                                    %
% -----------                                                             %
%   pos.C   : current electrodes position for forward problem             %
%   pos.P   : potential electrodes position used as current electrodes    %
%             for the calculation of reciprocal potential                 %
%             (sensit. calculation)                                       %
%   ind.C1C2: index of current electrodes (pos.C) in meas_posC1 and       %
%             meas_posC2                                                  %
%   L       : Contains the number of measures identified in 'ind' for     %
%             each 'pos' position                                         %
%   sign    : +1 for a C1 electrode and -1 for a C2 electrode             %
%   flag    : 1 for borehole electrodes and 2 for surface electrodes      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2007 Abderrezak BOUCHEDDA                                 %
%=====oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo====%
%          contact: ---------->\\\////                                    %
%                               |_ _|                                     %
%                               (@ @)                                     %
%               **********oooO***(_)***Oooo**********                     %
%               * -----> Abderrezak BOUCHEDDA<----- *                     %
%               *      bouchedda@geo.polymtl.ca     *                     %
%               *  Ecole polytechnique de Montreal  *                     %
%               *  depart. de geophysique appliquee *                     %
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

%% Initialisation

    eps=.0001;
    [nh1,~]=size(elect_bh1);
    [nh2,~]=size(elect_bh2);
    [nsr,~]=size(elect_sur);

    ind=[]; 
    L=[]; 
    sign=[]; 
    I=[];
    pos1=[]; 
    flag=[];
    
%% Calculation

    % find first borehole electrodes that are used in the measurements
    for i=1:nh1
        ind11 = find(meas_posC1(:,1) == elect_bh1(i,1) & meas_posC1(:,2) == elect_bh1(i,2));
        l1=length(ind11);
        ind12 = find(meas_posC2(:,1) == elect_bh1(i,1) & meas_posC2(:,2) == elect_bh1(i,2));
        l2=length(ind12);     

        if (l1==0 && l2==0)
            I=[I,i]; % indices to throw electrodes that do not contribute to the measurements
        end

        L11=length(ind11); sign1= 1*ones(size(ind11)); % times that electrode elect_bh1(i,1) is used in meas.PosC1
        L12=length(ind12); sign2=-1*ones(size(ind12)); % times that electrode elect_bh1(i,1) is used in meas.PosC2

        L1=L11+L12; 
        if L1==0
            L1=[]; 
        end

        ind=[ind;ind11;ind12];
        L=[L;L1];
        sign=[sign;sign1;sign2];
    end

    pos1=elect_bh1(I,:);flagg=ones(length(pos1),1);
    elect_bh1(I,:)=[]; I=[];
    flag=ones(length(elect_bh1),1);
    
    % find second borehole electrodes that are used in the measurement (C1 ou C2)
    for i=1:nh2

        ind11 = find(meas_posC1(:,1) == elect_bh2(i,1) & meas_posC1(:,2) == elect_bh2(i,2));
        l1=length(ind11);
        ind12 = find(meas_posC2(:,1) == elect_bh2(i,1) & meas_posC2(:,2) == elect_bh2(i,2));
        l2=length(ind12);     

        if (l1==0 && l2==0)
            I=[I,i];
        end

        L11=length(ind11); sign1=1*ones(size(ind11));
        L12=length(ind12); sign2=-1*ones(size(ind12));

        L1=L11+L12; 
        if L1==0
            L1=[]; 
        end

        ind=[ind;ind11;ind12];
        L=[L;L1];
        sign=[sign;sign1;sign2];

    end
    pos1=[pos1;elect_bh2(I,:)];
    flagg=[flagg;ones(length(elect_bh2(I,1)),1)];
    elect_bh2(I,:)=[]; 
    I=[];
    flag=[flag;ones(length(elect_bh2),1)];
    
    % find surface electrodes that are used in the measurement (C1 ou C2)
    for i=1:nsr

        ind11 = find(elect_sur(i,1) == meas_posC1(:,1) & elect_sur(i,2) == meas_posC1(:,2));
        l1=length(ind11);
        ind12 = find(elect_sur(i,1) == meas_posC2(:,1) & elect_sur(i,2) == meas_posC2(:,2));
        l2=length(ind12);

        if (l1==0 && l2==0)
            I=[I,i];
        end

        L11=length(ind11); 
        sign1=1*ones(size(ind11));
        L12=length(ind12);
        sign2=-1*ones(size(ind12));

        L1=L11+L12; 
        if L1==0
            L1=[]; 
        end

        ind=[ind;ind11;ind12];
        L=[L;L1];
        sign=[sign;sign1;sign2];

    end

    pos1=[pos1;elect_sur(I,:)];
    flagg=[flagg;2*ones(length(elect_sur(I,1)),1)];
    elect_sur(I,:)=[];
    I=[];
    flag=[flag;2*ones(length(elect_sur),1)];
    pos=[elect_bh1;elect_bh2;elect_sur];

    warning off
    ind_tmp = ind;
    clear ind
    ind.C1C2 = ind_tmp;

    % Search for nodes corresponding to potential electrodes (P1 & P2)
    n=length(ind.C1C2);
    ind_P1=[]; 
    ind_P2=[];

    if isnan(meas_posP2)
        for i=1:n
            k=ind.C1C2(i);
            ind11 = find( abs(grille(:,1)- meas_posP1(k,1)) < eps & abs(grille(:,2)- meas_posP1(k,2))< eps);    
            ind_P1=[ind_P1;ind11];
            ind_P2=[];
        end

    else
        for i=1:n
            k=ind.C1C2(i);
            ind11 = find( abs(grille(:,1) - meas_posP1(k,1)) < eps & abs(grille(:,2) - meas_posP1(k,2)) < eps);
            ind12 = find( abs(grille(:,1) - meas_posP2(k,1)) < eps & abs(grille(:,2) - meas_posP2(k,2)) < eps);    

            ind_P1=[ind_P1;ind11];
            ind_P2=[ind_P2;ind12];
        end
    end

    ind.P1 = ind_P1;
    if isempty(ind_P2)
        ind.P2=ones(size(ind_P1)); % P2 dummy values for pole-pole 
    else
        ind.P2 = ind_P2;
    end

    % Search for P1 & P2 electrodes position for sensitivity calculation
    n_ele=size(pos1);
    n_ele=n_ele(1,1);

    if n_ele ~=0 
        for i=1:n_ele
            ind11 = find(meas_posP1(:,1) == pos1(i,1) & meas_posP1(:,2) == pos1(i,2));l1=length(ind11);
            ind12 = find(meas_posP2(:,1) == pos1(i,1) & meas_posP2(:,2) == pos1(i,2));l2=length(ind12);    
            if (l1==0 && l2==0)
                I=[I,i];
            end
        end
        pos1(I,:) = [];
        flagg(I) = [];   
    end

    nb_meas=size(meas_posC1,1);
    position = [pos;pos1];

    pos_tmp = pos; 
    clear pos
    pos.C = pos_tmp;
    pos.P = pos1;
    pos.flag = [flag;flagg];
    warning on

    % Take electrodes indices in position touse them in order to rebuild
    % the data (for the sensitivity calculation)
    if isnan(meas_posP2(1,1)) && isnan(meas_posC2(1,1)) % pole-pole
        ind_meas = zeros(nb_meas,2);
        for i=1:nb_meas
            ind1 = find(position(:,1) == meas_posC1(i,1)  & position(:,2)== meas_posC1(i,2));
            ind3 = find(position(:,1) == meas_posP1(i,1)  & position(:,2)== meas_posP1(i,2));
            ind_meas(i,:) = [ind1 ind3];
        end

    elseif isnan(meas_posC2(1,1)) && ~isnan(meas_posP2(1,1)) % pole-dipole
        ind_meas=zeros(nb_meas,3);
        for i=1:nb_meas
            ind1 = find(position(:,1) == meas_posC1(i,1) & position(:,2) == meas_posC1(i,2));   
            ind3 = find(position(:,1) == meas_posP1(i,1) & position(:,2) == meas_posP1(i,2));
            ind4 = find(position(:,1) == meas_posP2(i,1) & position(:,2) == meas_posP2(i,2));   
            ind_meas(i,:)=[ind1 ind3 ind4];
        end

    else % dipole-dipole (default)
        ind_meas=zeros(nb_meas,4);
        for i=1:nb_meas
            ind1 = find(position(:,1) == meas_posC1(i,1) & position(:,2) == meas_posC1(i,2));
            ind2 = find(position(:,1) == meas_posC2(i,1) & position(:,2) == meas_posC2(i,2));   
            ind3 = find(position(:,1) == meas_posP1(i,1) & position(:,2) == meas_posP1(i,2));
            ind4 = find(position(:,1) == meas_posP2(i,1) & position(:,2) == meas_posP2(i,2));   
            ind_meas(i,:) = [ind1 ind2 ind3 ind4];
        end

    end

    ind.meas = ind_meas;

end