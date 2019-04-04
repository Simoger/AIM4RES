function K = geometric_factor(XYZ,C1C2_config)

%%========================================================================%
%                                                                         %
% Geometric factor calculation                                            %
%                                                                         %
%%========================================================================%
%                                                                         %
% In:                                                                     %
% -----------                                                             %
%   XYZ: electrodes position structure                                    %
%   C1C2_config:                                                          %
%     C1C2_config==1 (C1,C2) borehole electrodes                          %
%     C1C2_config==2 (C1,C2) surface electrodes                           %
%     C1C2_config==3 (C1) borehole electrodes ; (C2) surface electrodes   %
%     C1C2_config==4 (C2) borehole electrodes ; (C2) surface electrodes   %
%                                                                         %
% Out:                                                                    %
% -----------                                                             %
%   K: vector of the geometric factors corresponding to each quadrupole   %
%      (Guo et al. 2014)                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2009 Abderrezak BOUCHEDDA
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

    ind1 = find(C1C2_config == 1);
    ind2 = find(C1C2_config == 2);
    ind3 = find(C1C2_config == 3);
    ind4 = find(C1C2_config == 4);

    K = zeros(size(C1C2_config));
    

    if ~isempty(ind1)
        if isnan(XYZ.MEAS.C2(1,1)) && isnan(XYZ.MEAS.P2(1,1)) % pole-pole
            r2=0; 
            r3=0; 
            r4=0; 
            r21=0;
            r31=0;
            r41=0;  
        elseif isnan(XYZ.MEAS.C2(1,1)) && ~isnan(XYZ.MEAS.P2(1,1)) % pole-dipole
            r2=0; 
            r4=0; 
            r21=0;
            r41=0;
            r3=1./((XYZ.MEAS.C1(ind1,1)-XYZ.MEAS.P2(ind1,1)).^2 +(XYZ.MEAS.C1(ind1,2)-XYZ.MEAS.P2(ind1,2)).^2).^.5; % C1-P2
            r31=1./((XYZ.MEAS.C1(ind1,1)-XYZ.MEAS.P2(ind1,1)).^2 +(-XYZ.MEAS.C1(ind1,2)-XYZ.MEAS.P2(ind1,2)).^2).^.5; %C1 mirror -P2
        elseif ~isnan(XYZ.MEAS.C2(1,1)) && ~isnan(XYZ.MEAS.P2(1,1)) % dipole-dipole     
            r2=1./((XYZ.MEAS.C2(ind1,1)-XYZ.MEAS.P1(ind1,1)).^2 +(XYZ.MEAS.C2(ind1,2)-XYZ.MEAS.P1(ind1,2)).^2).^.5; % C2-P1
            r3=1./((XYZ.MEAS.C1(ind1,1)-XYZ.MEAS.P2(ind1,1)).^2 +(XYZ.MEAS.C1(ind1,2)-XYZ.MEAS.P2(ind1,2)).^2).^.5; % C1-P2
            r4=1./((XYZ.MEAS.C2(ind1,1)-XYZ.MEAS.P2(ind1,1)).^2 +(XYZ.MEAS.C2(ind1,2)-XYZ.MEAS.P2(ind1,2)).^2).^.5; % C2-P2
            r21=1./((XYZ.MEAS.C2(ind1,1)-XYZ.MEAS.P1(ind1,1)).^2 +(-XYZ.MEAS.C2(ind1,2)-XYZ.MEAS.P1(ind1,2)).^2).^.5; %C2 mirror -P1
            r31=1./((XYZ.MEAS.C1(ind1,1)-XYZ.MEAS.P2(ind1,1)).^2 +(-XYZ.MEAS.C1(ind1,2)-XYZ.MEAS.P2(ind1,2)).^2).^.5; %C1 mirror -P2
            r41=1./((XYZ.MEAS.C2(ind1,1)-XYZ.MEAS.P2(ind1,1)).^2 +(-XYZ.MEAS.C2(ind1,2)-XYZ.MEAS.P2(ind1,2)).^2).^.5; %C2 mirror -P2 
        else 
            errordlg('geometric factor function : array type is not defined ');
        end

        r1=1./((XYZ.MEAS.C1(ind1,1)-XYZ.MEAS.P1(ind1,1)).^2 +(XYZ.MEAS.C1(ind1,2)-XYZ.MEAS.P1(ind1,2)).^2).^.5; % C1-P1
        r11=1./((XYZ.MEAS.C1(ind1,1)-XYZ.MEAS.P1(ind1,1)).^2 +(-XYZ.MEAS.C1(ind1,2)-XYZ.MEAS.P1(ind1,2)).^2).^.5; % C1 mirror -P1

        K(ind1)=(4*pi)./(r1+r11-r2-r21-r3-r31+r4+r41);
    end

    if ~isempty(ind2)
        if isnan(XYZ.MEAS.C2(1,1)) && isnan(XYZ.MEAS.P2(1,1)) % pole-pole
            r2=0; 
            r3=0; 
            r4=0;
        elseif isnan(XYZ.MEAS.C2(1,1)) && ~isnan(XYZ.MEAS.P2(1,1)) % pole-dipole
            r2=0; 
            r4=0; 
            r3=1./((XYZ.MEAS.C1(ind2,1)-XYZ.MEAS.P2(ind2,1)).^2 +(XYZ.MEAS.C1(ind2,2)-XYZ.MEAS.P2(ind2,2)).^2).^.5; % C1-P2
        elseif ~isnan(XYZ.MEAS.C2(1,1)) && ~isnan(XYZ.MEAS.P2(1,1)) % dipole-dipole     
            r2=1./((XYZ.MEAS.C2(ind2,1)-XYZ.MEAS.P1(ind2,1)).^2 +(XYZ.MEAS.C2(ind2,2)-XYZ.MEAS.P1(ind2,2)).^2).^.5; % C2-P1
            r3=1./((XYZ.MEAS.C1(ind2,1)-XYZ.MEAS.P2(ind2,1)).^2 +(XYZ.MEAS.C1(ind2,2)-XYZ.MEAS.P2(ind2,2)).^2).^.5; % C1-P2
            r4=1./((XYZ.MEAS.C2(ind2,1)-XYZ.MEAS.P2(ind2,1)).^2 +(XYZ.MEAS.C2(ind2,2)-XYZ.MEAS.P2(ind2,2)).^2).^.5; % C2-P2
        else 
            errordlg('geometric_factor function : array type is not defined ');
        end

        r1=1./((XYZ.MEAS.C1(ind2,1)-XYZ.MEAS.P1(ind2,1)).^2 +(XYZ.MEAS.C1(ind2,2)-XYZ.MEAS.P1(ind2,2)).^2).^.5; % C1-P1

        K(ind2)=(2*pi)./(r1-r2-r3+r4);    
    end

    if ~isempty(ind3)
        if ~isnan(XYZ.MEAS.C2(1,1)) && ~isnan(XYZ.MEAS.P2(1,1)) % dipole-dipole     
            r2=1./((XYZ.MEAS.C2(ind3,1)-XYZ.MEAS.P1(ind3,1)).^2 +(XYZ.MEAS.C2(ind3,2)-XYZ.MEAS.P1(ind3,2)).^2).^.5; % C2-P1
            r3=1./((XYZ.MEAS.C1(ind3,1)-XYZ.MEAS.P2(ind3,1)).^2 +(XYZ.MEAS.C1(ind3,2)-XYZ.MEAS.P2(ind3,2)).^2).^.5; % C1-P2
            r4=1./((XYZ.MEAS.C2(ind3,1)-XYZ.MEAS.P2(ind3,1)).^2 +(XYZ.MEAS.C2(ind3,2)-XYZ.MEAS.P2(ind3,2)).^2).^.5; % C2-P2
            r21=0; %C2 mirror -P1
            r31=1./((XYZ.MEAS.C1(ind3,1)-XYZ.MEAS.P2(ind3,1)).^2 +(-XYZ.MEAS.C1(ind3,2)-XYZ.MEAS.P2(ind3,2)).^2).^.5; %C1 mirror -P2
            r41=0; %C2 mirror -P2 
        else 
            errordlg('geometric factor function : array type is not defined ');
        end

        r1=1./((XYZ.MEAS.C1(ind3,1)-XYZ.MEAS.P1(ind3,1)).^2 +(XYZ.MEAS.C1(ind3,2)-XYZ.MEAS.P1(ind3,2)).^2).^.5; % C1-P1
        r11=1./((XYZ.MEAS.C1(ind3,1)-XYZ.MEAS.P1(ind3,1)).^2 +(-XYZ.MEAS.C1(ind3,2)-XYZ.MEAS.P1(ind3,2)).^2).^.5; % C1 mirror -P1

        K(ind3)=(2*pi)./(.5*(r1+r11)-r2-r21-.5*(r3+r31)+r4+r41);
    end

    if ~isempty(ind4)

        if ~isnan(XYZ.MEAS.C2(1,1)) && ~isnan(XYZ.MEAS.P2(1,1)) % dipole-dipole     
            r2=1./((XYZ.MEAS.C2(ind4,1)-XYZ.MEAS.P1(ind4,1)).^2 +(XYZ.MEAS.C2(ind4,2)-XYZ.MEAS.P1(ind4,2)).^2).^.5; % C2-P1
            r3=1./((XYZ.MEAS.C1(ind4,1)-XYZ.MEAS.P2(ind4,1)).^2 +(XYZ.MEAS.C1(ind4,2)-XYZ.MEAS.P2(ind4,2)).^2).^.5; % C1-P2
            r4=1./((XYZ.MEAS.C2(ind4,1)-XYZ.MEAS.P2(ind4,1)).^2 +(XYZ.MEAS.C2(ind4,2)-XYZ.MEAS.P2(ind4,2)).^2).^.5; % C2-P2

            r21=1./((XYZ.MEAS.C2(ind4,1)-XYZ.MEAS.P1(ind4,1)).^2 +(-XYZ.MEAS.C2(ind4,2)-XYZ.MEAS.P1(ind4,2)).^2).^.5;%C2 mirror -P1
            r31=0; %C1 mirror -P2
            r41=1./((XYZ.MEAS.C2(ind4,1)-XYZ.MEAS.P2(ind4,1)).^2 +(-XYZ.MEAS.C2(ind4,2)-XYZ.MEAS.P2(ind4,2)).^2).^.5; %C2 mirror -P2 
        else 
            errordlg('geometric factor function : array type is not defined ');
        end

        r1=1./((XYZ.MEAS.C1(ind4,1)-XYZ.MEAS.P1(ind4,1)).^2 +(XYZ.MEAS.C1(ind4,2)-XYZ.MEAS.P1(ind4,2)).^2).^.5; % C1-P1
        r11=0; % C1 mirror -P1

        K(ind4)=(2*pi)./((r1+r11)-.5*(r2+r21)-(r3+r31)+.5*(r4+r41));
    end

end