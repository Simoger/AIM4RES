%% Homogeneous model analytical solution (Li and Uren)

clearvars -except u_num param XYZ

t0 = cputime;

% Resistivites

angle_degres = 30;

angle = angle_degres*pi/180; 
R = [cos(angle) 0 -sin(angle) ; 0 1 0 ; sin(angle) 0 cos(angle)];
RT = R';
res = RT * [100 , 0 , 0
            0 , 100 , 0
            0 , 0 , 400] * R;
                
% Initialisation and parameters definition

r1 = (res(2,2)*res(1,3) - res(1,2)*res(2,3))/(res(1,1)*res(2,2) - res(1,2)^2);
r2 = (res(1,1)*res(2,3) - res(1,2)*res(1,3))/(res(1,1)*res(2,2) - res(1,2)^2);

I = sqrt(det(res));

etaP = 0;
etaR = 0;
pot = zeros(1,2);

for pos = 1:length(XYZ.MEAS.C1)
    
    % Electrodes positions 

        % Source point 1 coordinates
        xP(1) = XYZ.MEAS.C1(pos,1);
        yP(1) = 0;
        zP(1) = XYZ.MEAS.C1(pos,2); 
        I1 = I;
        Int(1) = I;

        % Source point 2 coordinates
        xP(2) = XYZ.MEAS.C2(pos,1);
        yP(2) = 0;
        zP(2) = XYZ.MEAS.C2(pos,2); 
        Int(2) = -I;

        % Measurement point 1 coordinates
        x_pot(1) = XYZ.MEAS.P1(pos,1);
        y_pot(1) = 0;
        z_pot(1) = XYZ.MEAS.P1(pos,2);

        % Measurement point 2 coordinates
        x_pot(2) = XYZ.MEAS.P2(pos,1);
        y_pot(2) = 0;
        z_pot(2) = XYZ.MEAS.P2(pos,2);
 
    for position_potentiel = 1:2
    
        for position_courant = 1:2

            x = x_pot(position_potentiel);
            y = y_pot(position_potentiel);
            z = z_pot(position_potentiel);

            % P0
            XP = x - xP(position_courant);
            YP = y - yP(position_courant);
            ZP = z - zP(position_courant);

            % R
            xR = 2*r1*zP(position_courant) + xP(position_courant);
                XR = x - xR;
            yR = 2*r2*zP(position_courant) + yP(position_courant);
                YR = y - yR;
            zR = - zP(position_courant);
                ZR = z - zR;
                
            % etaP
            etaP = sqrt(res(1,1)*XP^2 + res(2,2)*YP^2 + res(3,3)*ZP^2 + ...
                2*res(1,2)*XP*YP + 2*res(1,3)*XP*ZP + 2*res(2,3)*YP*ZP);

            % etaR
            etaR = sqrt(res(1,1)*XR^2 + res(2,2)*YR^2 + res(3,3)*ZR^2 + ...
                2*res(1,2)*XR*YR + 2*res(1,3)*XR*ZR + 2*res(2,3)*YR*ZR);

            pot_temp = Int(position_courant)/(4*pi)*(1/abs(etaP) + 1/abs(etaR));
            
            % Superposition principle
            pot(position_potentiel) = pot_temp + pot(position_potentiel);
            
            if position_potentiel == 1 & position_courant == 1
                u_P1C1 = pot_temp;
            elseif position_potentiel == 1 & position_courant == 2
                u_P1C2 = pot_temp;
            elseif position_potentiel == 2 & position_courant == 1
                u_P2C1 = pot_temp;
            elseif position_potentiel == 2 & position_courant == 2
                u_P2C2 = pot_temp;
            end

            etaP = 0;
            etaR = 0;
    
        end
        
        end
    
            u(pos) = (u_P1C1-u_P1C2)-(u_P2C1-u_P2C2);
        
    pot_tot(pos,:) = pot;
    pot = zeros(1,2);
    
end

diff_pot = pot_tot(:,1)-pot_tot(:,2);

K = param.K;

u_anal = diff_pot(:);
u_num = u_num(:);

u = u(:);
app_res = K.*u;

%% Plots

ind = true(size(u_anal));
u_numind = u_num(ind);
u_analind = u_anal(ind);
Kind = K(ind);

figure
    subplot(4,4,[1 5 9])
        imagesc(0:25:150,0:5:15,zeros(4,7))
        colormap jet
        hold on
        xm = -12.5:12.5:162.5;
        ym = -2.5*ones(size(xm));
        plot(xm,ym,'ok','MarkerFaceColor','w','MarkerSize',10);
        ym = -2.5:2.5:17.5;
        xm = 75*ones(size(ym));
        plot(xm,ym,'ok','MarkerFaceColor','w','MarkerSize',10);
        set(gca,'FontSize',20,'XTick',[],'YTick',[])
        title('Resistivity model')
        text(-2,2.5,'\rho_1 = 100 \Omega.m','FontSize',16,'FontWeight','bold','Color','k')
        text(-2,7.5,'\rho_3 = 400 \Omega.m','FontSize',16,'FontWeight','bold','Color','k')
        text(-2,12.5,'\theta = 30°','FontSize',16,'FontWeight','bold','Color','k')
        xlabel('x (m)')
        ylabel('y (m)')
    ha(1) = subplot(4,4,[2 3 6 7]);
        plot(u_numind.*Kind)
        hold on
        plot(u_analind.*Kind,'.')
        ylabel('\rho_a (\Omega.m)')
        ylim([-1000 1000])
        xlim(([1 666]))
        set(gca,'FontSize',20,'YTick',-1000:500:1000,'XTick',0:100:666)
        xticklabels({})
        legend('Numerical','Analytical')
    ha(2) = subplot(4,4,[10 11]);
        err_form = 2*(u_numind.*Kind-u_analind.*Kind)./(abs(u_numind.*Kind)+abs(u_analind.*Kind))*100;
        plot(err_form,'.-')
        set(gca,'FontSize',20,'YTick',0:10:20,'XTick',0:100:666)
        xlabel('Measures')
        ylabel('error (%)')
        ylim([-6 20])
        xlim(([1 666]))
    linkaxes(ha, 'x');
    
    subplot(4,4,[4 8 12]);
        histogram(err_form,'BinWidth',.5)
        set(gca,'FontSize',20,'XTick',-5:5:25)
        xlim([-5 25])
        title('Relative error')
        xlabel('error (%)')
        ylabel('# measures')