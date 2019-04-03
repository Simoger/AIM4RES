%% 2 layer model analytical solution (Telford)

rhoX = param.rho.xx;
rho_calc = param.K(:).*u_num(:);

nb_points = 50;
nb_images = 100;

rho_h1 = 100; 
rho_v1 = 400; 
rho_n1 = sqrt(rho_v1*rho_h1);
rho_h2 = 10; 
rho_v2 = 40; 
rho_n2 = sqrt(rho_v2*rho_h2);
kn = (rho_n2-rho_n1)/(rho_n2+rho_n1);
f = sqrt(rho_v1/rho_h1);
z = 4;
n = 2; % dipole-dipole separation

rho_aa = zeros(nb_points,1);
summ = zeros(nb_points,1);
p = zeros(nb_points,1);
r1 = zeros(nb_points,1);
r2 = zeros(nb_points,1);
r3 = zeros(nb_points,1);
r4 = zeros(nb_points,1);

for i = 1:nb_points
    % Wenner
    j = i;
    r1(i) = j; 
    r2(i) = 2*j; 
    r3(i) = 2*j; 
    r4(i) = j;
    ABs2(i) = 3*j/2;
    
    p(i) = (1/r1(i)-1/r2(i)-1/r3(i)+1/r4(i))^(-1);
    summ(i) = 0;
    operande = 0;
    for m = 1:nb_images
        operande = 1/sqrt(r1(i)^2+4*m^2*f^2*z^2) - 1/sqrt(r2(i)^2+4*m^2*f^2*z^2)...
            -1/sqrt(r3(i)^2+4*m^2*f^2*z^2) + 1/sqrt(r4(i)^2+4*m^2*f^2*z^2);
        operande = kn^m*operande;
        summ(i) = summ(i) + operande;
    end
    rho_aa(i) = rho_n1*(1+2*p(i)*summ(i));
    
end

ABs2 = ABs2(2:end)';
rho_aa = rho_aa(2:end);
rho_calc = rho_calc(2:end);

%% Plots

figure
    subplot(4,3,[1 4 7])
        imagesc(0:25:150,0:5:15,rhoX(1:60,:))
        colormap jet
        hold on
        xm = 0:10:150;
        ym = -.1*ones(1,16);
        plot(xm,ym,'ok','MarkerFaceColor','w','MarkerSize',10);
        xlabel('x (m)')
        ylabel('z (m)')
        set(gca,'FontSize',20)
        title('Resistivity model')
        text(10,2,'\rho_H = 100 \Omega.m','FontSize',20,'FontWeight','bold','Color','w')
        text(80,2,'\rho_V = 400 \Omega.m','FontSize',20,'FontWeight','bold','Color','w')
        text(10,10,'\rho_H = 10 \Omega.m','FontSize',20,'FontWeight','bold','Color','w')
        text(80,10,'\rho_V = 40 \Omega.m','FontSize',20,'FontWeight','bold','Color','w')
    subplot(4,3,[2 3 5 6])
        semilogx(ABs2,rho_aa,'ok')
        set(gca,'FontSize',20,'YTick',0:50:220)
        xticklabels({})
        ylabel('\rho_a (\Omega.m)')
        ylim([0 220])
        hold on
        x = (6:3:150)/2;
        semilogx(x,rho_calc,'LineWidth',4)
        legend('Analytical','Numerical')
    subplot(4,3,[8 9])
        semilogx(ABs2,(rho_calc-rho_aa)./rho_aa*100,'.-')
        set(gca,'FontSize',20,'YTick',-1:2)
        xlabel('AB/2 (m)')
        ylabel('error (%)')
        ylim([-1 2])