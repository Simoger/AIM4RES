%% Inversion plots for anisotropic models

graphFontSize = 20;

param.cell_size = [4 4]; 

disp('Ki2 and rms values :')
Ki2_rms = [(1:length(Inv.rho))' Inv.Ki2 cell2mat(Inv.rms)']

n = length(Inv.rho); 

x = param.nb_raff(1)*param.cell_size(1); % cell number in 1m
z = param.nb_raff(2)*param.cell_size(2); % cell number in 1m
dist_from_xborders = 15;
dist_from_surface1 = 0; % should be = 0 to run the whole script, else you'll have an error
dist_from_surface2 = 15;
margin_h = (dist_from_xborders+param.nb_surr./param.cell_size(1))*x; 
depth1 = (dist_from_surface1)*z+1;
depth2 = (dist_from_surface2)*z;        
inv_xx = Inv.rho{1,n}.rho_1(param.h_z==1/param.cell_size(2), param.h_x==1/param.cell_size(1));
inv_zz = Inv.rho{1,n}.rho_2(param.h_z==1/param.cell_size(2), param.h_x==1/param.cell_size(1));
    inv_xx = inv_xx(depth1:depth2, margin_h+1:(size(inv_xx,2)-margin_h+1));
    inv_zz = inv_zz(depth1:depth2, margin_h+1:(size(inv_zz,2)-margin_h+1));
coef_anis = sqrt(inv_zz./inv_xx);

    
% figure
figure('units','normalized','outerposition',[0 0 1 1]);
            
    subplot(2,3,1)
%     subplot(222)
        imagesc(0:10:50, 0:15, param.rho.xx(1:(15*param.cell_size(1)),1:(50*param.cell_size(2))))
            xticks(0:10:50)
            title('Synthetic anisotropic model')
            colormap jet(64)
            set(gca,'fontsize',graphFontSize)
            ylabel('z (m)')
            xlabel('x (m)')
            text(2,2,'\rho_H = 400 \Omega.m','FontSize',graphFontSize,'FontWeight','bold','Color','w')
            text(28,2,'\rho_V = 100 \Omega.m','FontSize',graphFontSize,'FontWeight','bold','Color','w')
            text(2,10,'\rho_H = 40 \Omega.m','FontSize',graphFontSize,'FontWeight','bold','Color','w')
            text(28,10,'\rho_V = 10 \Omega.m','FontSize',graphFontSize,'FontWeight','bold','Color','w')
            xm = [25*ones(1,15) 0:2:50];
            ym = [1:15 zeros(1,length(0:2:50))];
            hold on
            plot(xm,ym,'ok','MarkerFaceColor','w','MarkerSize',7);
            axis tight
            
    subplot(2,3,2)
        minv = log(ceil(min(inv_xx(:))));
        maxv = log(floor(max(inv_xx(:))));
        imagesc(dist_from_xborders:size(inv_xx,2)/x+dist_from_xborders, 0:dist_from_surface2, log(inv_xx))
            title('\rho_H')
            colormap jet(128)
            h = colorbar('FontSize',16,'YTick',linspace(minv,maxv,8)...
                ,'YTickLabel',round(exp(linspace(minv,maxv,8))));
            title(h, '\Omega/m')
            set(gca,'fontsize',graphFontSize)
            ylabel('z (m)')
            xlabel('x (m)')
            axis tight
            
    subplot(2,3,3)
        minv = log(ceil(min(inv_zz(:))));
        maxv = log(floor(max(inv_zz(:))));
        imagesc(dist_from_xborders:size(inv_zz,2)/x+dist_from_xborders, 0:dist_from_surface2, log(inv_zz))
            title('\rho_V')
            colormap jet(128)
            h = colorbar('FontSize',16,'YTick',linspace(minv,maxv,8),'YTickLabel',round(exp(linspace(minv,maxv,8))));
            title(h, '\Omega/m')
            set(gca,'fontsize',graphFontSize)
            ylabel('z (m)')
            xlabel('x (m)')
            axis tight

    subplot(2,3,4)
        imagesc(dist_from_xborders:size(inv_zz,2)/x+dist_from_xborders, 0:.25:(dist_from_surface2-0.25), (coef_anis))
            colorbar
            colormap jet(128)
            title('\lambda')
            set(gca, 'fontsize',graphFontSize)
            h = colorbar;
            title(h, '[-]')
            xlabel('x (m)')
            ylabel('z (m)')
            caxis([0 4])
            axis tight
            
    subplot(2,3,5)
        his = histogram(coef_anis,10);
            title('\lambda distribution')
            xlabel('\lambda')
            ylabel('# cells')
            set(gca,'fontsize',graphFontSize)
            axis([0 max(his.Data(:))+.5 0 max(his.Values)+10])
            box off

%% Data fit plots

[K] = geometric_factor(XYZ,param.flag.geo_factor);
param.K = K;
res_num = Inv.d_cal{n}; 
res_obs = log(param.K.*param.MEAS.Res);
res_obs = res_obs(Inv.rho_app_pos_index{n});

subplot(2,3,6)
    his = histogram(2*(res_obs-res_num)./(abs(res_obs)+abs(res_num))*100,10);
    title('Relative error')
    xlabel('error (%)')
    ylabel('# cells')
    box off
    ylim([0 max(his.Values)+10])
    set(gca,'FontSize',graphFontSize)