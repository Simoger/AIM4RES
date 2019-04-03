%% Inversion plots for isotropic models

graphFontSize = 20;

param.cell_size = [4 4]; 

disp('Ki2 and rms values :')
Ki2_rms = [(1:length(Inv.rho))' Inv.Ki2 cell2mat(Inv.rms)']

n = length(Inv.rho); 

res_num = Inv.d_cal{n}; 
res_obs = log(param.K.*param.MEAS.Res);

x = param.nb_raff(1)*param.cell_size(1); 
z = param.nb_raff(2)*param.cell_size(2); 
dist_from_xborders = 15;
dist_from_surface1 = 0; % should be = 0 to run the whole script, else you'll have an error
dist_from_surface2 = 15;
margin_h = (dist_from_xborders+param.nb_surr./param.cell_size(1))*x; 
depth1 = (dist_from_surface1)*z+1;
depth2 = (dist_from_surface2)*z;        
inv = Inv.rho{1,n}(param.h_z==1/param.cell_size(2), param.h_x==1/param.cell_size(1));
    inv = inv(depth1:depth2, margin_h+1:(size(inv,2)-margin_h+1));

% figure
figure('units','normalized','outerposition',[0 0 1 1]);
            
    subplot(1,3,1)
        imagesc(0:10:50, 0:15, param.rho.xx(1:(15*param.cell_size(1)),1:(50*param.cell_size(2))))
            xticks(0:10:50)
            title('Synthetic anisotropic model')
            colormap jet(64)
            set(gca,'fontsize',graphFontSize)
            ylabel('z (m)')
            xlabel('x (m)')
            text(2,2,'\rho = 200 \Omega.m','FontSize',graphFontSize,'FontWeight','bold','Color','w')
            text(2,10,'\rho = 20 \Omega.m','FontSize',graphFontSize,'FontWeight','bold','Color','w')
            xm = [25*ones(1,15) 0:2:50];
            ym = [1:15 zeros(1,length(0:2:50))];
            hold on
            plot(xm,ym,'ok','MarkerFaceColor','w','MarkerSize',7);
            axis square
            
    subplot(1,3,2)
        minv = log(10);
        maxv = log(300);
        imagesc(dist_from_xborders:size(inv_xx,2)/x+dist_from_xborders, 0:dist_from_surface2, log(inv))
            title('\rho_H')
            colormap jet(128)
            h = colorbar('FontSize',16,'YTick',linspace(minv,maxv,8)...
                ,'YTickLabel',round(exp(linspace(minv,maxv,8))));
            title(h, '\Omega/m')
            set(gca,'fontsize',graphFontSize)
            ylabel('z (m)')
            xlabel('x (m)')
            caxis([minv maxv])
            axis square
            ax = gca;
            axpos = ax.Position;
            ax.Position = axpos;

    subplot(1,3,3)
        his = histogram(2*(res_obs-res_num)./(abs(res_obs)+abs(res_num))*100);
        title('Relative error')
        xlabel('error (%)')
        ylabel('# cells')
        box off
        xlim([-10 10])
        ylim([0 max(his.Values)+10])
        set(gca,'FontSize',graphFontSize)