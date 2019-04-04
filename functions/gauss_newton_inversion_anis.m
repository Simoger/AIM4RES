function [param, XYZ, Inv] = gauss_newton_inversion_anis(param,XYZ)

%%========================================================================%
%                                                                         %
% Gauss Newton inversion                                                  %
%                                                                         %
%%========================================================================%
%                                                                         %
% In :                                                                    %
% -----------                                                             %
%   param: parameters structure                                           %
%   XYZ: electrodes position structure                                    %
%   draw_plots: 0: no figure ; 1: dray figures                            %
%                                                                         %
% Out :                                                                   %
% -----------                                                             %
%   These structures are saved after each iteration in a iter_#.mat file: %
%   param: parameters structure                                           %
%   XYZ: electrodes position structure                                    %
%   Inv: inversion results at each iteration:                             %
%        Inv.rho: different rho components sections                       %
%        Inv.D: data weighting matrix (from Calc_data_weight.m)           %
%        Inv.d_cal: apparent resistivities calculated at the ith iteration%
%        Inv.rho_app_pos_index: data considered. Negative apparent        %
%                               resistivities are ignored at each iter.   %
%        Inv.rms: rms value                                               %
%        Inv.beta_weighting: weigthing coefficients                       %
%        Inv.sens_weighting: weigthing section                            %
%        Inv.Ki2: Ki^2 value                                              %
%        Inv.CTC: model regularization matrix (from CtC_anis.m)           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%     Abderrezak.Bouchedda@ete.inrs.ca                                    %
%     Institut National de la Recherche Scientifique                      %
%     Centre Eau-Terre-Environnement                                      %
%     http://www.ete.inrs.ca/                                             %
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

%%========================================================================%

    %% ------ Initialize -----------------------------------------
    
    tic

    X = reshape(param.grille(:,1),param.nb_ligne,param.nb_col);
    Z = -reshape(param.grille(:,2),param.nb_ligne,param.nb_col);

    if strcmp(param.inv.invparam,'log resistivity') 

        if isempty(param.const.rho_init) % if no initial model
            [K] = geometric_factor(XYZ,param.flag.geo_factor);
            nb_meas = length(param.MEAS.Res);
            mm = median(K.*param.MEAS.Res);
            
            if param.flag.inv.p == 2
                % solution_vector : [rho1_init ; rho2_init]
                solution_vector = log(mm)*ones(param.nb_ligne*param.nb_col,1);
                solution_vector = [solution_vector ; log(param.anis_init*mm)*ones(param.nb_ligne*param.nb_col,1)]; % + rho_2 initial
            elseif param.flag.inv.p == 3
                % solution_vector : [rho1_init ; rho2_init ; angles_init]
                solution_vector = log(mm)*ones(param.nb_ligne*param.nb_col,1); % rho_1 initial
                solution_vector = [solution_vector ; log(param.anis_init*mm)*ones(param.nb_ligne*param.nb_col,1)]; % + rho_2 initial
                solution_vector = [solution_vector ; log(60)*ones(param.nb_ligne*param.nb_col,1)]; % + angle initial
            end
            
            [rho]= vec2model(exp(solution_vector),param.nb_ligne,param.nb_col, param.flag.inv.p);
            
            % constraints
            if param.const_TrueFalse == true
                param.const.vect = log(param.const_vect);
            end
            
        else % if initial model (has to be a structure)
            nb_meas = length(param.MEAS.Res);
            
            if param.flag.inv.p == 2
                param.const.rho_init.xx = param.const.rho_init.xx';
                param.const.rho_init.zz = param.const.rho_init.zz';
                solution_vector = [param.const.rho_init.xx(:) ; param.const.rho_init.zz(:)];
            end
            
            solution_vector = log(solution_vector(:));
            [rho] = vec2model(exp(solution_vector),param.nb_ligne,param.nb_col, param.flag.inv.p);
            
            if param.const_TrueFalse == true
                param.const.vect = log(param.const_vect);
            end
        end

        if ~isempty(param.const.rho_ref) && param.const.rho_ref ~= 0
            param.const.rho_ref = log(param.const.rho_ref);
        end

    else
        errordlg('gauss_newton_inversion_anis function : inversion parameter type is not defined or not suitable',...
                 'inversion parameter Error');
             
    end

    % define reference model
    if isempty(param.const.rho_ref)
        param.const.rho_ref=0;
    end
    
    %  Set the constraints values
    if isfield(param,'const_ind') && param.const_TrueFalse == true
        param.const.P = ones(length(solution_vector),1);
        param.const.P(param.const_ind) = 0;
        solution_vector(param.const_ind) = param.const.vect;
        cglscd_init = ones(length(solution_vector),1);
        cglscd_init(param.const_ind) = param.const.vect;
        disp('Constraints added')
    else
        param.const.P = [];
        cglscd_init = ones(length(solution_vector),1);
    end

    % weight for matrix C
    wt = ones(size(solution_vector)); 

    % allocate some numbers
    gc = 1; 
    normg0 = 1;
    Ki2 = [];
    
    % data weight matrix calculation
    [D,dtw_pre] = Calc_data_weight(param);
    D_orig = D;
    dtw_pre_orig = dtw_pre;

    % flatness or smoothness or covariance matrix calculation
    [CTC,~,~] = CtC_anis(wt,param,param.flag.inv.p); 

    %% Start Gauss-Newton loop

    itc = 0; 

    while(norm(gc)/normg0 > param.inv.tol && itc < param.inv.maxit && norm(gc) > 1e-20)
        %% iteration count

        itc = itc+1;
        disp(['Iteration #',num2str(itc)])
        toc

        if itc > 1
            param.inv.BETA = max(param.inv.BETA/2,10^-6);
        end

        % potential calculation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        disp(['Sensitivity calculus iteration #',num2str(itc)]);
        [d,S] = calcul_u_S_anis(param,XYZ,rho,1);

        K = geometric_factor(XYZ,param.flag.geo_factor);
        if strcmp(param.inv.invparam,'log resistivity')
            dobs = log((K.*param.MEAS.Res));
        else
            dobs = (K.*param.MEAS.Res);
        end
        D = D_orig;
        dtw_pre = dtw_pre_orig;

        % We ignore for this each iteration the quadripoles measuring a negative apparent resistivity
        rho_app = (K'.*d);
        if strcmp(param.inv.invparam,'log resistivity')
            rho_app_pos_index = (rho_app > .1);
        end

        d = d(rho_app_pos_index);                                                                                               
        K = K(rho_app_pos_index);
        dobs = dobs(rho_app_pos_index);
        D = D(rho_app_pos_index,rho_app_pos_index);
        dtw_pre = dtw_pre(rho_app_pos_index);

        if param.flag.inv.p == 2
            J = [S.xx(rho_app_pos_index,:) S.zz(rho_app_pos_index,:)];
        elseif param.flag.inv.p == 3
            J = [S.xx(rho_app_pos_index,:) S.zz(rho_app_pos_index,:) S.xz(rho_app_pos_index,:)];
        end

        %% Weigthing ------------------------------------------------------
      
        if param.inv.weight == true && itc >= param.inv.weightMinit && itc <= param.inv.weightMaxit
        % -----------------------------------------------------------------

            % Distance weighting - recommended
            if strcmp(param.inv.weightFun,'distance') 
                % beta_weight = .3;
                    beta_weight_1 = .06;
                    beta_weight_2 = .06;
                message_weight = strcat('Weighting based weighting function applied (iter #',num2str(itc),')');
                disp(message_weight)
                wt1 = distance_weighting(param,beta_weight_1);
                wt2 = distance_weighting(param,beta_weight_2);
                wt1 = wt1';
                wt2 = wt2';
                %     figure
                %     imagesc(0:40,0:9,wt1(param.h_x==1/x, param.h_z==1/z));colormap jet;colorbar
                % wt = [wt1 , ones(size(wt2))];
                % wt = [ones(size(wt1)) , wt2];
                wt = [wt1 , wt2];

                wt = wt(:);
            
            % Sensitivity weighting
            elseif strcmp(param.inv.weightFun,'sensitivity') 
                message_weight = strcat('Sensitivity based weighting function applied (iter #',num2str(itc),', Ki = ',num2str(Ki),')');
                disp(message_weight)

                wt = (sum(J.*J));

                beta_weight_1 = 0.1;
                beta_weight_2 = 0.1;

                wt1 = wt(1:param.nb_col*param.nb_ligne).^(beta_weight_1/4);
                wt1 = reshape(wt1,param.nb_col,param.nb_ligne)';
                wt1(:,[1:param.nb_pad_bloc+20+1 , end-(param.nb_pad_bloc+16):end]) = 1;
                % wt1(end-(param.nb_pad_bloc+1):end,:) = 1;
                % wt1(:,[1:param.nb_pad_bloc+param.nb_raff(1)*(param.nb_surr) , end-(param.nb_pad_bloc+param.nb_raff(1)*(param.nb_surr)-1):end]) = 1;
                wt1(end-(param.nb_pad_bloc+param.nb_raff(2)*param.nb_surr):end,:) = 1;

                wt2 = wt(param.nb_col*param.nb_ligne+1:end).^(beta_weight_2/4);
                wt2 = reshape(wt2,param.nb_col,param.nb_ligne)';
                wt2(:,[1:param.nb_pad_bloc+20+1 , end-(param.nb_pad_bloc+16):end]) = 1;
                % wt2(end-(param.nb_pad_bloc+1):end,:) = 1;
                % wt2(:,[1:param.nb_pad_bloc+param.nb_raff(1)*(param.nb_surr) , end-(param.nb_pad_bloc+param.nb_raff(1)*(param.nb_surr)-1):end]) = 1;
                wt2(end-(param.nb_pad_bloc+param.nb_raff(2)*param.nb_surr):end,:) = 1;

                wt = [wt1' , wt2'];
                % wt = [wt1' , ones(size(wt2'))];

                fig1 = figure('visible','off');
                    subplot(121)
                        imagesc(wt1);
                        title('\omega_1')
                        colormap jet;
                        colorbar

                    subplot(122)
                        imagesc(wt2);
                        title('\omega_2')
                        colormap jet;
                        colorbar

                    set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
                    iterName = num2str(itc);
                    pathName = strcat('./savefig/weightings',iterName,'.fig');
                    savefig(fig1,pathName)
                    close

            [CTC,~,~] = CtC_anis(wt,param,param.flag.inv.p); %
            
            else
                error('Weighting function has to be param.inv.weightFun=''distance'' or param.inv.weightFun=''sensitivity''')
                
            end

            
        % Unweight the model after max_itc_weights iterations
        elseif param.inv.weight == true && itc > param.inv.weightMaxit
            beta_weight_1 = nan; 
            beta_weight_2 = nan; 
            wt = ones(size(solution_vector));

            [CTC,~,~] = CtC_anis(wt,param,param.flag.inv.p); %


        % Conserve last applied weighting vector after max_itc_weights iterations 
        else 
            beta_weight_1 = nan; 
            beta_weight_2 = nan; 

        end

        if ~all(K'.*d > 0)
            disp('Negative apparent resistivity')
            disp(' ')
        end

        d = log(K'.*d)';

        clear S

        %% Evaluate model perturbation

        [s,error_cglscd,~] = cglscd(J,cglscd_init, (dobs-d),param.inv.BETA,CTC,solution_vector-param.const.rho_ref,D,param.const.P, 1000, 10^-6,param.inv.appl_fct_reg);
        error_cglscd % indication of the iteration efficiency

        if itc == 1 
            normg0 = norm(gc); 
        end

        %% Test for convergence

        if max(abs(s)) < 1e-6
            errordlg('gauss_newton_inversion_anis function : STEP size too small CONVERGE ',...
            'inversion stop'); 
        end;

        % Try the step
        mu1 = 1; 
        
        % Calculate the value of data objective function 
        fd = 0.5*(d-dobs)'*(D.^2)*(d-dobs);

        %%  Armijo Line search
        
        % Calculate the value of model objective function
        fm = 0.5*( (solution_vector-param.const.rho_ref)'*param.inv.BETA*CTC*(solution_vector-param.const.rho_ref) );

        % Add them for the total Objective function
        fc = fd + fm;

        % Evaluate the gradient
        gc = (((D.^2)*(d-dobs))'*J)' + param.inv.BETA*CTC*(solution_vector-param.const.rho_ref);
        clear J

        for ils = 1 : 6

            disp(['Line search calculation #',num2str(ils)]);

            xt = solution_vector + mu1*s;

            %% Evaluate the new objective function

            if strcmp(param.inv.invparam,'log resistivity') == 1 % log parameters
                if ~isempty(param.const.rho_min)
                    xt(exp(xt) < param.const.rho_min) = param.const.rho_min; 
                end
                if ~isempty(param.const.rho_max)
                    xt(exp(xt) > param.const.rho_max) = param.const.rho_max; 
                end
                rho = vec2model(exp(xt),param.nb_ligne,param.nb_col,param.flag.inv.p);
            end

            [d] = calcul_u_S_anis(param,XYZ,rho,0);

            K = geometric_factor(XYZ,param.flag.geo_factor);
            dobs = log(K.*param.MEAS.Res);
            D = D_orig;
            dtw_pre = dtw_pre_orig;

            % We ignore for each iteration the quadripoles producing a negative apparent resistivity
            rho_app = (K'.*d);
            rho_app_pos_index = (rho_app > .1);

            d = d(rho_app_pos_index);
            K = K(rho_app_pos_index);
            dobs = dobs(rho_app_pos_index);
            D = D(rho_app_pos_index,rho_app_pos_index);
            dtw_pre = dtw_pre(rho_app_pos_index);

            if ~all(K'.*d > 0)
                disp('Negative apparent resistivity')
                disp(' ')
            end

            d = log(K'.*d)';

            fd = 0.5*(d-dobs)'*(D.^2)*(d-dobs);
            fm = 0.5*((xt-param.const.rho_ref)'*param.inv.BETA*CTC*(xt-param.const.rho_ref));
            ft = fd+fm;

            fgoal = fc - param.inv.alp*mu1*(s'*gc);

            if ft < fgoal
                break
            else
                mu1 = mu1/2;    
            end;  

        end % end line search

        rms_model = 100*sum(sqrt((((solution_vector-xt).^2)./(1e-10+solution_vector.^2))/(param.nb_ligne*param.nb_col)));

        if rms_model < param.inv.rms_model
            errordlg('gauss_newton_inversion function : rms_model < desired rms_model',...
            'inversion convergence');
            break; 
        end

        %% Update model

        solution_vector = xt;

        e = (D*(d-dobs));
        Ki = sum(e.^2) / nb_meas;

        % data reweigthing ( L1 or W-estimator)
        if ~strcmp(param.inv.dataweight.type,'const.weight')
            [D,~] = Calc_data_weight(param,e,dtw_pre);
        end

        Ki2 = [Ki2;Ki];
        rms_data = 100*sqrt(sum((d-dobs).^2))/(norm(dobs)); % misfit in percentages
        
        % Stock inv. infos into Inv structure
        Inv.rho{itc} = rho;
        Inv.D{itc} = D;
        Inv.d_cal{itc} = d;
        Inv.rho_app_pos_index{itc} = rho_app_pos_index;
        Inv.rms{itc} = rms_data;
        Inv.beta_weighting{itc} = [beta_weight_1 ; beta_weight_2];
        Inv.sens_weighting{itc} = wt;
        Inv.Ki2 = Ki2;
        Inv.CTC{itc} = CTC;
        
        % Display convergence info : # iter, Ki2, rmse
        [(1:itc)' , Ki2 , cell2mat(Inv.rms)'] %#ok<NOPRT>

        str_save = strcat('iter_',num2str(itc),'.mat');
        save(str_save,'param','XYZ','Inv')

    end
    
end