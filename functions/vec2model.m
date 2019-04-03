function [rho]= vec2model(solution_vector,nb_ligne,nb_col,flag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Convert the resistivity solution vector to a structure which elements are 
% s.xx, s.zz and s.xz
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% flag == 2 rhoxx rhozz
% flag == 3 rhoxx rhozz rhoxz
    
    if flag == 2
       n = length(solution_vector)/2;
       rho.xx = reshape(solution_vector(1:n),nb_col,nb_ligne)';
%        rho.xx(rho.xx < 1)=1;
       rho.yy = rho.xx;
       rho.zz = reshape(solution_vector(1+n:n*2),nb_col,nb_ligne)';
%        rho.zz(rho.zz < 1)=1;
       rho.xz = zeros(nb_ligne,nb_col);
       rho.angle = zeros(nb_ligne,nb_col);
       rho.rho_1 = rho.xx;
       rho.rho_2 = rho.zz;
    elseif flag == 3
        n = length(solution_vector)/3;
        rho.rho_1 = reshape(solution_vector(1:n),nb_col,nb_ligne)';
        rho.rho_2 = reshape(solution_vector(1+n:n*2),nb_col,nb_ligne)';
        rho.rho_1(rho.rho_1 < 1) = 1;
        rho.rho_2(rho.rho_2 < 1) = 1;
        rho.angle = reshape(solution_vector(1+2*n:end),nb_col,nb_ligne)';
		rho.xx = rho.rho_1.*cosd(rho.angle).^2 + rho.rho_2.*sind(rho.angle).^2;
        rho.zz = rho.rho_1.*sind(rho.angle).^2 + rho.rho_2.*cosd(rho.angle).^2;
        rho.xz = (rho.rho_2-rho.rho_1).*cosd(rho.angle).*sind(rho.angle);
        rho.yy = rho.rho_1;
        
    end
end