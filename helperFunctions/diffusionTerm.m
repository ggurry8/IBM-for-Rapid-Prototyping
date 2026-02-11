% applies sutherlands law for finding viscosity: https://www.cfd-online.com/Wiki/Sutherland's_law
% 

function [mu_du_dx, mu_dv_dy] = diffusionTerm(UL_x, UR_x, UL_y, UR_y, dx, dy, R)

% derive temperature face centered values
[rhoL_x, ~, ~, pL_x] = deal(UL_x(:,:,1), UL_x(:,:,2), UL_x(:,:,3), UL_x(:,:,4));
[rhoL_y, ~, ~, pL_y] = deal(UL_y(:,:,1), UL_y(:,:,2), UL_y(:,:,3), UL_y(:,:,4));
[rhoR_x, ~, ~, pR_x] = deal(UR_x(:,:,1), UR_x(:,:,2), UR_x(:,:,3), UR_x(:,:,4));
[rhoR_y, ~, ~, pR_y] = deal(UR_y(:,:,1), UR_y(:,:,2), UR_y(:,:,3), UR_y(:,:,4));

% 
T_x = (pR_x ./ (rhoR_x .* R) + pL_x ./ (rhoL_x .* R))/2; % Temperature face centered value (along x)
T_y = (pR_y ./ (rhoR_y .* R) + pL_y ./ (rhoL_y .* R))/2; % Temperature face centered value (along x)

% coefficients
% mu_ref = 1.716 * 10 ^ -5;
% T_ref = 273.15;
S = 110.4;
C1 = 1.458 * 10 ^ -6; 

mu_x = (C1 * T_x .^ 1.5) ./ (T_x + S);
mu_y = (C1 * T_y .^ 1.5) ./ (T_y + S);

du_dx = (UR_x(:,:,2) - UL_x(:,:,2))/dx;
dv_dy = (UR_y(:,:,3) - UL_y(:,:,3))/dy;

mu_du_dx = mu_x .* du_dx;
mu_dv_dy = mu_y .* dv_dy;

end