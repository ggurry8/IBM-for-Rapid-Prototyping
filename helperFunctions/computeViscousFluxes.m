function [Fv, Gv] = computeViscousFluxes(U, dx, dy, mu, R, gamma)

rho = U(:,:,1);
u   = U(:,:,2);
v   = U(:,:,3);
p   = U(:,:,4);

[Ny, Nx] = size(rho);

% calculate thermodynamics values
T  = p ./ (rho * R);
Cp = gamma * R / (gamma - 1);
Pr = 0.72;
k  = mu * Cp / Pr;

%% 1. gradients via finite difference central differencing 

du_dx = zeros(Ny,Nx);
du_dy = zeros(Ny,Nx);
dv_dx = zeros(Ny,Nx);
dv_dy = zeros(Ny,Nx);
dT_dx = zeros(Ny,Nx);
dT_dy = zeros(Ny,Nx);

% interior points (simple central difference form)
du_dx(:,2:Nx-1) = (u(:,3:Nx) - u(:,1:Nx-2)) / (2*dx);
dv_dx(:,2:Nx-1) = (v(:,3:Nx) - v(:,1:Nx-2)) / (2*dx);
dT_dx(:,2:Nx-1) = (T(:,3:Nx) - T(:,1:Nx-2)) / (2*dx);

du_dy(2:Ny-1,:) = (u(3:Ny,:) - u(1:Ny-2,:)) / (2*dy);
dv_dy(2:Ny-1,:) = (v(3:Ny,:) - v(1:Ny-2,:)) / (2*dy);
dT_dy(2:Ny-1,:) = (T(3:Ny,:) - T(1:Ny-2,:)) / (2*dy);

% divergence
divU = du_dx + dv_dy;

%% 2. CELL-CENTER STRESSES
tau_xx = 2*mu*du_dx - (2/3)*mu*divU; % normal viscous stresses w compressibility
tau_yy = 2*mu*dv_dy - (2/3)*mu*divU; 
tau_xy = mu*(du_dy + dv_dx); % shear stress

%% 3. FACE AVERAGING
% x-faces (Ny x Nx-1)
u_fx = 0.5 * (u(:,1:Nx-1) + u(:,2:Nx));
v_fx = 0.5 * (v(:,1:Nx-1) + v(:,2:Nx));

tau_xx_f = 0.5 * (tau_xx(:,1:Nx-1) + tau_xx(:,2:Nx));
tau_xy_f = 0.5 * (tau_xy(:,1:Nx-1) + tau_xy(:,2:Nx));
qx_f     = -k * 0.5 * (dT_dx(:,1:Nx-1) + dT_dx(:,2:Nx)); % heat flux

% y-faces (Ny-1 x Nx)
u_fy = 0.5 * (u(1:Ny-1,:) + u(2:Ny,:));
v_fy = 0.5 * (v(1:Ny-1,:) + v(2:Ny,:));

tau_yy_f = 0.5 * (tau_yy(1:Ny-1,:) + tau_yy(2:Ny,:));
tau_xy_f2 = 0.5 * (tau_xy(1:Ny-1,:) + tau_xy(2:Ny,:));
qy_f      = -k * 0.5 * (dT_dy(1:Ny-1,:) + dT_dy(2:Ny,:)); 


%% 4. VISCOUS FLUXES
% X-direction flux
Fv = zeros(Ny, Nx-1, 4); % preallocate
Fv(:,:,1) = 0;
Fv(:,:,2) = tau_xx_f;
Fv(:,:,3) = tau_xy_f;
Fv(:,:,4) = u_fx .* tau_xx_f + v_fx .* tau_xy_f - qx_f;

% Y-direction flux 
Gv = zeros(Ny-1, Nx, 4); % preallocate
Gv(:,:,1) = 0;
Gv(:,:,2) = tau_xy_f2;
Gv(:,:,3) = tau_yy_f;
Gv(:,:,4) = u_fy .* tau_xy_f2 + v_fy .* tau_yy_f - qy_f;

end
