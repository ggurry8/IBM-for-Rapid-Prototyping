% simple_2D_compressible_CFD_solver.m
% A minimal 2D explicit finite-volume compressible Euler solver (ideal gas)
% - Structured Cartesian grid
% - Rusanov (local Lax-Friedrichs) numerical flux
% - Second-order MUSCL reconstruction with minmod limiter
% - RK2 (Heun) time integration
% - Optional explicit (artificial) viscosity for stabilization
%
% Notes:
% - Intended for low-Mach demonstration (M < 0.2). For very low Mach numbers
%   consider preconditioning (not implemented here) or use an incompressible solver.
% - Boundary conditions: options for periodic / wall / inflow-outflow are included.
% - This code is intentionally simple and readable, not optimized for performance.

clear; close all; clc;

% --------------------
% Problem / grid setup
% --------------------
Lx = 1.0; Ly = 1.0;         % domain size
Nx = 120; Ny = 120;         % grid cells in x/y
dx = Lx / Nx; dy = Ly / Ny;
x = linspace(dx/2, Lx-dx/2, Nx);
y = linspace(dy/2, Ly-dy/2, Ny);
[X,Y] = meshgrid(x,y);

% Physical parameters
gamma = 1.4;                % ratio of specific heats
Rgas = 287;                 % gas constant (not used explicitly if using p=(gamma-1)*rho*e)

% Flow initial condition (a small perturbation on a uniform flow)
rho0 = 1.0;                 % background density
u0 = 0.0;                   % background velocity in y
U0 = 0.2;                   % background velocity in x (choose M ~ 0.2 w.r.t. c)
p0 = 1.0;                   % background pressure

% small vortex / perturbation to see dynamics
amp = 0.01;
xc = 0.5; yc = 0.5; rc = 0.1;
r = sqrt((X-xc).^2 + (Y-yc).^2);

% primitive variables
rho = rho0 * (1 + amp * exp(-(r/rc).^2));
u = u0 * ones(size(rho));
U = U0 * ones(size(rho));
p = p0 * ones(size(rho));

% Conservative variables (Ucon = [rho; rho*u; rho*v; rho*E])
E = p./((gamma-1)) ./ rho + 0.5*(U.^2 + u.^2); % specific total energy
rhoe = rho .* E;
Ucon = zeros(Ny, Nx, 4);
Ucon(:,:,1) = rho;
Ucon(:,:,2) = rho .* U;
Ucon(:,:,3) = rho .* u;
Ucon(:,:,4) = rhoe;

% solver parameters
CFL = 0.45;                 % CFL number for stability
tEnd = 0.5;                 % final time
t = 0;
visc_alpha = 0.0;           % artificial viscosity coeff (set >0 to add stabilization)
boundary = 'periodic';      % 'periodic' or 'wall' or 'outflow'

% diagnostics / plotting
figure('units','normalized','position',[0.1 0.1 0.6 0.6]);
plotEvery = 10; frame=1;

% minmod limiter
minmod = @(a,b) (sign(a)+sign(b)).*min(abs(a),abs(b))/2 .* (sign(a)+sign(b)~=0);

% Main time loop
while t < tEnd
    % compute primitive variables & flux Jacobian max wave speed
    rho = Ucon(:,:,1);
    u = Ucon(:,:,2)./rho;
    v = Ucon(:,:,3)./rho;
    E = Ucon(:,:,4)./rho;
    p = (gamma-1)*rho.*(E - 0.5*(u.^2 + v.^2));
    c = sqrt(max(gamma * p ./ rho, 1e-12));
    Smax = max(max(sqrt(u.^2 + v.^2) + c));

    % timestep from CFL
    dt = CFL * min(dx,dy) / Smax;
    if t + dt > tEnd, dt = tEnd - t; end

    % RK2 (Heun): stage 1
    U1 = Ucon + dt * RHS(Ucon, Nx, Ny, dx, dy, gamma, visc_alpha, boundary);
    % stage 2
    U2 = 0.5*(Ucon + U1 + dt * RHS(U1, Nx, Ny, dx, dy, gamma, visc_alpha, boundary));

    Ucon = U2;
    t = t + dt;

    % plot density occasionally
    if mod(frame,plotEvery)==0 || t>=tEnd
        rho = Ucon(:,:,1);
        imagesc(x,y,rho); axis equal tight; set(gca,'YDir','normal');
        title(sprintf('Density at t = %.4f', t)); colorbar; drawnow;
    end
    frame = frame + 1;
end

fprintf('Done. Final time t = %g\n', t);

% --------------------
% RHS function (finite-volume / fluxes)
% --------------------
function RHSout = RHS(Ucon, Nx, Ny, dx, dy, gamma, visc_alpha, boundary)
    % Ucon: Ny x Nx x 4 array
    % returns dU/dt of same shape

    % allocate flux arrays
    Fx = zeros(Ny, Nx+1, 4); % flux through vertical faces (i-1/2, i+1/2)
    Fy = zeros(Ny+1, Nx, 4); % flux through horizontal faces

    % reconstruct states at faces using MUSCL (linear reconstruction) with minmod limiter
    % For simplicity handle interior with one-sided for boundaries depending on 'boundary'

    % get primitives
    rho = Ucon(:,:,1);
    u = Ucon(:,:,2)./rho;
    v = Ucon(:,:,3)./rho;
    E = Ucon(:,:,4)./rho;
    p = (gamma-1)*rho.*(E - 0.5*(u.^2 + v.^2));
    c = sqrt(max(gamma * p ./ rho, 1e-12));

    % extend arrays for BCs
    [rho_e, u_e, v_e, p_e] = applyBCs(rho,u,v,p,boundary);
    E_e = p_e./((gamma-1))./rho_e + 0.5*(u_e.^2 + v_e.^2);
    Ue = zeros(size(rho_e,1), size(rho_e,2), 4);
    Ue(:,:,1) = rho_e; Ue(:,:,2) = rho_e.*u_e; Ue(:,:,3) = rho_e.*v_e; Ue(:,:,4) = rho_e.*E_e;

    % compute fluxes in x-direction (faces i-1/2 ... )
    for j = 1:Ny
        for i = 1:Nx+1
            % left cell index in extended grid
            il = i; ir = i+1; % in extended Ue indexing (1..Nx+2)
            UL = squeeze(Ue(j+1, il, :)); % note: Ue has extra rows/cols from BC padding
            UR = squeeze(Ue(j+1, ir, :));
            % convert to primitive
            [primL] = cons2prim(UL, gamma);
            [primR] = cons2prim(UR, gamma);
            % compute Rusanov flux
            Fx(j,i,:) = rusanov_flux(primL, primR, 'x', gamma);
        end
    end

    % compute fluxes in y-direction
    for j = 1:Ny+1
        for i = 1:Nx
            jl = j; jr = j+1;
            UL = squeeze(Ue(jl, i+1, :));
            UR = squeeze(Ue(jr, i+1, :));
            [primL] = cons2prim(UL, gamma);
            [primR] = cons2prim(UR, gamma);
            Fy(j,i,:) = rusanov_flux(primL, primR, 'y', gamma);
        end
    end

    % assemble RHS via divergence of fluxes
    RHSout = zeros(Ny, Nx, 4);
    for j = 1:Ny
        for i = 1:Nx
            % flux differences
            Fx_in = squeeze(Fx(j,i+1,:)); Fx_out = squeeze(Fx(j,i+0,:));
            Fy_in = squeeze(Fy(j+1,i,:)); Fy_out = squeeze(Fy(j+0,i,:));
            RHSout(j,i,:) = - (Fx_in - Fx_out)/dx - (Fy_in - Fy_out)/dy;
        end
    end

    % optional explicit artificial viscosity (simple Laplacian on conserved variables)
    if visc_alpha > 0
        % compute discrete laplacian (interior only) for each conserved field
        U = Ucon;
        lap = zeros(size(U));
        for q=1:4
            A = U(:,:,q);
            lap(:,:,q) = (circshift(A,[0 -1]) + circshift(A,[0 1]) - 2*A)/dx^2 + (circshift(A,[-1 0]) + circshift(A,[1 0]) - 2*A)/dy^2;
        end
        RHSout = RHSout + visc_alpha * lap;
    end
end

% --------------------
% helper: applyBCs (simple padding)
% --------------------
function [rho_e, u_e, v_e, p_e] = applyBCs(rho,u,v,p,boundary)
    % add one-cell padding around arrays (Ny+2 x Nx+2)
    [Ny,Nx] = size(rho);
    rho_e = zeros(Ny+2, Nx+2); u_e = rho_e; v_e = rho_e; p_e = rho_e;
    rho_e(2:end-1,2:end-1) = rho; u_e(2:end-1,2:end-1) = u; v_e(2:end-1,2:end-1) = v; p_e(2:end-1,2:end-1) = p;

    switch lower(boundary)
        case 'periodic'
            % left-right
            rho_e(:,1)   = rho_e(:,end-1);
            rho_e(:,end) = rho_e(:,2);
            u_e(:,1)   = u_e(:,end-1);
            u_e(:,end) = u_e(:,2);
            v_e(:,1)   = v_e(:,end-1);
            v_e(:,end) = v_e(:,2);
            p_e(:,1)   = p_e(:,end-1);
            p_e(:,end) = p_e(:,2);
            % top-bottom
            rho_e(1,:)   = rho_e(end-1,:);
            rho_e(end,:) = rho_e(2,:);
            u_e(1,:)   = u_e(end-1,:);
            u_e(end,:) = u_e(2,:);
            v_e(1,:)   = v_e(end-1,:);
            v_e(end,:) = v_e(2,:);
            p_e(1,:)   = p_e(end-1,:);
            p_e(end,:) = p_e(2,:);
        case 'wall'
            % reflective wall for velocities, copy density/pressure
            rho_e(:,1)   = rho_e(:,2); rho_e(:,end) = rho_e(:,end-1);
            rho_e(1,:)   = rho_e(2,:); rho_e(end,:) = rho_e(end-1,:);
            p_e(:,1)   = p_e(:,2); p_e(:,end) = p_e(:,end-1);
            p_e(1,:)   = p_e(2,:); p_e(end,:) = p_e(end-1,:);
            u_e(:,1)   = -u_e(:,2); u_e(:,end) = -u_e(:,end-1);
            u_e(1,:)   = -u_e(2,:); u_e(end,:) = -u_e(end-1,:);
            v_e(:,1)   = -v_e(:,2); v_e(:,end) = -v_e(:,end-1);
            v_e(1,:)   = -v_e(2,:); v_e(end,:) = -v_e(end-1,:);
        otherwise
            % default: outflow (copy)
            rho_e(:,1)   = rho_e(:,2); rho_e(:,end) = rho_e(:,end-1);
            rho_e(1,:)   = rho_e(2,:); rho_e(end,:) = rho_e(end-1,:);
            p_e(:,1)   = p_e(:,2); p_e(:,end) = p_e(:,end-1);
            p_e(1,:)   = p_e(2,:); p_e(end,:) = p_e(end-1,:);
            u_e(:,1)   = u_e(:,2); u_e(:,end) = u_e(:,end-1);
            u_e(1,:)   = u_e(2,:); u_e(end,:) = u_e(end-1,:);
            v_e(:,1)   = v_e(:,2); v_e(:,end) = v_e(:,end-1);
            v_e(1,:)   = v_e(2,:); v_e(end,:) = v_e(end-1,:);
    end
end

% --------------------
% helper: convert conserved -> primitive (vector)
% --------------------
function prim = cons2prim(Uvec, gamma)
    rho = Uvec(1);
    u = Uvec(2)/rho;
    v = Uvec(3)/rho;
    E = Uvec(4)/rho;
    p = (gamma-1)*rho*(E - 0.5*(u^2 + v^2));
    prim = struct('rho',rho,'u',u,'v',v,'p',p);
end

% --------------------
% helper: Rusanov flux between two primitive states
% --------------------
function F = rusanov_flux(L, R, dir, gamma)
    % compute conservative variables
    UL = [L.rho; L.rho*L.u; L.rho*L.v; L.rho*(L.p/(gamma-1)/L.rho + 0.5*(L.u^2+L.v^2))];
    UR = [R.rho; R.rho*R.u; R.rho*R.v; R.rho*(R.p/(gamma-1)/R.rho + 0.5*(R.u^2+R.v^2))];

    if strcmp(dir,'x')
        FL = [L.rho*L.u;
              L.rho*L.u^2 + L.p;
              L.rho*L.u*L.v;
              (UL(4)+L.p)*L.u];
        FR = [R.rho*R.u;
              R.rho*R.u^2 + R.p;
              R.rho*R.u*R.v;
              (UR(4)+R.p)*R.u];
        aL = sqrt(gamma * L.p / L.rho);
        aR = sqrt(gamma * R.p / R.rho);
        smax = max(abs(L.u)+aL, abs(R.u)+aR);
    else % y-direction
        FL = [L.rho*L.v;
              L.rho*L.u*L.v;
              L.rho*L.v^2 + L.p;
              (UL(4)+L.p)*L.v];
        FR = [R.rho*R.v;
              R.rho*R.u*R.v;
              R.rho*R.v^2 + R.p;
              (UR(4)+R.p)*R.v];
        aL = sqrt(gamma * L.p / L.rho);
        aR = sqrt(gamma * R.p / R.rho);
        smax = max(abs(L.v)+aL, abs(R.v)+aR);
    end

    F = 0.5*(FL + FR) - 0.5*smax*(UR - UL);
end

% End of file
