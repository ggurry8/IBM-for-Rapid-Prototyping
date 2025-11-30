% NOTE: THIS CODE WAS INITIALLY CREATED WITH ChatGPT - however it has been modified significantly from the original version 
% Original GPT chat - https://chatgpt.com/share/68dfd53a-d1d4-8011-9519-0565c66807f9
% The code from this chat is located within the same directory as this file and is called "GPT_original.m"

%% simple_2D_compressible_CFD_solver.m
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
%
% Notations (WIP)
% FLOW:
%   U - cell variable averages array [rho; u; v; p]
%   use of con, ext suffixes on U denotes conservative and 
%   extended (for application of wall BCs) arrays
%   rho - density
%   u,v - x,y velocity
%   p - pressure
%   
% GRID:
%   Lx,Ly - Domain physical size (x,y)
%   Nx,Ny - Domain grid size (x,y)
%   

function main_loop(app)

% --------------------
% Problem / grid setup
% --------------------
Lx = app.XDomainEditField.Value; Ly = app.YDomainEditField.Value;         % domain size
Nx = app.XGridEditField.Value; Ny = app.YGridEditField.Value;         % grid cells in x/y
dx = Lx / Nx; dy = Ly / Ny; % grid spacing (x,y) 
x = linspace(dx/2, Lx-dx/2, Nx);
y = linspace(dy/2, Ly-dy/2, Ny);
[X,Y] = meshgrid(x,y);

% Import Geometry
[xGeometry, yGeometry] = aerofoil_generation(app.NACA4DigitCodeEditField.Value, app.ChordLengthEditField.Value,...
    app.AngleofAttackEditField.Value, app.NumberofPointsEditField.Value);
xGeometry = xGeometry + Lx * app.CoordinatesofLeadingEdgeX.Value/100;
yGeometry = yGeometry + Ly * app.CoordinatesofLeadingEdgeY.Value/100; % adjust position of leading edge
geometryObj = geometry(X, Y, xGeometry, yGeometry);
if ~any(geometryObj.withinGeometry,"all")
    uialert(app, "Geometry may be too large relative to grid size",'Warning: ','Icon','Warning')
end
geometryObj = geometryObj.geoInterpret;

% Physical parameters
gamma = 1.4;                % ratio of specific heats
R = 287;                 % gas constant (not used explicitly if using p=(gamma-1)*rho*e)

% Flow initial condition (uniform flow)
rho0 = 1.0;                 % background density
U0 = app.InletVelocityEditField.Value;                   % background velocity in x
V0 = 0.0;                   % background velocity in y
p0 = 100000;                   % background pressure in Pa

% primitive variables
rho = rho0;
U = U0 * ones(size(rho));
V = V0 * ones(size(rho));
p = p0 * ones(size(rho));

% Conservative variables (Ucon = [rho; rho*u; rho*v; rho*E])
E = p./((gamma-1)) ./ rho + 0.5*(U.^2 + V.^2); % specific total energy
rhoe = rho .* E;
Ucon = zeros(Ny, Nx, 4); % initialise variable
Ucon(:,:,1) = rho;       % conservative variable flow properties
Ucon(:,:,2) = rho .* U;
Ucon(:,:,3) = rho .* V;
Ucon(:,:,4) = rhoe;

% solver parameters
CFL = app.CFLEditField.Value;                 % CFL number for stability
tEnd = app.EndTimesEditField.Value;           % final time
t = 0;                                          % initialise variable
visc_alpha = app.ArtificialViscEditField.Value;  % artificial viscosity coeff (set >0 to add stabilization)

boundaryTypes = app.boundaryTypes;           % 4x1 struct of bcs for each wall
boundaryVelocities = app.boundaryVelocities;

residuals = []; % initialise variable

% diagnostics / plotting
% figure('units','normalized','position',[0.1 0.1 0.6 0.6]);
plotEvery = app.UpdateFreqEditField.Value;
iteration = 1;

% minmod limiter
% minmod = @(a,b) (sign(a)+sign(b)).*min(abs(a),abs(b))/2 .* (sign(a)+sign(b)~=0); % TODO this isn't used

% ---------------------
% Main time loop
% ---------------------
while t < tEnd
    % compute primitive variables & flux Jacobian max wave speed
    Uinitial = cons2prim(Ucon, gamma);
    [rho, u, v, p] = deal(Uinitial(:,:,1), Uinitial(:,:,2), Uinitial(:,:,3), Uinitial(:,:,4)); % assign the primitives to individual variables
    c = sqrt(max(gamma * p ./ rho, 1e-12));
    Smax = max(max(sqrt(u.^2 + v.^2) + c));

    % timestep from CFL
    dt = CFL * min(dx,dy) / Smax;
    if t + dt > tEnd, dt = tEnd - t; end

    % RK2 (Heun): stage 1
    [RHSout, iterationResiduals] = RHS(Ucon, R, dx, dy, gamma, visc_alpha, boundaryTypes, boundaryVelocities);
    residuals = [residuals; iterationResiduals]; %#ok - residuals vector is relatively small memory wise
    U1 = Ucon + dt * RHSout;

    % stage 2
    U2 = 0.5 * (Ucon + U1 + dt * RHS(U1, R, dx, dy, gamma, visc_alpha, boundaryTypes, boundaryVelocities));

    Ucon = U2; % update Ucon with new values
    Ucon = applyGeometryBCs(Ucon, geometryObj);
    t = t + dt;

    Ufinal = cons2prim(Ucon, gamma);
    app.simData = simDataStore(Ufinal, X, Y, gamma, R, t, geometryObj, residuals);

    % update flow visualiser occasionally
    if mod(iteration,plotEvery)==0 || t>=tEnd
        varToDisplay = app.DisplayVarDropDown.Value;
        flow_visualiser(app, varToDisplay);
        plotResiduals(app.ResidualsPlot, iteration, residuals);
    end

    % perform end of loop checks
    app.SimulationTimeField.Value = t;
    iteration = iteration + 1;
    app.IterationNumberField.Value = iteration;
    if app.TerminateSimProperty
        app.TerminateSimProperty = false; % reset property
        % pass data out of app
        break
    elseif app.PauseSimulationButton.Value
        % pass sim data out of app
        while app.PauseSimulationButton.Value
            if app.TerminateSimProperty
                app.TerminateSimProperty = false; % reset property
            break
            end
            pause(1)
        end
    end
end
end

% --------------------
% RHS function (finite-volume / fluxes)
% --------------------
function [RHSout, residuals] = RHS(Ucon, R, dx, dy, gamma, visc_alpha, boundary, velocities)
    % Ucon: Ny x Nx x 4 array
    % returns dU/dt of same shape

    % reconstruct states at faces using MUSCL (linear reconstruction) with minmod limiter
    % For simplicity handle interior with one-sided for boundaries depending on 'boundary'

    % get primitives
    U = cons2prim(Ucon, gamma);

    % extend arrays for BCs
    Uext = applyWallBCs(U, boundary, velocities); % adds one-cell padding around arrays (Ny+2 x Nx+2), and applies BCs there
    Uext_con = prim2cons(Uext, gamma);

    % compute Rusanov flux, x-direction (faces perpendicular to x)
    [F, UL_x, UR_x] = rusanov_flux(Uext_con, 'x', gamma);

    % compute Rusanov flux, y-direction (faces perpendicular to y)
    [G, UL_y, UR_y] = rusanov_flux(Uext_con, 'y', gamma);

    % compute viscosity
    [mu_du_dx, mu_dv_dy] = diffusionTerm(UL_x, UR_x, UL_y, UR_y, dx, dy, R);

    % fluxes including diffusion/viscosity
    F = F - mu_du_dx;
    G = G - mu_dv_dy;
    
    % residuals
    residuals = residualsEuler(F, G, dx,dy);

    % assemble RHS via divergence of fluxes
    RHSout = flow_divergence(F, G, dx, dy);

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
% helper: applyGeometryBCs (simple padding)
% --------------------
function Ucon = applyGeometryBCs(Ucon, geometryObj)

% preallocate Ughost
Ughost = Ucon; % initialise variable

% loop through fluid cells within or bordering the geometry
for j = 2:size(Ucon,1)-1
    for i = 2:size(Ucon,2)-1
        if geometryObj.geometryBorder(j,i) || geometryObj.withinGeometry(j,i)
            Ughost(j,i,2) = 0; % x velocity
            Ughost(j,i,3) = 0; % y velocity
        end
    end
end
% update
Ucon = Ughost;
end

% --------------------
% helper: Divergence of fluxes 
% --------------------
function RHSout = flow_divergence(Fx, Fy, dx, dy)

% remove BCs
Fx_out = Fx(2:end-1,1:end-1,:);
Fx_in = Fx(2:end-1,2:end,:);
Fy_out = Fy(1:end-1,2:end-1,:);
Fy_in = Fy(2:end,2:end-1,:);

RHSout = - (Fx_in - Fx_out)/dx - (Fy_in - Fy_out)/dy;  % this is essentially Fxdx/dV - Fdy/dV, but cancelled out

end

