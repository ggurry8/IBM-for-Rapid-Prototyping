% NOTE: THIS CODE WAS INITIALLY CREATED WITH ChatGPT - however it has been modified significantly from the original version 
% Original GPT chat - https://chatgpt.com/share/68dfd53a-d1d4-8011-9519-0565c66807f9
% The code from this chat is located within the same directory as this file and is called "GPT_original.m"

%% main_loop.m
% Core function for controlling the 2D explicit finite-volume compressible (ideal gas) N-S solver for testing IBM approaches.
% - Uniform Structured Cartesian grid
% - HLLC numerical flux
% - Second-order MUSCL reconstruction with minmod limiter
% - RK2 (Heun) time integration
% - Immersed Boundary Methods
% - Optional explicit (artificial) viscosity for stabilization
%
% Notes:
% This code requires the Image Processing and App Designer Toolboxes to function correctly.
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
minGridSpacing = min(dx,dy);

% Import Geometry
NACA4Code = app.NACA4DigitCodeEditField.Value;
chordLength = app.ChordLengthEditField.Value;
thickPcent = str2double([NACA4Code(3:4)])/100;
refA = thickPcent * chordLength;
[xGeometry, yGeometry] = aerofoil_generation(NACA4Code, chordLength,...
    app.AngleofAttackEditField.Value, app.NumberofPointsEditField.Value);
xGeometry = xGeometry + Lx * app.CoordinatesofLeadingEdgeX.Value/100;
yGeometry = yGeometry + Ly * app.CoordinatesofLeadingEdgeY.Value/100; % adjust position of leading edge
geometryObj = geometry(X, Y, xGeometry, yGeometry);
if ~any(geometryObj.withinGeometry,"all")
    uialert(app.UIFigure, "Geometry may be too large relative to grid size",'Warning: ','Icon','Warning')
end
geometryObj = geometryObj.geoInterpret;

% Physical parameters
gamma = 1.4;                % ratio of specific heats
R = 287;                 % gas constant (not used explicitly if using p=(gamma-1)*rho*e)

% Flow initial condition (uniform flow)
rho0 = 1.225;                 % background density
U0 = app.InletVelocityEditField.Value;                   % background velocity in x
V0 = 0.0;                   % background velocity in y
p0 = 100000;                   % background pressure in Pa
freestreamConditions = [rho0, U0, V0, p0];

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

% determine initial Reynolds and Mach Numbers (WARNING: may be incorrect if inlet x velocity not set correctly or crossflow)
velFS = app.InletVelocityEditField.Value;
L_char = app.ChordLengthEditField.Value;
T0 = (p0 ./ (rho0 .* R));
Re0 = determineRe(velFS, L_char, rho0, T0);
app.ReynoldsNumberEditField.Value = Re0;
M0 = velFS / sqrt(gamma * R * T0);
app.MachNoEditField.Value = M0;

% solver parameters
CFL = app.CFLEditField.Value;                 % CFL number for stability
holdCFL = 0;                                  % initialise variable
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

% ---------------------
% Main time loop
% ---------------------
while t < tEnd
    % compute primitive variables & wave speeds (HLLC)
    Uinitial = cons2prim(Ucon, gamma);
    [rho, u, v, p] = deal(Uinitial(:,:,1), Uinitial(:,:,2), Uinitial(:,:,3), Uinitial(:,:,4)); % assign the primitives to individual variables
    c = sqrt(max(gamma * p ./ rho, 1e-12));
    Smax = max(max(sqrt(u.^2 + v.^2) + c));

    % set CFL -> timestep from CFL
    [CFL, holdCFL] = CFL_Controller(CFL, holdCFL, residuals, iteration-1);
    app.CurrentCFLEditField.Value = CFL;
    dt = CFL * minGridSpacing / Smax;
    if t + dt > tEnd, dt = tEnd - t; end

    % RK2 (Heun): stage 1 (predictor)
    RKstage = 1;
    [RHSout, ~, ~] = RHS(Ucon, R, dx, dy, gamma, visc_alpha, boundaryTypes, boundaryVelocities, freestreamConditions, RKstage, geometryObj, rho0, velFS, refA, dt);
    U1 = Ucon + dt * RHSout;

    % stage 2 
    RKstage = 2;
    [RHSoutS2, iterationResiduals, iterationCoefficients] = RHS(U1, R, dx, dy, gamma, visc_alpha, boundaryTypes, boundaryVelocities,  freestreamConditions, RKstage, geometryObj, rho0, velFS, refA, dt);
    residuals = [residuals; iterationResiduals]; %#ok - residuals vector is relatively small memory wise
    app.C_DEditField.Value = iterationCoefficients(1);

    U2 = 0.5 * (Ucon + U1 + dt * RHSoutS2);
    Ucon = U2; % update Ucon with new values
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
    drawnow limitrate
end
end

% --------------------
% RHS function (finite-volume / fluxes)
% --------------------
function [RHSout, currentResiduals, forceCoefficients] = RHS(Ucon, R, dx, dy, gamma, visc_alpha, boundary, velocities, freestreamConditions, stage, geoObj, rho0, velFS, refA, dt)
    % Ucon: Ny x Nx x 4 array
    % returns dU/dt of same shape

    % reconstruct states at faces using MUSCL (linear reconstruction) with minmod limiter
    % For simplicity handle interior with one-sided for boundaries depending on 'boundary'

    % get primitives
    U = cons2prim(Ucon, gamma);

    % extend arrays for BCs
    Uext = applyWallBCs(U, boundary, velocities, freestreamConditions); % adds one-cell padding around arrays (Ny+2 x Nx+2), and applies BCs there
    Uext_con = prim2cons(Uext, gamma); % get conservative vars for padded arrays

    forcingMask = geoObj.geometryBorder + geoObj.withinGeometry;
    % compute HLLC flux, x-direction (faces perpendicular to x)
    [F, UL_x, UR_x] = HLLC(Uext_con, 'x', gamma, dt, forcingMask);

    % compute HLLC flux, y-direction (faces perpendicular to y)
    [G, UL_y, UR_y] = HLLC(Uext_con, 'y', gamma, dt, forcingMask);

    % compute viscosity
    [mu_du_dx, mu_dv_dy] = diffusionTerm(UL_x, UR_x, UL_y, UR_y, dx, dy, R);

    % fluxes including diffusion/viscosity
    F = F - mu_du_dx;
    G = G - mu_dv_dy;

    % % Apply IBM boundary conditions and assemble RHS via divergence of fluxes including IBM source terms
    % RHSout = applyIBMForcing(Ucon, geoObj, gamma, dt, F, G, dx, dy);

    %% Appply momentum forcing and update velocities
    % remove ghost BC cells
    F_out = F(2:end-1,1:end-1,:);
    F_in = F(2:end-1,2:end,:);
    G_out = G(1:end-1,2:end-1,:);
    G_in = G(2:end,2:end-1,:);

    % determine momentum and density fluxes
    RHSout = - (F_out - F_in)/dx - (G_out - G_in)/dy;  % this is essentially Fxdx/dV - Fdy/dV, but cancelled out

    % Force Coeffs
    if stage == 2
        contVolFaces = boundBox(geoObj);
        forces = calculateForces(contVolFaces, F, G, dx, dy);
        dragCoeff = 2 * forces(1) / (rho0 * velFS ^ 2 * refA);
        forceCoefficients = dragCoeff; % TODO - add other coeffs
    else 
        forceCoefficients = [];
    end

    % residuals
    for idx = 1:4
        currentResiduals(idx) = max(abs(RHSout(:,:,idx)), [], "all"); %#ok - 4 elements is tiny memory wise
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