% NOTE: THIS CODE WAS INITIALLY CREATED WITH ChatGPT - however it has been modified significantly from the original version 
% Original GPT chat - https://chatgpt.com/share/68dfd53a-d1d4-8011-9519-0565c66807f9
% The code from this chat is located within the same directory as this file and is called "GPT_original.m"

%% simple_2D_compressible_CFD_solver.m
% A 2D explicit finite-volume compressible N-S solver (ideal gas)
% - Structured Cartesian grid
% - HLL numerical flux
% - Second-order MUSCL reconstruction with minmod limiter
% - RK2 (Heun) time integration
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
numGeoPoints = app.NumberofPointsEditField.Value;
geoAngle = app.AngleEditField.Value;
geoLength = app.LengthEditField.Value;
geoCentre = [Lx * app.CoordinatesofLeadingEdgeX.Value/100, Ly * app.CoordinatesofLeadingEdgeY.Value/100];
switch app.SwitchType.Value
    case "Aerofoil"
        NACA4Code = app.NACA4DigitCodeEditField.Value;
        thickPcent = str2double([NACA4Code(3:4)])/100;
        refA = thickPcent * geoLength;
        [xGeometry, yGeometry] = aerofoil_generation(NACA4Code, geoLength, geoAngle, numGeoPoints);
        xGeometry = xGeometry + geoCentre(1);
        yGeometry = yGeometry + geoCentre(2); % adjust position of leading edge
    case "Shape"
        geoShape = app.ShapeDropDown.Value;
        [xGeometry, yGeometry] = createGeometry(geoShape, geoCentre, geoLength, geoAngle, numGeoPoints);
        xGeometry = xGeometry';
        yGeometry = yGeometry';
        refA = geoLength*2;
end
geometryObj = geometry(X, Y, xGeometry, yGeometry);
% if ~any(geometryObj.withinGeometry,"all")
%     uialert(app.UIFigure, "Geometry may be too large relative to grid size",'Warning: ','Icon','Warning')
% end
geometryObj = geometryObj.geoInterpret;

% Physical parameters
gamma = 1.4;            % ratio of specific heats
R = 287;                % gas constant (not used explicitly if using p=(gamma-1)*rho*e)
mu = 0.0084035;         % enforce viscosity for Re = 1e-4

% Flow initial condition (uniform flow)
rho0 = 1.225;                 % background density
U0 = app.InletVelocityEditField.Value;                   % background velocity in x
V0 = 0.0;                     % background velocity in y
p0 = 101325;                  % background pressure in Pa
freestreamConditions = [rho0, U0, V0, p0];

% determine initial Reynolds and Mach Numbers (WARNING: may be incorrect if inlet x velocity not set correctly or crossflow)
velFS = app.InletVelocityEditField.Value;
L_char = 2 * geoLength;
T0 = (p0 ./ (rho0 .* R));
Re0 = rho0 * U0 * L_char / mu;
app.ReynoldsNumberEditField.Value = Re0;
M0 = velFS / sqrt(gamma * R * T0);
app.MachNoEditField.Value = M0;

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
%CFL = app.CFLEditField.Value;                 % CFL number for stability
tEnd = app.EndTimesEditField.Value;           % final time
t = 0;                                          % initialise variable
visc_alpha = app.ArtificialViscEditField.Value;  % artificial viscosity coeff (set >0 to add stabilization)

boundaryTypes = app.boundaryTypes;           % 4x1 struct of bcs for each wall
boundaryVelocities = app.boundaryVelocities;

residuals = []; % initialise variable
times = []; % initialise variable
coefficientHistory = []; % initialise variable

% diagnostics / plotting
% figure('units','normalized','position',[0.1 0.1 0.6 0.6]);
plotEvery = app.UpdateFreqEditField.Value;
iteration = 1;

% ---------------------
% Main time loop
% ---------------------
while t < tEnd
    % compute primitive variables & flux Jacobian max wave speed
    Uinitial = cons2prim(Ucon, gamma);
    [rho, u, v, p] = deal(Uinitial(:,:,1), Uinitial(:,:,2), Uinitial(:,:,3), Uinitial(:,:,4)); % assign the primitives to individual variables
    c = sqrt(max(gamma * p ./ rho, 1e-12));
    Smax = max(max(sqrt(u.^2 + v.^2) + c));

    % set CFL -> timestep from CFL
    % CFL = CFL_Controller(CFL, iteration-1);
    CFL = app.CFLEditField.Value;                 % CFL number for stability
    app.CurrentCFLEditField.Value = CFL;
    tEnd = app.EndTimesEditField.Value;           % check final time input
    dt = CFL * minGridSpacing / Smax;
    if t + dt > tEnd, dt = tEnd - t; end

    % RK2 (Heun): stage 1 (predictor)
    RKstage = 1;
    [RHSout, ~, ~] = RHS(Ucon, dx, dy, gamma, visc_alpha, boundaryTypes, boundaryVelocities, freestreamConditions, RKstage, geometryObj, rho0, velFS, refA, dt);
    U1 = Ucon + dt * RHSout;

    % stage 2 
    RKstage = 2;
    [RHSoutS2, iterationResiduals, iterationCoefficients] = RHS(U1, dx, dy, gamma, visc_alpha, boundaryTypes, boundaryVelocities,  freestreamConditions, RKstage, geometryObj, rho0, velFS, refA, dt);
    residuals = [residuals; iterationResiduals]; %#ok - residuals vector is relatively small memory wise
    times = [times t]; %#ok - times vector is relatively small memory wise
    
    coefficientHistory = [coefficientHistory; iterationCoefficients(1), iterationCoefficients(2)]; %#ok - coefficients vector is relatively small memory wise
    app.C_DEditField.Value = mean(coefficientHistory(:,1));
    app.C_LEditField.Value = mean(coefficientHistory(:,2));

    U2 = 0.5 * (Ucon + U1 + dt * RHSoutS2);
    Ucon = U2; % update Ucon with new values
    t = t + dt;

    Ufinal = cons2prim(Ucon, gamma);
    app.simData = simDataStore(Ufinal, X, Y, gamma, R, t, geometryObj, residuals, coefficientHistory);

    % update flow visualiser occasionally
    if mod(iteration,plotEvery)==0 || t>=tEnd
        varToDisplay = app.DisplayVarDropDown.Value;
        flow_visualiser(app, varToDisplay);
        plotResiduals(app.ResidualsPlot, times, residuals);
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


if t == tEnd % only export results if end time reached
    % Export results
    savepath = pwd;
    if savepath == 0
        return
    end
    timeAndDateChar = char(datetime('now','Format','yyyy-MM-dd_HH-mm-ss'));
    numberOfPlots = length(app.DisplayVarDropDown.Items);

    % Flow Field Property Plots
    for nPlot = 1:numberOfPlots
        % loop through output variables
        currentDisplayVar = app.DisplayVarDropDown.Items(nPlot);
        flow_visualiser(app, currentDisplayVar);
        % create filename and export
        currentFileNameAndPath = [savepath '\' char(app.DisplayVarDropDown.Items(nPlot)) '_' timeAndDateChar '.png'];
        exportgraphics(app.FlowVisualiser, currentFileNameAndPath);
    end
    % Residuals Plot
    currentFileNameAndPath = [savepath '\Residuals_' timeAndDateChar '.png'];
    exportgraphics(app.ResidualsPlot, currentFileNameAndPath); % TODO may want to allow choice of colourmap from GUI

    % Simulation Data
    currentFileNameAndPath = [savepath '\SimulationData_' timeAndDateChar];
    simdataVariable = app.simData;
    save(currentFileNameAndPath, 'simdataVariable');

    currentFileNameAndPath = [savepath '\coefficientHistory_' timeAndDateChar];
    save(currentFileNameAndPath, 'coefficientHistory');
end

end

% --------------------
% RHS function (finite-volume / fluxes)
% --------------------
function [RHSout, currentResiduals, forceCoefficients] = RHS(Ucon, dx, dy, gamma, visc_alpha, boundary, velocities, freestreamConditions, stage, geoObj, rho0, velFS, refA, dt)

    % reconstruct states at faces using MUSCL (linear reconstruction) with minmod limiter
    % For simplicity handle interior with one-sided for boundaries depending on 'boundary'

    % get primitives
    U = cons2prim(Ucon, gamma);

    % extend arrays for BCs
    Uext = applyWallBCs(U, boundary, velocities, freestreamConditions); % adds one-cell padding around arrays (Ny+2 x Nx+2), and applies BCs there

    % Enforce adiabatic walls for geometry boundaries
    % Uext = applyGeometryAdiabaticBCs(Uext, geoObj, gamma, dt, dx, dy);
    
    % get conservative vars for padded arrays
    Uext_con = prim2cons(Uext, gamma); 

    % Enforce viscosity to enforce desired Reynolds number
    mu = 0.0084035;

    % compute HLL flux, x-direction (faces perpendicular to x)
    F = HLL_flux(Uext_con, 'x', gamma, mu, dx);

    % compute HLL flux, y-direction (faces perpendicular to y)
    G = HLL_flux(Uext_con, 'y', gamma, mu, dy);

    % assemble RHS via divergence of fluxes
    RHS = flow_divergence(F, G, dx, dy);

    % Determine and apply IBM source terms
    S = applyGeometryBCs(RHS, Ucon, geoObj, gamma, dt);
    RHSout = RHS + S;

    % Force Coeffs
    if stage == 2
        contVolFaces = boundBox(geoObj);
        forces = calculateForces(contVolFaces, F, G, dx, dy);
        forces(2) = -sum(S(:,:,2),"all")*(dx*dy); % TODO - fix alternate method
        dragCoeff = 2 * forces(1) / (rho0 * velFS ^ 2 * refA);
        dragCoeffAlt = 2 * forces(2) / (rho0 * velFS ^ 2 * refA);
        forceCoefficients = [dragCoeff dragCoeffAlt] ; % TODO - add lift coeff
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