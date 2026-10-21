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
profile on
%% Problem / grid setup
% Physical parameters
gamma = 1.4;            % ratio of specific heats
R = 287;                % gas constant (not used explicitly if using p=(gamma-1)*rho*e)
mu = 0.0084035;         % enforce viscosity

% Flow initial condition (uniform flow)
rho0 = 1.225;                 % background density
U0 = app.InletVelocityEditField.Value;                   % background velocity in x
V0 = 0.0;                     % background velocity in y
p0 = 101325;                  % background pressure in Pa
freestreamConditions = [rho0, U0, V0, p0];

if app.loadFromFile % if loading from save file use initial conditions from save
    % load domain information
    Lx = max(app.simData.geometryData.XGrid,[],"all"); Ly = max(app.simData.geometryData.YGrid,[],"all");         % domain size
    [Ny, Nx] = size(app.simData.geometryData.XGrid);         % grid cells in x/y
    dx = Lx / Nx; dy = Ly / Ny; % grid spacing (x,y)
    [X,Y] = deal(app.simData.geometryData.XGrid,app.simData.geometryData.YGrid);
    minGridSpacing = min(dx,dy);
    physicalTime = app.simData.time(end);

    % load primitive variables
    rho = app.simData.density;
    u = app.simData.xVelocity;
    v = app.simData.yVelocity;
    p = app.simData.pressure;
    U = cat(3, rho, u, v, p);

    % load geometry data
    geometryObj = app.simData.geometryData;
    [geoLength, refA] = deal(geometryObj.loadUnloadInfo(1), geometryObj.loadUnloadInfo(2));
else
    % initialise domain variables
    Lx = app.XDomainEditField.Value; Ly = app.YDomainEditField.Value;         % domain size
    Nx = app.XGridEditField.Value; Ny = app.YGridEditField.Value;         % grid cells in x/y
    dx = Lx / Nx; dy = Ly / Ny; % grid spacing (x,y)
    x = linspace(dx/2, Lx-dx/2, Nx);
    y = linspace(dy/2, Ly-dy/2, Ny);
    [X,Y] = meshgrid(x,y);
    minGridSpacing = min(dx,dy);
    physicalTime = 0;

    % initialise primitive variables
    rho = rho0 * ones([Ny Nx]);
    u = U0 * ones(size(rho));
    v = V0 * ones(size(rho));
    p = p0 * ones(size(rho));
    U = cat(3, rho, u, v, p);

    % define geometry
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
    app.NewSimulationButton.Text = "Interpreting Geometry...";
    pause(0.1) % allow time for updating GUI
    geoInterpMethod = app.GeoInterpDD;
    geometryObj = geometry(X, Y, xGeometry, yGeometry);
    geometryObj = geometryObj.geoInterpret(geoInterpMethod);
    geometryObj.loadUnloadInfo = [geoLength; refA];
end

% determine initial Reynolds and Mach Numbers (WARNING: may be incorrect if inlet x velocity not set correctly or crossflow)
velFS = app.InletVelocityEditField.Value;
L_char = 2 * geoLength;
T0 = (p0 ./ (rho0 .* R));
Re0 = rho0 * U0 * L_char / mu;
app.ReynoldsNumberEditField.Value = Re0;
M0 = velFS / sqrt(gamma * R * T0);
app.MachNoEditField.Value = M0;

% Initialise conservative variables (Ucon = [rho; rho*u; rho*v; rho*E])
E = p./((gamma-1)) ./ rho + 0.5*(u.^2 + v.^2); % specific total energy
rhoe = rho .* E;
Ucon = zeros(Ny, Nx, 4); % initialise variable
Ucon(:,:,1) = rho;       % conservative variable flow properties
Ucon(:,:,2) = rho .* u;
Ucon(:,:,3) = rho .* v;
Ucon(:,:,4) = rhoe;

% solver parameters
tEnd = app.EndTimesEditField.Value;           % final time
visc_alpha = 0; %app.ArtificialViscEditField.Value;  % artificial viscosity coeff (set >0 to add stabilization)

% external boundaries
boundaryTypes = app.boundaryTypes;           % 4x1 struct of bcs for each wall
boundaryVelocities = app.boundaryVelocities;

residuals = []; % initialise variable
times = []; % initialise variable
coefficientHistory = []; % initialise variable

% diagnostics / plotting / autosave
% figure('units','normalized','position',[0.1 0.1 0.6 0.6]);
autosaveFrequency = app.AutosaveFrequencyField.Value;
autosaveFreqCurrent = autosaveFrequency;
plotEvery = app.UpdateFreqEditField.Value;
iteration = 1;

% store initial field
fieldAvg = U;
fieldAvgPostTransient = zeros([size(X) 4]); % preallocate for update later
app.simData = simDataStore(U, X, Y, gamma, R, times, geometryObj, residuals, coefficientHistory, fieldAvg, fieldAvgPostTransient, freestreamConditions);

app.NewSimulationButton.Text = "Simulation Running...";
app.TerminateSimulationButton.Enable = "on";          
app.PauseSimulationButton.Enable = "on";
pause(0.1) % allow time for updating GUI

%% Main time loop
while physicalTime < tEnd
    drawnow limitrate % update GUI with iteration, time, etc

    % compute primitive variables & flux Jacobian max wave speed
    [rho, u, v, p] = deal(U(:,:,1), U(:,:,2), U(:,:,3), U(:,:,4)); % assign the primitives from the previous final U to individual variables
    c = sqrt(max(gamma * p ./ rho, 1e-12));
    Smax = max(sqrt(u.^2 + v.^2) + c, [], "all");

    % CFL control parameters
    iterationResiduals = []; % initialise variable
    settlingTime = 1; % how long the initial transient lasts (seconds) - also used for starting the running average post transient
    if max(iterationResiduals) < 0.1
        if physicalTime > settlingTime && CFL < 0.29
            app.CFLEditField.Value = 0.3;
        elseif physicalTime > settlingTime && CFL < 0.39
            app.CFLEditField.Value = 0.4;
        end
    elseif max(iterationResiduals) > 0.5 % residual spike
        app.CFLEditField.Value = 0.2; % bring CFL down to "safe" value
    end
    % set CFL -> timestep from CFL
    CFL =  app.CFLEditField.Value;                % CFL number for stability

    % Determine timestep size based on CFL
    tEnd = app.EndTimesEditField.Value;           % check final time input
    acouDT = CFL * minGridSpacing / Smax; % acoustic time step
    viscDT = CFL * minGridSpacing .^2 / (mu/rho0); % viscous time step
    iterationTimestep = min(acouDT,viscDT);
    app.CurrentdtEditField.Value = iterationTimestep;
    if physicalTime + iterationTimestep > tEnd, iterationTimestep = tEnd - physicalTime; end

    % RK2 (Heun): stage 1 (predictor)
    RKstage = 1;
    [RHSout, ~, ~] = RHS(Ucon, dx, dy, gamma, R, visc_alpha, boundaryTypes, boundaryVelocities, freestreamConditions, RKstage, geometryObj, rho0, velFS, refA, iterationTimestep, mu);
    Uconint = Ucon + iterationTimestep * RHSout; % intermediate field approximation (after predictor step)

    % RK2 (Heun): stage 2 (corrector)
    RKstage = 2;
    [RHSoutS2, iterationResiduals, iterationCoefficients] = RHS(Uconint, dx, dy, gamma, R, visc_alpha, boundaryTypes, boundaryVelocities,  freestreamConditions, RKstage, geometryObj, rho0, velFS, refA, iterationTimestep, mu);
    residuals = [residuals; iterationResiduals]; %#ok - residuals vector is relatively small memory wise
    times = [times physicalTime]; %#ok - times vector is relatively small memory wise

    Ucon = 0.5 * (Ucon + Uconint + iterationTimestep * RHSoutS2); % final field approximation (after corrector step)
    physicalTime = physicalTime + iterationTimestep;

    U = cons2prim(Ucon, gamma); % update U
    
    coefficientHistory = [coefficientHistory; iterationCoefficients(1), iterationCoefficients(2)]; %#ok - coefficients vector is relatively small memory wise
    app.C_DEditField.Value = mean(coefficientHistory(:,1));
    app.C_LEditField.Value = mean(coefficientHistory(:,2));

    fieldAvg = fieldAvg + U;
    if physicalTime > settlingTime % begin updating running avg after settling time reached
        fieldAvgPostTransient = fieldAvgPostTransient + U;
    end

    % update flow visualiser occasionally
    if mod(iteration,plotEvery)==0 || physicalTime>=tEnd
        app.simData = simDataStore(U, X, Y, gamma, R, times, geometryObj, residuals, coefficientHistory, fieldAvg, fieldAvgPostTransient, freestreamConditions);
        app.UpdateDisplayButton.Enable = "on"; % can be activated once sim data is available
        varToDisplay = app.DisplayVarDropDown.Value;
        flow_visualiser(app, varToDisplay);
        plotResiduals(app.ResidualsPlot, times, residuals);
    end

    if app.TerminateSimProperty
        % pass sim data out to app
        app.simData = simDataStore(U, X, Y, gamma, R, times, geometryObj, residuals, coefficientHistory, fieldAvg, fieldAvgPostTransient, freestreamConditions);
        app.UpdateDisplayButton.Enable = "on"; % can be activated once sim data is available
        app.TerminateSimProperty = false; % reset property
        break
    elseif app.PauseSimulationButton.Value
        % pass sim data out to app
        app.simData = simDataStore(U, X, Y, gamma, R, times, geometryObj, residuals, coefficientHistory, fieldAvg, fieldAvgPostTransient, freestreamConditions);
        app.UpdateDisplayButton.Enable = "on"; % can be activated once sim data is available
        while app.PauseSimulationButton.Value
            if app.TerminateSimProperty
                app.TerminateSimProperty = false; % reset property
            break
            end
            pause(1)
        end
    end

    % perform end of loop checks
    app.SimulationTimeField.Value = physicalTime;
    iteration = iteration + 1;
    app.IterationNumberField.Value = iteration;

    if autosaveFrequency == 0
        % do nothing
    elseif autosaveFreqCurrent == 0
        app.simData = simDataStore(U, X, Y, gamma, R, times, geometryObj, residuals, coefficientHistory, fieldAvg, fieldAvgPostTransient, freestreamConditions);
        app.UpdateDisplayButton.Enable = "on"; % can be activated once sim data is available
        exportResults(app)
        autosaveFreqCurrent = autosaveFrequency;
    else
        autosaveFreqCurrent = autosaveFreqCurrent - 1;
    end
end

% p = profile('info');
% save mainprofiledata p

if physicalTime == tEnd % only export results if end time reached
    app.simData = simDataStore(U, X, Y, gamma, R, times, geometryObj, residuals, coefficientHistory, fieldAvg, fieldAvgPostTransient, freestreamConditions);
    app.UpdateDisplayButton.Enable = "on"; % can be activated once sim data is available
    exportResults(app)
end
p = profile('info');
save('profileRun','p');
end

%% RHS function (finite-volume / fluxes)
% compute fluxes via HLL flux with 2nd order MUSCL reconstruction + minmod limiter
function [RHS, currentResiduals, forceCoefficients] = RHS(Ucon, dx, dy, gamma, R, visc_alpha, boundary, velocities, freestreamConditions, stage, geoObj, rho0, velFS, refA, dt, mu)
    % recalculate primitives from Ucon
    U = cons2prim(Ucon, gamma);

    % Determine IBM source terms
    S = IBMForcing(U, geoObj, dt);

    % extend arrays for BCs
    Uext = applyWallBCs(U, boundary, velocities, freestreamConditions); % adds one-cell padding around arrays (Ny+2 x Nx+2), and applies BCs there

    % get conservative vars for padded arrays
    Uext_con = prim2cons(Uext, gamma); 

    % compute HLL flux, x-direction (faces perpendicular to x)
    F = HLLC(Uext_con, 'x', gamma, mu, dx, geoObj);

    % compute HLL flux, y-direction (faces perpendicular to y)
    G = HLLC(Uext_con, 'y', gamma, mu, dy, geoObj);

    % Assemble RHS inviscid via divergence of fluxes
    RHSinv = flow_divergence(F, G, dx, dy);

    % compute viscous fluxes and assemble
    [Fv, Gv] = computeViscousFluxes(Uext, dx, dy, mu, R, gamma);
    RHSvisc = flow_divergence(Fv, Gv, dx, dy);
    
    % Assemble full RHS
    RHS = RHSinv - RHSvisc + S;

    % Force Coeffs
    if stage == 2
        forces(1) = -sum(S(:,:,2),"all")*(dx*dy); % drag force is x-momentum forcing
        forces(2) = -sum(S(:,:,3),"all")*(dx*dy); % lift force is y-momentum forcing
        dragCoeff = 2 * forces(1) / (rho0 * velFS ^ 2 * refA);
        liftCoeff = 2 * forces(2) / (rho0 * velFS ^ 2 * refA);
        forceCoefficients = [dragCoeff liftCoeff] ;
    else 
        forceCoefficients = [];
    end

    % residuals
    for idx = 1:4
        currentResiduals(idx) = max(abs(RHSinv(:,:,idx)), [], "all"); %#ok - 4 elements is tiny memory wise
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
        RHS = RHS + visc_alpha * lap;
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
