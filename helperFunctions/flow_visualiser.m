function flowFieldPlot = flow_visualiser(app, param) % simData should be simDataStore class
   
    % extract variables from simData
    param = string(param);
    simData = app.simData;
    rho = simData.density;
    u = simData.xVelocity; 
    v = simData.yVelocity;    
    p = simData.pressure;
    x = simData.meshGridXY(1,:,1);
    y = simData.meshGridXY(:,1,2);
    t = simData.time;
    geometryData = simData.geometryData;

    switch param % TODO - check all the units are correct here
        case "Density"
            dataToDisplay = rho;
            colourbarLabel = "kg/m^3";
        case "X Velocity"
            dataToDisplay = u;
            colourbarLabel = "m/s";
        case "Y Velocity"
            dataToDisplay = v;
            colourbarLabel = "m/s";
        case "Pressure"
            dataToDisplay = p;
            colourbarLabel = "Pa";
        case {"Velocity_Magnitude", "Mach", "Temperature"} % derived variables
            % derive parameters
            velocityVector = sqrt(simData.xVelocity.^2 + simData.yVelocity.^2);
            gamma = simData.gamma;
            R = simData.R;
            c = sqrt(max((gamma .* simData.pressure) ./ simData.density, 1e-12));
            if param == "Velocity_Magnitude"
                dataToDisplay = velocityVector;
                colourbarLabel = "m/s";
            elseif param == "Mach"
                dataToDisplay = velocityVector./c;
                colourbarLabel = "";
            else % temperature
                dataToDisplay = p./(rho.*R);
                colourbarLabel = "K";
            end
    end
    
    imageSizeX = [x(1) x(end)];
    imageSizeY = [y(1) y(end)];
    flowFieldPlot = imagesc(app.FlowVisualiser, imageSizeX, imageSizeY, dataToDisplay);

    % plot geometry shape if requested
    if app.OverlayGeometryCheckBox.Value
        xGeometry = geometryData.XYGeometry(:,1);
        yGeometry = geometryData.XYGeometry(:,2);
        hold(app.FlowVisualiser, "on");
        plot(app.FlowVisualiser, xGeometry, yGeometry, "r");
        hold(app.FlowVisualiser, "off");
    end

    set(app.FlowVisualiser, 'YDir', 'normal');
    title(app.FlowVisualiser, sprintf('%s', param));
    %title(app.FlowVisualiser, sprintf('%s at t = %.4f', param, t));
    colourbar = colorbar(app.FlowVisualiser);
    colourbar.Label.String = colourbarLabel;
    app.SimulationTimeField.Value = t;
    drawnow;
end