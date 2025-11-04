function flow_visualiser(app, param) % simData should be simDataStore class
   
    simData = app.simData;

    x = simData.meshGridXY(1,:,1);
    y = simData.meshGridXY(:,1,2);
    t = simData.time;

    switch param
        case "Density"
            dataToDisplay = simData.density;
        case "X Velocity"
            dataToDisplay = simData.xVelocity;
        case "Y Velocity"
            dataToDisplay = simData.yVelocity;
        case "Energy"
            dataToDisplay = simData.energy;
    end
    
    imageSizeX = [x(1) x(end)];
    imageSizeY = [y(1) y(end)];
    imagesc(app.FlowVisualiser, imageSizeX, imageSizeY, dataToDisplay);
    set(app.FlowVisualiser, 'YDir', 'normal');
    title(app.FlowVisualiser, sprintf('%s at t = %.4f', param, t));
    colorbar(app.FlowVisualiser);
    app.SimulationTimeField.Value = t;
    drawnow;
end