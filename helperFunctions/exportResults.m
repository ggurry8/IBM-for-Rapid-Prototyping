function exportResults(app)

    % Export results
    currentDateTime = datetime("now");
    saveDir = fullfile(pwd, "save_" + string(currentDateTime, "yyyy-MM-dd_HH-mm-ss"));
    mkdir(saveDir);
    timeAndDateChar = char(datetime('now','Format','yyyy-MM-dd_HH-mm-ss'));
    numberOfPlots = length(app.DisplayVarDropDown.Items);

    % Flow Field Property Plots
    for nPlot = 1:numberOfPlots
        % loop through output variables
        currentDisplayVar = app.DisplayVarDropDown.Items(nPlot);
        flow_visualiser(app, currentDisplayVar);
        % create filename and export
        currentFileNameAndPath = [saveDir '\' char(app.DisplayVarDropDown.Items(nPlot)) '_' timeAndDateChar '.png'];
        exportgraphics(app.FlowVisualiser, join(currentFileNameAndPath, ''));
    end

    % Residuals Plot
    currentFileNameAndPath = [saveDir '\Residuals_' timeAndDateChar '.png'];
    exportgraphics(app.ResidualsPlot, join(currentFileNameAndPath, '')); 

    % Simulation Data
    currentFileNameAndPath = [saveDir '\SimulationData_' timeAndDateChar];
    simdataVariable = app.simData;
    save(join(currentFileNameAndPath, ''), 'simdataVariable');

end
