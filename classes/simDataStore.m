classdef simDataStore
    % simData - store sim data

    properties
        density double
        xVelocity double
        yVelocity double
        pressure double
        fieldAvg double
        fieldAvgPostTransient double
        gamma double
        R double
        time double
        meshGridXY 
        geometryData geometry
        residuals
        coefficients
        freestreamConditions
    end

    methods
        function simDataStore = simDataStore(U, X, Y, gamma, R, t, geometryObj, residuals, coefficients, fieldAvg, fieldAvgPostTransient, freestreamConditions)
            % simData - Construct an instance of this class
            simDataStore.density = U(:,:,1);
            simDataStore.xVelocity = U(:,:,2);
            simDataStore.yVelocity = U(:,:,3);
            simDataStore.pressure = U(:,:,4);
            simDataStore.fieldAvg = fieldAvg; % store the initial field - running average will be updated each iteration
            simDataStore.fieldAvgPostTransient = fieldAvgPostTransient; % store zeros - running average will be updated each iteration after settling time reached
            simDataStore.gamma = gamma;
            simDataStore.R = R;
            simDataStore.time = t;
            simDataStore.geometryData = geometryObj;
            simDataStore.meshGridXY = cat(3, X, Y);
            simDataStore.residuals = residuals;
            simDataStore.coefficients = coefficients;
            simDataStore.freestreamConditions = freestreamConditions;
        end
    end
end