classdef simDataStore
    %simData - store sim data for app purposes when main loop paused
    %   

    properties
        density double
        xVelocity double
        yVelocity double
        pressure double
        gamma double
        R double
        time double
        meshGridXY 
        geometryData geometry
        residuals
    end

    methods
        function simDataStore = simDataStore(U, X, Y, gamma, R, t, geometryObj, residuals)
            % simData - Construct an instance of this class
            simDataStore.density = U(:,:,1);
            simDataStore.xVelocity = U(:,:,2);
            simDataStore.yVelocity = U(:,:,3);
            simDataStore.pressure = U(:,:,4);
            simDataStore.gamma = gamma;
            simDataStore.R = R;
            simDataStore.time = t;
            simDataStore.geometryData = geometryObj;
            simDataStore.meshGridXY = cat(3, X, Y);
            simDataStore.residuals = residuals;
        end
    end
end