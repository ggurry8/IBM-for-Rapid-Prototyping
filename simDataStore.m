classdef simDataStore
    %simData - store sim data for app purposes when main loop paused
    %   

    properties
        density
        xVelocity
        yVelocity
        energy
        time
        meshGridXY
    end

    methods
        function simDataStore = simDataStore(Ucon, t)
            %simData - Construct an instance of this class
            %   
            simDataStore.density = Ucon(:,:,1);
            simDataStore.xVelocity = Ucon(:,:,2);
            simDataStore.yVelocity = Ucon(:,:,3);
            simDataStore.energy = Ucon(:,:,4);
            simDataStore.time = t;
        end
        function simDataStore = MeshInfo(simDataStore, X, Y)
            %MeshInfo - store info about the mesh
            simDataStore.meshGridXY = cat(3, X, Y);
        end
    end
end