classdef geometry
    % geometry class for creating a geometry object containing all the
    % relevant information to the object within the flow
    %  

    properties
        XYGeometry
        geoCentre
        XGrid
        YGrid
        geometryCells
        geometryBorder
        withinGeometry
        geoGradient
        forcingMask
        SDF
        loadUnloadInfo
    end

    methods
        function geoObj = geometry(XGrid, YGrid, xGeometry, yGeometry)
            % geometry: Construct geometry class
            %   create the XYGeometry property which contains the X,Y
            %   coordinates of the polygon that represents the geometry
            geoObj.XYGeometry = [xGeometry, yGeometry];
            geoObj.XGrid = XGrid;
            geoObj.YGrid = YGrid;
            geoObj.geometryCells = inpolygon(XGrid, YGrid, xGeometry, yGeometry); % logical mask for cells which contain the geometry
            geoObj.geometryBorder = bwperim(geoObj.geometryCells); % find the cells which contain the geometry edge
            geoObj.withinGeometry = logical(geoObj.geometryCells - geoObj.geometryBorder); % find the cells entirely within the geometry edge
            geoObj.geoCentre = [mean(XGrid(geoObj.geometryCells)) mean(YGrid(geoObj.geometryCells))]; % find the geometric centre
        end

        function geoObj = geoInterpret(geoObj, geoInterpMethod)
            % geoInterpret 
            switch geoInterpMethod.Value
                case "None (<1)"
                geoObj.forcingMask = geoObj.geometryCells;
                case "Volume Fraction Weighting (~1)"
                geoObj = volumeFractions(geoObj.XYGeometry(:,1), geoObj.XYGeometry(:,2), geoObj.XGrid(1,:), geoObj.YGrid(:,1), geoObj);
                case "Velocity Interpolation (~2)"
                geoObj = gridInterpret(geoObj.XYGeometry(:,1), geoObj.XYGeometry(:,2), geoObj.XGrid(1,:), geoObj.YGrid(:,1), geoObj);
            end
        end
    end
end