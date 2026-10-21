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
            if nnz(geoObj.geometryCells) == 0 % most likely flat plate
                x1 = xGeometry(1); y1 = yGeometry(1);
                x2 = xGeometry(end); y2 = yGeometry(end);
                dx = x2 - x1;
                dy = y2 - y1;

                dist = abs(dy.*XGrid - dx.*YGrid + x2*y1 - y2*x1) ...
                    / sqrt(dx^2 + dy^2);

                dx_grid = abs(XGrid(1,2) - XGrid(1,1));
                dy_grid = abs(YGrid(2,1) - YGrid(1,1));

                tol = 1 * min(dx_grid, dy_grid);

                geoObj.geometryCells = dist < tol;

                % parameter along the line (0 → 1 means within segment)
                t = ((XGrid - x1)*dx + (YGrid - y1)*dy) / (dx^2 + dy^2);

                geoObj.geometryCells = (dist < tol) & (t >= 0) & (t <= 1);
            end
            geoObj.geometryBorder = bwperim(geoObj.geometryCells); % find the cells which contain the geometry edge
            geoObj.withinGeometry = logical(geoObj.geometryCells - geoObj.geometryBorder); % find the cells entirely within the geometry edge
            if nnz(geoObj.withinGeometry) == 0 
                geoObj.withinGeometry = geoObj.geometryCells & ~geoObj.geometryBorder;
            end
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