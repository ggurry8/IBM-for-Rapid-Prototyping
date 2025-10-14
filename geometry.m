classdef geometry
    % geometry class for creating a geometry object containing all the
    % relevant information to the object within the flow
    %  

    properties
        XYGeometry
        geometryCells
        geometryBorder
        withinGeometry
        geoGradient
    end

    methods
        function geoObj = geometry(XGrid, YGrid, xGeometry,yGeometry)
            % geometry: Construct geometry class
            %   create the XYGeometry property which contains the X,Y
            %   coordinates of the polygon that represents the geometry
            geoObj.XYGeometry = [xGeometry, yGeometry];
            geoObj.geometryCells = inpolygon(XGrid, YGrid, xGeometry, yGeometry); % logical mask for cells which contain the geometry
            geoObj.geometryBorder = bwperim(geoObj.geometryCells); % find the cells which contain the geometry edge
            geoObj.withinGeometry = logical(geoObj.geometryCells - geoObj.geometryBorder); % find the cells entirely within the geometry edge
        end

        function geoObj = geoInterpret(geoObj)
            % geoInterpret 
            %   
            [gradx, grady] = gradient(double(geoObj.geometryCells));
            lengthOfNormalVector = sqrt(gradx.^2 + grady.^2);
            gradx = gradx ./ (lengthOfNormalVector + eps);
            grady = grady ./ (lengthOfNormalVector + eps);
            geoObj.geoGradient = cat(3, gradx, grady);
        end
    end
end