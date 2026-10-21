function [x, y] = createGeometry(shape, centre, radius, orientation, points)

orientRads = orientation * pi/180; % convert to radians

switch shape
    case "Square"
        spacing = 2 * radius/(points-1);
        x1 = (-radius:spacing:radius);
        x2 = ones(1,points) * radius;
        x3 = (radius:-spacing:-radius);
        x4 = ones(1,points) * -radius;
        y1 = ones(1,points) * radius;
        y2 = (radius:-spacing:-radius);
        y3 = ones(1,points) * -radius;
        y4 = (-radius:spacing:radius);

        x = [x1 x2 x3 x4];
        y = [y1 y2 y3 y4];

    case "Cylinder"
        spacing = 2 * pi/points;
        x = zeros(1,points);
        y = zeros(1,points);
        for idx = 1:points+1
            x(idx) = centre(1) + radius * cos(spacing * idx);
            y(idx) = centre(2) + radius * sin(spacing * idx);
        end
        return
    case "Triangle"
        
    case "Plate"
        spacing = 2 * radius/(points-1);
        x = (-radius:spacing:radius);
        y = zeros(1,points);
end

% rotate and translate
coords = [x; y];
R = [cos(orientRads) -sin(orientRads); sin(orientRads) cos(orientRads)]; % rotation matrix
rotatedCoords = R * coords;

x = rotatedCoords(1,:) + centre(1);
y = rotatedCoords(2,:) + centre(2);

end