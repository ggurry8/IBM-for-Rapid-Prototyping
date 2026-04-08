function geoObj = volumeFractions(x_poly, y_poly,x ,y, geoObj)
% compute volume fractions for a polygon over a mesh grid
% x_poly, y_poly - x and y coordinates of a polygon
% x, y - x, y coordinates of grid points (length N vectors defining grid points in each direction)

[X, Y] = meshgrid(x, y);

dx = x(2) - x(1);
dy = y(2) - y(1);
cellArea = dx * dy;

% Cell corner offsets
ox = [-0.5 0.5 0.5 -0.5] * dx;
oy = [-0.5 -0.5 0.5 0.5] * dy;

% Expand to all cells
Xc = X(:);
Yc = Y(:);

% Each column = one vertex
rectX = Xc + ox;
rectY = Yc + oy;

[in, ~] = inpolygon(X, Y, x_poly, y_poly);

areaFrac = double(in);   % start with 0/1

boundary = ~in & ...
    (inpolygon(X+dx/2, Y, x_poly, y_poly) | ...
     inpolygon(X-dx/2, Y, x_poly, y_poly) | ...
     inpolygon(X, Y+dy/2, x_poly, y_poly) | ...
     inpolygon(X, Y-dy/2, x_poly, y_poly)); % find cells with centres not within boundary but cut by it
boundary = boundary + geoObj.geometryBorder; % add boundary cells with centres fully within as identified by bwperim

idx = find(boundary); % get indices of boundary cells

areaFrac(idx) = arrayfun(@(i) ...
    polyarea_wrapper(rectX(i,:), rectY(i,:), x_poly, y_poly) / cellArea, idx); 
% for each of the boundary cell indices, find the area cut by the boundary

geoObj.forcingMask = reshape(areaFrac, size(X));

end

function A = polyarea_wrapper(rx, ry, px, py)
    p_cell = polyshape(rx, ry);
    p_poly = polyshape(px, py);
    p_int = intersect(p_cell, p_poly);
    if isempty(p_int.Vertices)
        A = 0;
    else
        A = area(p_int);  % area of polyshape object
    end
end