function geoObj = gridInterpret(xGeoCoords, yGeoCoords, x, y, geoObj)

nCurvePoints = length(xGeoCoords);

[Xmesh, Ymesh] = meshgrid(x, y);

cellCentresX = Xmesh(:);
cellCentresY = Ymesh(:);

cellCentresXunique = unique(cellCentresX);
cellCentresYunique = unique(cellCentresY);

[Ncells_y, Ncells_x] = size(Xmesh);

%% Find boundary intersections with grid nodes
% add NaN at the end of every grid line
Xmeshlines = [Xmesh; NaN([1,Ncells_x])];
Ymeshlines = [Ymesh; NaN([1,Ncells_x])];
Xmeshlines = [Xmeshlines NaN([Ncells_y+1,1])];
Ymeshlines = [Ymeshlines NaN([Ncells_y+1,1])];

% vertical gridlines, separated with NaNs
Xvgrid = Xmeshlines(:);
Yvgrid = Ymeshlines(:);

% find vertical intersections
[xv_i,yv_i] = polyxpoly(xGeoCoords, yGeoCoords, Xvgrid, Yvgrid);

% horizontal gridlines
Xtmesh = Xmeshlines';
Ytmesh = Ymeshlines';

Xhgrid = Xtmesh(:);
Yhgrid = Ytmesh(:);

% find horizontal intersections
[xh_i,yh_i] = polyxpoly(xGeoCoords, yGeoCoords, Xhgrid, Yhgrid);

intersectX = [xv_i; xh_i];
intersectY = [yv_i; yh_i];

%% Find local gradients at intersection points
[~, cellNeighborXidx] = min(abs(intersectX - cellCentresXunique'), [], 2);    % find nearest X cell centre to the intersection point (for forcing)
[~, cellNeighborYidx] = min(abs(intersectY - cellCentresYunique'), [], 2);    % find nearest Y cell centre to the intersection point (for forcing)

[~, uniqueIndices]= unique([cellCentresXunique(cellNeighborXidx) cellCentresXunique(cellNeighborYidx)], 'rows'); % get unique cell centre indices
nUniqueCellCentres = length(uniqueIndices);

uniqueIntersectX = intersectX(uniqueIndices);
uniqueIntersectY = intersectY(uniqueIndices);

num_neighbors = 4; % number of points before and after to use for gradient estimation
localNormIsect = zeros(nUniqueCellCentres, 2);
globalNormIsect = zeros(nUniqueCellCentres, 2);
interpCellIdxX = zeros(nUniqueCellCentres, 1);
interpCellIdxY = zeros(nUniqueCellCentres, 1);

% find nearest geometry neighbor indices for each unique cell centre neighbor of the intersection points
[~, idx_nearest] = min(sqrt((uniqueIntersectX - xGeoCoords') .^ 2 + (uniqueIntersectY - yGeoCoords') .^ 2),[],2); % Distances from each cell centre neighbor point to every geometry point

for idxUniqueCCNeighbor = 1:nUniqueCellCentres

    idxCCNeighbor = uniqueIndices(idxUniqueCCNeighbor); % get index of cell 

    % find extended neighbor coordinates (relies on geometry points being ordered around the curve loop!!!)
    start_idx = idx_nearest(idxUniqueCCNeighbor) - num_neighbors;
    end_idx = idx_nearest(idxUniqueCCNeighbor) + num_neighbors;

    % shift window if number of neighbors exceeds geometry coordinates
    if start_idx < 1 
        end_idx = end_idx + (1 - start_idx);
        start_idx = 1;
    end
    if end_idx > nCurvePoints
        start_idx = start_idx - (end_idx - nCurvePoints);
        end_idx = nCurvePoints;
    end

    % use neighbor indices to determine local gradient
    neighbor_indices = start_idx:end_idx;
    currentNeighbors = [xGeoCoords(neighbor_indices), yGeoCoords(neighbor_indices)];
    [~, currentNormal] = PCA(currentNeighbors); % use PCA for full rotational invariance in determining local normal
    
    % check normal direction to ensure normal is outward
    testPoint = [uniqueIntersectX(idxUniqueCCNeighbor) + 1e-10*currentNormal(1), uniqueIntersectY(idxUniqueCCNeighbor) + 1e-10*currentNormal(2)]; % step a tiny amount in the normal direction
    if inpolygon(testPoint(1), testPoint(2), xGeoCoords, yGeoCoords) % check if point lies within geometry
        currentNormal = -currentNormal; % flip normal if point lies within geometry
    end
     localNormIsect(idxUniqueCCNeighbor,:) = currentNormal';

    % calculate global vector to use for velocity interpolation
    theta(idxUniqueCCNeighbor) = atan2d(localNormIsect(idxUniqueCCNeighbor,2), localNormIsect(idxUniqueCCNeighbor,1)); %#ok - TODO maybe no need to store theta for each point
    if theta(idxUniqueCCNeighbor) <= 45 && theta(idxUniqueCCNeighbor) > -45 % use global x
        globalNormIsect(idxUniqueCCNeighbor,:) = [1 0];
        dirX = 1; dirY = 0;
    elseif theta(idxUniqueCCNeighbor) <= 135 && theta(idxUniqueCCNeighbor) > 45 % use global y
        globalNormIsect(idxUniqueCCNeighbor,:) = [0 1];
        dirX = 0; dirY = 1;
    elseif theta(idxUniqueCCNeighbor) <= -45 && theta(idxUniqueCCNeighbor) > -135 % use global -y
        globalNormIsect(idxUniqueCCNeighbor,:) = [0 -1];
        dirX = 0; dirY = -1;
    else % use global -x
        globalNormIsect(idxUniqueCCNeighbor,:) = [-1 0];
        dirX = -1; dirY = 0;
    end 
    
    interpCellIdxX(idxUniqueCCNeighbor) = cellNeighborXidx(idxCCNeighbor) + dirX; % x index of cell to use for interpolation
    interpCellIdxY(idxUniqueCCNeighbor) = cellNeighborYidx(idxCCNeighbor) + dirY; % x index of cell to use for interpolation

    % interpCellCoordX(idxUniqueCCNeighbor) = cellCentresXunique(interpCellX);
    % interpCellCoordY(idxUniqueCCNeighbor) = cellCentresYunique(interpCellY);

end

% forcingMask is 3 dim with (:,:,1) containing intersection coordinates, (:,:,2) containing forcing cell indices and (:,:,3) containing the indices for respective cells to use for interpolation
geoObj.forcingMask = cat(3, [uniqueIntersectX, uniqueIntersectY],[cellNeighborXidx(uniqueIndices), cellNeighborYidx(uniqueIndices)], [interpCellIdxX, interpCellIdxY], globalNormIsect); 

% Plot all lines and intersections (for debug)
% figure
% hold on
% mapshow(xGeoCoords, yGeoCoords) % plot geometry
% mapshow(Xvgrid, Yvgrid, "LineWidth", 0.05) % plot vertical lines
% mapshow(Xhgrid, Yhgrid, "LineWidth", 0.05) % plot horizontal lines
% mapshow(uniqueIntersectX,uniqueIntersectY,'DisplayType','point','Marker','x', 'MarkerEdgeColor','g', 'LineWidth', 0.2) % plot intersection points
% mapshow(CCNeighborsUnique(:,1), CCNeighborsUnique(:,2), 'DisplayType','point','Marker','s', 'MarkerEdgeColor','k', 'LineWidth', 0.2) % plot cell centre neighbors to each intersection point
% mapshow(interpCellCoordX, interpCellCoordY, 'DisplayType','point','Marker','s', 'MarkerEdgeColor','m', 'LineWidth', 0.2) % plot cell centre neighbors to each intersection point
% % quiver(CCNeighborsUnique(:,1), CCNeighborsUnique(:,2), 0.1*localNormIsect(1,:)', 0.1*localNormIsect(2,:)', 0, 'r', 'LineWidth', 1.5); % plot local normals from each intersection point on cell centre neighbors
% quiver(CCNeighborsUnique(:,1), CCNeighborsUnique(:,2), 0.1*globalNormIsect(1,:)', 0.1*globalNormIsect(2,:)', 0, 'r', 'LineWidth', 1.5); % plot global normals from each intersection point on cell centre neighbors
% hold off

%% Create Discretised SDF for cell centres
distancesCell = min(sqrt((xGeoCoords' - cellCentresX).^2 + (yGeoCoords' - cellCentresY).^2),[],2);  % distances from each geometric point to each cell centre
% TODO this is dependent on the number of points used to define the geometry, would be better to find the distance to the curve itself
inside = geoObj.withinGeometry(:);
signMask = (inside * -2) + 1; % convert all 1's to -1 and 0's to 1
SDF = distancesCell .* signMask; % apply sign to distance
SDF = reshape(SDF, size(Xmesh));
geoObj.SDF = SDF;

% % plotting stuff for debug
% 
% imageSizeX = [Xmesh(1) Xmesh(end)];
% imageSizeY = [Ymesh(1) Ymesh(end)];
% imagesc(imageSizeX, imageSizeY, SDF)
% colormap hsv
% colorbar
% hold on
% plot(xGeoCoords,yGeoCoords,"g");

end

function [tangent, normal] = PCA(points)
% Principal Component Analysis to find tangent and normal of a set of points
pts_centered = points - mean(points,1);
[~,~,V] = svd(pts_centered, 'econ'); % Singular Value Decomposition
tangent = V(:,1); % first principal component
normal  = V(:,2); % second principal component

end