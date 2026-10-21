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

[~, cellNeighborXidx] = min(abs(intersectX - cellCentresXunique'), [], 2);    % find nearest X cell centre to the intersection point (for forcing)
[~, cellNeighborYidx] = min(abs(intersectY - cellCentresYunique'), [], 2);    % find nearest Y cell centre to the intersection point (for forcing)

%% When 2 intersections match to same gridpoint, find nearest grid node neighbor and filter to remove duplicate
% Step 1: Convert grid-node subscripts into unique row identifiers
neighborPairs = [cellNeighborXidx, cellNeighborYidx];

% Step 2: Find duplicate grid-node assignments
[uniquePairs, ~, groupIdx] = unique(neighborPairs, 'rows');

% Step 3: Preallocate logical mask for intersections to keep
keepMask = false(size(intersectX));

% Step 4: Loop through each group
for k = 1:size(uniquePairs,1)
    currentMembers = find(groupIdx == k);
    
    % If only one intersection assigned to this node, keep it
    if numel(currentMembers) == 1
        keepMask(currentMembers) = true;
    else
        % Grid node coordinates
        gx = cellCentresXunique(uniquePairs(k,1));
        gy = cellCentresYunique(uniquePairs(k,2));
        
        % Compute distances of all duplicate intersections to grid node
        dx = intersectX(currentMembers) - gx;
        dy = intersectY(currentMembers) - gy;
        distSq = dx.^2 + dy.^2;   % squared distance
        
        % Keep only the closest one
        [~, minIdx] = min(distSq);
        keepMask(currentMembers(minIdx)) = true;
    end
end

% Step 5: Filter results
uniqueIntersectX = intersectX(keepMask);
uniqueIntersectY = intersectY(keepMask);
cellNeighborXidx_filtered = cellNeighborXidx(keepMask);
cellNeighborYidx_filtered = cellNeighborYidx(keepMask);

%% Find local gradients at intersection points
nUniqueCellCentres = length(cellNeighborXidx_filtered);

num_neighbors = 4; % number of points before and after to use for gradient estimation
% preallocations
localNormIsect = zeros(nUniqueCellCentres, 2);
globalNormIsect = zeros(nUniqueCellCentres, 2);
cellNeighborXCoord = zeros(nUniqueCellCentres, 1); 
cellNeighborYCoord = zeros(nUniqueCellCentres, 1);
interpCellIdxX = zeros(nUniqueCellCentres, 1); 
interpCellIdxY = zeros(nUniqueCellCentres, 1);
interpCellCoordX = zeros(nUniqueCellCentres, 1); 
interpCellCoordY = zeros(nUniqueCellCentres, 1);
borderMask = zeros(size(Xmesh));

% find nearest geometry neighbor indices for each unique cell centre neighbor of the intersection points
[~, idx_nearest] = min(sqrt((uniqueIntersectX - xGeoCoords') .^ 2 + (uniqueIntersectY - yGeoCoords') .^ 2),[],2); % Distances from each cell centre neighbor point to every geometry point

for idxLoopUniqueCCNeighbor = 1:nUniqueCellCentres

    % find extended neighbor coordinates (relies on geometry points being ordered around the curve loop!!!)
    start_idx = idx_nearest(idxLoopUniqueCCNeighbor) - num_neighbors;
    end_idx = idx_nearest(idxLoopUniqueCCNeighbor) + num_neighbors;

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
    testPoint = [uniqueIntersectX(idxLoopUniqueCCNeighbor) + 1e-10*currentNormal(1), uniqueIntersectY(idxLoopUniqueCCNeighbor) + 1e-10*currentNormal(2)]; % step a tiny amount in the normal direction
    if inpolygon(testPoint(1), testPoint(2), xGeoCoords, yGeoCoords) % check if point lies within geometry
        currentNormal = -currentNormal; % flip normal if point lies within geometry
    end
     localNormIsect(idxLoopUniqueCCNeighbor,:) = currentNormal';

    % calculate global vector to use for velocity interpolation
    theta(idxLoopUniqueCCNeighbor) = atan2d(localNormIsect(idxLoopUniqueCCNeighbor,2), localNormIsect(idxLoopUniqueCCNeighbor,1)); %#ok - TODO maybe no need to store theta for each point
    if theta(idxLoopUniqueCCNeighbor) <= 45 && theta(idxLoopUniqueCCNeighbor) > -45 % use global x
        globalNormIsect(idxLoopUniqueCCNeighbor,:) = [1 0];
        dirX = 1; dirY = 0;
    elseif theta(idxLoopUniqueCCNeighbor) <= 135 && theta(idxLoopUniqueCCNeighbor) > 45 % use global y
        globalNormIsect(idxLoopUniqueCCNeighbor,:) = [0 1];
        dirX = 0; dirY = 1;
    elseif theta(idxLoopUniqueCCNeighbor) <= -45 && theta(idxLoopUniqueCCNeighbor) > -135 % use global -y
        globalNormIsect(idxLoopUniqueCCNeighbor,:) = [0 -1];
        dirX = 0; dirY = -1;
    else % use global -x
        globalNormIsect(idxLoopUniqueCCNeighbor,:) = [-1 0];
        dirX = -1; dirY = 0;
    end 
    
    cellNeighborXCoord(idxLoopUniqueCCNeighbor) = cellCentresXunique(cellNeighborXidx_filtered(idxLoopUniqueCCNeighbor));
    cellNeighborYCoord(idxLoopUniqueCCNeighbor) = cellCentresYunique(cellNeighborYidx_filtered(idxLoopUniqueCCNeighbor));

    interpCellIdxX(idxLoopUniqueCCNeighbor) = cellNeighborXidx_filtered(idxLoopUniqueCCNeighbor) + dirX; % x index of cell to use for interpolation
    interpCellIdxY(idxLoopUniqueCCNeighbor) = cellNeighborYidx_filtered(idxLoopUniqueCCNeighbor) + dirY; % x index of cell to use for interpolation

    interpCellCoordX(idxLoopUniqueCCNeighbor) = cellCentresXunique(interpCellIdxX(idxLoopUniqueCCNeighbor),1);
    interpCellCoordY(idxLoopUniqueCCNeighbor) = cellCentresYunique(interpCellIdxY(idxLoopUniqueCCNeighbor),1);

    borderMask(cellNeighborYidx_filtered(idxLoopUniqueCCNeighbor), cellNeighborXidx_filtered(idxLoopUniqueCCNeighbor)) = 1;
end

% forcingMask is 3 dimensional with:
% (:,:,1) containing intersection coordinates, 
% (:,:,2) containing forcing cell indices, 
% (:,:,3) containing the indices for respective cells to use for interpolation, 
% (:,:,4) containing global normal directions to use for interpolation, and
geoObj.forcingMask = cat(3, [uniqueIntersectX, uniqueIntersectY],[cellNeighborXidx_filtered, cellNeighborYidx_filtered], [interpCellIdxX, interpCellIdxY], globalNormIsect); 
geoObj.geometryBorder = borderMask;
geoObj.geometryCells = imfill(borderMask,"holes");
geoObj.withinGeometry = geoObj.geometryCells - geoObj.geometryBorder;

% Plot all lines and intersections (for debug)
figure
hold on
mapshow(xGeoCoords, yGeoCoords, "Color",'magenta') % plot geometry
mapshow(Xvgrid, Yvgrid, "LineWidth", 0.05) % plot vertical lines
mapshow(Xhgrid, Yhgrid, "LineWidth", 0.05) % plot horizontal lines
mapshow(intersectX,intersectY,'DisplayType','point','Marker','x', 'MarkerEdgeColor','r', 'LineWidth', 0.2) % plot all intersection points - these should be the futher away points in case of duplicate cell centre neighbor
mapshow(uniqueIntersectX,uniqueIntersectY,'DisplayType','point','Marker','x', 'MarkerEdgeColor','g', 'LineWidth', 0.2) % plot filtered (unique) intersection points - these should be the closest points in case of duplicate cell centre neighbor
mapshow(cellNeighborXCoord,cellNeighborYCoord,'DisplayType','point','Marker','o', 'MarkerEdgeColor','k', 'LineWidth', 0.2) % plot all cell centre neighbors
mapshow(interpCellCoordX,interpCellCoordY,'DisplayType','point','Marker','s', 'MarkerEdgeColor','k', 'LineWidth', 0.2) % plot all velocity interpolation cell centres
legend('Geometry','','',"Intersections","Unique Intersections","Forcing Points", "Interpolation Points",'Location','southeast')
hold off

%% Create Discretised SDF for cell centres
distancesCell = min(sqrt((xGeoCoords' - cellCentresX).^2 + (yGeoCoords' - cellCentresY).^2),[],2);  % distances from each geometric point to each cell centre
% TODO this is dependent on the number of points used to define the geometry, would be better to find the distance to the curve itself
inside = geoObj.withinGeometry(:);
signMask = (inside * -2) + 1; % convert all 1's to -1 and 0's to 1
SDF = distancesCell .* signMask; % apply sign to distance
SDF = reshape(SDF, size(Xmesh));
geoObj.SDF = SDF;

% % plotting stuff for SDF debug
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