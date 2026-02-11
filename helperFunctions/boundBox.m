function boundBoxFaces = boundBox(geometryData)
%% Bounding Box
% find dimensions for bounding box
[xBoundCells, yBoundCells] = find(geometryData.geometryBorder);

domainDims = size(geometryData.geometryBorder);

xmin = min(xBoundCells);
xmax = max(xBoundCells);
ymin = min(yBoundCells);
ymax = max(yBoundCells);

xChar = max([xmax - xmin, ymax - ymin]); % length scale for buffer region

% Max buffer sizes
xBufferDS = round(10 * xChar); % greater buffer downstream for wake effects
xBufferUS = round(5 * xChar); 
yBuffer = round(5 * xChar);

% Ensure bounding box doesn't exceed domain size
if xBufferDS > domainDims(1) - xmax
    xBufferDS = round((domainDims(1) - xmax)/2);
end
if xBufferUS > xmin
    xBufferUS = round(xmin/2);
end
if yBuffer > ymin || yBuffer > domainDims(2) - ymax
    yBuffer = round(min(ymin, domainDims(2) - ymax)/2);
end

% create box
boundBox = false(domainDims); % initialise variable
boundBox(xmin - xBufferDS:xmax + xBufferUS, ymin - yBuffer:ymax + yBuffer) = true; % fill control volume bounding box with true
boundBoxFaces = bwperim(boundBox); % get perimeter of control volume

% cut out aerofoil and show the bounding box (for debugging only)
% boundBox = boundBox - simData.geometryData.geometryBorder;
% imshow(boundBox)


end
