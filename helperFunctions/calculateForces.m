function forces = calculateForces(contVolFaces, F, G, dx, dy)
[xContIndices, yContIndices] = find(contVolFaces);

% find left and right (y-normal) faces
yContIndMin = min(yContIndices);
yContIndMax = max(yContIndices);

% Find top and bottom (x-normal) faces
xContIndMin = min(xContIndices);
xContIndMax = max(xContIndices);

% Find fluxes
xMomFluxXdir = F(:,:,2);
xMomFluxMaskXdir = zeros(size(xMomFluxXdir));
xMomFluxMaskXdir(xContIndMin:xContIndMax, yContIndMin) = 1; % left face
xMomFluxMaskXdir(xContIndMin:xContIndMax, yContIndMax) = -1; % right face
xMomFluxXdirMasked = xMomFluxXdir .* xMomFluxMaskXdir; % apply mask

xMomFluxYdir = G(:,:,2);
xMomFluxMaskYdir = zeros(size(xMomFluxYdir));
xMomFluxMaskYdir(xContIndMin, yContIndMin:yContIndMax) = 1; % top face - TODO might be unnecessary for drag
xMomFluxMaskYdir(xContIndMax, yContIndMin:yContIndMax) = -1; % bottom face - TODO might be unnecessary
xMomFluxYdirMasked = xMomFluxYdir .* xMomFluxMaskYdir; % apply mask

drag = sum(xMomFluxXdirMasked, 'all') * dy + sum(xMomFluxYdirMasked, 'all') * dx;

forces = drag;
end