% --------------------
% helper: velocityInterpolation
% --------------------
% (For further information see: Roberto Verzicco, An Introduction to Immersed Boundary Methods. Cambridge University Press, 2025.)
% 
% Take forcing/interpolation pairs and calculate the desired forcing velocities
% Uses a linear interpolation to find required velocity at A given known velocities at W and B
%
% Notation: W - intersection point of boundary with CFD grid (nodes at cell centres)
%           A - forcing point defined as closest cell centre to W
%           B - interpolation point used to define the forcing required at A
% 
function [u_A, v_A] = velocityInterpolation(forcingMask, U, geoObj)

    XGridUnique = unique(geoObj.XGrid);
    YGridUnique = unique(geoObj.YGrid);

    u = U(:,:,2); % x velocity field
    v = U(:,:,3); % y velocity field

    [isectCoords, forcingCellIdxs, interpCellIdxs] = deal(forcingMask(:,:,1), forcingMask(:,:,2), forcingMask(:,:,3));

    nForcedCells = length(isectCoords);

    internalForcingMask = geoObj.withinGeometry;

    u_W = 0; % boundary X Velocity
    v_W = 0; % boundary Y Velocity

    int_u = 0; % internal velocities
    int_v = 0;

    u_A = zeros(size(u));
    v_A = zeros(size(v));

    for idx = 1:nForcedCells

        forceXIdx = forcingCellIdxs(idx,1);
        forceYIdx = forcingCellIdxs(idx,2);
        interpXIdx = interpCellIdxs(idx,1);
        interpYIdx = interpCellIdxs(idx,2);

        % determine direction of interpolation
        if forceYIdx == interpYIdx % x dir
            isectPoint_W = isectCoords(idx, 1);
            forcingPoint_A = XGridUnique(forceXIdx);
            interpPoint_B = XGridUnique(interpXIdx);
        elseif forceXIdx == interpXIdx % y dir
            isectPoint_W = isectCoords(idx, 2);
            forcingPoint_A = YGridUnique(forceYIdx);
            interpPoint_B = YGridUnique(interpYIdx);
        end

        distance_WA = forcingPoint_A - isectPoint_W;
        distance_WB = interpPoint_B - isectPoint_W;

        phi_B = distance_WA / distance_WB;

        phi_W = 0; % zeroed temporarily % TODO how to determine coefficient ϕw when boundary velocity isn't 0

        u_B = u(interpYIdx, interpXIdx);
        v_B = v(interpYIdx, interpXIdx);
        u_A(forceYIdx, forceXIdx) = phi_B * u_B + phi_W * u_W;
        v_A(forceYIdx, forceXIdx) = phi_B * v_B + phi_W * v_W;

    end
    % add internal (within boundary) to desired velocity matrix
    u_A = u_A + int_u * internalForcingMask;
    v_A = v_A + int_v * internalForcingMask;
end