% --------------------
% helper: applyGeometryBCs
% --------------------
function S = applyGeometryBCs(RHSout, Ucon, geometryObj, gamma, dt)

arguments
    RHSout
    Ucon
    geometryObj
    gamma
    dt (1,1)
end

U = cons2prim(Ucon, gamma);
[rho, u, v] = deal(U(:,:,1), U(:,:,2), U(:,:,3));

forcingMask = geometryObj.forcingMask; % + geometryObj.withinGeometry;

if ndims(forcingMask) == 3 % reconstruction method is velocity interpolation % MUST check this doesn't fail if new reconstruction method is introduced
    [desiredxVelocity, desiredyVelocity] = velocityInterpolation(forcingMask, U, geometryObj);
    forcingMask = geometryObj.geometryCells;
else 
    desiredxVelocity = 0;
    desiredyVelocity = 0;
end

xMomForcing = (-RHSout(:,:,2) + (desiredxVelocity - u) / dt) .* rho .* forcingMask;
yMomForcing = (-RHSout(:,:,3) + (desiredyVelocity - v) / dt) .* rho .* forcingMask;
energyForcing = zeros(size(rho)); % let pressure evolve within the boundary
% energyForcing = -RHSout(:,:,4) .* forcingMask; % dont evolve pressure inside the boundary

S = cat(3, zeros(size(rho)), xMomForcing, yMomForcing, energyForcing); % forcing 3D array

end
