% --------------------
% helper: applyGeometryBCs (simple padding)
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
forcingMask = geometryObj.volFractions; % + geometryObj.withinGeometry;

desiredxVelocity = 0;
desiredyVelocity = 0;

xForcing = (-RHSout(:,:,2) + (desiredxVelocity - u) / dt) .* forcingMask;
yForcing = (-RHSout(:,:,3) + (desiredyVelocity - v) / dt) .* forcingMask;

S = cat(3, zeros(size(rho)), rho .* xForcing, rho .* yForcing, zeros(size(rho))); % forcing 3D array

end
