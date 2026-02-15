% --------------------
% helper: applyGeometryBCs (simple padding)
% --------------------
function S = applyGeometryBCs(Ucon, geometryObj, gamma, dt)

arguments
    Ucon
    geometryObj
    gamma
    dt (1,1)
end

U = cons2prim(Ucon, gamma);
[rho, u, v] = deal(U(:,:,1), U(:,:,2), U(:,:,3));
forcingMask = geometryObj.geometryBorder + geometryObj.withinGeometry;

desiredxVelocity = 0;
desiredyVelocity = 0;

xForcing = (desiredxVelocity - u .* forcingMask) / dt;
yForcing = (desiredyVelocity - v .* forcingMask) / dt;
EForcing = (xForcing .* u + yForcing .* v);

S = cat(3, zeros(size(rho)), rho .* xForcing, rho .* yForcing, rho .* EForcing); % forcing 3D array

end