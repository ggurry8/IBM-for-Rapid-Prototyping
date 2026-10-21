% --------------------
% helper: applyGeometryBCs
% --------------------
function S = IBMForcing(U, geometryObj, dt)

arguments
    U
    geometryObj
    dt (1,1)
end

[rho, u, v] = deal(U(:,:,1), U(:,:,2), U(:,:,3));

if ndims(geometryObj.forcingMask) == 3 % reconstruction method is velocity interpolation % MUST check this doesn't fail if new reconstruction method is introduced
    [u_target, v_target] = velocityInterpolation(geometryObj.forcingMask, U, geometryObj);
    forcingMask = geometryObj.geometryCells;
else 
    forcingMask = geometryObj.forcingMask; % + geometryObj.withinGeometry;
    u_target = 0;
    v_target = 0;
end

% vel_target = sqrt(u_target.^2 + v_target.^2)

Sx = ((u_target - u) / dt) .* rho .* forcingMask;
Sy = ((v_target - v) / dt) .* rho .* forcingMask;
SE = (u_target .* Sx + v_target .* Sy) .* forcingMask;

S = cat(3, zeros(size(rho)), Sx, Sy, SE); %zeros(size(rho))); % forcing 3D array

end
