% --------------------
% helper: applyGeometryBCs (simple padding)
% --------------------
function RHSout = applyIBMForcing(UconArray, geometryObj, gamma, dt, F, G, dx, dy)
%% Calculate initial forcing
% Ucon = struct('rho',UconArray(:,:,1), 'rho_u',UconArray(:,:,2), 'rho_v',UconArray(:,:,3), 'rho_E',UconArray(:,:,4));

U = cons2prim(UconArray, gamma);
[rho, u, v] = deal(U(:,:,1), U(:,:,2), U(:,:,3));







S = cat(3, zeros(size(rho)), rho .* xForcing, rho .* yForcing, zeros(size(rho))); % forcing 3D array

end