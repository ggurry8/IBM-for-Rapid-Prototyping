% --------------------
% helper: convert conserved -> primitive (vector)
% --------------------
function prim = cons2prim(cons, gamma)

rho = cons(:,:,1);
u = cons(:,:,2)./rho;
v = cons(:,:,3)./rho;
E = cons(:,:,4)./rho;
p = (gamma-1).*rho.*(E - 0.5*(u.^2 + v.^2));
prim = cat(3, rho, u, v, p);

% BETTER EFFICIENCY??
% rho = cons(:,:,1);
% u = cons(:,:,2) .* invRho;
% v = cons(:,:,3) .* invRho;
% E = cons(:,:,4) .* invRho;
% KE = 0.5 *(u.*u + v.*v);
% p = (gamma-1).*rho.*(E - KE);
% prim = cat(3, rho, u, v, p);
%
% % prim = zeros(size(cons), 'like', cons);
% % invRho = 1 ./ rho; % prefer multiplication over division for speed
% % prim(:,:,1) = rho
% % prim(:,:,2) = u
% % prim(:,:,3) = v
% % prim(:,:,4) = p
end