% --------------------
% helper: convert conserved -> primitive (vector)
% --------------------
function prim = cons2prim(Uvec, gamma)
    rho = Uvec(:,:,1);
    u = Uvec(:,:,2)./rho;
    v = Uvec(:,:,3)./rho;
    E = Uvec(:,:,4)./rho;
    p = (gamma-1).*rho.*(E - 0.5*(u.^2 + v.^2));
    prim = cat(3, rho, u, v, p);
end