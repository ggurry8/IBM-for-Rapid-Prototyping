% --------------------
% helper: convert conserved -> primitive (vector)
% --------------------
function prim = cons2prim(cons, gamma)
    rho = cons(:,:,1); 
    u = cons(:,:,2)./rho; 
    v = cons(:,:,3)./rho; 
    E = cons(:,:,4)./rho; 
    p = (gamma-1).*rho.*(E - 0.5 * (u.*u + v.*v)); 
    prim = cat(3, rho, u, v, p);
end