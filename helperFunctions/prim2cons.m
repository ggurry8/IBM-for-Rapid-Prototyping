% --------------------
% helper: convert primitive -> conserved (vector)
% --------------------
function cons = prim2cons(prim, gamma)
    rho = prim(:,:,1);
    rho_u = prim(:,:,2).*rho;
    rho_v = prim(:,:,3).*rho;
    E = prim(:,:,4)./((gamma-1)*rho) + 0.5*(rho_u.^2 + rho_v.^2);
    rho_E = E.*rho;
    cons = cat(3, rho, rho_u, rho_v, rho_E);
end