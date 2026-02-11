function Ucon = applyReflectiveBC(Ucon, solid)
% Apply reflective wall boundary conditions on arbitrary solid geometry.
%
% Inputs:
%   Ucon  - Conservative variables [rho, rho*u, rho*v, rho*E]
%   solid - Logical mask (1 = inside solid geometry, 0 = fluid)
%
% Output:
%   Ucon  - Updated conservative variables with wall BCs applied

% get geometry gradients -> approximate wall normals
[nx, ny] = gradient(double(solid));
nlen = sqrt(nx.^2 + ny.^2) + eps;
nx = nx ./ nlen;
ny = ny ./ nlen;

% find boundary cells (fluid cells next to solid)
wallMask = bwperim(solid);

% preallocate ghost arrays
Ughost = Ucon;

% loop through boundary-adjacent fluid cells
for j = 2:size(Ucon,1)-1
    for i = 2:size(Ucon,2)-1
        if wallMask(j,i) && ~solid(j,i)
            % extract primitive vars
            rho = Ucon(j,i,1);
            u = Ucon(j,i,2)/rho;
            v = Ucon(j,i,3)/rho;
            E = Ucon(j,i,4)/rho;

            % local wall normal
            nx_loc = nx(j,i);
            ny_loc = ny(j,i);
            if nx_loc == 0 && ny_loc == 0
                continue
            end

            % decompose velocity into normal/tangential components
            un = u*nx_loc + v*ny_loc;
            ut = -u*ny_loc + v*nx_loc;

            % reflect normal velocity (no penetration)
            un_g = -un;
            ut_g = ut;

            % reconstruct reflected velocity
            u_g = un_g*nx_loc - ut_g*ny_loc;
            v_g = un_g*ny_loc + ut_g*nx_loc;

            % mirror rho and pressure
            rho_g = rho;
            p = (1.4-1)*rho*(E - 0.5*(u^2 + v^2));
            E_g = p/((1.4-1)*rho_g) + 0.5*(u_g^2 + v_g^2);

            % assign ghost state back to conservative vars
            Ughost(j,i,1) = rho_g;
            Ughost(j,i,2) = rho_g * u_g;
            Ughost(j,i,3) = rho_g * v_g;
            Ughost(j,i,4) = rho_g * E_g;
        end
    end
end

% update
Ucon = Ughost;
end
