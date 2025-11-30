% --------------------
% helper: Rusanov flux between two primitive states
% --------------------
% inputs: Uext_con - flow variables array, extended conservative
%         dir - x or y direction of flux computation
%         gamma - ratio of specific heats
% L & R denote the left and right faces of the cells relative to direction of flux computation (dir)
function [F, UL, UR] = rusanov_flux(Uext_con, dir, gamma)

    % 1st order upwind differencing
    switch dir
        case 'x' % create left and right side of cell face arrays by removing the outer columns
            ULcon = Uext_con(:,1:end-1,:);
            URcon = Uext_con(:,2:end,:);
        case 'y' % create left and right side (along y-direction) of cell face arrays by removing the outer rows
            ULcon = Uext_con(1:end-1,:,:);
            URcon = Uext_con(2:end,:,:);
    end

    UL = cons2prim(ULcon, gamma);
    UR = cons2prim(URcon, gamma);
    % create L & R structs containing the cell values either side of the face for the given direction
    L = struct('rho',ULcon(:,:,1), 'rho_u',ULcon(:,:,2), 'rho_v',ULcon(:,:,3), 'rho_E',ULcon(:,:,4), 'p',UL(:,:,4));
    R = struct('rho',URcon(:,:,1), 'rho_u',URcon(:,:,2), 'rho_v',URcon(:,:,3), 'rho_E',URcon(:,:,4), 'p',UR(:,:,4));

    [FL, FR,  smax] = cellFaceFluxes(L, R, UL, UR, gamma, dir);
    
    F = 0.5 .* (FL + FR) - 0.5 .* smax .* (URcon - ULcon);
end

% --------------------
% helper: 
% --------------------
function [FL, FR, smax] = cellFaceFluxes(L, R, UL, UR, gamma, dir)

% initialise pressure terms, set to 0 to ensure unassigned terms are ignored and don't
% break this function later on
pxL = 0;
pxR = 0;
pyL = 0;
pyR = 0;

switch dir
    case 'x'
        % use u to calculate fluxes
        velocityL = UL(:,:,2);
        velocityR = UR(:,:,2);
        % pressure term required in second element of flux array only for x dir
        pxL = UL(:,:,4);    
        pxR = UR(:,:,4);

    case 'y'
        % use v to calculate fluxes
        velocityL = UL(:,:,3);
        velocityR = UR(:,:,3);
        % pressure term required in third element of flux array only for y dir
        pyL = UR(:,:,4);    
        pyR = UR(:,:,4);
end

% left face flux
FL = L.rho .* velocityL;
FL = cat(3, FL, (L.rho_u .* velocityL) + pxL);
FL = cat(3, FL, (L.rho_v .* velocityL) + pyL);
FL = cat(3, FL, ((L.rho_E + L.p) .* velocityL));
% right face flux
FR = R.rho .* velocityR;
FR = cat(3, FR, (R.rho_u .* velocityR) + pxR);
FR = cat(3, FR, (R.rho_v .* velocityR) + pyR);
FR = cat(3, FR, ((R.rho_E + R.p) .* velocityR));

% other values
aL = sqrt(gamma .* L.p ./ L.rho); % left cell speed of sound
aR = sqrt(gamma .* R.p ./ R.rho); % right cell speed of sound
smax = max(abs(velocityL) + aL, abs(velocityR) + aR); % smax is finding the largest speed at the cell face
% TODO - see if 'max' is doing what it is supposed to here i.e. comparing the 2 adjacent cell values and finding the max

end
