% --------------------
% helper: HLL flux between two primitive states
% --------------------
% inputs: Uext_con - flow variables array, extended conservative
%         dir - x or y direction of flux computation
%         gamma - ratio of specific heats
% L & R denote the left and right faces of the cells relative to direction of flux computation (dir)
function F = HLL_flux(Uext_con, dir, gamma, mu, spacing)
% define minmod limiter anonymous function
minmod = @(a,b) 0.5 * (sign(a) + sign(b)) .* min(abs(a), abs(b));

% 2nd Order MUSCL reconstruction with minmod limiter
    switch dir
        case 'x' % create left and right side of cell face arrays by removing the outer columns
            diffL = Uext_con(:,2:end-1,:) - Uext_con(:,1:end-2,:);
            diffR = Uext_con(:,3:end,:) - Uext_con(:,2:end-1,:);
    
            dU_i = zeros(size(Uext_con));
            dU_i(:,2:end-1,:) = minmod(diffL, diffR); % Ui-1, Ui, Ui+1
    
            ULcon = Uext_con(:,1:end-1,:) + 0.5 .* dU_i(:,1:end-1,:);
            URcon = Uext_con(:,2:end,:) - 0.5 .* dU_i(:,2:end,:);
   
            UL = cons2prim(ULcon, gamma);
            UR = cons2prim(URcon, gamma);
            velL = UL(:,:,2);
            velR = UR(:,:,2);
        case 'y' % create 'left' and 'right' side of cell face arrays by removing the outer rows (take y-direction aligned left -> right)
            diffL = Uext_con(2:end-1,:,:) - Uext_con(1:end-2,:,:);
            diffR = Uext_con(3:end,:,:) - Uext_con(2:end-1,:,:);
    
            dU_i = zeros(size(Uext_con));
            dU_i(2:end-1,:,:) = minmod(diffL, diffR); % Ui-1, Ui, Ui+1
    
            ULcon = Uext_con(1:end-1,:,:) + 0.5 .* dU_i(1:end-1,:,:);
            URcon = Uext_con(2:end,:,:) - 0.5 .* dU_i(2:end,:,:);
    
            UL = cons2prim(ULcon, gamma);
            UR = cons2prim(URcon, gamma);
            velL = UL(:,:,3);
            velR = UR(:,:,3);
    end

    Fvisc = mu * (velR - velL)/spacing; % mu * du/dx
    
    % create L & R structs containing the cell values either side of the face for the given direction
    ULstruct = struct('rho',ULcon(:,:,1), 'rho_u',ULcon(:,:,2), 'rho_v',ULcon(:,:,3), 'rho_E',ULcon(:,:,4));
    URStruct = struct('rho',URcon(:,:,1), 'rho_u',URcon(:,:,2), 'rho_v',URcon(:,:,3), 'rho_E',URcon(:,:,4));

    % Calculate fluxes
    Finv = cellFaceFluxes(ULstruct, URStruct, UL, UR, gamma, dir);

    F = Finv - Fvisc;
end

% --------------------
% helper: cellFaceFluxes
% --------------------
function F = cellFaceFluxes(ULcon, URcon, UL, UR, gamma, dir)

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
        pyL = UL(:,:,4);
        pyR = UR(:,:,4);
end

% wave speeds
aL = sqrt(gamma .* (pxL + pyL) ./ ULcon.rho); % left cell speed of sound
aR = sqrt(gamma .* (pxL + pyL) ./ URcon.rho); % right cell speed of sound
SL = min(velocityL - aL, velocityR - aR); % left wave 
SR = max(velocityL + aL, velocityR + aR); % right wave 

% logic to determine correct wave to use
FL_mask = SL > 0;
Fhll_mask = (SL <= 0) & (SR > 0);
FR_mask = SR <= 0;

% preallocate FL, FR
FL = zeros([size(ULcon.rho), 4]);
FR = zeros([size(URcon.rho), 4]);

% left face flux
FL(:,:,1) = ULcon.rho .* velocityL; % Continuity
FL(:,:,2) = (ULcon.rho_u .* velocityL) + pxL; % x-momentum
FL(:,:,3) = (ULcon.rho_v .* velocityL) + pyL; % y-momentum
FL(:,:,4) = (ULcon.rho_E + pxL + pyL) .* velocityL; % Energy % note: in (pxL + pyL) either pxL or pyL will always be 0

% right face flux
FR(:,:,1) = URcon.rho .* velocityR; % Continuity
FR(:,:,2) = (URcon.rho_u .* velocityR) + pxR; % x-momentum
FR(:,:,3) = (URcon.rho_v .* velocityR) + pyR; % y-momentum
FR(:,:,4) = (URcon.rho_E + (pxL + pyL)) .* velocityR; % Energy % note: in (pxL + pyL) either pxL or pyL will always be 0

Fhll = (SR .* FL - SL .* FR + (SR .* SL) .* (UR - UL)) ./ (SR - SL); % HLL flux

F = FL_mask .* FL + Fhll_mask .* Fhll + FR_mask .* FR;

end
