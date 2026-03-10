% --------------------
% helper: HLL flux between two primitive states
% --------------------
% inputs: Uext_con - flow variables array, extended conservative
%         dir - x or y direction of flux computation
%         gamma - ratio of specific heats
% L & R denote the left and right faces of the cells relative to direction of flux computation (dir)
function F = HLL_flux(Uext_con, dir, gamma)

% define minmod limiter anonymous function
minmod = @(a,b) 0.5 * (sign(a) + sign(b)) .* min(abs(a), abs(b));

    % 2nd Order MUSCL reconstruction with minmod limiter
    switch dir
        case 'x' % create left and right side of cell face arrays by removing the outer columns
            U_im_con = Uext_con(:,1:end-2,:); % U_i-1
            U_i_con = Uext_con(:,2:end-1,:);   % U_i
            U_ip_con = Uext_con(:,3:end,:);   % U_i+1

            ULcon = Uext_con(:,1:end-1,:); 
            URcon = Uext_con(:,2:end,:); 

            dU_i = zeros(size(Uext_con));

            dU_i(:,2:end-1,:) = minmod(U_i_con - U_im_con, U_ip_con - U_i_con); % Ui-1, Ui, Ui+1

            ULcon = ULcon + 0.5 .* dU_i(:,1:end-1,:);
            URcon = URcon - 0.5 .* dU_i(:,2:end,:);

        case 'y' % create 'left' and 'right' side of cell face arrays by removing the outer rows (take y-direction aligned left -> right)
            U_im_con = Uext_con(1:end-2,:,:); % U_i-1
            U_i_con = Uext_con(2:end-1,:,:);   % U_i
            U_ip_con = Uext_con(3:end,:,:);   % U_i+1

            ULcon = Uext_con(1:end-1,:,:); 
            URcon = Uext_con(2:end,:,:); 

            dU_i = zeros(size(Uext_con));

            dU_i(2:end-1,:,:) = minmod(U_i_con - U_im_con, U_ip_con - U_i_con); % Ui-1, Ui, Ui+1

            ULcon = ULcon + 0.5 .* dU_i(1:end-1,:,:);
            URcon = URcon - 0.5 .* dU_i(2:end,:,:);
    end

    UL = cons2prim(ULcon, gamma);
    UR = cons2prim(URcon, gamma);
    
    % create L & R structs containing the cell values either side of the face for the given direction
    ULstruct = struct('rho',ULcon(:,:,1), 'rho_u',ULcon(:,:,2), 'rho_v',ULcon(:,:,3), 'rho_E',ULcon(:,:,4), 'p',UL(:,:,4));
    URStruct = struct('rho',URcon(:,:,1), 'rho_u',URcon(:,:,2), 'rho_v',URcon(:,:,3), 'rho_E',URcon(:,:,4), 'p',UR(:,:,4));

    % Calculate fluxes
    F = cellFaceFluxes(ULstruct, URStruct, UL, UR, gamma, dir);
end

% --------------------
% helper: 
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
        pyL = UR(:,:,4);    
        pyR = UR(:,:,4);
end

% left face flux
FL = ULcon.rho .* velocityL; % Continuity
FL = cat(3, FL, (ULcon.rho_u .* velocityL) + pxL); % x-momentum
FL = cat(3, FL, (ULcon.rho_v .* velocityL) + pyL); % y-momentum
FL = cat(3, FL, ((ULcon.rho_E + ULcon.p) .* velocityL)); % Energy

% right face flux
FR = URcon.rho .* velocityR;
FR = cat(3, FR, (URcon.rho_u .* velocityR) + pxR);
FR = cat(3, FR, (URcon.rho_v .* velocityR) + pyR);
FR = cat(3, FR, ((URcon.rho_E + URcon.p) .* velocityR));

% wave speeds
aL = sqrt(gamma .* ULcon.p ./ ULcon.rho); % left cell speed of sound
aR = sqrt(gamma .* URcon.p ./ URcon.rho); % right cell speed of sound
SL = min(velocityL - aL, velocityR - aR); % left wave 
SR = max(velocityL + aL, velocityR + aR); % right wave 

Fhll = (SR .* FL - SL .* FR + (SR .* SL) .* (UR - UL)) ./ (SR - SL);

FL_mask = SL > 0;
Fhll_mask = (SL <= 0) & (SR > 0);
FR_mask = SR <= 0;

F = FL_mask .* FL + Fhll_mask .* Fhll + FR_mask .* FR;

end
