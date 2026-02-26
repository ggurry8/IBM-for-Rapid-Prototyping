% --------------------
% helper: HLLC between two primitive states
% --------------------
% inputs: Uext_con - flow variables array, extended conservative
%         dir - x or y direction of flux computation
%         gamma - ratio of specific heats
% L & R denote the left and right faces of the cells relative to direction of flux computation (dir)
function [F, UL, UR] = HLLC(Uext_con, dir, gamma, dt, forcingMask)

% define minmod limiter anonymous function
minmod = @(a,b) 0.5 * (sign(a) + sign(b)) .* min(abs(a), abs(b));

    % 2nd Order MUSCL reconstruction with minmod limiter
    switch dir
        case 'x' % create left and right side of cell face arrays by removing the outer columns
            U_im_con = Uext_con(:,1:end-2,:); % U_i-1
            U_i_con = Uext_con(:,2:end-1,:);   % U_i
            U_ip_con = Uext_con(:,3:end,:);   % U_i+1

            UL_i = Uext_con(:,1:end-1,:); 
            UR_i = Uext_con(:,2:end,:); 

            dU_i = zeros(size(Uext_con));
            dU_i(:,2:end-1,:) = minmod(U_i_con - U_im_con, U_ip_con - U_i_con); % Ui-1, Ui, Ui+1

            ULcon = UL_i + 0.5 .* dU_i(:,1:end-1,:);
            URcon = UR_i - 0.5 .* dU_i(:,2:end,:);

        case 'y' % create 'left' and 'right' side of cell face arrays by removing the outer rows (take y-direction aligned left -> right)
            U_im_con = Uext_con(1:end-2,:,:); % U_i-1
            U_i_con = Uext_con(2:end-1,:,:);   % U_i
            U_ip_con = Uext_con(3:end,:,:);   % U_i+1

            UL_i = Uext_con(1:end-1,:,:); 
            UR_i = Uext_con(2:end,:,:); 

            dU_i = zeros(size(Uext_con));
            dU_i(2:end-1,:,:) = minmod(U_i_con - U_im_con, U_ip_con - U_i_con); % Ui-1, Ui, Ui+1

            ULcon = UL_i + 0.5 .* dU_i(1:end-1,:,:);
            URcon = UR_i - 0.5 .* dU_i(2:end,:,:);
    end

    UL = cons2prim(ULcon, gamma);
    UR = cons2prim(URcon, gamma);
    
    % create L & R structs containing the cell values either side of the face for the given direction
    ULconStruct = struct('rho',ULcon(:,:,1), 'rho_u',ULcon(:,:,2), 'rho_v',ULcon(:,:,3), 'rho_E',ULcon(:,:,4));
    URconStruct = struct('rho',URcon(:,:,1), 'rho_u',URcon(:,:,2), 'rho_v',URcon(:,:,3), 'rho_E',URcon(:,:,4));
    
    % Calculate fluxes
    F = cellFaceFluxes(ULconStruct, URconStruct, UL, UR, gamma, dir, dt, forcingMask);
end

% --------------------
% helper: 
% --------------------
function F = cellFaceFluxes(ULcon, URcon, ULarray, URarray, gamma, dir, dt, forcingMask)

UL = struct('rho',ULarray(:,:,1), 'u',ULarray(:,:,2), 'v',ULarray(:,:,3), 'p',ULarray(:,:,4));
UR = struct('rho',URarray(:,:,1), 'u',URarray(:,:,2), 'v',URarray(:,:,3), 'p',URarray(:,:,4));

forcingMaskL = zeros(size(UL.u));
forcingMaskR = zeros(size(UL.u));

switch dir % choose correct velocity for given flux and apply forcing
    case 'x'
        desiredxVelocity = 0;
        forcingMaskL(2:end-1,2:end-1) = forcingMask(:,1:end-1,:);
        forcingMaskR(2:end-1,2:end-1) = forcingMask(:,2:end,:);
        xForcingL = (desiredxVelocity - UL.u .* forcingMaskL) ./ dt;
        xForcingR = (desiredxVelocity - UR.u .* forcingMaskR) ./ dt;

        % use u to calculate fluxes and wave speeds
        momentumL = ULcon.rho_u + 0.5 * dt * xForcingL;
        momentumR = URcon.rho_u + 0.5 * dt * xForcingR;
        velocityL = momentumL ./ UL.rho;
        velocityR = momentumR ./ UR.rho;
        
        % pressure term required in second element of flux array only for x dir
        pxL = UL.p;     
        pxR = UR.p;
        pyL = 0;
        pyR = 0;
    case 'y'
        desiredyVelocity = 0;
        forcingMaskL(2:end-1,2:end-1) = forcingMask(1:end-1,:,:);
        forcingMaskR(2:end-1,2:end-1) = forcingMask(2:end,:,:);
        yForcingL = (desiredyVelocity - UL.v .* forcingMaskL) ./ dt;
        yForcingR = (desiredyVelocity - UR.v .* forcingMaskR) ./ dt;

        % use u to calculate fluxes and wave speeds
        momentumL = ULcon.rho_v + 0.5 * dt * yForcingL;
        momentumR = URcon.rho_v + 0.5 * dt * yForcingR;
        velocityL = momentumL ./ UL.rho;
        velocityR = momentumR ./ UR.rho;

        % pressure term required in third element of flux array only for y dir
        pxL = 0;
        pxR = 0;
        pyL = UL.p;    
        pyR = UR.p;
end

% left face flux
FL = ULcon.rho .* velocityL; % Continuity
FL = cat(3, FL, (ULcon.rho_u .* velocityL) + pxL); % x-momentum
FL = cat(3, FL, (ULcon.rho_v .* velocityL) + pyL); % y-momentum
FL = cat(3, FL, ((ULcon.rho_E + UL.p) .* velocityL)); % Energy

% right face flux
FR = URcon.rho .* velocityR;
FR = cat(3, FR, (URcon.rho_u .* velocityR) + pxR);
FR = cat(3, FR, (URcon.rho_v .* velocityR) + pyR);
FR = cat(3, FR, ((URcon.rho_E + UR.p) .* velocityR));

% wave speeds
aL = sqrt(gamma .* UL.p ./ ULcon.rho); % left cell speed of sound
aR = sqrt(gamma .* UR.p ./ URcon.rho); % right cell speed of sound
SL = min(velocityL - aL, velocityR - aR); % left wave 
SR = max(velocityL + aL, velocityR + aR); % right wave 

% construct U* array
SC = (UR.p - UL.p + momentumL .* (SL - velocityL) - momentumR .* (SR - velocityR))./((UL.rho .* (SL - velocityL)) - (UR.rho .* (SR - velocityR))); % contact wave

% ρ* 
rho_star_L = UL.rho .* (SL - velocityL) ./ (SL - SC);
rho_star_R = UR.rho .* (SR - velocityR) ./ (SR - SC);

% u* & v*
switch dir
    case 'x'
        xmom_star_L = rho_star_L .* SC;
        ymom_star_L = rho_star_L .* UL.v;
        xmom_star_R = rho_star_R .* SC;
        ymom_star_R = rho_star_R .* UR.v;
    case 'y'
        xmom_star_L = rho_star_L .* UL.u;
        ymom_star_L = rho_star_L .* SC;
        xmom_star_R = rho_star_R .* UR.u;
        ymom_star_R = rho_star_R .* SC;
end

p_star = UL.p + UL.rho .* (SL - velocityL) .* (SC - velocityL);
% p_star_R = UR.p + UR.rho .* (SR - velocityR) .* (SC - velocityR);

energy_star_L = ((SL - velocityL) .* ULcon.rho_E - UL.p .* velocityL + p_star .* SC) ./ (SL - SC);
energy_star_R = ((SR - velocityR) .* URcon.rho_E - UR.p .* velocityR + p_star .* SC) ./ (SR - SC);

UL_star = cat(3, rho_star_L, xmom_star_L, ymom_star_L, energy_star_L);
UR_star = cat(3, rho_star_R, xmom_star_R, ymom_star_R, energy_star_R);

FL_star = FL + SL .* (UL_star - ULarray);
FR_star = FR + SR .* (UR_star - URarray);

maskFL = SL > 0;
maskFL_star = (SL <= 0) & (SC > 0);
maskFR_star = (SC <= 0) & (SR > 0);
maskFR = SR <= 0;

F = FL .* maskFL + FL_star .* maskFL_star + FR_star .* maskFR_star + FR .* maskFR;

end
