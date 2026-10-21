% --------------------
% helper: HLLC between two primitive states
% --------------------
% inputs: Uext_con - flow variables array, extended conservative
%         dir - x or y direction of flux computation
%         gamma - ratio of specific heats
% L & R denote the left and right faces of the cells relative to direction of flux computation (dir)
function [F, UL, UR] = HLLC(Uext_con, dir, gamma, mu, spacing, geoObj)

% define minmod limiter anonymous function
% minmod = @(a,b) 0.5 * (sign(a) + sign(b)) .* min(abs(a), abs(b));

% Van Albada slope limiter
% vanAlbada = @(a,b) ((a.*b + abs(a.*b)) .* (a + b)) ./ (a.^2 + b.^2 + eps);

% Monotonised Central (MC) slope limiter
% MCLimiter = @(a,b) 0.5*(sign(a)+sign(b)) .* min([2*abs(a); 2*abs(b); 0.5*abs(a+b)], [], 1);

% Venkatakrishnan limiter
epsV = 1e-10;   % small number to avoid division issues
K = 1;        % limiter aggressiveness parameter (typical: 0.3–1)

venkat = @(a,b) ( ...
    ((a + b).^2 + epsV) .* (2*a) .* (2*b) ) ./ ...
    ( (a.^2 + b.^2 + epsV) .* ( (2*a).^2 + (2*b).^2 + epsV ) + K^3 );

    % 2nd Order MUSCL reconstruction with MC limiter
    switch dir
        case 'x' % create left and right side of cell face arrays by removing the outer columns
            U_im_con = Uext_con(:,1:end-2,:); % U_i-1
            U_i_con = Uext_con(:,2:end-1,:);   % U_i
            U_ip_con = Uext_con(:,3:end,:);   % U_i+1

            UL_i = Uext_con(:,1:end-1,:); 
            UR_i = Uext_con(:,2:end,:); 

            dU_i = zeros(size(Uext_con));
            dU_i(:,2:end-1,:) = venkat(U_i_con - U_im_con, U_ip_con - U_i_con); % Ui-1, Ui, Ui+1

            ULcon = UL_i + 0.5 .* dU_i(:,1:end-1,:);
            URcon = UR_i - 0.5 .* dU_i(:,2:end,:);

            UL = cons2prim(ULcon, gamma);
            UR = cons2prim(URcon, gamma);
        case 'y' % create 'left' and 'right' side of cell face arrays by removing the outer rows (take y-direction aligned left -> right)
            U_im_con = Uext_con(1:end-2,:,:); % U_i-1
            U_i_con = Uext_con(2:end-1,:,:);   % U_i
            U_ip_con = Uext_con(3:end,:,:);   % U_i+1

            UL_i = Uext_con(1:end-1,:,:);
            UR_i = Uext_con(2:end,:,:);

            dU_i = zeros(size(Uext_con));
            dU_i(2:end-1,:,:) = venkat(U_i_con - U_im_con, U_ip_con - U_i_con); % Ui-1, Ui, Ui+1

            ULcon = UL_i + 0.5 .* dU_i(1:end-1,:,:);
            URcon = UR_i - 0.5 .* dU_i(2:end,:,:);

            UL = cons2prim(ULcon, gamma);
            UR = cons2prim(URcon, gamma);
    end
    
    % create L & R structs containing the cell values either side of the face for the given direction
    ULconStruct = struct('rho',ULcon(:,:,1), 'rho_u',ULcon(:,:,2), 'rho_v',ULcon(:,:,3), 'rho_E',ULcon(:,:,4));
    URconStruct = struct('rho',URcon(:,:,1), 'rho_u',URcon(:,:,2), 'rho_v',URcon(:,:,3), 'rho_E',URcon(:,:,4));
    
    % Calculate fluxes
    F = cellFaceFluxes(ULconStruct, URconStruct, UL, UR, ULcon, URcon, gamma, dir);
end

% --------------------
% helper: cellFaceFluxes
% --------------------
% calculate HLLC flux across faces
function F = cellFaceFluxes(ULcon_struct, URcon_struct, ULarray, URarray, ULcon, URcon, gamma, dir)

UL = struct('rho',ULarray(:,:,1), 'u',ULarray(:,:,2), 'v',ULarray(:,:,3), 'p',ULarray(:,:,4));
UR = struct('rho',URarray(:,:,1), 'u',URarray(:,:,2), 'v',URarray(:,:,3), 'p',URarray(:,:,4));

% initialise pressure terms, set to 0 to ensure unassigned terms are ignored and don't
% break this function later on
pxL = 0;
pxR = 0;
pyL = 0;
pyR = 0;

switch dir
    case 'x'
        % use u to calculate fluxes
        velocityL = UL.u;
        velocityR = UR.u;
        % pressure term required in second element of flux array only for x dir
        pxL = UL.p;
        pxR = UR.p;
    case 'y'
        % use v to calculate fluxes
        velocityL = UL.v;
        velocityR = UR.v;
        % pressure term required in third element of flux array only for y dir
        pyL = UL.p;
        pyR = UR.p;
end

% left face flux
FL = ULcon_struct.rho .* velocityL; % Continuity
FL = cat(3, FL, (ULcon_struct.rho_u .* velocityL) + pxL); % x-momentum
FL = cat(3, FL, (ULcon_struct.rho_v .* velocityL) + pyL); % y-momentum
FL = cat(3, FL, ((ULcon_struct.rho_E + UL.p) .* velocityL)); % Energy

% right face flux
FR = URcon_struct.rho .* velocityR;
FR = cat(3, FR, (URcon_struct.rho_u .* velocityR) + pxR);
FR = cat(3, FR, (URcon_struct.rho_v .* velocityR) + pyR);
FR = cat(3, FR, ((URcon_struct.rho_E + UR.p) .* velocityR));

% wave speeds
aL = sqrt(gamma .* UL.p ./ ULcon_struct.rho); % left cell speed of sound
aR = sqrt(gamma .* UR.p ./ URcon_struct.rho); % right cell speed of sound
SL = min(velocityL - aL, velocityR - aR); % left wave 
SR = max(velocityL + aL, velocityR + aR); % right wave 

% construct U* array
switch dir
    case 'x'
        SC = (UR.p - UL.p + ULcon_struct.rho_u .* (SL - UL.u) - URcon_struct.rho_u.* (SR - UR.u)) ./ ...
            ((UL.rho .* (SL - UL.u)) - (UR.rho .* (SR - UR.u))); % contact wave
    case 'y'
        SC = (UR.p - UL.p + ULcon_struct.rho_v .* (SL - UL.v) - URcon_struct.rho_v.* (SR - UR.v)) ./ ...
            ((UL.rho .* (SL - UL.v)) - (UR.rho .* (SR - UR.v))); % contact wave
end

ULstar = ones(size(ULarray)) .* UL.rho .* ((SL - velocityL)./(SL - SC));
URstar = ones(size(URarray)) .* UR.rho .* ((SR - velocityR)./(SR - SC));

switch dir
    case 'x'
        ULstar(:,:,2) = ULstar(:,:,2) .* SC;
        ULstar(:,:,3) = ULstar(:,:,3) .* UL.v;
        ULstar(:,:,4) = ULstar(:,:,4) .* (ULcon_struct.rho_E./UL.rho + (SC - velocityL) .* (SC + UL.p ./ (UL.rho .* (SL - velocityL))));

        URstar(:,:,2) = URstar(:,:,2) .* SC;
        URstar(:,:,3) = URstar(:,:,3) .* UR.v;
        URstar(:,:,4) = URstar(:,:,4) .* (URcon_struct.rho_E ./ UR.rho + (SC - velocityR) .* (SC + UR.p ./ (UR.rho .* (SR - velocityR))));
  
    case 'y'
        ULstar(:,:,2) = ULstar(:,:,2) .* UL.u;
        ULstar(:,:,3) = ULstar(:,:,3) .* SC;
        ULstar(:,:,4) = ULstar(:,:,4) .* (ULcon_struct.rho_E./UL.rho + (SC - velocityL) .* (SC + UL.p ./ (UL.rho .* (SL - velocityL))));

        URstar(:,:,2) = URstar(:,:,2) .* UR.u;
        URstar(:,:,3) = URstar(:,:,3) .* SC;
        URstar(:,:,4) = URstar(:,:,4) .* (URcon_struct.rho_E ./ UR.rho + (SC - velocityR) .* (SC + UR.p ./ (UR.rho .* (SR - velocityR))));
end

FL_star = FL + SL .* (ULstar - ULcon);
FR_star = FR + SR .* (URstar - URcon);

maskFL      = SL > 0;
maskFL_star = (SL <= 0) & (SC > 0);
maskFR_star = (SC <= 0) & (SR > 0);
maskFR      = SR <= 0;

F = FL .* maskFL + FL_star .* maskFL_star + FR_star .* maskFR_star + FR .* maskFR;

end
