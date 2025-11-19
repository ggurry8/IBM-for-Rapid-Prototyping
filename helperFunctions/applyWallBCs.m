function Uext = applyWallBCs(U, boundary, velocities)
% --------------------
% helper: applyWallBCs (simple padding)
% --------------------
arguments 
    U (:,:,4)
    boundary (4,1) % 4x1 string array containing boundary type data, order is [T; R; B; L]
    velocities % 4x1 string array containing velocities for bcs, order is [T; R; B; L], unused if no velocity inlets are used
end

    % Sanity check bcs
    velocityInletCheck = any(boundary == 'velocity inlet');
    if velocityInletCheck
        outletCheck = any(boundary == 'outflow');
        if ~outletCheck
            error('Velocity Inlet BC requires at least 1 outflow BC')
        end
    end

    % add one-cell padding around arrays (Ny+2 x Nx+2)
    [Ny,Nx] = size(U(:,:,1));
    % allocate extended arrays
    rho_e = zeros(Ny+2, Nx+2); 
    u_e = zeros(Ny+2, Nx+2); 
    v_e = zeros(Ny+2, Nx+2); 
    p_e = zeros(Ny+2, Nx+2);

    % fill middle with domain values
    rho_e(2:end-1,2:end-1) = U(:,:,1); 
    u_e(2:end-1,2:end-1) = U(:,:,2); 
    v_e(2:end-1,2:end-1) = U(:,:,3);  
    p_e(2:end-1,2:end-1) = U(:,:,4); 

    domainWallList = ["top"; "right"; "bottom"; "left"];
    % apply bcs to each wall
    for wallIdx = 1:4 % [T; R; B; L]
        domainWall = domainWallList(wallIdx);
        switch lower(boundary(wallIdx))
            case 'periodic'
                [rho_e, u_e, v_e, p_e] = periodicBC(rho_e, u_e, v_e, p_e, domainWall);
            case 'wall'
                [rho_e, u_e, v_e, p_e] = wallBC(rho_e, u_e, v_e, p_e, domainWall);
            case 'velocity inlet'
                [rho_e, u_e, v_e, p_e] = velocityInletBC(rho_e, u_e, v_e, p_e, domainWall, velocities);
            otherwise % outflow is default
                [rho_e, u_e, v_e, p_e] = outflowBC(rho_e, u_e, v_e, p_e, domainWall);
        end
    end

    Uext = cat(3, rho_e, u_e, v_e, p_e);
end

% ------------- Periodic BC -------------
function [rho_e, u_e, v_e, p_e] = periodicBC(rho_e, u_e, v_e, p_e, domainWall)
    switch domainWall
        case 'left'
            rho_e(:,1)   = rho_e(:,end-1);
            u_e(:,1)   = u_e(:,end-1);
            v_e(:,1)   = v_e(:,end-1);
            p_e(:,1)   = p_e(:,end-1);
        case 'right'
            rho_e(:,end) = rho_e(:,2);
            u_e(:,end) = u_e(:,2);
            v_e(:,end) = v_e(:,2);
            p_e(:,end) = p_e(:,2);
        case 'top'
            rho_e(1,:)   = rho_e(end-1,:);
            u_e(1,:)   = u_e(end-1,:);
            v_e(1,:)   = v_e(end-1,:);
            p_e(1,:)   = p_e(end-1,:);
        case 'bottom'
            rho_e(end,:) = rho_e(2,:);
            u_e(end,:) = u_e(2,:);
            v_e(end,:) = v_e(2,:);
            p_e(end,:) = p_e(2,:);
    end
end

% ------------- Wall BC -------------
function [rho_e, u_e, v_e, p_e] = wallBC(rho_e, u_e, v_e, p_e, domainWall)
    % reflective wall for velocities, copy density/pressure
    switch domainWall
        case 'left'
            rho_e(:,1)   = rho_e(:,2);
            p_e(:,1)   = p_e(:,2);
            u_e(:,1)   = -u_e(:,2);
            v_e(:,1)   = -v_e(:,2);
        case 'right'
            rho_e(:,end) = rho_e(:,end-1);
            p_e(:,end) = p_e(:,end-1);
            u_e(:,end) = -u_e(:,end-1);
            v_e(:,end) = -v_e(:,end-1);
        case 'top'
            rho_e(1,:)   = rho_e(2,:);
            p_e(1,:)   = p_e(2,:);
            u_e(1,:)   = -u_e(2,:);
            v_e(1,:)   = -v_e(2,:);
        case 'bottom'
            rho_e(end,:) = rho_e(end-1,:);
            p_e(end,:) = p_e(end-1,:);
            u_e(end,:) = -u_e(end-1,:);
            v_e(end,:) = -v_e(end-1,:);
    end
end

% ------------- Outflow BC -------------
function [rho_e, u_e, v_e, p_e] = outflowBC(rho_e, u_e, v_e, p_e, domainWall)
    % outflow (copy)
    switch domainWall
        case 'left'
            rho_e(:,1)   = rho_e(:,2);
            p_e(:,1)   = p_e(:,2);
            u_e(:,1)   = u_e(:,2);
            v_e(:,1)   = v_e(:,2);

        case 'right'
            rho_e(:,end) = rho_e(:,end-1);
            p_e(:,end) = p_e(:,end-1);
            u_e(:,end) = u_e(:,end-1);
            v_e(:,end) = v_e(:,end-1);
        case 'top'
            rho_e(1,:)   = rho_e(2,:);
            p_e(1,:)   = p_e(2,:);
            u_e(1,:)   = u_e(2,:);
            v_e(1,:)   = v_e(2,:);
        case 'bottom'
            rho_e(end,:) = rho_e(end-1,:);
            p_e(end,:) = p_e(end-1,:);
            u_e(end,:) = u_e(end-1,:);
            v_e(end,:) = v_e(end-1,:);
    end
end

function [rho_e, u_e, v_e, p_e] = velocityInletBC(rho_e, u_e, v_e, p_e, domainWall, velocities)
    % all velocity inlet BCs must be specified in normal direction to the wall
    switch domainWall
        case 'left' 
            rho_e(:,1) = rho_e(:,2);
            p_e(:,1)   = p_e(:,2);
            u_e(:,1)   = velocities(4);
            v_e(:,1)   = v_e(:,2);

        case 'right'
            rho_e(:,end) = rho_e(:,end-1);
            p_e(:,end) = p_e(:,end-1);
            u_e(:,end) = -velocities(2);
            v_e(:,end) = v_e(:,end-1);
        case 'top'
            rho_e(1,:)   = rho_e(2,:);
            p_e(1,:)   = p_e(2,:);
            u_e(1,:)   = u_e(2,:);
            v_e(1,:)   = -velocities(1);
        case 'bottom'
            rho_e(end,:) = rho_e(end-1,:);
            p_e(end,:) = p_e(end-1,:);
            u_e(end,:) = u_e(end-1,:);
            v_e(end,:) = velocities(3);
    end
end