function Uext = applyWallBCs(U, boundary, velocities, freestreamConditions)
% --------------------
% helper: applyWallBCs (simple padding)
% --------------------
arguments 
    U (:,:,4)
    boundary (4,1) % 4x1 string array containing boundary type data, order is [T; R; B; L]
    velocities % 4x1 string array containing velocities for bcs, order is [T; R; B; L], unused if no velocity inlets are used
    freestreamConditions % 1 x 4 array containing initial: [rho, u, v, pressure]
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
                [rho_e, u_e, v_e, p_e] = velocityInletBC(rho_e, u_e, v_e, p_e, domainWall, velocities, freestreamConditions); % , freestreamConditions
            case 'pressure outlet'
                [rho_e, u_e, v_e, p_e] = pressureOutBC(rho_e, u_e, v_e, p_e, domainWall, freestreamConditions);
            otherwise % outflow is default
                [rho_e, u_e, v_e, p_e] = outflowBC(rho_e, u_e, v_e, p_e, domainWall, freestreamConditions);
        end
    end
    Uext = cat(3, rho_e, u_e, v_e, p_e);
end

% ------------- Outflow/Farfield BC -------------
function [rho_e, u_e, v_e, p_e] = outflowBC(rho_e, u_e, v_e, p_e, domainWall, freestreamConditions)
    % outflow (copy)
    switch domainWall
        case 'left'
            rho_e(:,1)   = freestreamConditions(1);
            p_e(:,1)   = freestreamConditions(4);
            u_e(:,1)   = freestreamConditions(2);
            v_e(:,1)   = freestreamConditions(3);
        case 'right'
            rho_e(:,end) = freestreamConditions(1);
            p_e(:,end) = freestreamConditions(4);
            u_e(:,end) = freestreamConditions(2);
            v_e(:,end) = freestreamConditions(3);
        case 'top'
            rho_e(1,:)   = freestreamConditions(1);
            p_e(1,:)   = freestreamConditions(4);
            u_e(1,:)   = freestreamConditions(2);
            v_e(1,:)   = freestreamConditions(3);
        case 'bottom'
            rho_e(end,:) = freestreamConditions(1);
            p_e(end,:) = freestreamConditions(4);
            u_e(end,:) = freestreamConditions(2);
            v_e(end,:) = freestreamConditions(3);
    end
end

% ------------- Velocity Inlet BC -------------
function [rho_e, u_e, v_e, p_e] = velocityInletBC(rho_e, u_e, v_e, p_e, domainWall, velocities, freestreamConditions) % , freestreamConditions
    % Dirichlet condition - set velocities, use FS density, allow rest to float
    switch domainWall
        case 'left' 
            rho_e(:,1) = freestreamConditions(1); % rho_e(:,2);
            p_e(:,1)   = p_e(:,2);
            u_e(:,1)   = velocities(4);
            v_e(:,1)   = 0; % v_e(:,2);
        case 'right'
            rho_e(:,end) = freestreamConditions(1); % rho_e(:,end-1);
            p_e(:,end) = p_e(:,end-1);
            u_e(:,end) = -velocities(2);
            v_e(:,end) = 0; % v_e(:,end-1);
        case 'top'
            rho_e(1,:)   = freestreamConditions(1); % rho_e(2,:);
            p_e(1,:)   = p_e(2,:);
            u_e(1,:)   = 0; % u_e(2,:);
            v_e(1,:)   = -velocities(1);
        case 'bottom'
            rho_e(end,:) = freestreamConditions(1); % rho_e(end-1,:);
            p_e(end,:) = p_e(end-1,:);
            u_e(end,:) = 0; % u_e(end-1,:);
            v_e(end,:) = velocities(3);
    end
end

% ------------- Pressure Outlet BC -------------
function [rho_e, u_e, v_e, p_e] = pressureOutBC(rho_e, u_e, v_e, p_e, domainWall, freestreamConditions)
% Dirichlet condition - allow velocities, density to float, set pressure
    switch domainWall
        case 'left'
            rho_e(:,1)   = rho_e(:,2);
            p_e(:,1)   = freestreamConditions(4);
            u_e(:,1)   = u_e(:,2);
            v_e(:,1)   = v_e(:,2);
        case 'right'
            % rho_e(:,end) = rho_e(:,end-1);
            % p_e(:,end) = freestreamConditions(4);
            % u_e(:,end) = u_e(:,end-1);
            % v_e(:,end) = v_e(:,end-1);

            % relaxed version
            alpha = 0.05;

            rho_e(:,end) = rho_e(:,end-1);
            u_e(:,end)   = u_e(:,end-1);
            v_e(:,end)   = v_e(:,end-1);
            p_e(:,end) = (1-alpha)*p_e(:,end-1) + alpha*freestreamConditions(4);
        case 'top'
            rho_e(1,:)   = rho_e(2,:);
            p_e(1,:)   = freestreamConditions(4);
            u_e(1,:)   = u_e(2,:);
            v_e(1,:)   = v_e(2,:);
        case 'bottom'
            rho_e(end,:) = rho_e(end-1,:);
            p_e(end,:) = freestreamConditions(4);
            u_e(end,:) = u_e(end-1,:);
            v_e(end,:) = v_e(end-1,:);
    end
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