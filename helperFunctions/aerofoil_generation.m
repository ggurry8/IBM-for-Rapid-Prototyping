% plotting function for NACA 4 digit aerofoils
% see http://airfoiltools.com/airfoil/naca4digit
function [xCoords, yCoords, centroid] = aerofoil_generation(MPXX, chord, aoaDegrees, nPoints)
    arguments
        MPXX (1,4) char
        chord (1,1) double = 1
        aoaDegrees (1,1) double = 0
        nPoints (1,1) double = 100
    end

    nPointsUL = nPoints/2; % number of points for each upper and lower surface

    % ------------ aerofoil geometry -----------------
    % a chord of 1 is assumed initally, and actual chord is applied later 
    B = linspace(0,pi,nPointsUL); % range of points in cosine space
    xLinear = linspace(0,1,nPointsUL); % x coordinates for linear spacing - necessary for determining point of max camber in the cosine data
    xCosine = (1-cos(B))/2; % x coordinates for cosine spacing
    aoaRads = aoaDegrees * pi/180; % convert to radians
    M = str2double(MPXX(1))/100; 
    P = str2double(MPXX(2))/10;
    T = str2double(MPXX(3:4))/100 * chord;
    
    % ----------- calculate camber line ------------------
    yCamber = zeros(1,nPointsUL);
    yCamberGradient = zeros(1,nPointsUL);
    if M ~= 0 && P ~= 0
        idx = 1;
        idxMaxCamber = ceil(nPointsUL * P); % x index of maximum camber
        xLinearMaxCamber = xLinear(idxMaxCamber);
        [~,xCosineMaxCamberidx] = min(abs(xCosine-xLinearMaxCamber));
        while idx <= xCosineMaxCamberidx
            yCamber(idx) = M/P^2 * (2 * P * xCosine(idx) - xCosine(idx) ^ 2); % y coordinates of camber line
            yCamberGradient(idx) = (2*M)/P^2 * (P - xCosine(idx));
            idx = idx + 1;
        end
        for idx2 = idx:nPointsUL
            yCamber(idx2) = M/(1-P)^2 * (1 - 2*P + 2*P * xCosine(idx2) - xCosine(idx2) ^ 2);
            yCamberGradient(idx) = (2*M)/(1-P)^2 * (P - xCosine(idx));
        end
    elseif M == 0 && P ~= 0
        error("Invalid NACA format!")
    end

    % --------------- coordinate creation ------------------------
    
    thickness = zeros(1,nPointsUL);
    yupper = zeros(1,nPointsUL); % y coordinates of upper surface
    ylower = zeros(1,nPointsUL); % y coordinates of lower surface
    xupper = zeros(1,nPointsUL); % y coordinates of upper surface
    xlower = zeros(1,nPointsUL); % y coordinates of lower surface
    for idx = 1:nPointsUL
        thickness(idx) = T/0.2 * (0.2969 * xCosine(idx) ^ 0.5 + -0.126 * xCosine(idx) + -0.3516 * xCosine(idx) ^ 2 ...
            + 0.2843 * xCosine(idx) ^ 3 + -0.1036 * xCosine(idx) ^ 4); % calculate aerofoil thickness at a given x coordinate
        theta = atan(yCamberGradient(idx));
        yupper(idx) = yCamber(idx)  * chord + thickness(idx) * cos(theta);
        ylower(idx) = yCamber(idx)  * chord - thickness(idx) * cos(theta);
        xupper(idx) = (xCosine(idx) * chord - thickness(idx) * sin(theta));
        xlower(idx) = (xCosine(idx) * chord + thickness(idx) * sin(theta));
    end

    % ----------- rotation --------------
    upper = [xupper' yupper'];
    lower = [xlower' ylower'];
    R = [cos(aoaRads) -sin(aoaRads); sin(aoaRads) cos(aoaRads)]; % rotation matrix
    xyRotatedUpper = upper*R;
    xyRotatedLower = lower*R;
    xCoords = [xyRotatedUpper(:,1); flip(xyRotatedLower(:,1))];
    yCoords = [xyRotatedUpper(:,2); flip(xyRotatedLower(:,2))];
    
    % ----------- plotting --------------
    % figure
    % plot(xCoords,yCoords)
    % hold on
    % xlim([0 chord])
    % ylim([-chord/2 chord/2])
    % grid on

end