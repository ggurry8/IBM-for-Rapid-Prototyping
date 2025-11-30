function residuals = residualsEuler(F, G, dx, dy)

% residuals - 4 x 1 column vector of continuityResidual, xMomentumResidual,  yMomentumResidual, energyResidual

% F, G - fluxes in west-east direction and south-north direction
% dx, dy - grid spacing in x, y ----- THESE MUST BE CONSTANTS -----


FL = F(:,1:end-1,:);
FR = F(:,2:end,:);
GL = G(1:end-1,:,:);
GR = G(2:end,:,:);

% Continuity
RFC = max((FR(:,:,1) - FL(:,:,1))/dx, [], "all");
RGC = max((GR(:,:,1) - GL(:,:,1))/dy, [], 'all');
continuityResidual = max(RFC, RGC);


% X-Momentum
RFX = max((FR(:,:,2) - FL(:,:,2))/dx, [], "all");
RGX = max((GR(:,:,2) - GL(:,:,2))/dy, [], 'all');
xMomentumResidual = max(RFX, RGX);

% Y-Momentum
RFY = max((FR(:,:,3) - FL(:,:,3))/dx, [], "all");
RGY = max((GR(:,:,3) - GL(:,:,3))/dy, [], 'all');
yMomentumResidual = max(RFY, RGY);

% Energy
RFE = max((FR(:,:,4) - FL(:,:,4))/dx, [], "all");
RGE = max((GR(:,:,4) - GL(:,:,4))/dy, [], 'all');
energyResidual = max(RFE, RGE);

residuals = [continuityResidual, xMomentumResidual,  yMomentumResidual, energyResidual];

end