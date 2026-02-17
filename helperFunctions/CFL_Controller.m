function [CFL, holdCFL] = CFL_Controller(CFL0, holdCFL, residuals, iteration)
% CFL0 -> initial CFL
% holdCFL -> number of remaining iterations to hold at current CFL (for stability)
%
% u = U(:,:,2); 
% v = U(:,:,3);
% 
% lambda = (abs(u) + c)/dx + (abs(v) + c)/dy;
    
    % Parameters
    CFLMin = 0.1;
    CFLMax = 0.55;
    growthFactor = 1.05;
    shrinkFactor = 0.5;
    
    % Update CFL
    resNorm = normalize(residuals,'range');
    if iteration > 25 && holdCFL == 0 
        if max(resNorm(iteration, :),[], "all") < 1e-3
            CFL = min(CFL0 * growthFactor, CFLMax);
            holdCFL = 10;
        elseif max(resNorm(iteration, :),[], "all") > 1e-1
           CFL = max(CFL0 * shrinkFactor, CFLMin);
           holdCFL = 10;
        else
           CFL = CFL0;
        end
    else
        CFL = CFL0;
        if holdCFL > 0 
            holdCFL = holdCFL - 1;
        else
            holdCFL = 0;
        end
    end
end