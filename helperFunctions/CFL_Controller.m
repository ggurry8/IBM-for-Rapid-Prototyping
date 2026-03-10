function CFL = CFL_Controller(CFL0, iteration)
% CFL0 -> current CFL
% holdCFL -> number of remaining iterations to hold at current CFL (for stability)
% CFL -> new CFL
%
% u = U(:,:,2); 
% v = U(:,:,3);
% 
% lambda = (abs(u) + c)/dx + (abs(v) + c)/dy;
    
    % Parameters
    % CFLMin = 0.2;
    CFLMid = 0.3;
    CFLMax = 0.4;

    if iteration == 300
        CFL = CFLMid;
    elseif iteration == 500
        CFL = CFLMax;
    else
        CFL = CFL0;
    end
    % growthFactor = 1.05;
    % shrinkFactor = 0.5;
    
    % Update CFL
    % resNorm = normalize(residuals,'range');
    % if iteration > 25 && holdCFL == 0 
    %     if max(resNorm(iteration, :),[], "all") < 1e-2
    %         CFL = min(CFL0 * growthFactor, CFLMax);
    %         holdCFL = 10;
    %     elseif max(resNorm(iteration, :),[], "all") > 1e-1
    %        CFL = max(CFL0 * shrinkFactor, CFLMin);
    %        holdCFL = 10;
    %     else
    %        CFL = CFL0;
    %     end
    % else
    %     CFL = CFL0;
    %     if holdCFL > 0 
    %         holdCFL = holdCFL - 1;
    %     else
    %         holdCFL = 0;
    %     end
    % end
end