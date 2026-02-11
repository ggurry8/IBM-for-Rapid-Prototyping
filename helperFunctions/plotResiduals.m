function plotResiduals(targetAxes, iteration, residuals)

semilogy(targetAxes, 1:iteration, normalize(residuals, "range")); % for non normalised residuals - just use residuals
title(targetAxes, "Residuals")
xlabel(targetAxes, "Iteration")
ylabel(targetAxes, "Normalised Residual")
legend(targetAxes, "Continuity", "X-Momentum", "Y-Momentum", "Energy")
grid(targetAxes, "on")

end