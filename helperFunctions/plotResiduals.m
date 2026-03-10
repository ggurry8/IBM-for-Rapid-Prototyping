function plotResiduals(targetAxes, time, residuals)

semilogy(targetAxes, time, normalize(residuals, "range")); % for non normalised residuals - just use residuals
title(targetAxes, "Residuals")
xlabel(targetAxes, "Physical Time (s)")
ylabel(targetAxes, "Normalised Residual")
legend(targetAxes, "Continuity", "X-Momentum", "Y-Momentum", "Energy")
grid(targetAxes, "on")

end