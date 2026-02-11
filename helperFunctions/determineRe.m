function Re = determineRe(V, L, rho, T)

S = 110.4;
C1 = 1.458 * 10 ^ -6; 
mu = (C1 * T .^ 1.5) ./ (T + S);

Re = rho .* V .* L ./ mu;

end