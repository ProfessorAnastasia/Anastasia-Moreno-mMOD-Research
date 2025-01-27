function N = shapeFunctionsHexahedral(xi, eta, zeta)
    N = zeros(1,8);
    N(1) = 0.125 * (1 - xi) * (1 - eta) * (1 - zeta); % N1
    N(2) = 0.125 * (1 + xi) * (1 - eta) * (1 - zeta); % N2
    N(3) = 0.125 * (1 + xi) * (1 + eta) * (1 - zeta); % N3
    N(4) = 0.125 * (1 - xi) * (1 + eta) * (1 - zeta); % N4
    N(5) = 0.125 * (1 - xi) * (1 - eta) * (1 + zeta); % N5
    N(6) = 0.125 * (1 + xi) * (1 - eta) * (1 + zeta); % N6
    N(7) = 0.125 * (1 + xi) * (1 + eta) * (1 + zeta); % N7
    N(8) = 0.125 * (1 - xi) * (1 + eta) * (1 + zeta); % N8
end