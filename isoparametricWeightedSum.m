function weightedSum = isoparametricWeightedSum(param_coordinates, node_coordinates)
    xi = param_coordinates(1);
    eta = param_coordinates(2);
    zeta = param_coordinates(3);

    N=shapeFunctionsHexahedral(xi, eta, zeta);
    weightedSum = N * node_coordinates;
end

