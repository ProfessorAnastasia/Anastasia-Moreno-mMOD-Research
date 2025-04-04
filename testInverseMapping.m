%% testInverseMapping.m
% This script sets up a 4-node quadrilateral element (with node coordinates in x and y),
% computes a physical point from known parent coordinates (xi, eta) using the forward mapping,
% and then applies both the Newton-Raphson inverse mapping and the analytic inverse mapping
% (based on Hua's method) to recover the parent coordinates.

clc; clear;

%% Define node coordinates for the quadrilateral element (4x2 matrix)
% Node numbering (assumed as in Hua's paper):
%   Node 1: bottom left, Node 2: bottom right, Node 3: top right, Node 4: top left
node_coords = [
    0.0, 0.0;    % Node 1
    2.0, 0.1;    % Node 2
    1.8, 1.2;    % Node 3
    0.2, 1.0     % Node 4
];

%% Choose known parent coordinates (xi, eta) and compute the corresponding physical point
xi_true  = 0.3;
eta_true = -0.2;

% Compute the bilinear shape functions at (xi, eta)
N = zeros(1,4);
N(1) = 0.25 * (1 - xi_true) * (1 - eta_true);
N(2) = 0.25 * (1 + xi_true) * (1 - eta_true);
N(3) = 0.25 * (1 + xi_true) * (1 + eta_true);
N(4) = 0.25 * (1 - xi_true) * (1 + eta_true);

% Forward mapping: compute the physical point (x, y)
phys_pt = N * node_coords;
fprintf('Forward Mapping:\n');
fprintf('  True parent (xi, eta) = (%.4f, %.4f) \n', xi_true, eta_true);
fprintf('  Computed physical (x, y) = (%.4f, %.4f) \n\n', phys_pt(1), phys_pt(2));

%% Test the Newton-Raphson inverse mapping
try
    param_coords_NR = MapPhysicalToParentQuad(node_coords, phys_pt);
    fprintf('Newton-Raphson Inverse Mapping:\n');
    fprintf('  Recovered parent (xi, eta) = (%.4f, %.4f) \n\n', param_coords_NR(1), param_coords_NR(2));
catch ME
    fprintf('Newton-Raphson Inverse Mapping Error: %s\n\n', ME.message);
end

%% Test the Analytic Inverse Mapping (Hua method)
try
    param_coords_analytic = analyticalInverseMappingQuad(node_coords, phys_pt);
    fprintf('Analytic Inverse Mapping (Hua method):\n');
    fprintf('  Recovered parent (xi, eta) = (%.4f, %.4f) \n\n', param_coords_analytic(1), param_coords_analytic(2));
catch ME
    fprintf('Analytic Inverse Mapping Error: %s\n', ME.message);
end
