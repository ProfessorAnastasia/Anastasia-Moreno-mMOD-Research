function param_coords = MapPhysicalToParent(node_coordinates, phys_point)
    % Convert a point from the physical domain to the parent domain of a hexahedral element.
    %
    % Inputs:
    %   node_coordinates : 8x3 matrix containing the coordinates of the element's nodes.
    %   phys_point       : 1x3 vector representing the point in the physical domain.
    %
    % Output:
    %   param_coords     : 1x3 vector of parent domain coordinates [xi, eta, zeta].
    
    % Validate inputs
    if size(node_coordinates, 1) ~= 8 || size(node_coordinates, 2) ~= 3
        error('node_coordinates must be an 8x3 matrix.');
    end
    if numel(phys_point) ~= 3
        error('phys_point must be a 3-element vector.');
    end
    phys_point = phys_point(:)'; % Ensure it's a row vector
    
    % --- Check if physical point is within axis-aligned bounding box of the element ---
    x_coords = node_coordinates(:,1);
    y_coords = node_coordinates(:,2);
    z_coords = node_coordinates(:,3);
    x_min = min(x_coords); x_max = max(x_coords);
    y_min = min(y_coords); y_max = max(y_coords);
    z_min = min(z_coords); z_max = max(z_coords);
    
    tolerance_bbox = 1e-6; % Small tolerance for floating-point errors
    if phys_point(1) < x_min - tolerance_bbox || phys_point(1) > x_max + tolerance_bbox || ...
       phys_point(2) < y_min - tolerance_bbox || phys_point(2) > y_max + tolerance_bbox || ...
       phys_point(3) < z_min - tolerance_bbox || phys_point(3) > z_max + tolerance_bbox
        error('Physical point is outside the axis-aligned bounding box of the element.');
    end
    
    % --- Perform Newton-Raphson inversion ---
    param_coords = localNewtonRaphson(phys_point, node_coordinates);
    
    % --- Check if computed parent coordinates are within valid range [-1,1]^3 ---
    tolerance_parent = 1e-6;
    if any(abs(param_coords) > 1 + tolerance_parent)
        error('Physical point maps outside the parent domain [-1,1]^3.');
    end
end

function param_coords = localNewtonRaphson(phys_pt, node_coords)
    % Newton-Raphson method to find parent coordinates from physical point.
    tol = 1e-6;    % Convergence tolerance
    max_iter = 20; % Maximum iterations
    xi = 0; eta = 0; zeta = 0; % Initial guess
    
    converged = false;
    for iter = 1:max_iter
        % Evaluate shape functions
        N = shapeFunctionsHexahedral(xi, eta, zeta);
        current_pt = N * node_coords; % Current physical point
        residual = (phys_pt - current_pt)';
        
        if norm(residual) < tol
            converged = true;
            break; % Converged
        end
        
        % Compute Jacobian matrix
        [dN_dxi, dN_deta, dN_dzeta] = ShapeFunctionDerivatives(xi, eta, zeta);
        J = zeros(3, 3);
        for i = 1:8
            J(:,1) = J(:,1) + dN_dxi(i) * node_coords(i,:)';
            J(:,2) = J(:,2) + dN_deta(i) * node_coords(i,:)';
            J(:,3) = J(:,3) + dN_dzeta(i) * node_coords(i,:)';
        end
        
        % Solve for parameter update
        delta = J \ residual;
        xi = xi + delta(1);
        eta = eta + delta(2);
        zeta = zeta + delta(3);
    end
    
    if ~converged
        error('Newton-Raphson failed to converge within %d iterations.', max_iter);
    end
    
    param_coords = [xi, eta, zeta];
end

function N = shapeFunctionsHexahedral(xi, eta, zeta)
    % Trilinear shape functions for an 8-node hexahedral element.
    N = zeros(1, 8);
    N(1) = 0.125 * (1 - xi) * (1 - eta) * (1 - zeta);
    N(2) = 0.125 * (1 + xi) * (1 - eta) * (1 - zeta);
    N(3) = 0.125 * (1 + xi) * (1 + eta) * (1 - zeta);
    N(4) = 0.125 * (1 - xi) * (1 + eta) * (1 - zeta);
    N(5) = 0.125 * (1 - xi) * (1 - eta) * (1 + zeta);
    N(6) = 0.125 * (1 + xi) * (1 - eta) * (1 + zeta);
    N(7) = 0.125 * (1 + xi) * (1 + eta) * (1 + zeta);
    N(8) = 0.125 * (1 - xi) * (1 + eta) * (1 + zeta);
end

function [dN_dxi, dN_deta, dN_dzeta] = ShapeFunctionDerivatives(xi, eta, zeta)
    % Derivatives of the shape functions with respect to parent coordinates.
    dN_dxi = zeros(1, 8);
    dN_deta = zeros(1, 8);
    dN_dzeta = zeros(1, 8);
    
    % Node 1
    dN_dxi(1) = -0.125 * (1 - eta) * (1 - zeta);
    dN_deta(1) = -0.125 * (1 - xi) * (1 - zeta);
    dN_dzeta(1) = -0.125 * (1 - xi) * (1 - eta);
    
    % Node 2
    dN_dxi(2) = 0.125 * (1 - eta) * (1 - zeta);
    dN_deta(2) = -0.125 * (1 + xi) * (1 - zeta);
    dN_dzeta(2) = -0.125 * (1 + xi) * (1 - eta);
    
    % Node 3
    dN_dxi(3) = 0.125 * (1 + eta) * (1 - zeta);
    dN_deta(3) = 0.125 * (1 + xi) * (1 - zeta);
    dN_dzeta(3) = -0.125 * (1 + xi) * (1 + eta);
    
    % Node 4
    dN_dxi(4) = -0.125 * (1 + eta) * (1 - zeta);
    dN_deta(4) = 0.125 * (1 - xi) * (1 - zeta);
    dN_dzeta(4) = -0.125 * (1 - xi) * (1 + eta);
    
    % Node 5
    dN_dxi(5) = -0.125 * (1 - eta) * (1 + zeta);
    dN_deta(5) = -0.125 * (1 - xi) * (1 + zeta);
    dN_dzeta(5) = 0.125 * (1 - xi) * (1 - eta);
    
    % Node 6
    dN_dxi(6) = 0.125 * (1 - eta) * (1 + zeta);
    dN_deta(6) = -0.125 * (1 + xi) * (1 + zeta);
    dN_dzeta(6) = 0.125 * (1 + xi) * (1 - eta);
    
    % Node 7
    dN_dxi(7) = 0.125 * (1 + eta) * (1 + zeta);
    dN_deta(7) = 0.125 * (1 + xi) * (1 + zeta);
    dN_dzeta(7) = 0.125 * (1 + xi) * (1 + eta);
    
    % Node 8
    dN_dxi(8) = -0.125 * (1 + eta) * (1 + zeta);
    dN_deta(8) = 0.125 * (1 - xi) * (1 + zeta);
    dN_dzeta(8) = 0.125 * (1 - xi) * (1 + eta);
end