function param_coords = MapPhysicalToParentQuad(node_coordinates, phys_point)
    % MapPhysicalToParentQuad converts a point from the physical domain
    % to the parent (natural) domain of a quadrilateral element.
    %
    % Inputs:
    %   node_coordinates : 4x2 matrix containing the coordinates of the element's nodes.
    %   phys_point       : 1x2 vector representing the point in the physical domain.
    %
    % Output:
    %   param_coords     : 1x2 vector of parent domain coordinates [xi, eta].

    % Validate inputs
    if size(node_coordinates, 1) ~= 4 || size(node_coordinates, 2) ~= 2
        error('node_coordinates must be a 4x2 matrix.');
    end
    if numel(phys_point) ~= 2
        error('phys_point must be a 2-element vector.');
    end
    phys_point = phys_point(:)'; % Ensure it's a row vector

    % --- Check if physical point is within axis-aligned bounding box of the element ---
    x_coords = node_coordinates(:, 1);
    y_coords = node_coordinates(:, 2);
    x_min = min(x_coords); x_max = max(x_coords);
    y_min = min(y_coords); y_max = max(y_coords);
    
    tol_bbox = 1e-6; % Tolerance for floating-point errors
    if phys_point(1) < x_min - tol_bbox || phys_point(1) > x_max + tol_bbox || ...
       phys_point(2) < y_min - tol_bbox || phys_point(2) > y_max + tol_bbox
        error('Physical point is outside the axis-aligned bounding box of the element.');
    end

    % --- Perform Newton-Raphson inversion ---
    param_coords = localNewtonRaphsonQuad(phys_point, node_coordinates);
    
    % --- Check if computed parent coordinates are within valid range [-1,1]^2 ---
    tol_parent = 1e-6;
    if any(abs(param_coords) > 1 + tol_parent)
        error('Physical point maps outside the parent domain [-1,1]^2.');
    end
end

function param_coords = localNewtonRaphsonQuad(phys_pt, node_coords)
    % localNewtonRaphsonQuad performs the Newton-Raphson iteration to invert
    % the mapping from physical to parent coordinates for a quadrilateral element.
    %
    % Inputs:
    %   phys_pt   : 1x2 physical point.
    %   node_coords : 4x2 matrix of node coordinates.
    %
    % Output:
    %   param_coords : 1x2 parent domain coordinates [xi, eta].

    tol = 1e-6;    % Convergence tolerance
    max_iter = 20; % Maximum iterations
    xi = 0; eta = 0; % Initial guess

    for iter = 1:max_iter
        % Evaluate shape functions at current (xi, eta)
        N = shapeFunctionsQuad(xi, eta); % 1x4 vector
        current_pt = N * node_coords;     % Map to physical coordinates (1x2 vector)
        residual = (phys_pt - current_pt)'; % 2x1 residual vector

        if norm(residual) < tol
            param_coords = [xi, eta];
            return; % Convergence achieved
        end

        % Compute Jacobian matrix (2x2)
        [dN_dxi, dN_deta] = shapeFunctionDerivativesQuad(xi, eta); % 1x4 vectors
        J = zeros(2, 2);
        for i = 1:4
            % Accumulate contributions for each node
            J(1, 1) = J(1, 1) + dN_dxi(i) * node_coords(i, 1);
            J(1, 2) = J(1, 2) + dN_deta(i) * node_coords(i, 1);
            J(2, 1) = J(2, 1) + dN_dxi(i) * node_coords(i, 2);
            J(2, 2) = J(2, 2) + dN_deta(i) * node_coords(i, 2);
        end

        % Solve for update in parent coordinates
        delta = J \ residual;
        xi  = xi  + delta(1);
        eta = eta + delta(2);
    end

    error('Newton-Raphson failed to converge within %d iterations.', max_iter);
end

function N = shapeFunctionsQuad(xi, eta)
    % shapeFunctionsQuad computes the bilinear shape functions for a 4-node
    % quadrilateral element in the natural coordinate system.
    %
    % Inputs:
    %   xi, eta : Natural coordinates (range: -1 to 1)
    %
    % Output:
    %   N : 1x4 vector of shape functions
    N = zeros(1, 4);
    N(1) = 0.25 * (1 - xi) * (1 - eta);
    N(2) = 0.25 * (1 + xi) * (1 - eta);
    N(3) = 0.25 * (1 + xi) * (1 + eta);
    N(4) = 0.25 * (1 - xi) * (1 + eta);
end

function [dN_dxi, dN_deta] = shapeFunctionDerivativesQuad(xi, eta)
    % shapeFunctionDerivativesQuad computes the derivatives of the bilinear shape
    % functions for a 4-node quadrilateral element with respect to xi and eta.
    %
    % Inputs:
    %   xi, eta : Natural coordinates
    %
    % Outputs:
    %   dN_dxi  : 1x4 vector of derivatives with respect to xi
    %   dN_deta : 1x4 vector of derivatives with respect to eta
    dN_dxi = zeros(1, 4);
    dN_deta = zeros(1, 4);

    % Node 1
    dN_dxi(1)  = -0.25 * (1 - eta);
    dN_deta(1) = -0.25 * (1 - xi);
    % Node 2
    dN_dxi(2)  =  0.25 * (1 - eta);
    dN_deta(2) = -0.25 * (1 + xi);
    % Node 3
    dN_dxi(3)  =  0.25 * (1 + eta);
    dN_deta(3) =  0.25 * (1 + xi);
    % Node 4
    dN_dxi(4)  = -0.25 * (1 + eta);
    dN_deta(4) =  0.25 * (1 - xi);
end
