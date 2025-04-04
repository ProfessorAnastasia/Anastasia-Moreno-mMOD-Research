function param_coords = analyticalInverseMappingQuad(node_coordinates, phys_point)
% analyticalInverseMappingQuad analytically computes the parent coordinates
% (xi, eta) from the physical coordinates (x,y) for a 4-node quadrilateral
% isoparametric element based on Hua's inverse transformation method.
%
% Inputs:
%   node_coordinates - 4x2 matrix with rows [xi, yi] for nodes 1 to 4. The node
%                 numbering is assumed as:
%                 Node 1: (x1,y1)
%                 Node 2: (x2,y2)
%                 Node 3: (x3,y3)
%                 Node 4: (x4,y4)
%
%   phys_point     - 1x2 vector [x, y] representing the physical coordinates.
%
% Output:
%   param_coords - 1x2 vector [xi, eta] in the parent domain (should lie in [-1,1]^2).

    % Validate inputs
    if size(node_coordinates,1) ~= 4 || size(node_coordinates,2) ~= 2
        error('node_coordinates must be a 4x2 matrix.');
    end
    if numel(phys_point) ~= 2
        error('phys_point must be a 2-element vector.');
    end

    % Extract nodal coordinates
    x1 = node_coordinates(1,1);  y1 = node_coordinates(1,2);
    x2 = node_coordinates(2,1);  y2 = node_coordinates(2,2);
    x3 = node_coordinates(3,1);  y3 = node_coordinates(3,2);
    x4 = node_coordinates(4,1);  y4 = node_coordinates(4,2);
    
    % Compute d1 and d2
    d1 = 4*phys_point(1) - (x1 + x2 + x3 + x4);
    d2 = 4*phys_point(2) - (y1 + y2 + y3 + y4);
    
    % Compute coefficients (using local numbering as in Fig.2(a))
    a1 = x2 - x1;
    b1 = y2 - y1;
    a2 = x3 - x4;
    b2 = y3 - y4;
    c1 = x3 - x2;
    c2 = y3 - y2;
    
    % Check for the degenerate case (c1 nearly zero)
    tol = 1e-12;
    if abs(c1) < tol
        error('Case with c1 == 0 is not implemented in this analytic solver.');
    end

    % Formulate the quadratic equation in eta:
    % A*eta^2 + B*eta + C = 0
    A = a2 * b1;
    B = - (c2 * b1 - a2 * d1 - b2 * c1);

    C = -(d2 * c1 - c2 * d1);
    
    % Solve quadratic equation
    discriminant = B^2 - 4*A*C;
    if discriminant < 0
        error('No real solution exists for the inverse mapping (discriminant < 0).');
    end
    sqrt_disc = sqrt(discriminant);
    eta_candidates = [(-B + sqrt_disc) / (2*A), (-B - sqrt_disc) / (2*A)];
    
    % Choose the eta candidate that lies in [-1, 1]
    valid_eta = [];
    for eta_val = eta_candidates
        if eta_val >= -1-1e-6 && eta_val <= 1+1e-6
            valid_eta(end+1) = eta_val; %#ok<AGROW>
        end
    end
    if isempty(valid_eta)
        error('No valid eta solution in [-1, 1] found.');
    end
    
    % For each valid eta, compute xi and check that it lies in [-1,1].
    param_coords = [];
    for eta = valid_eta
        xi = (d1 - b1 * eta) / c1;
        if xi >= -1-1e-6 && xi <= 1+1e-6
            param_coords = [xi, eta];
            break;
        end
    end
    if isempty(param_coords)
        error('No valid (xi, eta) solution in [-1,1]^2 found.');
    end
end
