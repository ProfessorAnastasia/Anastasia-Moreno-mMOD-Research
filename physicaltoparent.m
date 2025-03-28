function userInputSlidersPhysicalToParentMapping()
    clc; clear;
    % --- Define physical node coordinates for our hexahedron ---
    node_coordinates = [...
        0.00 0.00 0.00;  % Node 1
        3.00 0.00 0.00;  % Node 2
        2.50 1.00 0.00;  % Node 3
        0.00 3.00 0.00;  % Node 4
        0.00 0.00 1.00;  % Node 5
        2.00 0.00 2.20;  % Node 6
        2.00 1.00 1.20;  % Node 7
        0.00 1.00 1.00;  % Node 8
    ];
    
    % --- Get axis-aligned bounding box from node coordinates (no extra margin) ---
    xmin = min(node_coordinates(:,1));
    xmax = max(node_coordinates(:,1));
    ymin = min(node_coordinates(:,2));
    ymax = max(node_coordinates(:,2));
    zmin = min(node_coordinates(:,3));
    zmax = max(node_coordinates(:,3));
    
    % --- Create figure and subplots ---
    fig = figure('Name','Physical -> Parent Domain Mapping via Newton-Raphson (Bounded)',...
                 'Color','w',...
                 'Units','normalized',...
                 'Position',[0.1 0.1 0.8 0.7],...
                 'NumberTitle','off');
    
    % Parametric (parent) domain subplot (left)
    ax1 = subplot(1,2,1);
    hold(ax1, 'on'); axis(ax1, 'equal'); grid(ax1, 'on');
    xlabel(ax1, '\xi'); ylabel(ax1, '\eta'); zlabel(ax1, '\zeta');
    title(ax1, 'Parent Domain [-1,1]^3');
    % Draw the parent domain as a box
    param_corners = [-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1];
    plotHexahedralEdges(param_corners, 'b', 1.5)
    xlim(ax1, [-1.2, 1.2]); ylim(ax1, [-1.2, 1.2]); zlim(ax1, [-1.2, 1.2]);
    view(ax1, 30, 20);
    
    % Physical domain subplot (right)
    ax2 = subplot(1,2,2);
    hold(ax2, 'on'); axis(ax2, 'equal'); grid(ax2, 'on');
    xlabel(ax2, 'X'); ylabel(ax2, 'Y'); zlabel(ax2, 'Z');
    title(ax2, 'Physical Hexahedron');
    plotHexahedralEdges(node_coordinates, 'k', 1.5)
    xlim(ax2, [xmin, xmax]);
    ylim(ax2, [ymin, ymax]);
    zlim(ax2, [zmin, zmax]);
    view(ax2, 30, 20);
    
    % --- Create slider panel for physical coordinates ---
    panelHeight = 0.15;
    sliderPanel = uipanel('Parent', fig,...
                          'Position',[0 0 1 panelHeight],...
                          'BackgroundColor',[0.95 0.95 0.95]);
    
    sliderWidth = 0.25;
    labelWidth = 0.08;
    
    % X slider
    uicontrol('Parent', sliderPanel, 'Style','text',...
              'Units','normalized',...
              'Position',[0.05 0.6 labelWidth 0.3],...
              'String','X:',...
              'BackgroundColor',[0.95 0.95 0.95]);
    xSlider = uicontrol('Parent', sliderPanel, 'Style','slider',...
                        'Units','normalized',...
                        'Position',[0.05+labelWidth 0.6 sliderWidth 0.3],...
                        'Min',xmin, 'Max',xmax, 'Value',(xmin+xmax)/2,...
                        'Tag','xSlider');
    xValue = uicontrol('Parent', sliderPanel, 'Style','text',...
                       'Units','normalized',...
                       'Position',[0.05+labelWidth+sliderWidth+0.01 0.6 0.1 0.3],...
                       'String',sprintf('%.2f',(xmin+xmax)/2),...
                       'BackgroundColor',[0.95 0.95 0.95]);
    
    % Y slider
    uicontrol('Parent', sliderPanel, 'Style','text',...
              'Units','normalized',...
              'Position',[0.05 0.3 labelWidth 0.3],...
              'String','Y:',...
              'BackgroundColor',[0.95 0.95 0.95]);
    ySlider = uicontrol('Parent', sliderPanel, 'Style','slider',...
                        'Units','normalized',...
                        'Position',[0.05+labelWidth 0.3 sliderWidth 0.3],...
                        'Min',ymin, 'Max',ymax, 'Value',(ymin+ymax)/2,...
                        'Tag','ySlider');
    yValue = uicontrol('Parent', sliderPanel, 'Style','text',...
                       'Units','normalized',...
                       'Position',[0.05+labelWidth+sliderWidth+0.01 0.3 0.1 0.3],...
                       'String',sprintf('%.2f',(ymin+ymax)/2),...
                       'BackgroundColor',[0.95 0.95 0.95]);
    
    % Z slider
    uicontrol('Parent', sliderPanel, 'Style','text',...
              'Units','normalized',...
              'Position',[0.05 0.0 labelWidth 0.3],...
              'String','Z:',...
              'BackgroundColor',[0.95 0.95 0.95]);
    zSlider = uicontrol('Parent', sliderPanel, 'Style','slider',...
                        'Units','normalized',...
                        'Position',[0.05+labelWidth 0.0 sliderWidth 0.3],...
                        'Min',zmin, 'Max',zmax, 'Value',(zmin+zmax)/2,...
                        'Tag','zSlider');
    zValue = uicontrol('Parent', sliderPanel, 'Style','text',...
                       'Units','normalized',...
                       'Position',[0.05+labelWidth+sliderWidth+0.01 0.0 0.1 0.3],...
                       'String',sprintf('%.2f',(zmin+zmax)/2),...
                       'BackgroundColor',[0.95 0.95 0.95]);
    
    % --- Markers ---
    % Red dot for physical point from slider values
    phys_point  = plot3(ax2, 0, 0, 0, 'ro', 'MarkerSize',8, 'LineWidth',2);
    % Green dot for computed parent coordinates
    param_point = plot3(ax1, 0, 0, 0, 'go', 'MarkerSize',8, 'LineWidth',2);
    
    % --- Update callback for sliders ---
    function updatePlot(~,~)
        % Read slider values (physical coordinates)
        X = xSlider.Value;
        Y = ySlider.Value;
        Z = zSlider.Value;
        
        % Update the slider labels
        xValue.String = sprintf('%.2f', X);
        yValue.String = sprintf('%.2f', Y);
        zValue.String = sprintf('%.2f', Z);
        
        % Physical point from slider values
        phys_pt = [X, Y, Z];
        
        % Use Newton-Raphson inverse mapping to compute parent coords
        computed_params = NewtonRaphsonInverseMapping(phys_pt, node_coordinates);
        
        % If any parent coordinate is out of [-1,1], then we clamp them,
        % recalc the physical point via forward mapping, and update the sliders.
        if any(abs(computed_params) > 1)
            % Clamp parent coords to [-1,1]
            clamped_params = max(min(computed_params, 1), -1);
            % Recompute physical point with clamped parent coords
            corrected_phys_pt = isoparametricWeightedSum(clamped_params, node_coordinates);
            % Update slider values to reflect the corrected physical point
            X = corrected_phys_pt(1);
            Y = corrected_phys_pt(2);
            Z = corrected_phys_pt(3);
            xSlider.Value = X;
            ySlider.Value = Y;
            zSlider.Value = Z;
            xValue.String = sprintf('%.2f', X);
            yValue.String = sprintf('%.2f', Y);
            zValue.String = sprintf('%.2f', Z);
            % Use the clamped parent parameters from now on
            computed_params = clamped_params;
        end
        
        % Update markers on the plots
        phys_point.XData = X;
        phys_point.YData = Y;
        phys_point.ZData = Z;
        param_point.XData = computed_params(1);
        param_point.YData = computed_params(2);
        param_point.ZData = computed_params(3);
        
        % Update subplot titles with computed values
        ax1.Title.String = sprintf('Computed Parent: [%.3f, %.3f, %.3f]', ...
                                   computed_params(1), computed_params(2), computed_params(3));
        ax2.Title.String = sprintf('Physical Hexahedron (Input): [%.3f, %.3f, %.3f]', X, Y, Z);
        drawnow;
    end

    % Set slider callbacks
    xSlider.Callback = @updatePlot;
    ySlider.Callback = @updatePlot;
    zSlider.Callback = @updatePlot;
    
    % Do an initial update
    updatePlot();
end

%% --- Supporting Functions ---

% Quick shape functions for a trilinear hexahedral element (nothing fancy)
function N = shapeFunctionsHexahedral(xi, eta, zeta)
    N = zeros(1,8);
    N(1) = 0.125 * (1 - xi) * (1 - eta) * (1 - zeta);
    N(2) = 0.125 * (1 + xi) * (1 - eta) * (1 - zeta);
    N(3) = 0.125 * (1 + xi) * (1 + eta) * (1 - zeta);
    N(4) = 0.125 * (1 - xi) * (1 + eta) * (1 - zeta);
    N(5) = 0.125 * (1 - xi) * (1 - eta) * (1 + zeta);
    N(6) = 0.125 * (1 + xi) * (1 - eta) * (1 + zeta);
    N(7) = 0.125 * (1 + xi) * (1 + eta) * (1 + zeta);
    N(8) = 0.125 * (1 - xi) * (1 + eta) * (1 + zeta);
end

% Forward mapping: compute physical point from parent coordinates
function weightedSum = isoparametricWeightedSum(param_coordinates, node_coordinates)
    xi   = param_coordinates(1);
    eta  = param_coordinates(2);
    zeta = param_coordinates(3);
    N = shapeFunctionsHexahedral(xi, eta, zeta);
    weightedSum = N * node_coordinates;
end

% Quick function to plot hexahedral edges given corner points
function plotHexahedralEdges(corners, colorSpec, lineWidth)
    edges = [1 2; 2 3; 3 4; 4 1; ...  % bottom face
             5 6; 6 7; 7 8; 8 5; ...  % top face
             1 5; 2 6; 3 7; 4 8 ];    % vertical edges
    for e = 1:size(edges,1)
        i1 = edges(e,1);
        i2 = edges(e,2);
        pts = [corners(i1,:); corners(i2,:)];
        plot3(pts(:,1), pts(:,2), pts(:,3), '-', 'Color',colorSpec, 'LineWidth',lineWidth);
    end
end

% Newton-Raphson based inverse mapping function (yep, that's its name!)
function param_coords = NewtonRaphsonInverseMapping(phys_point, node_coords)
    tol = 1e-6;    % tolerance for convergence (not super strict)
    max_iter = 20; % maximum iterations (should be enough)
    % Start with a guess right at the center of the parent domain
    xi = 0; eta = 0; zeta = 0;
    
    for iter = 1:max_iter
        % Get shape functions at current guess
        N = shapeFunctionsHexahedral(xi, eta, zeta);
        % Compute current physical position from our guess
        current_point = N * node_coords;
        % Compute error/residual as a column vector (3x1)
        r = (phys_point - current_point)'; 
        % Quick break if we're close enough
        if norm(r) < tol
            break;
        end
        % Get derivatives of shape functions 
        [dN_dxi, dN_deta, dN_dzeta] = ShapeFunctionDerivatives(xi, eta, zeta);
        % Assemble Jacobian matrix (3x3) by summing contributions from each node
        J = zeros(3,3);
        for i = 1:8
            J(:,1) = J(:,1) + dN_dxi(i) * node_coords(i,:)';
            J(:,2) = J(:,2) + dN_deta(i) * node_coords(i,:)';
            J(:,3) = J(:,3) + dN_dzeta(i) * node_coords(i,:)';
        end
        % Solve for parameter update delta (3x1 vector)
        delta = J \ r;
        % Update our guess
        xi   = xi   + delta(1);
        eta  = eta  + delta(2);
        zeta = zeta + delta(3);
    end
    % Return computed parent coordinates
    param_coords = [xi, eta, zeta];
end

% Informal derivative calculator for our shape functions 
function [dN_dxi, dN_deta, dN_dzeta] = ShapeFunctionDerivatives(xi, eta, zeta)
    dN_dxi   = zeros(1,8);
    dN_deta  = zeros(1,8);
    dN_dzeta = zeros(1,8);
    % Node 1:
    dN_dxi(1)   = -0.125 * (1 - eta) * (1 - zeta);
    dN_deta(1)  = -0.125 * (1 - xi)  * (1 - zeta);
    dN_dzeta(1) = -0.125 * (1 - xi)  * (1 - eta);
    % Node 2:
    dN_dxi(2)   =  0.125 * (1 - eta) * (1 - zeta);
    dN_deta(2)  = -0.125 * (1 + xi)  * (1 - zeta);
    dN_dzeta(2) = -0.125 * (1 + xi)  * (1 - eta);
    % Node 3:
    dN_dxi(3)   =  0.125 * (1 + eta) * (1 - zeta);
    dN_deta(3)  =  0.125 * (1 + xi)  * (1 - zeta);
    dN_dzeta(3) = -0.125 * (1 + xi)  * (1 + eta);
    % Node 4:
    dN_dxi(4)   = -0.125 * (1 + eta) * (1 - zeta);
    dN_deta(4)  =  0.125 * (1 - xi)  * (1 - zeta);
    dN_dzeta(4) = -0.125 * (1 - xi)  * (1 + eta);
    % Node 5:
    dN_dxi(5)   = -0.125 * (1 - eta) * (1 + zeta);
    dN_deta(5)  = -0.125 * (1 - xi)  * (1 + zeta);
    dN_dzeta(5) =  0.125 * (1 - xi)  * (1 - eta);
    % Node 6:
    dN_dxi(6)   =  0.125 * (1 - eta) * (1 + zeta);
    dN_deta(6)  = -0.125 * (1 + xi)  * (1 + zeta);
    dN_dzeta(6) =  0.125 * (1 + xi)  * (1 - eta);
    % Node 7:
    dN_dxi(7)   =  0.125 * (1 + eta) * (1 + zeta);
    dN_deta(7)  =  0.125 * (1 + xi)  * (1 + zeta);
    dN_dzeta(7) =  0.125 * (1 + xi)  * (1 + eta);
    % Node 8:
    dN_dxi(8)   = -0.125 * (1 + eta) * (1 + zeta);
    dN_deta(8)  =  0.125 * (1 - xi)  * (1 + zeta);
    dN_dzeta(8) =  0.125 * (1 - xi)  * (1 + eta);
end
