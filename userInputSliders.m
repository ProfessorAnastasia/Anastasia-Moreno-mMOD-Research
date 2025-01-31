function userInputSliders()
clc; clear;
    % Node coordinates (same as before)
    node_coordinates = [...
        0.00 0.00 0.00;  % Node 1
        2.00 0.00 0.00;  % Node 2
        2.00 1.00 0.00;  % Node 3
        0.00 1.00 0.00;  % Node 4
        0.00 0.00 1.00;  % Node 5
        2.00 0.00 1.20;  % Node 6
        2.00 1.00 1.20;  % Node 7
        0.00 1.00 1.00;  % Node 8
    ];

    % Create main figure with adjusted layout for sliders
    fig = figure('Name','Interactive Isoparametric Mapping',...
                 'Color','w',...
                 'Units','normalized',...
                 'Position',[0.1 0.1 0.8 0.7],...
                 'NumberTitle','off');
    
    % Create subplots (same as before)
    ax1 = subplot(1,2,1);
    hold(ax1, 'on'); axis(ax1, 'equal'); grid(ax1, 'on');
    xlabel(ax1, '\xi'); ylabel(ax1, '\eta'); zlabel(ax1, '\zeta');
    title(ax1, 'Parametric Domain [-1,1]^3');
    param_corners = [-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1];
    plotHexahedralEdges(param_corners, 'b', 1.5)
    xlim(ax1, [-1.2, 1.2]); ylim(ax1, [-1.2, 1.2]); zlim(ax1, [-1.2, 1.2]);
    view(ax1, 30, 20);
    
    ax2 = subplot(1,2,2);
    hold(ax2, 'on'); axis(ax2, 'equal'); grid(ax2, 'on');
    xlabel(ax2, 'X'); ylabel(ax2, 'Y'); zlabel(ax2, 'Z');
    title(ax2, 'Physical Hexahedron');
    plotHexahedralEdges(node_coordinates, 'k', 1.5)
    xlim(ax2, [min(node_coordinates(:,1))-0.5, max(node_coordinates(:,1))+0.5]);
    ylim(ax2, [min(node_coordinates(:,2))-0.5, max(node_coordinates(:,2))+0.5]);
    zlim(ax2, [min(node_coordinates(:,3))-0.5, max(node_coordinates(:,3))+0.5]);
    view(ax2, 30, 20);

    % Create slider panel at the bottom
    panelHeight = 0.15;
    sliderPanel = uipanel('Parent', fig,...
                          'Position',[0 0 1 panelHeight],...
                          'BackgroundColor',[0.95 0.95 0.95]);
    
    % Create sliders and labels
    sliderWidth = 0.25;
    sliderSpacing = 0.05;
    labelWidth = 0.08;
    
    % Xi slider
    uicontrol('Parent', sliderPanel, 'Style','text',...
              'Units','normalized',...
              'Position',[0.05 0.6 labelWidth 0.3],...
              'String','ξ:',...
              'BackgroundColor',[0.95 0.95 0.95]);
    
    xiSlider = uicontrol('Parent', sliderPanel, 'Style','slider',...
                        'Units','normalized',...
                        'Position',[0.05+labelWidth 0.6 sliderWidth 0.3],...
                        'Min',-1, 'Max',1, 'Value',0,...
                        'Tag','xiSlider');
    
    xiValue = uicontrol('Parent', sliderPanel, 'Style','text',...
                       'Units','normalized',...
                       'Position',[0.05+labelWidth+sliderWidth+0.01 0.6 0.1 0.3],...
                       'String','0.00',...
                       'BackgroundColor',[0.95 0.95 0.95]);
    
    % Eta slider (same structure as xi)
    uicontrol('Parent', sliderPanel, 'Style','text',...
              'Units','normalized',...
              'Position',[0.05 0.3 labelWidth 0.3],...
              'String','η:',...
              'BackgroundColor',[0.95 0.95 0.95]);
    
    etaSlider = uicontrol('Parent', sliderPanel, 'Style','slider',...
                         'Units','normalized',...
                         'Position',[0.05+labelWidth 0.3 sliderWidth 0.3],...
                         'Min',-1, 'Max',1, 'Value',0,...
                         'Tag','etaSlider');
    
    etaValue = uicontrol('Parent', sliderPanel, 'Style','text',...
                        'Units','normalized',...
                        'Position',[0.05+labelWidth+sliderWidth+0.01 0.3 0.1 0.3],...
                        'String','0.00',...
                        'BackgroundColor',[0.95 0.95 0.95]);
    
    % Zeta slider (same structure)
    uicontrol('Parent', sliderPanel, 'Style','text',...
              'Units','normalized',...
              'Position',[0.05 0.0 labelWidth 0.3],...
              'String','ζ:',...
              'BackgroundColor',[0.95 0.95 0.95]);
    
    zetaSlider = uicontrol('Parent', sliderPanel, 'Style','slider',...
                          'Units','normalized',...
                          'Position',[0.05+labelWidth 0.0 sliderWidth 0.3],...
                          'Min',-1, 'Max',1, 'Value',0,...
                          'Tag','zetaSlider');
    
    zetaValue = uicontrol('Parent', sliderPanel, 'Style','text',...
                         'Units','normalized',...
                         'Position',[0.05+labelWidth+sliderWidth+0.01 0.0 0.1 0.3],...
                         'String','0.00',...
                         'BackgroundColor',[0.95 0.95 0.95]);

    % Create initial plot points
    param_point = plot3(ax1, 0, 0, 0, 'ro', 'MarkerSize',8, 'LineWidth',2);
    phys_point = plot3(ax2, 0, 0, 0, 'ro', 'MarkerSize',8, 'LineWidth',2);
    
    % Update function for all sliders
    function updatePlot(~,~)
        % Get current slider values
        xi = xiSlider.Value;
        eta = etaSlider.Value;
        zeta = zetaSlider.Value;
        
        % Update value displays
        xiValue.String = sprintf('%.2f', xi);
        etaValue.String = sprintf('%.2f', eta);
        zetaValue.String = sprintf('%.2f', zeta);
        
        % Calculate mapped point
        mapped_point = isoparametricWeightedSum([xi, eta, zeta], node_coordinates);
        
        % Update plot points
        param_point.XData = xi;
        param_point.YData = eta;
        param_point.ZData = zeta;
        
        phys_point.XData = mapped_point(1);
        phys_point.YData = mapped_point(2);
        phys_point.ZData = mapped_point(3);
        
        drawnow;
    end

    % Set callback for all sliders
    xiSlider.Callback = @updatePlot;
    etaSlider.Callback = @updatePlot;
    zetaSlider.Callback = @updatePlot;
    
    % Initial update
    updatePlot();
end