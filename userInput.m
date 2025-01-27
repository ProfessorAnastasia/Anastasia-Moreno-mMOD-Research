function userInput()

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

    defaultAnswer = {'0 0 0'};
    answer = inputdlg('Enter [xi eta zeta] in the range [-1,1]:', ...
        'Parametric Coordinates', [1 50], defaultAnswer);
    if isempty(answer)
        disp('Exiting Program...');
        return;
    end

    param_coordinates = str2num(answer{1});
    if length(param_coordinates) ~= 3
        error('Expected 3 numeric values (xi, eta, zeta).')
    end

    mapped_point = isoparametricWeightedSum(param_coordinates, node_coordinates);
    disp('Parametric coordinates:'), disp(param_coordinates);
    disp('Mapped physical point:'), disp(mapped_point);

    figure('Name','Trilinear Hexahedral Isoparametric Mapping','Color','w','Units','normalized','Position',[0.1 0.2 0.8 0.5]);
    
    % (a) Left subplot: Parametric domain
    subplot(1,2,1)
    hold on; axis equal; grid on;
    xlabel('\xi'); ylabel('\eta'); zlabel('\zeta');
    title('Parametric Domain [-1,1]^3');
    % Draw edges of the parametric cube
    param_corners = [ ...
        -1 -1 -1;
         1 -1 -1;
         1  1 -1;
        -1  1 -1;
        -1 -1  1;
         1 -1  1;
         1  1  1;
        -1  1  1 ];
    plotHexahedralEdges(param_corners, 'b', 1.5)
    
    % Plot user parametric point
    plot3(param_coordinates(1), param_coordinates(2), param_coordinates(3), 'ro','MarkerSize',8,'LineWidth',2)
    xlim([-1.2, 1.2]); ylim([-1.2, 1.2]); zlim([-1.2, 1.2]);
    view(30,20);
    
    % (b) Right subplot: Physical domain
    subplot(1,2,2)
    hold on; axis equal; grid on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Physical Hexahedron');
    % Draw edges of the actual hexahedron
    plotHexahedralEdges(node_coordinates, 'k', 1.5)
    
    % Plot the mapped physical point
    plot3(mapped_point(1), mapped_point(2), mapped_point(3), 'ro','MarkerSize',8,'LineWidth',2)
    % Adjust axes so everything is nicely visible
    xlim([min(node_coordinates(:,1))-0.5, max(node_coordinates(:,1))+0.5])
    ylim([min(node_coordinates(:,2))-0.5, max(node_coordinates(:,2))+0.5])
    zlim([min(node_coordinates(:,3))-0.5, max(node_coordinates(:,3))+0.5])
    view(30,20);
    
    hold off;

end

