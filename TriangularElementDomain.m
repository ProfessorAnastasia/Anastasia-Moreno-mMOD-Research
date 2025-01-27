function TriangularElementDomain(coords)
% TRIANGULARREF2PHYSICAL  Demonstrate the forward isoparametric mapping 
% from the reference triangle (r,s) to a physical triangle (x,y).
%
% Usage:
%   TriangularRef2Physical(coords)
%
%   where coords is a 3x2 matrix:
%     coords(1,:) = [x1, y1]
%     coords(2,:) = [x2, y2]
%     coords(3,:) = [x3, y3]
%
% The function:
% 1) Creates a fine grid of (r,s) points in the reference triangle.
% 2) Maps them to the physical triangle using linear shape functions.
% 3) Plots both the reference triangle and the physical triangle side by side.

    %---------------------------------------------------------
    % 1. Define the reference triangle and create a grid of (r,s) points.
    %    The reference triangle is:
    %         (0,0), (1,0), (0,1) 
    %    with r>=0, s>=0, r+s<=1.
    %---------------------------------------------------------
    Ngrid = 20;            % number of divisions along each edge for plotting
    r_vals = linspace(0,1,Ngrid);
    s_vals = linspace(0,1,Ngrid);

    % We'll collect valid (r,s) points (where r+s <= 1)
    rs_points = [];
    for i = 1:Ngrid
        for j = 1:Ngrid
            r = r_vals(i);
            s = s_vals(j);
            if (r + s <= 1.0)
                rs_points = [rs_points; r, s]; %#ok<AGROW>
            end
        end
    end

    % Also, store the corners of the reference triangle
    refCorners = [0,0; 1,0; 0,1];

    %---------------------------------------------------------
    % 2. Extract (x1, y1), (x2, y2), (x3, y3) from coords
    %---------------------------------------------------------
    x1 = coords(1,1);  y1 = coords(1,2);
    x2 = coords(2,1);  y2 = coords(2,2);
    x3 = coords(3,1);  y3 = coords(3,2);

    %---------------------------------------------------------
    % 3. Map each (r,s) -> (x,y) using linear shape functions:
    %    x(r,s) = x1 * (1 - r - s) + x2*r + x3*s
    %    y(r,s) = y1 * (1 - r - s) + y2*r + y3*s
    %---------------------------------------------------------
    Npoints = size(rs_points,1);
    xy_points = zeros(Npoints, 2);

    for k = 1:Npoints
        r = rs_points(k,1);
        s = rs_points(k,2);
        N1 = 1 - r - s;
        N2 = r;
        N3 = s;
        xk = x1*N1 + x2*N2 + x3*N3;
        yk = y1*N1 + y2*N2 + y3*N3;
        xy_points(k,:) = [xk, yk];
    end

    % Map the reference corners as well (for plotting edges):
    refCornerXY = zeros(3, 2);
    for i = 1:3
        r = refCorners(i,1);
        s = refCorners(i,2);
        N1 = 1 - r - s;
        N2 = r;
        N3 = s;
        refCornerXY(i,1) = x1*N1 + x2*N2 + x3*N3;
        refCornerXY(i,2) = y1*N1 + y2*N2 + y3*N3;
    end

    %---------------------------------------------------------
    % 4. Plot results
    %---------------------------------------------------------
    figure;

    %--- Subplot (1): Reference domain
    subplot(1,2,1);
    % Plot all interior points
    scatter(rs_points(:,1), rs_points(:,2), 20, 'filled', 'MarkerFaceColor',[0.2 0.6 1]);
    hold on;
    % Plot the boundary of the reference triangle
    fill(refCorners(:,1), refCorners(:,2), 'r','FaceAlpha',0.2);
    plot(refCorners([1,2,3,1],1), refCorners([1,2,3,1],2), 'k-', 'LineWidth',1.5);
    hold off;
    axis equal;
    xlabel('r'); ylabel('s');
    title('Reference Triangle (r,s)');

    %--- Subplot (2): Physical domain
    subplot(1,2,2);
    % Plot mapped points
    scatter(xy_points(:,1), xy_points(:,2), 20, 'filled', 'MarkerFaceColor',[0.2 0.6 1]);
    hold on;
    % Plot the edges using the mapped corners
    fill(refCornerXY(:,1), refCornerXY(:,2), 'g','FaceAlpha',0.2);
    plot(refCornerXY([1,2,3,1],1), refCornerXY([1,2,3,1],2), 'k-', 'LineWidth',1.5);
    % Plot the actual physical corners as black dots
    plot(coords(:,1), coords(:,2), 'ko','MarkerFaceColor','k');
    text(x1, y1, ' (x_1,y_1)', 'Color','k','FontSize',8);
    text(x2, y2, ' (x_2,y_2)', 'Color','k','FontSize',8);
    text(x3, y3, ' (x_3,y_3)', 'Color','k','FontSize',8);
    hold off;
    axis equal;
    xlabel('x'); ylabel('y');
    title('Physical Triangle (x,y)');
end