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