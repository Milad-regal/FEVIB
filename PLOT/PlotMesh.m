% Load the mesh data
load('Ex_34_Slot_60x30_quartLoad.mat'); % Ensure the correct filename

% Extract node data
node_numbers = mesh.X(:,1);  % Node numbers
X_coords = mesh.X(:,2);      % X coordinates
Y_coords = mesh.X(:,3);      % Y coordinates

% Extract element connectivity
element_numbers = mesh.IX(:,1);  % Element numbers
elements = mesh.IX(:,2:5);       % Node connectivity (assuming quadrilateral elements)

% Plot the mesh
figure; hold on;
for i = 1:size(elements,1)
    % Get node indices for this element
    node_ids = elements(i, :);
    
    % Get corresponding coordinates
    x = X_coords(node_ids);
    y = Y_coords(node_ids);
    
    % Close the element shape
    x = [x; x(1)];
    y = [y; y(1)];
    
    % Plot the element
    plot(x, y, 'k-', 'LineWidth', 2);
end

% Plot nodes
scatter(X_coords, Y_coords, 50, 'filled', 'r');

% Labels and formatting
% xlabel('X Coordinate');
% ylabel('Y Coordinate');
% title('Mesh Visualization');
axis equal;
axis off;
grid on;

% Annotate node numbers (Optional)
for i = 1:length(node_numbers)
    text(X_coords(i), Y_coords(i), sprintf('%d', node_numbers(i)), ...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
         'FontSize', 20, 'Color', 'blue');
end

% Annotate element numbers (Optional)
for i = 1:length(element_numbers)
    % Compute element centroid
    node_ids = elements(i, :);
    x_centroid = mean(X_coords(node_ids));
    y_centroid = mean(Y_coords(node_ids));
    
    text(x_centroid, y_centroid, sprintf('%d', element_numbers(i)), ...
         'FontSize', 20, 'Color', 'black', 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center');
end

hold off;