
clear all;
close all;
clc;

%% Define output directory for figures
outputDir = 'New Folder 7';
% Create the directory if it doesn't exist
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Add paths to required directories
addpath('FEM')
addpath('MESH')
addpath('PLOT')

%% Define parameters
% Plate parameters
a = 6.0;  % plate length [m]
b = 2.0;  % plate width [m]
h = 0.015; % plate thickness [m]
E = 80e9; % Young's modulus [Pa]
nu = 0.33;   % Poisson's ratio
rho = 2700; % density [kg/m³]

% Mesh parameters
nx = 25; % elements in x-direction
ny = 10; % elements in y-direction
support = 'simple'; % boundary conditions

% Create mesh
disp('Creating mesh and assembling matrices...')
mesh = StructMeshGenerator(a, b, nx, ny, E, nu, rho, h, support);

% Add point load at (a/4, b/4)
nodes = mesh.X(:,2:3);
dist = sqrt((nodes(:,1) - a/4).^2 + (nodes(:,2) - b/4).^2);
[~, load_node] = min(dist);
mesh.PointLoads = [load_node 3 1000]; % 1000 N point load
mesh = CompressNumbering(mesh);

% Assemble the system
study = struct();
study.element = 'mindlin';
study.assembly = 'standard';
study.intMethod = 'Full';
opt = struct();
opt = AssemblyMindlin(mesh, study, opt);

% We need to create the mass matrix M explicitly
% Since it's not provided by AssemblyMindlin
% Create empty mass matrix
M = sparse(opt.neqn, opt.neqn);

% Assemble mass matrix from element matrices
for e = 1:size(mesh.IX, 1)
    % Get element nodes
    nen = mesh.IX(e, 2:5);
    
    % Get coordinates
    xy = mesh.X(nen, 2:3);
    
    % Get element DOFs
    edof = zeros(1, 12);
    for i = 1:4
        edof(3*i-2) = 3*nen(i)-2;
        edof(3*i-1) = 3*nen(i)-1;
        edof(3*i-0) = 3*nen(i)-0;
    end
    
    % Get material parameters
    matID = mesh.IX(e, 6);
    thk = mesh.Material(matID, 1);
    E_val = mesh.Material(matID, 2);
    nu_val = mesh.Material(matID, 3);
    rho_val = mesh.Material(matID, 4);
    
    % Get element matrices
    [ke, me, ~] = elementMatrixMindlin(xy(:,1)', xy(:,2)', E_val, thk, nu_val, rho_val, study.intMethod);
    
    % Assemble mass matrix
    for i = 1:length(edof)
        for j = 1:length(edof)
            M(edof(i), edof(j)) = M(edof(i), edof(j)) + me(i, j);
        end
    end
end

% Fix the boundary conditions more explicitly
% Modify stiffness matrix for BCs
K_modified = opt.K;
for i = 1:opt.neqn
    if opt.Null(i,i) == 0  % This is a constrained DOF
        K_modified(i,:) = 0;
        K_modified(:,i) = 0;
        K_modified(i,i) = 1;
    end
end

% Modify mass matrix similarly
M_modified = M;
for i = 1:opt.neqn
    if opt.Null(i,i) == 0  % This is a constrained DOF
        M_modified(i,:) = 0;
        M_modified(:,i) = 0;
        M_modified(i,i) = 1e-6;  % Small mass for numerical stability
    end
end

% Modify force vector
F = opt.P(:,1);
F_modified = F;
for i = 1:opt.neqn
    if opt.Null(i,i) == 0  % This is a constrained DOF
        F_modified(i) = 0;
    end
end

%% Exercise 2.1: Eigenvalue analysis
disp('EXERCISE 2.1: EIGENVALUE ANALYSIS')
disp('Running eigenvalue analysis...')

% Solve eigenvalue problem
nModes = 20;
[eigenvectors, eigenvalues_diag] = eigs(K_modified, M_modified, nModes, 'sm', 'Tolerance', 1e-5, 'MaxIterations', 1000);
eigenvalues = diag(eigenvalues_diag);

% Sort eigenvalues and eigenvectors
[eigenvalues, idx] = sort(eigenvalues);
eigenvectors = eigenvectors(:, idx);

% Calculate natural frequencies in Hz
nat_freqs = sqrt(eigenvalues) / (2*pi);
disp(['First 5 natural frequencies (Hz): ' num2str(real(nat_freqs(1:5))', '%.2f ')])

% Create frequency range around natural frequencies
omega_min = 0.5 * sqrt(real(eigenvalues(1)));
omega_max = 1.5 * sqrt(real(eigenvalues(end)));
n_points = 500;
omega_range = linspace(omega_min, omega_max, n_points);

%% Exercise 2.1: Direct frequency response
disp('EXERCISE 2.1: DIRECT FREQUENCY RESPONSE')
disp('Running direct frequency response analysis...')

% Initialize response matrix
nDOF = size(K_modified, 1);
nFreq = length(omega_range);
U_freq_direct = zeros(nDOF, nFreq);

% Loop through frequencies and solve
tic;
for i = 1:nFreq
    omega = omega_range(i);
    
    % For undamped forced vibration: (K - omega^2*M)*U = F
    A = K_modified - omega^2 * M_modified;
    
    % Solve the system
    U_freq_direct(:, i) = A \ F_modified;
end
direct_time = toc;
disp(['Direct method completed in ' num2str(direct_time) ' seconds'])

% Get response at loaded node
response_dof = 3*load_node; % Vertical DOF at loaded node
u_direct = abs(U_freq_direct(response_dof, :));

% Plot frequency response
figure;
semilogx(omega_range/(2*pi), 20*log10(u_direct), 'b-', 'LineWidth', 1.5);
hold on;
grid on;
xlabel('Frequency [Hz]');
ylabel('Displacement [dB re 1m]');
title('Vertical Response at Point Load (a/4, b/4) - Direct Method');

% Mark natural frequencies
for i = 1:length(nat_freqs)
    if real(nat_freqs(i)) <= omega_max/(2*pi) && isreal(nat_freqs(i))
        xline(real(nat_freqs(i)), 'k:', ['Mode ' num2str(i)]);
    end
end

%% Exercise 2.2: Forced response - modal superposition
disp('EXERCISE 2.2: MODAL SUPERPOSITION')

% Test with different numbers of modes
mode_counts = [3, 6, 9, 15, 20];
u_modal_results = zeros(length(mode_counts), length(omega_range));
modal_times = zeros(size(mode_counts));

% New figure for comparing modal results
figure;
semilogx(omega_range/(2*pi), 20*log10(u_direct), 'k-', 'LineWidth', 2, 'DisplayName', 'Direct Method');
hold on;
grid on;
xlabel('Frequency [Hz]');
ylabel('Displacement [dB re 1m]');
title('Modal Superposition vs Direct Method - Effect of Number of Modes');

% Loop through different mode counts
for i = 1:length(mode_counts)
    nModes_i = mode_counts(i);
    disp(['Running modal superposition with ' num2str(nModes_i) ' modes...'])
    
    % Initialize response matrix
    U_freq_modal = zeros(nDOF, nFreq);
    
    % Get reduced set of eigenvectors and eigenvalues
    eigenvectors_i = eigenvectors(:, 1:nModes_i);
    eigenvalues_i = eigenvalues(1:nModes_i);
    
    % Calculate modal force vector
    F_modal = eigenvectors_i' * F_modified;
    
    % Start timing
    tic;
    
    % Loop through frequencies
    for j = 1:nFreq
        omega = omega_range(j);
        
        % Initialize modal response
        U_modal = zeros(nModes_i, 1);
        
        % Calculate modal response for each mode
        for k = 1:nModes_i
            % Modal response: q_j = F_j / (omega_j^2 - omega^2)
            % Add small damping term to avoid division by zero at resonance
            omega_n = sqrt(eigenvalues_i(k));
            damping_term = 1e-4 * omega_n^2;
            U_modal(k) = F_modal(k) / (omega_n^2 - omega^2 + 1i*damping_term);
        end
        
        % Transform back to physical coordinates: U = Phi * q
        U_freq_modal(:, j) = eigenvectors_i * U_modal;
    end
    
    % Stop timing
    modal_times(i) = toc;
    
    % Extract response
    u_modal = abs(U_freq_modal(response_dof, :));
    u_modal_results(i, :) = u_modal;
    
    % Plot the result
    semilogx(omega_range/(2*pi), 20*log10(u_modal), 'LineWidth', 1.2, 'DisplayName', sprintf('%d Modes', nModes_i));
    
    disp(['Modal method with ' num2str(nModes_i) ' modes completed in ' num2str(modal_times(i)) ' seconds']);
    disp(['   Modal method took ' num2str(100*modal_times(i)/direct_time, '%.1f') '% of direct method time']);
end

legend('show');




%% Exercise 2.3: Forced response - modal acceleration
disp('EXERCISE 2.3: MODAL ACCELERATION')

% Test with different numbers of modes for modal acceleration
figure;
semilogx(omega_range/(2*pi), 20*log10(u_direct), 'k-', 'LineWidth', 2, 'DisplayName', 'Direct Method');
hold on;
grid on;
xlabel('Frequency [Hz]');
ylabel('Displacement [dB re 1m]');
title('Modal Acceleration vs Direct Method - Effect of Number of Modes');

% Calculate static response
U_static = K_modified \ F_modified;

% Loop through different mode counts for modal acceleration
modal_acc_times = zeros(size(mode_counts));
for i = 1:length(mode_counts)
    nModes_i = mode_counts(i);
    disp(['Running modal acceleration with ' num2str(nModes_i) ' modes...'])
    
    % Initialize response matrix
    U_freq_modal_acc = zeros(nDOF, nFreq);
    
    % Get reduced set of eigenvectors and eigenvalues
    eigenvectors_i = eigenvectors(:, 1:nModes_i);
    eigenvalues_i = eigenvalues(1:nModes_i);
    
    % Calculate modal force vector
    F_modal = eigenvectors_i' * F_modified;
    
    % Start timing
    tic;
    
    % Loop through frequencies
    for j = 1:nFreq
        omega = omega_range(j);
        
        % Initialize modal correction
        U_correction = zeros(nDOF, 1);
        
        % Calculate dynamic correction for each mode
        for k = 1:nModes_i
            % Modal acceleration formula: 
            % u = u_static - sum[ (omega^2/omega_j^2) * F_j * phi_j / (omega_j^2 - omega^2)]
            omega_n = sqrt(eigenvalues_i(k));
            damping_term = 1e-4 * omega_n^2;
            modal_factor = (omega^2 / omega_n^2) * F_modal(k) / (omega_n^2 - omega^2 + 1i*damping_term);
            U_correction = U_correction + modal_factor * eigenvectors_i(:, k);
        end
        
        % Combine static solution and modal correction
        U_freq_modal_acc(:, j) = U_static - U_correction;
    end
    
    % Stop timing
    modal_acc_times(i) = toc;
    
    % Extract response
    u_modal_acc = abs(U_freq_modal_acc(response_dof, :));
    
    % Plot the result
    semilogx(omega_range/(2*pi), 20*log10(u_modal_acc), 'LineWidth', 1.2, 'DisplayName', sprintf('%d Modes', nModes_i));
    
    disp(['Modal acceleration with ' num2str(nModes_i) ' modes completed in ' num2str(modal_acc_times(i)) ' seconds']);
    disp(['   Modal acceleration took ' num2str(100*modal_acc_times(i)/direct_time, '%.1f') '% of direct method time']);
end

legend('show');

%% Exercise 2.4: Different excitation - Compare load locations
disp('EXERCISE 2.4: DIFFERENT EXCITATION LOCATIONS')

% Define different load positions to test
load_positions = [
    a/2, b/2;    % Center
    a/3, b/3;    % Interior
    a/10, b/10   % Near boundary
];

position_names = {'Center (a/2, b/2)', 'Interior (a/3, b/3)', 'Near boundary (a/10, b/10)'};

% Create a new figure
figure;
subplot(3, 1, 1);
ylabel('Direct [dB]');
title('Effect of Load Position on Frequency Response');
hold on; grid on;

subplot(3, 1, 2);
ylabel('Modal (10 modes) [dB]');
hold on; grid on;

subplot(3, 1, 3);
xlabel('Frequency [Hz]');
ylabel('Modal acc (10 modes) [dB]');
hold on; grid on;

% Use different colors for different positions
colors = {'b', 'r', 'g'};

% Loop through different load positions
for i = 1:size(load_positions, 1)
    disp(['Analyzing load position: ' position_names{i}])
    
    % Create a new mesh
    mesh_pos = StructMeshGenerator(a, b, nx, ny, E, nu, rho, h, support);
    
    % Add point load at current position
    load_x = load_positions(i, 1);
    load_y = load_positions(i, 2);
    nodes = mesh_pos.X(:,2:3);
    dist = sqrt((nodes(:,1) - load_x).^2 + (nodes(:,2) - load_y).^2);
    [~, load_node] = min(dist);
    
    % Clear existing point loads and add new one
    mesh_pos.PointLoads = [load_node 3 1000];
    
    % Compress the node and element numbering
    mesh_pos = CompressNumbering(mesh_pos);
    
    % Vertical response DOF
    response_dof = 3*load_node;
    
    % Assemble the system
    opt_pos = struct();
    opt_pos = AssemblyMindlin(mesh_pos, study, opt_pos);
    
    % Create mass matrix
    M_pos = sparse(opt_pos.neqn, opt_pos.neqn);
    
    % Assemble mass matrix
    for e = 1:size(mesh_pos.IX, 1)
        % Get element nodes
        nen = mesh_pos.IX(e, 2:5);
        
        % Get coordinates
        xy = mesh_pos.X(nen, 2:3);
        
        % Get element DOFs
        edof = zeros(1, 12);
        for j = 1:4
            edof(3*j-2) = 3*nen(j)-2;
            edof(3*j-1) = 3*nen(j)-1;
            edof(3*j-0) = 3*nen(j)-0;
        end
        
        % Get material parameters
        matID = mesh_pos.IX(e, 6);
        thk = mesh_pos.Material(matID, 1);
        E_val = mesh_pos.Material(matID, 2);
        nu_val = mesh_pos.Material(matID, 3);
        rho_val = mesh_pos.Material(matID, 4);
        
        % Get element matrices
        [ke, me, ~] = elementMatrixMindlin(xy(:,1)', xy(:,2)', E_val, thk, nu_val, rho_val, study.intMethod);
        
        % Assemble mass matrix
        for j = 1:length(edof)
            for k = 1:length(edof)
                M_pos(edof(j), edof(k)) = M_pos(edof(j), edof(k)) + me(j, k);
            end
        end
    end
    
    % Apply boundary conditions
    K_pos = opt_pos.K;
    for j = 1:opt_pos.neqn
        if opt_pos.Null(j,j) == 0  % This is a constrained DOF
            K_pos(j,:) = 0;
            K_pos(:,j) = 0;
            K_pos(j,j) = 1;
            
            M_pos(j,:) = 0;
            M_pos(:,j) = 0;
            M_pos(j,j) = 1e-6;
        end
    end
    
    % Modify force vector
    F_pos = opt_pos.P(:,1);
    for j = 1:opt_pos.neqn
        if opt_pos.Null(j,j) == 0
            F_pos(j) = 0;
        end
    end
    
    % Eigenvalue analysis
    [eigenvectors_pos, eigenvalues_diag_pos] = eigs(K_pos, M_pos, nModes, 'sm', 'Tolerance', 1e-5, 'MaxIterations', 1000);
    eigenvalues_pos = diag(eigenvalues_diag_pos);
    [eigenvalues_pos, idx] = sort(eigenvalues_pos);
    eigenvectors_pos = eigenvectors_pos(:, idx);
    
    % Direct method
    U_pos_direct = zeros(opt_pos.neqn, nFreq);
    for j = 1:nFreq
        omega = omega_range(j);
        A = K_pos - omega^2 * M_pos;
        U_pos_direct(:, j) = A \ F_pos;
    end
    
    % Modal method with 10 modes
    U_pos_modal = zeros(opt_pos.neqn, nFreq);
    nModes_use = 10;
    eigenvectors_use = eigenvectors_pos(:, 1:nModes_use);
    eigenvalues_use = eigenvalues_pos(1:nModes_use);
    F_modal = eigenvectors_use' * F_pos;
    
    for j = 1:nFreq
        omega = omega_range(j);
        U_modal = zeros(nModes_use, 1);
        
        for k = 1:nModes_use
            omega_n = sqrt(eigenvalues_use(k));
            damping_term = 1e-4 * omega_n^2;
            U_modal(k) = F_modal(k) / (omega_n^2 - omega^2 + 1i*damping_term);
        end
        
        U_pos_modal(:, j) = eigenvectors_use * U_modal;
    end
    
    % Modal acceleration with 10 modes
    U_pos_modal_acc = zeros(opt_pos.neqn, nFreq);
    U_static_pos = K_pos \ F_pos;
    
    for j = 1:nFreq
        omega = omega_range(j);
        U_correction = zeros(opt_pos.neqn, 1);
        
        for k = 1:nModes_use
            omega_n = sqrt(eigenvalues_use(k));
            damping_term = 1e-4 * omega_n^2;
            modal_factor = (omega^2 / omega_n^2) * F_modal(k) / (omega_n^2 - omega^2 + 1i*damping_term);
            U_correction = U_correction + modal_factor * eigenvectors_use(:, k);
        end
        
        U_pos_modal_acc(:, j) = U_static_pos - U_correction;
    end
    
    % Extract responses
    u_pos_direct = abs(U_pos_direct(response_dof, :));
    u_pos_modal = abs(U_pos_modal(response_dof, :));
    u_pos_modal_acc = abs(U_pos_modal_acc(response_dof, :));
    
    % Plot results
    subplot(3, 1, 1);
    semilogx(omega_range/(2*pi), 20*log10(u_pos_direct), colors{i}, 'LineWidth', 1.2, 'DisplayName', position_names{i});
    
    subplot(3, 1, 2);
    semilogx(omega_range/(2*pi), 20*log10(u_pos_modal), colors{i}, 'LineWidth', 1.2, 'DisplayName', position_names{i});
    
    subplot(3, 1, 3);
    semilogx(omega_range/(2*pi), 20*log10(u_pos_modal_acc), colors{i}, 'LineWidth', 1.2, 'DisplayName', position_names{i});
end

% Add legends
subplot(3, 1, 1); legend('show', 'Location', 'best');
subplot(3, 1, 2); legend('show', 'Location', 'best');
subplot(3, 1, 3); legend('show', 'Location', 'best');

%% Exercise 2.5: Performance comparison
disp('EXERCISE 2.5: PERFORMANCE COMPARISON')

% Summarize timing results
fprintf('\n===== Performance comparison =====\n');
fprintf('Direct method: %.2f seconds\n', direct_time);

fprintf('\nModal superposition method:\n');
for i = 1:length(mode_counts)
    fprintf('  %d modes: %.2f seconds (%.1f%% of direct time)\n', ...
        mode_counts(i), modal_times(i), 100*modal_times(i)/direct_time);
end

fprintf('\nModal acceleration method:\n');
for i = 1:length(mode_counts)
    fprintf('  %d modes: %.2f seconds (%.1f%% of direct time)\n', ...
        mode_counts(i), modal_acc_times(i), 100*modal_acc_times(i)/direct_time);
end

% Comprehensive timing plot
figure;
bar([direct_time, modal_times, modal_acc_times]);
set(gca, 'XTickLabel', {'Direct', ...
    sprintf('Modal-%d', mode_counts(1)), ...
    sprintf('Modal-%d', mode_counts(2)), ...
    sprintf('Modal-%d', mode_counts(3)), ...
    sprintf('Modal-%d', mode_counts(4)), ...
    sprintf('Modal-Acc-%d', mode_counts(1)), ...
    sprintf('Modal-Acc-%d', mode_counts(2)), ...
    sprintf('Modal-Acc-%d', mode_counts(3)), ...
    sprintf('Modal-Acc-%d', mode_counts(4))});
ylabel('Computation time [s]');
title('Performance Comparison');
xtickangle(45);
grid on;

%% Exercise 2.6: Frequency range refinement
disp('EXERCISE 2.6: FREQUENCY RANGE REFINEMENT')

discretization_points = [50, 100, 500, 1000];
figure;
title('Effect of Frequency Resolution');
xlabel('Frequency [Hz]');
ylabel('Displacement [dB re 1m]');
grid on;
hold on;

for i = 1:length(discretization_points)
    n_pts = discretization_points(i);
    disp(['Testing with ' num2str(n_pts) ' frequency points...'])
    omega_disc = linspace(omega_min, omega_max, n_pts);
    
    % Run direct analysis with this discretization
    U_disc = zeros(nDOF, n_pts);
    
    for j = 1:n_pts
        omega = omega_disc(j);
        A = K_modified - omega^2 * M_modified;
        U_disc(:, j) = A \ F_modified;
    end
    
    % Extract response
    u_disc = abs(U_disc(response_dof, :));
    
    % Plot result
    semilogx(omega_disc/(2*pi), 20*log10(u_disc), 'LineWidth', 1.2, 'DisplayName', sprintf('%d points', n_pts));
end

legend('show', 'Location', 'best');

%% Exercise 2.7: Anti-resonances
disp('EXERCISE 2.7: ANTI-RESONANCES')

% Define different load positions to test
anti_positions = [
    a/2, b/2;    % Center
    a/4, b/4;    % Quarter point
    a/6, b/6;    % Near corner
];

position_names = {'Center (a/2, b/2)', 'Quarter point (a/4, b/4)', 'Near corner (a/6, b/6)'};

figure;
title('Anti-resonances with Different Load Positions');
xlabel('Frequency [Hz]');
ylabel('Displacement [dB re 1m]');
grid on;
hold on;

% Use different colors for different positions
colors = {'b', 'r', 'g'};

% Use finer frequency resolution for anti-resonance investigation
finer_omega = linspace(omega_min, omega_max, 1000);

% Loop through different load positions
for i = 1:size(anti_positions, 1)
    disp(['Analyzing anti-resonances for position: ' position_names{i}])
    
    % Create a new mesh
    mesh_anti = StructMeshGenerator(a, b, nx, ny, E, nu, rho, h, support);
    
    % Add point load at current position
    load_x = anti_positions(i, 1);
    load_y = anti_positions(i, 2);
    nodes = mesh_anti.X(:,2:3);
    dist = sqrt((nodes(:,1) - load_x).^2 + (nodes(:,2) - load_y).^2);
    [~, load_node] = min(dist);
    
    % Clear existing point loads and add new one
    mesh_anti.PointLoads = [load_node 3 1000];
    
    % Compress the node and element numbering
    mesh_anti = CompressNumbering(mesh_anti);
    
    % Vertical response DOF
    response_dof = 3*load_node;
    
    % Assemble the system
    opt_anti = struct();
    opt_anti = AssemblyMindlin(mesh_anti, study, opt_anti);
    
    % Create mass matrix
    M_anti = sparse(opt_anti.neqn, opt_anti.neqn);
    
    % Assemble mass matrix
    for e = 1:size(mesh_anti.IX, 1)
        % Get element nodes
        nen = mesh_anti.IX(e, 2:5);
        
        % Get coordinates
        xy = mesh_anti.X(nen, 2:3);
        
        % Get element DOFs
        edof = zeros(1, 12);
        for j = 1:4
            edof(3*j-2) = 3*nen(j)-2;
            edof(3*j-1) = 3*nen(j)-1;
            edof(3*j-0) = 3*nen(j)-0;
        end
        
        % Get material parameters
        matID = mesh_anti.IX(e, 6);
        thk = mesh_anti.Material(matID, 1);
        E_val = mesh_anti.Material(matID, 2);
        nu_val = mesh_anti.Material(matID, 3);
        rho_val = mesh_anti.Material(matID, 4);
        
        % Get element matrices
        [ke, me, ~] = elementMatrixMindlin(xy(:,1)', xy(:,2)', E_val, thk, nu_val, rho_val, study.intMethod);
        
        % Assemble mass matrix
        for j = 1:length(edof)
            for k = 1:length(edof)
                M_anti(edof(j), edof(k)) = M_anti(edof(j), edof(k)) + me(j, k);
            end
        end
    end
    
    % Apply boundary conditions
    K_anti = opt_anti.K;
    for j = 1:opt_anti.neqn
        if opt_anti.Null(j,j) == 0  % This is a constrained DOF
            K_anti(j,:) = 0;
            K_anti(:,j) = 0;
            K_anti(j,j) = 1;
            
            M_anti(j,:) = 0;
            M_anti(:,j) = 0;
            M_anti(j,j) = 1e-6;
        end
    end
    
    % Modify force vector
    F_anti = opt_anti.P(:,1);
    for j = 1:opt_anti.neqn
        if opt_anti.Null(j,j) == 0
            F_anti(j) = 0;
        end
    end
    
    % Direct method with finer resolution
    U_anti = zeros(opt_anti.neqn, length(finer_omega));
    
    for j = 1:length(finer_omega)
        omega = finer_omega(j);
        A = K_anti - omega^2 * M_anti;
        U_anti(:, j) = A \ F_anti;
    end
    
    % Extract response
    u_anti = abs(U_anti(response_dof, :));
    
    % Plot result
    semilogx(finer_omega/(2*pi), 20*log10(u_anti), colors{i}, 'LineWidth', 1.2, 'DisplayName', position_names{i});
end

legend('show', 'Location', 'best');

%% Save results for report
disp('Saving results...')
save('Exercise2_Results.mat', 'omega_range', 'u_direct', 'u_modal_results', ...
    'mode_counts', 'direct_time', 'modal_times', 'modal_acc_times');

% Plot mode shapes for visualization
figure;
for i = 1:min(3, size(eigenvectors, 2))
    subplot(1, 3, i);
    
    % Get mode shape for this eigenmode
    mode_shape = eigenvectors(:, i);
    
    % Create visualization data for plotting
    nno = size(mesh.X, 1);
    Vp = zeros(nno, 3);
    Vp(:, 1) = mode_shape(1:3:end);  % w
    Vp(:, 2) = mode_shape(2:3:end);  % theta_x
    Vp(:, 3) = mode_shape(3:3:end);  % theta_y
    
    % Plot deformation using patched quad mesh
    nel = size(mesh.IX, 1);
    xdata = zeros(nel, 4);
    ydata = zeros(nel, 4);
    zdata = zeros(nel, 4);
    
    for e = 1:nel
        % Get node indices for this element
        nen = mesh.IX(e, 2:5);
        
        % Get coordinates
        xdata(e, :) = mesh.X(nen, 2);
        ydata(e, :) = mesh.X(nen, 3);
        zdata(e, :) = real(Vp(nen, 1));  % Use w component (vertical)
    end
    
    % Plot the mode shape
    patch(xdata', ydata', zdata', 'edgecolor', 'none');
    title(['Mode ' num2str(i) ' - ' num2str(real(nat_freqs(i)), '%.2f') ' Hz']);
    axis equal tight;
    colormap(jet);
    colorbar;
    view(3);
end

%% Write field data for ParaView visualization
disp('Creating visualization files...')
% Extract first few modes for visualization
nvis = min(5, size(eigenvectors, 2));
Vp = zeros(size(mesh.X, 1), nvis);

% Extract vertical component for each mode
for i = 1:nvis
    Vp(:, i) = real(eigenvectors(1:3:end, i));
end

% Create names for visualization fields
VpNames = cell(nvis, 1);
for i = 1:nvis
    VpNames{i} = sprintf('Mode%d_%.2fHz', i, real(nat_freqs(i)));
end

% Write mode shapes to VTU file
VtuWriter(Vp, VpNames, [], [], mesh.X(:,2:3), mesh.IX(:,2:5), 4, 'Quadrilateral', 2, 'Exercise2_ModeShapes');

% Extract frequency response data for visualization
nvis_freq = min(5, length(omega_range));
vis_indices = round(linspace(1, length(omega_range), nvis_freq));
Vp_freq = zeros(size(mesh.X, 1), nvis_freq);

% Extract vertical displacement at selected frequencies
for i = 1:nvis_freq
    idx = vis_indices(i);
    Vp_freq(:, i) = abs(U_freq_direct(1:3:end, idx));
end

% Create names for frequency response fields
VpNames_freq = cell(nvis_freq, 1);
for i = 1:nvis_freq
    freq_hz = omega_range(vis_indices(i))/(2*pi);
    VpNames_freq{i} = sprintf('Resp_%.2fHz', freq_hz);
end

% Write frequency response to VTU file
VtuWriter(Vp_freq, VpNames_freq, [], [], mesh.X(:,2:3), mesh.IX(:,2:5), 4, 'Quadrilateral', 2, 'Exercise2_FreqResponse');

disp('All exercises completed successfully!')

%% Save all figures
disp('Saving all figures to output directory...')

% Get all open figure handles
figHandles = findall(0, 'Type', 'figure');

% Define descriptive names for each figure
figNames = {
    'Fig_Ex2_1_DirectMethod',
    'Fig_Ex2_2_ModalSuperposition',
    'Fig_Ex2_3_ModalAcceleration',
    'Fig_Ex2_4_LoadPositions',
    'Fig_Ex2_5_PerformanceComparison',
    'Fig_Ex2_6_FrequencyResolution',
    'Fig_Ex2_7_Antiresonances',
    'Fig_Ex2_ModeShapes'
};

% If we have more figures than names, add generic names
if length(figHandles) > length(figNames)
    for i = (length(figNames)+1):length(figHandles)
        figNames{i} = ['Fig_Ex2_Extra_' num2str(i-length(figNames))];
    end
end

% Save each figure
for i = 1:length(figHandles)
    % Make sure we don't try to access names that don't exist
    if i <= length(figNames)
        figName = figNames{i};
    else
        figName = ['Fig_Ex2_' num2str(i)];
    end
    
    % Activate the figure
    figure(figHandles(i));
    
    % Save as PNG and FIG
    saveas(figHandles(i), fullfile(outputDir, [figName '.png']));
    saveas(figHandles(i), fullfile(outputDir, [figName '.fig']));
    
    disp(['Saved figure ' num2str(i) ' as ' figName]);
end

% Also update these save commands
save(fullfile(outputDir, 'Exercise2_Results.mat'), 'omega_range', 'u_direct', 'u_modal_results', ...
    'mode_counts', 'direct_time', 'modal_times', 'modal_acc_times');

VtuWriter(Vp, VpNames, [], [], mesh.X(:,2:3), mesh.IX(:,2:5), 4, 'Quadrilateral', 2, fullfile(outputDir, 'Exercise2_ModeShapes'));
VtuWriter(Vp_freq, VpNames_freq, [], [], mesh.X(:,2:3), mesh.IX(:,2:5), 4, 'Quadrilateral', 2, fullfile(outputDir, 'Exercise2_FreqResponse'));

disp(['All figures and data saved to ' outputDir]);