clear
clc
close all

%% Main Code 
%% Add paths to FEA, MESH, and VISUALIZATION (ParaView)
addpath('FEM')
addpath('MESH')
addpath('PLOT')

%% Set up the study
study.element = 'mindlin'; 
study.assembly = 'standard'; 
study.analysis = 'freq_direct'; % Eigenvalue analysis (Change to 'freq_direct' for forced vibration or )
study.intMethod = 'Selective'; % Integration method
study.omega = 0:10:1000;
 % Frequency range for forced vibration (Modify as needed)

%% Create a structured mesh
Lx = 4;       % Length in x-direction
Ly = 2;       % Length in y-direction
nelx = 80;    % Number of elements in x-direction
nely = 40;    % Number of elements in y-direction
E = 70e9;     % Young's modulus
thk = 0.025;  % Thickness
nu = 0.33;    % Poisson's ratio
rho = 2700;   % Density
support = 'simple'; % Boundary condition type

mesh = StructMeshGenerator(Lx, Ly, nelx, nely, E, nu, rho, thk, support);

%% Perform finite element analysis
opt = Controller(mesh, study); 

%% If forced response analysis, plot the frequency response
if strcmp(study.analysis, 'freq_direct')
    %figure;
    %semilogy(study.omega, abs(opt.U_response(500,:))); % Log scale for better visualization


% Frequency response plot (log scale recommended for clarity)

% Convert omega (rad/s) to frequency (Hz), if needed:
f = study.omega;  % Optional — only if you prefer Hz instead of rad/s

% Magnitude of displacement at loaded node
w_mag = abs(opt.w);

% Plot
figure;
semilogy(f, w_mag, 'b-o', 'LineWidth', 1.5);
xlabel('Frequency (Rad/sec)');
ylabel('|w| [m]');
title('Frequency Response (Direct Method)');
grid on;









    hold on;

% Mark natural frequencies
natural_frequencies = [119.79, 191.35, 311.21, 408.51, 479.51, 479.51, 598.40, 696.41, 765.61, 891.63];
for i = 1:length(natural_frequencies)
    xline(natural_frequencies(i), '--r', ['\omega_' num2str(i)], 'LabelHorizontalAlignment', 'right');
end
    xlabel('Frequency (rad/s)');
    ylabel('Log |Displacement|');
    title('Frequency Response of the Loaded Node');
end


% %% Display eigenvalues
% disp('Eigenvalues:');
% disp(opt.D.^0.5);
% 
% %% Compare with analytical solution
% Dflex = (E * thk^3) / (12 * (1 - nu^2));
% dens = 2700;
% t = 0.025;
% a = 4;
% b = 2;
% Q = -1;
% analysisType = 'eigenvalue';
% 
% Analytical_solution_results = plateanalyt(Dflex, dens, t, a, b, Q, analysisType);
% sorted = sort(Analytical_solution_results(:));
% disp('Analytical_Eigenvalues')
% disp(sorted(1:6))
% [row, col] = find(ismember(Analytical_solution_results, sorted(1:6)));
% 
% %% Plot fields in Matlab:
% NumSol = 1;
% par.title = ['Static solution - point load, w'];
% PlotFields(mesh, par, opt.U(:, NumSol));

% %% Write fields for ParaView
% disp('VtuWriter - write for ParaView')
% Vp = [opt.U(1:3:end,1) opt.U(2:3:end,1) opt.U(3:3:end,1) opt.U(1:3:end,2) opt.U(2:3:end,2) opt.U(3:3:end,2)];
% VpNames = {'00 w (p)'; '01 rot_x (p)'; '02 rot y (p)';'03 w (d)'; '04 rot_x (d)'; '05 rot y (d)'};
% VtuWriter(Vp, VpNames, [], [], mesh.X(:,2:3), mesh.IX(:,2:5), 4,...
%     'Quadrilateral', 2, 'ParaViewTestFile');
