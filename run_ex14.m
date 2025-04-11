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
study.analysis = 'eigen_lin'; % Eigenvalue analysis
study.intMethod = 'Selective'; % Integration method
study.neig = 10; % Number of eigenvalues to compute

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

%% Display eigenvalues
disp('Eigenvalues:');
disp(opt.D.^0.5);

%% Plot eigenvectors (optional)
% Add code here to plot the eigenvectors if needed.
% plate specifications
Dflex=(E*thk^3)/(12*(1-nu^2));
dens=2700;
t=0.025;
a=4;
b=2;
Q=-1;
analysisType='eigenvalue';

Analytical_solution_results=plateanalyt(Dflex,dens,t,a,b,Q,analysisType);
sorted=sort(Analytical_solution_results(:));
disp('Analytical_Eigenvalues')
disp(sorted(1:10))
[row, col] = find(ismember(Analytical_solution_results, sorted(1:10)))

%% Plot fields in Matlab:
NumSol = 1;
par.title = ['Static solution - point load, w'];
PlotFields(mesh,par,opt.U(:,NumSol));

% %% Write fields for ParaView
% disp('VtuWriter - write for ParaView')
% Vp = [opt.U(1:3:end,1) opt.U(2:3:end,1) opt.U(3:3:end,1) opt.U(1:3:end,2) opt.U(2:3:end,2) opt.U(3:3:end,2)];
% VpNames = {'00 w (p)'; '01 rot_x (p)'; '02 rot y (p)';'03 w (d)'; '04 rot_x (d)'; '05 rot y (d)'};
% VtuWriter(Vp,VpNames,[],[], mesh.X(:,2:3), mesh.IX(:,2:5), 4,...
%     'Quadrilateral', 2, 'ParaViewTestFile');