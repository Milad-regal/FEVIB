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
study.analysis = 'freq_modal'; % Eigenvalue analysis (Change to 'freq_direct' for direct method or 'freq_modal' superposition method)
study.intMethod = 'Selective'; % Integration method
study.omega = linspace(0, 1000, 1000);
study.neig=10;
 % Frequency range for forced vibration (Modify as needed)

%% Create a structured mesh
Lx = 4;       % Length in x-direction
Ly = 2;       % Length in y-direction
nelx = 10;    % Number of elements in x-direction
nely = 5;    % Number of elements in y-direction
E = 70e9;     % Young's modulus
thk = 0.025;  % Thickness
nu = 0.33;    % Poisson's ratio
rho = 2700;   % Density
support = 'simple'; % Boundary condition type

mesh = StructMeshGenerator(Lx, Ly, nelx, nely, E, nu, rho, thk, support);

%% Perform finite element analysis
opt = Controller(mesh, study); % This should now handle 'freq_direct'

%% If forced response analysis, plot the frequency response
if strcmp(study.analysis, 'freq_modal')
    % figure;
    % semilogy(study.omega, abs(opt.U_response(8,:)),'linewidth',2.5); % Log scale for better visualization
    % 
    % hold on;

% Convert omega (rad/s) to frequency (Hz), if needed:
f = study.omega;  % Optional — only if you prefer Hz instead of rad/s

% Magnitude of displacement at loaded node
w_mag = abs(opt.w);


% Plot
figure;
semilogy(f, w_mag, 'b-', 'LineWidth', 1.5);
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
    legend('Modal Superposition', 'Natural Frequencies');

end
