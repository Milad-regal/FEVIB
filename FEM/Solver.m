function opt = Solver(mesh,study,opt)
% function opt = Solver(mesh,study,p)
%
% Input 
% study.analysis =  'eigen_lin'
%                   'eigen_quad'
%                   'freq_direct'
%                   'freq_modal'
%                   'static'
% study.X with X = neig (int), omega (float(:)) 
%
% mesh: mesh data (X,IX,bound,....)
% opt: Matrices (p.K,p.C,p.M,p.P,p.Null)
%
% Output
% opt.U , opt.w , (opt.pdof)

%% Eigenvalue - damped, undamped
if strcmp(study.analysis,'eigen_lin')==1 || strcmp(study.analysis,'eigen_quad')==1
    % Call a eigensolver for quadratic problems (and a Dirichlet matrix)
    [opt.U,opt.D] = SolverEigen(opt.K,opt.M,opt.C,opt.Null,study.neig,study.analysis);
	
%% STATIC
elseif strcmp(study.analysis,'static')==1
    % Modidy stiffness matrix for BCs
    opt.K = opt.Null'*opt.K*opt.Null - (opt.Null-speye(opt.neqn,opt.neqn));
    opt.P = opt.Null*opt.P;
    % Solve static problem
    tic;
    opt.U = opt.K \ opt.P;
    time = toc;
    % Assign output parameter
    opt.pdof = 3*mesh.wdof(1) - (3-mesh.wdof(2));
    %p.pdof = find(p.P~=0);
    opt.w = opt.U(opt.pdof,:);
    fprintf('Static solution found in %f s\n',time);
	fprintf('FE predicted displacement of the center node in plate %f (m)\n',opt.w );

%% undamped forced vibration
   elseif strcmp(study.analysis,'freq_direct')==1
    % Initialize response matrix
    numFreqs = length(study.omega);
    opt.U_response = zeros(size(opt.P,1), numFreqs);
    opt.P = opt.P(:, 1); % Keep only the first column

    opt.K = opt.Null'*opt.K*opt.Null - (opt.Null-speye(opt.neqn,opt.neqn));
    opt.M = opt.Null'*opt.M*opt.Null - (opt.Null-speye(opt.neqn,opt.neqn));
    opt.P = opt.Null*opt.P;
    % Solve for each frequency
    for i = 1:numFreqs
        omega_i = study.omega(i);
        % Solve (K - ω^2 M) U = P
       
        K_eff = opt.K - omega_i^2 * opt.M;
        opt.U_response(:, i) =  K_eff \ opt.P;
        %opt.U_response(:, i) = pinv(full(K_eff)) * opt.P;

    end

    % Find displacement of the loaded node
    opt.pdof = 3*mesh.wdof(1) - (3-mesh.wdof(2)); % Adjust as needed
    opt.w = opt.U_response(opt.pdof, :);
    
    fprintf('Forced response analysis completed for %d frequencies\n', numFreqs);

%% undamped forced vibration (modal superposition method)

elseif strcmp(study.analysis, 'freq_modal') == 1
  

      % Apply boundary conditions
    opt.K = opt.Null' * opt.K * opt.Null - (opt.Null - speye(opt.neqn));
    opt.M = opt.Null' * opt.M * opt.Null - (opt.Null - speye(opt.neqn));
    opt.P = opt.Null * opt.P(:,1);  % Use only first column of force

    % Eigenvalue analysis (linear undamped)
    [opt.U_modes, opt.D] = SolverEigen(opt.K, opt.M, opt.C, opt.Null, study.neig, 'eigen_lin');

    % Extract natural frequencies
    omega_n = sqrt(diag(opt.D));

    % Truncate modes
    num_modes = study.neig;
    omega_n = omega_n(1:num_modes);
    Phi = opt.U_modes(:, 1:num_modes);

    % Project force vector into modal space
    F_modal = Phi' * opt.P;
    % Allocate response
    numFreqs = length(study.omega);
    opt.U_response = zeros(size(opt.P,1), numFreqs);

    % Small artificial damping (optional)
    zeta = 0.002;

    % Solve modal response for each frequency
    % Solve modal response for each frequency
for i = 1:numFreqs
    omega_i = study.omega(i);

    % Calculate modal response for each mode
    U_modal = zeros(num_modes,1);
    for j = 1:num_modes
        denom = (omega_n(j)^2 - omega_i^2) + 1i * 2 * zeta * omega_i * omega_n(j);
        U_modal(j) = F_modal(j) / denom;
    end

    % Transform back to physical coordinates (reduced system)
    U_reduced = Phi * U_modal;

    % Expand to full DOF set
    U_full = opt.Null * U_reduced;

    % Store full response
    opt.U_response(:, i) = U_full;
end


    % Extract vertical response at loaded DOF
    opt.pdof = 3*mesh.wdof(1) - (50 - mesh.wdof(2));  % Adjust if needed
    opt.w = opt.U_response(opt.pdof, :);

    fprintf('Modal superposition solution computed with %d modes\n', num_modes);


%% Forced vib-modal acceleration
elseif strcmp(study.analysis, 'freq_modal_acc') == 1


 opt.K = opt.Null'*opt.K*opt.Null - (opt.Null-speye(opt.neqn,opt.neqn));
    opt.M = opt.Null'*opt.M*opt.Null - (opt.Null-speye(opt.neqn,opt.neqn));
    opt.P = opt.Null'*opt.P(:,1);
    [opt.U, opt.D] = SolverEigen(opt.K, opt.M, opt.C, opt.Null, study.neig, 'eigen_lin');



    % Extract natural frequencies
    omega_n = sqrt(diag(opt.D));
    
    % Project force vector into modal space
    F_modal = opt.U' * opt.P(:,1);
    
    % Truncate to selected number of modes
    num_modes = study.neig;
    omega_n = omega_n(1:num_modes);
    F_modal = F_modal(1:num_modes);
    opt.U = opt.U(:, 1:num_modes);
    
    % Static response: U_static = K \ F
    U_static = opt.K \ opt.P(:,1);
    
    % Preallocate response matrix
    opt.U_response = zeros(size(opt.U, 1), length(study.omega));
    
    % Small damping ratio
    zeta = 0.002;

    % Loop over excitation frequencies
    for i = 1:length(study.omega)
        omega_i = study.omega(i);
        
        % Initialize modal correction
        U_corr = zeros(size(opt.U, 1), 1);
        
        % Loop over each mode
        for j = 1:num_modes
            % Modal acceleration correction term
            denom = (omega_n(j)^2 - omega_i^2) + 1i * 2 * zeta * omega_i * omega_n(j);
            modal_factor = (omega_i^2 / omega_n(j)^2) * F_modal(j) / denom;
            
            % Add correction in physical coordinates
            U_corr = U_corr + modal_factor * opt.U(:, j);
        end
        
        % Final response = static - correction
        opt.U_response(:, i) = U_static - U_corr;
    end

  
end
	
end

