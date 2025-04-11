function [mesh, study] = Solver2(mesh, study)
% SOLVER Solves the FEM problem defined by mesh and study
%   [mesh, study] = Solver(mesh, study) solves the study defined by
%   mesh and study structures and updates both with solution results

% Extract solution settings
analysisType = study.type;

switch analysisType
    case 'static'
        % Static analysis
        K = mesh.K;
        F = mesh.F;
        
        % Apply boundary conditions
        if isfield(mesh, 'fixedDOFs') && ~isempty(mesh.fixedDOFs)
            freeDOFs = setdiff(1:size(K,1), mesh.fixedDOFs);
            K_free = K(freeDOFs, freeDOFs);
            F_free = F(freeDOFs);
            
            % Solve the reduced system
            u_free = K_free \ F_free;
            
            % Expand solution to include fixed DOFs
            u = zeros(size(K,1), 1);
            u(freeDOFs) = u_free;
        else
            % No boundary conditions
            u = K \ F;
        end
        
        % Store the solution
        study.u = u;
        
    case 'eigen'
        % Eigenvalue analysis
        K = mesh.K;
        M = mesh.M;
        
        % Apply boundary conditions
        if isfield(mesh, 'fixedDOFs') && ~isempty(mesh.fixedDOFs)
            freeDOFs = setdiff(1:size(K,1), mesh.fixedDOFs);
            K_free = K(freeDOFs, freeDOFs);
            M_free = M(freeDOFs, freeDOFs);
            
            % Solve eigenvalue problem
            nModes = min(study.nModes, length(freeDOFs));
            [eigenvectors_free, eigenvalues_free] = eigs(K_free, M_free, nModes, 'sm');
            
            % Extract diagonal eigenvalues
            eigenvalues = diag(eigenvalues_free);
            
            % Expand eigenvectors to include fixed DOFs
            eigenvectors = zeros(size(K,1), nModes);
            eigenvectors(freeDOFs, :) = eigenvectors_free;
        else
            % No boundary conditions
            nModes = min(study.nModes, size(K,1));
            [eigenvectors, eigenvalues_diag] = eigs(K, M, nModes, 'sm');
            eigenvalues = diag(eigenvalues_diag);
        end
        
        % Sort eigenvalues and eigenvectors
        [eigenvalues, idx] = sort(eigenvalues);
        eigenvectors = eigenvectors(:, idx);
        
        % Store results
        study.eigenvalues = eigenvalues;
        study.eigenvectors = eigenvectors;
        
    case 'freq_direct'
        % Exercise 2.1: Direct frequency response analysis
        K = mesh.K;
        M = mesh.M;
        F = mesh.F;
        
        % Apply boundary conditions
        if isfield(mesh, 'fixedDOFs') && ~isempty(mesh.fixedDOFs)
            freeDOFs = setdiff(1:size(K,1), mesh.fixedDOFs);
            K_free = K(freeDOFs, freeDOFs);
            M_free = M(freeDOFs, freeDOFs);
            F_free = F(freeDOFs);
        else
            % No boundary conditions
            K_free = K;
            M_free = M;
            F_free = F;
            freeDOFs = 1:size(K,1);
        end
        
        % Initialize response matrix
        nDOF = size(K, 1);
        nFreq = length(study.omega);
        U_freq = zeros(nDOF, nFreq);
        
        % Loop through frequencies and solve
        for i = 1:nFreq
            omega = study.omega(i);
            
            % For undamped forced vibration: (K - omega^2*M)*U = F
            A = K_free - omega^2 * M_free;
            
            % Solve the system
            U_omega_free = A \ F_free;
            
            % Expand solution to include fixed DOFs
            U_omega = zeros(nDOF, 1);
            U_omega(freeDOFs) = U_omega_free;
            
            % Store the solution for this frequency
            U_freq(:, i) = U_omega;
        end
        
        % Store frequency response in study structure
        study.U_freq = U_freq;
        
    case 'freq_modal'
        % Exercise 2.2: Modal superposition method
        % Make sure we have eigenvalues and eigenvectors
        if ~isfield(study, 'eigenvectors') || ~isfield(study, 'eigenvalues')
            error('Eigenvectors and eigenvalues must be provided for modal analysis');
        end
        
        % Get required data
        K = mesh.K;
        M = mesh.M;
        F = mesh.F;
        nModes = study.nModes;
        eigenvectors = study.eigenvectors(:, 1:nModes);
        eigenvalues = study.eigenvalues(1:nModes);
        omega_n = sqrt(eigenvalues);
        
        % Apply boundary conditions
        if isfield(mesh, 'fixedDOFs') && ~isempty(mesh.fixedDOFs)
            freeDOFs = setdiff(1:size(K,1), mesh.fixedDOFs);
            F_free = F(freeDOFs);
            eigenvectors_free = eigenvectors(freeDOFs, :);
        else
            % No boundary conditions
            F_free = F;
            eigenvectors_free = eigenvectors;
            freeDOFs = 1:size(K,1);
        end
        
        % Calculate modal force vector
        F_modal = eigenvectors_free' * F_free;
        
        % Initialize response matrix
        nDOF = size(K, 1);
        nFreq = length(study.omega);
        U_freq = zeros(nDOF, nFreq);
        
        % Loop through frequencies
        for i = 1:nFreq
            omega = study.omega(i);
            U_modal = zeros(nModes, 1);
            
            % Calculate modal response for each mode
            for j = 1:nModes
                % Modal response: q_j = F_j / (omega_j^2 - omega^2)
                % Add small damping term to avoid division by zero at resonance
                damping_term = 1e-4 * omega_n(j)^2;
                U_modal(j) = F_modal(j) / (omega_n(j)^2 - omega^2 + 1i*damping_term);
            end
            
            % Transform back to physical coordinates: U = Phi * q
            U_physical_free = eigenvectors_free * U_modal;
            
            % Expand to full DOF set
            U_physical = zeros(nDOF, 1);
            U_physical(freeDOFs) = U_physical_free;
            
            % Store the solution
            U_freq(:, i) = U_physical;
        end
        
        % Store frequency response in study structure
        study.U_freq = U_freq;
        
    case 'freq_modal_acc'
        % Exercise 2.3: Modal acceleration method
        % Make sure we have eigenvalues and eigenvectors
        if ~isfield(study, 'eigenvectors') || ~isfield(study, 'eigenvalues')
            error('Eigenvectors and eigenvalues must be provided for modal analysis');
        end
        
        % Get required data
        K = mesh.K;
        M = mesh.M;
        F = mesh.F;
        nModes = study.nModes;
        eigenvectors = study.eigenvectors(:, 1:nModes);
        eigenvalues = study.eigenvalues(1:nModes);
        omega_n = sqrt(eigenvalues);
        
        % Apply boundary conditions
        if isfield(mesh, 'fixedDOFs') && ~isempty(mesh.fixedDOFs)
            freeDOFs = setdiff(1:size(K,1), mesh.fixedDOFs);
            K_free = K(freeDOFs, freeDOFs);
            F_free = F(freeDOFs);
            eigenvectors_free = eigenvectors(freeDOFs, :);
        else
            % No boundary conditions
            K_free = K;
            F_free = F;
            eigenvectors_free = eigenvectors;
            freeDOFs = 1:size(K,1);
        end
        
        % Calculate modal force vector
        F_modal = eigenvectors_free' * F_free;
        
        % Calculate static response (K^-1 * F)
        U_static_free = K_free \ F_free;
        
        % Initialize response matrix
        nDOF = size(K, 1);
        nFreq = length(study.omega);
        U_freq = zeros(nDOF, nFreq);
        
        % Loop through frequencies
        for i = 1:nFreq
            omega = study.omega(i);
            U_correction = zeros(length(freeDOFs), 1);
            
            % Calculate dynamic correction for each mode
            for j = 1:nModes
                % Modal acceleration formula: 
                % u = u_static - sum[ (omega^2/omega_j^2) * F_j * phi_j / (omega_j^2 - omega^2)]
                % Add small damping term to avoid division by zero at resonance
                damping_term = 1e-4 * omega_n(j)^2;
                modal_factor = (omega^2 / omega_n(j)^2) * F_modal(j) / (omega_n(j)^2 - omega^2 + 1i*damping_term);
                U_correction = U_correction + modal_factor * eigenvectors_free(:, j);
            end
            
            % Combine static solution and modal correction
            U_full_free = U_static_free - U_correction;
            
            % Expand to full DOF set
            U_full = zeros(nDOF, 1);
            U_full(freeDOFs) = U_full_free;
            
            % Store the solution
            U_freq(:, i) = U_full;
        end
        
        % Store frequency response in study structure
        study.U_freq = U_freq;
        
    otherwise
        error('Unknown analysis type: %s', analysisType);
end

end