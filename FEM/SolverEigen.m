function [U, D] = SolverEigen(K, M, C, Null, neig, type)
% Function that solves linear eigenvalue problems
% with the possibility to equate some entries to zero.
%
% Inputs:
%   K    - Stiffness matrix
%   M    - Mass matrix
%   C    - Damping matrix (not used for linear EVP)
%   Null - Null space matrix for applying boundary conditions
%   neig - Number of eigenvalues to compute
%   type - Type of eigenvalue problem ('eigen_lin' for linear EVP)
%
% Outputs:
%   U    - Eigenvectors (mass-normalized)
%   D    - Eigenvalues (sorted in ascending order)

% Apply boundary conditions to K and M
E = eye(size(Null));
K=Null'*K*Null-(E-Null);
M=Null'*M*Null;
% Solve the generalized eigenvalue problem: K * v = lambda * M * v
[U, D] = eigs(K, M, neig, 'smallestabs');

% Extract eigenvalues from the diagonal of D
D = diag(D);

% Sort eigenvalues and eigenvectors in ascending order
[D, sortIdx] = sort(D, 'ascend');
U = U(:, sortIdx);

% Mass-normalize eigenvectors
for i = 1:neig
    U(:, i) = U(:, i) / sqrt(U(:, i)' * M * U(:, i));
end
end