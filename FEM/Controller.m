function opt = Controller(mesh,study)
% Main driver for fea and postprocessing
% opt = Controller(mesh,study)
%
% Input containers:
% mesh.X            = nodal coordinates
% mesh.IX           = topology table (and material index)
% mesh.Material:    = list of materials (e,[thk E nu dens])
% mesh.bound:       = list of constrained nodes (and their dof)
% mesh.PointLoads:  = list of loaded nodes (and their dof + load)
%
% study.element     = 'mindlin'
% study.assembly    = 'standard'
% study.analysis    = 'static',
% study.intMethod   = 'Selective', 'Full', 'Reduced'
% study.neig        = number of eigenvalus/vectors to be found
% study.omega       = list of freqs for frequency responses
%
% Output container
% opt.K,C,M,P,U(field),D(eigenval.), and more... check it out !

%% Print overall information
fprintf('\n####### New analysis ########');
fprintf('\nElement type: %s\nAssembly type: %s\nIntegration: %s\nAnalysis type: %s\nMesh file: %s\n\n',...
    study.element,study.assembly,study.intMethod,study.analysis,mesh.filename);

% Problem size
fprintf('Number of elements: %i \nNumber of dofs: %i \n', size(mesh.IX,1),3*size(mesh.X,1))

%% Initialize empty output 
opt = [];

%% Assemble matrices and vectors
opt = AssemblyMindlin(mesh,study,opt);

%% Solve the finite element problem
opt = Solver(mesh,study,opt);

%% Reconstruction (CMS):
opt = Reconstruct(mesh,study,opt);

%% Postproc for energies, stresses, strains
opt = Postprocess(mesh,study,opt);

end