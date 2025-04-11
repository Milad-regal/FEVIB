function opt = AssemblyMindlin(mesh,study,opt)
% function opt = AssembleMindlin(mesh,study,opt)
% 
% Method that returns assemblies of K,M,C and P
% output: opt.K, opt.M, opt.C, opt.Null, opt.P

%% Global numbers (for safe keeping)
opt.nel = size(mesh.IX,1);
opt.neqn = size(mesh.X,1)*3;

%% SUPPORTS (The N-Matrix)
N = ones(opt.neqn,1);
opt.g=zeros(opt.neqn,1);
for i=1:size(mesh.bound,1)
   N(3*mesh.bound(i,1)-(3-mesh.bound(i,2))) = 0;
   opt.g(3*mesh.bound(i,1)-(3-mesh.bound(i,2))) = mesh.bound(i,3);
end
opt.Null = spdiags(N,0,opt.neqn,opt.neqn);

%% VECTOR (LOADS)
opt.P = zeros(opt.neqn,2); % first column is for points, second for press.
for i=1:size(mesh.PointLoads,1)
   opt.P(3*mesh.PointLoads(i,1)-(3-mesh.PointLoads(i,2)),1) = mesh.PointLoads(i,3); 
end
%opt.P(:,2)=-1;
%% MATRICES (STIFFMESS, MASS and DAMPING)
ldof=12;

% Zero the arrays for the triplets
I = zeros(opt.nel*ldof*ldof,1);
J = zeros(opt.nel*ldof*ldof,1);
KE = zeros(opt.nel*ldof*ldof,1);
ME = zeros(opt.nel*ldof*ldof, 1); %new array for saving the values of mass matrix

ntriplets = 0;

if strcmp(study.assembly,'standard')==1
    % Loop over element and integrate
    for e=1:opt.nel
        % element nodes
        nen = mesh.IX(e,2:5);
        
        % Get coordinates
        xy = mesh.X(nen,2:3);
        
        % Get element dofs
        for i=1:4
            edof(3*i-2) = 3*nen(i)-2;
            edof(3*i-1) = 3*nen(i)-1;
            edof(3*i-0) = 3*nen(i)-0;
        end
        
        % Get material parameters
        matID = mesh.IX(e,6);
        thk = mesh.Material(matID,1);        E   = mesh.Material(matID,2);
        nu  = mesh.Material(matID,3);        rho = mesh.Material(matID,4);        
        
        % Get the integrated element matrices
        [ke, me, fe] = elementMatrixMindlin(xy(:,1)',xy(:,2)',E,thk,nu,rho,study.intMethod);
        
        % add to global system (I,J,[KE,ME,CE])
        for krow = 1:ldof
            for kcol = 1:ldof
                ntriplets = ntriplets+1;
                I(ntriplets) = edof(krow);
                J(ntriplets) = edof(kcol);
                KE(ntriplets) = ke(krow,kcol);
                ME(ntriplets) = me(krow,kcol); % mass values matrix saving  
            end
        end
        
    end
    %assembe global matrices
    ind = find(I>0);
    opt.K = sparse(I(ind),J(ind),KE(ind),opt.neqn,opt.neqn);
    opt.M = sparse(I(ind),J(ind),ME(ind),opt.neqn,opt.neqn);% global mass matrice
    opt.C = sparse(I(ind),J(ind),ME(ind),opt.neqn,opt.neqn);% global mass matrice
end
end