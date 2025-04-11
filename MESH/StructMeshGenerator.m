function [mesh] = StructMeshGenerator(Lx,Ly,nelx,nely,E,nu,rho,thk,support)
% Simple, structured mesh generator
%
% mesh = StructMeshGenerator(Lx,Ly,nelx,nely,E,nu,rho,thk,support)
%
% INPUTS:
% Lx, Ly:       dimensions
% nelx, nely:   number of elements in the grid
% E,nu,rho,thk: Material and thickness parameters
% support:      'simple' or 'clamped'
%
% OUTPUT:
% "mesh" structure with fields:
% X:            nodal coordinates
% IX:           topology table
% Material:     [thk E nu dens]
% bound:        list of BCs
% PointLoads:    a single in the center 
% PressureLoads: distributed pressure load 

% Derived data
nel = nelx*nely;
nn =  (1+nelx)*(1+nely);

% Containers
mesh.X = zeros(nn,3);
mesh.IX = zeros(nel,6);
mesh.Material = [thk E nu rho];
mesh.bound = []; % clamped or simple
mesh.PointLoads = []; % in the center (or as close as poss.)
mesh.PressureLoads = []; % evenly distributed
% Name the mesh - with elements
mesh.filename = ['RectMesh_' num2str(nelx) 'x' num2str(nely) ...
    '_E_' num2str(E) '_nu_' num2str(nu) '_rho_' num2str(rho) ...
    '_thk_' num2str(thk) '_' support ];
% Nodes
n=0;
dx=Lx/nelx;
dy=Ly/nely;
for nx=1:nelx+1
    for ny=1:nely+1
        n=n+1;
        mesh.X(n,1) = n;
        mesh.X(n,2) = (nx-1)*dx;
        mesh.X(n,3) = (ny-1)*dy;
    end
end

% Elements
e=0;
for elx=1:nelx
    for ely=1:nely
        e=e+1;
        n1 = (nely+1)*(elx-1)+ ely;
        n2 = (nely+1)* elx   + ely;
        n3 = (nely+1)* elx   + ely+1;
        n4 = (nely+1)*(elx-1)+ ely+1;
        edof = [n1 n2 n3 n4];
        mesh.IX(e,1)=e; 
        mesh.IX(e,2:5)=edof; 
        mesh.IX(e,6)=1; % fixed material ! only ONE !
    end
end

% support
if strcmp(support,'clamped')
    mesh.bound = zeros(3 * (2*(nely+1)+2*(nelx-1)),3);
    b = 0;
    for nx=1:nelx:nelx+1
        for ny=1:nely+1
            nnum = (nely+1)*(nx-1)+ny;
            %fprintf('nx: %i, ny: %i , node: %i \n',nx,ny,nnum)
            b=b+1;
            mesh.bound(b,:) = [nnum 1 0];
            b=b+1;
            mesh.bound(b,:) = [nnum 2 0];
            b=b+1;
            mesh.bound(b,:) = [nnum 3 0];            
        end
    end
    
    for nx=2:nelx
        for ny=1:nely:nely+1
            nnum = (nely+1)*(nx-1)+ny;
            %fprintf('nx: %i, ny: %i , node: %i \n',nx,ny,nnum)
            b=b+1;
            mesh.bound(b,:) = [nnum 1 0];
            b=b+1;
            mesh.bound(b,:) = [nnum 2 0];
            b=b+1;
            mesh.bound(b,:) = [nnum 3 0];          
        end
    end
    
elseif strcmp(support,'simple')
        mesh.bound = zeros((2*(nely+1)+2*(nelx-1)),3);
    b = 0;
    for nx=1:nelx:nelx+1
        for ny=1:nely+1
            nnum = (nely+1)*(nx-1)+ny;
            %fprintf('nx: %i, ny: %i , node: %i \n',nx,ny,nnum)
            b=b+1;
            mesh.bound(b,:) = [nnum 1 0];
        end
    end
    
    for nx=2:nelx
        for ny=1:nely:nely+1
            nnum = (nely+1)*(nx-1)+ny;
            %fprintf('nx: %i, ny: %i , node: %i \n',nx,ny,nnum)
            b=b+1;
            mesh.bound(b,:) = [nnum 1 0];            
        end
    end
end

% Pointload in the center node
   mesh.PointLoads = [floor((nn+1)/2) 1 -1];
   mesh.wdof =  [floor((nn+1)/2) 1 -1];

% %Find the node closest to (1/4a, 1/4b)
% target_x = Lx / 4;
% target_y = Ly / 4;
% distances = (mesh.X(:,2) - target_x).^2 + (mesh.X(:,3) - target_y).^2;
% [~, load_node] = min(distances); % Find the closest node
% 
% % Apply point load at the selected node
% mesh.PointLoads = [load_node 1 -1]; % Load in vertical (w) direction
% mesh.wdof = [load_node 1 -1];
% 
% fprintf('Point load applied at node %d (closest to (%.3f, %.3f))\n', load_node, mesh.X(load_node,2), mesh.X(load_node,3));

end