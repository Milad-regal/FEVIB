function mesh = CompressNumbering(mesh)
% Special for ANSYS input deck
% Renumber the nodes and elements to be contigouos

% max nodenumber and number of nodes
nn = size(mesh.X,1);
nnmax = max(mesh.X(:,1));
ne = size(mesh.IX,1);
nemax = max(mesh.IX(:,1));

% Check if all is in order !
if nn == nnmax && ne == nemax
    return
end


% help vectors
nlist = 1:nn;
nvec = zeros(nnmax,1);
nvec(mesh.X(:,1)) = nlist;

% correct the nodes
mesh.X(:,1) = nlist;

% Correct the elements
mesh.IX(:,1) = (1:ne)';
for e=1:ne
   mesh.IX(e,2:5) = nvec(mesh.IX(e,2:5))';
end

% Bound
if size(mesh.bound,1)>0
    mesh.bound(:,1) = nvec(mesh.bound(:,1));
end

% Pointloads
if size(mesh.PointLoads,1)>0
    mesh.PointLoads(:,1) = nvec(mesh.PointLoads(:,1));
end

% Pressure load


end