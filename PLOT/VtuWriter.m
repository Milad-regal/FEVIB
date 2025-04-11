function VtuWriter(Vp,VpNames,Vc,VcNames, X, IX, nno, elemType, ndim, filename)
% function VtuWriter(Vp,VpNames,Vc,VcNames, X, IX, nno, elemType, ndim, filename)
%
% Vp: array of pointdata vectors (len(Vp(i)) must equal size(X,1))
% VpNames: List of names for the pointdata
% Vc: array of celldata vectors (len(Vc(i)) must equal size(IX,1))
% VcNames: List of names for the celldata
% X: Nodal coordinate matrix
% IX: Element node matrix
% nno: Number of element nodes
% elemType: Quadrilateral, Triangle, Tetrahedra, Hexahedra
% ndim: Number of dimensions 2 or 3.
% filename: Filename of the vtk file.

%% Initial parameters
ne = size(IX,1);
nn = size(X,1);  % Total of nodes
cell_id = 1;  % Cell id indicate the (VTK) id of the type of polygon 

%% ***** ***** ***** Writing the vtk file ***** ***** *****
fid = fopen([filename '.vtu'],'w');

%% Writing the header from the vtk file
fprintf(fid,'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
fprintf(fid,'<UnstructuredGrid>\n');

%% The cell_id depending of the shape of the finite element
switch elemType
    case 'Quadrilateral'
        % VTK_QUADRATIC_QUAD
        if nno == 8
            cell_id = '23';
            % Renumbering - waiting for Audun test
            %IX = IX(:, [1 3 5 7 2 4 6 8]);
        end
        % VTK_QUAD
        if nno == 4
            cell_id = '9';
        end
    case 'Triangle'
        % VTK_TRIANGLE
        if nno == 3
            cell_id = '5';
        end
        % VTK_QUADRATIC_TRIANGLE
        if nno == 6
            cell_id = '22';
            % Test if renum is nec
            %IX = IX(:, [1 3 5 2 4 6]);
        end
    case 'Tetrahedra'
        % VTK_TETRA
        if nno == 4
            cell_id = '10';
        end
        % VTK_QUADRATIC_TETRA
        if nno == 10
            cell_id = '24';
        end
    case 'Hexahedra'
        % VTK_HEXAHEDRON
        if nno == 8
            cell_id = '12';
        end
        % VTK_QUADRATIC_HEXAHEDRON
        if nno == 20
            cell_id = '25';
            % Test if renum is nec
            % IX(:,[1 2 3 4 5 6 7 8 9 10 11 12 17 18 19 20 13 14 15 16]) = IX;
        end
end

%% Writing the coordinates of each node
fprintf(fid,'<Piece NumberOfPoints="%i" NumberOfCells="%i">\n',nn,ne);
fprintf(fid,'<Points>\n');
fprintf(fid,'<DataArray type="Float32" NumberOfComponents="%i" format="ascii">\n', 3);
if ndim == 2
    X = [X zeros(nn,1)];
end
for i=1:nn
    fprintf(fid, [repmat('%e ',1,3) '\n'], X(i,:));
end
fprintf(fid, '</DataArray>\n');
fprintf(fid, '</Points>\n');

%% Writing the cells or nodes
IX = IX - 1; % Use C indexation

fprintf(fid, '<Cells>\n');
fprintf(fid, '<DataArray type="Int32" Name="connectivity" format="ascii">\n');
for i=1:ne
    fprintf(fid, [repmat('%d ',1,nno) '\n'],  IX(i,:));
end
fprintf(fid, '</DataArray>\n');
%% Offset
fprintf(fid, '<DataArray type="Int32" Name="offsets" format="ascii">\n');
offset = 0;
for i=1:ne
    offset = offset + nno;
    fprintf(fid, '%d\n', offset);
end
fprintf(fid, '</DataArray>\n');

%% Element types
fprintf(fid, '<DataArray type="UInt8" Name="types" format="ascii">\n');
for i=1:ne
    fprintf(fid, '%s\n', cell_id);
end
fprintf(fid, '</DataArray>\n');
fprintf(fid, '</Cells>\n');
%% End the mesh information
%% Cell data: FIELDS
if size(Vc,1)>0
    fprintf(fid, '<CellData>\n');
    for j=1:size(Vc,2)
        if exist('VcNames')==0
            VcNames(j) = ['u_' num2str(j)];
        end
        
        fprintf(fid, '<DataArray type="Float32" NumberOfComponents="1" Name="%s" format="ascii">\n',VcNames{j});
        for i=1:ne
            fprintf(fid, '%e\n', Vc(i,j));
        end
        fprintf(fid, '</DataArray>\n');
    end
    fprintf(fid, '</CellData>\n');
end

%% Point data: FIELDS
if size(Vp,1)>0
    fprintf(fid, '<PointData>\n');
    for j=1:size(Vp,2)
        if exist('VpNames')==0
            VpNames(j) = ['u_' num2str(j)];
        end
        fprintf(fid, '<DataArray type="Float32" NumberOfComponents="1" Name="%s" format="ascii">\n',VpNames{j});
        for i=1:nn
            fprintf(fid, '%e\n', Vp(i,j));
        end
        fprintf(fid, '</DataArray>\n');
    end
    fprintf(fid, '</PointData>\n');
end
%% End the file
fprintf(fid, '</Piece>\n');
fprintf(fid, '</UnstructuredGrid>\n');
fprintf(fid, '</VTKFile>\n');
fclose(fid);
return;