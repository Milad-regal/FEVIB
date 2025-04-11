function PlotFields(mesh,param,V)
% function opt = Postprocess(mesh,param,V)
% plot the solution
%
% mesh:     the mesh struct
% param:    parameters for visualization
%           .title = 'string'
%           .name = 'string'
% V:        Solutions to be plotted

for sol=1:size(V,2);
    
    % Create the patch
    nel = size(mesh.IX,1);
    xdata = zeros(nel,4);
    ydata = zeros(nel,4);
    udata = zeros(nel,4);
    vdata = zeros(nel,4);
    wdata = zeros(nel,4);
    
    for e=1:nel
        nen =  mesh.IX(e,2:5);
        xy = mesh.X(nen,2:3);
        xdata(e,:) = xy(:,1);
        ydata(e,:) = xy(:,2);
        for i=1:4
            edof(3*i-2) = 3*nen(i)-2;
            edof(3*i-1) = 3*nen(i)-1;
            edof(3*i-0) = 3*nen(i)-0;
        end
        udata(e,:)=V(edof(1:3:end),sol);
        vdata(e,:)=V(edof(2:3:end),sol);
        wdata(e,:)=V(edof(3:3:end),sol);
    end
    
    figure;
    subplot(3,1,1)
    patch(xdata',ydata',real(udata)','edgecolor','none')
    if ~isempty(param)
        % check if title is there
        title(param.title)
    else  
        title('w')
    end
    axis equal
    colormap jet
    colorbar
    
    subplot(3,1,2)
    patch(xdata',ydata',real(vdata)','edgecolor','none')
    title(['\theta_x'])
    axis equal
    colormap jet
    colorbar
    
    subplot(3,1,3)
    patch(xdata',ydata',real(wdata)','edgecolor','none')
    title(['\theta_y'])
    axis equal
    colormap jet
    colorbar
end

end