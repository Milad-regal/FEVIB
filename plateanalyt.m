function [w] = plateanalyt(Dflex,dens,t,a,b,Q,analysisType,loadType)
% Analytic response for deflection and eigenvalues for a simply supported
% plate

if strcmp(analysisType,'eigenvalue')==1
    % EIGENVALUE
    w=zeros(6);
    for m=1:6
        for n=1:6
    w(m,n) = sqrt(Dflex/dens/t)*((m*pi/a)^2+(n*pi/b)^2);
        end
    end

elseif  strcmp(analysisType,'static')==1
  
    w = 0;
    for n = 1:1000
        for m = 1:1000
            if strcmp(loadType,'point_load')==1
                p = 4*Q/(a*b)*sin(m*pi*0.5)*sin(n*pi*0.5);
            elseif strcmp(loadType,'pressure_load')==1
                pres = Q/(a*b);
                p = 4*pres/(pi^2*m*n)*(cos(m*pi)-1)*(cos(n*pi)-1);
            else
                error('loadType not recognized')
            end
            den = pi^4*n^4*((m*b/(n*a))^2+1)^2;
            w  = w + b^4/den*(p/Dflex)*sin(0.5*pi*m)*sin(0.5*pi*n);
        end
    end
else
    error('analysisType not recognized')
end

end

