function [ke, me, fe] = elementMatrixMindlin(x,y,E,t,nu,dens,intMethod)

% FIXED PRESSURE FOR EACH ELEMENT !!!! SHOULD BE FIXED
Q=-1;

% Zero the two term element matrix
% number of element dofs for a mindlin Q4
ldof=12; 
kb = zeros(ldof,ldof);                
ks = zeros(ldof,ldof);
me = zeros(ldof,ldof);
fe = zeros(ldof,1);

% Get Gausspoints for given scheme
[kbGauss, ksGauss] = numGaussPoints(intMethod);
[kbgaussPoints,kbw] = gaussFunc(kbGauss);
[ksgaussPoints,ksw] = gaussFunc(ksGauss);

% Constitutive law
G = E/(2*(1+nu));
DD = E*t^3/(12*(1-nu^2));
D = [DD      nu*DD    0;
     nu*DD    DD       0;
     0      0      (1-nu)*DD/2];
Dm = zeros(5,5);
Dm(1:3,1:3) = D;
k = 5/6; % The shear correction factor !
Dm(4,4) = k*G*t;
Dm(5,5) = k*G*t;

% Integrate bending stiffness
for i = 1:kbGauss
    xi = kbgaussPoints(i);
    for j = 1:kbGauss
       eta = kbgaussPoints(j);
       Bb = zeros(5,ldof);
      [Bm,J,N] = shapeFuncQ4(x,y,eta,xi);
       Bb(1:3,:) = Bm(1:3,:);
       dJ = det(J);
       kb = kb+Bb'*Dm*Bb*kbw(i)*kbw(j)*dJ;
       fe = fe + N(1,:)'*Q*kbw(i)*kbw(j)*dJ;
    end
end

% Integrate shear stiffness
for i = 1:ksGauss
    xi = ksgaussPoints(i);
    for j = 1:ksGauss
       eta = ksgaussPoints(j);
       Bs = zeros(5,ldof);
       [Bm,J,~] = shapeFuncQ4(x,y,eta,xi);
       Bs(4:5,:) = Bm(4:5,:);
       dJ = det(J);
       ks = ks+Bs'*Dm*Bs*ksw(i)*ksw(j)*dJ;
    end
end
% Combine the two
ke = kb + ks;

% Integrate mass matrix
meGauss = 2;
[megaussPoints,mew] = gaussFunc(meGauss);
for i = 1:meGauss
    xi = megaussPoints(i);
    for j = 1:meGauss
       eta = megaussPoints(j);
       [~,J,N] = shapeFuncQ4(x,y,eta,xi);
       Mc = dens*[t 0 0;
             0  t^3/12 0;
             0 0 t^3/12];
       dJ = det(J);
       me = me+N'*Mc*N*mew(i)*mew(j)*dJ;
    end
end

end

% Shapefunction and strain displacement for a quadrilateral
function [Bm,J,N] = shapeFuncQ4(x,y,eta,xi)
% Zero the strain-displacement, jacobian, shapefunction matrices
Bm = zeros(5,12);
J = zeros(2,2);
N = zeros(3,12);

% The nodal shape functions
NN = zeros(4,1);
NN(1,1) = ((1 - xi) * (1 - eta)) / 4;
NN(2,1) = ((1 + xi) * (1 - eta)) / 4;
NN(3,1) = ((1 + xi) * (1 + eta)) / 4;
NN(4,1) = ((1 - xi) * (1 + eta)) / 4;

% The shape functions for the displ. and two rotations
for i = 1:4
    N(1,i*3-2) = NN(i,1);
    N(2,i*3-1) = NN(i,1);
    N(3,i*3) = NN(i,1);
end

% Jacobian
J(1,1) = 1/4*[(eta-1) (1-eta) (1+eta) (-1-eta)]*x.';
J(1,2) = 1/4*[(eta-1) (1-eta) (1+eta) (-1-eta)]*y.';
J(2,1) = 1/4*[(xi-1) (-1-xi) (1+xi) (1-xi)]*x.';
J(2,2) = 1/4*[(xi-1) (-1-xi) (1+xi) (1-xi)]*y.';
% inverse of J
gamma = inv(J);

% Derivative of shapefunction in local space
dN = 1/4*[(eta-1) (1-eta) (1+eta) (-1-eta)
          (xi-1) (-1-xi) (1+xi) (1-xi)];

% Construct the strain displacement matrix
for i = 1:4
   j = 3*i-2;
   Bm(1,j+1) = gamma(1,1)*dN(1,i)+gamma(1,2)*dN(2,i);
   Bm(2,j+2) = gamma(2,1)*dN(1,i)+gamma(2,2)*dN(2,i);
   Bm(3,j+1) = gamma(2,1)*dN(1,i)+gamma(2,2)*dN(2,i);
   Bm(3,j+2) = gamma(1,1)*dN(1,i)+gamma(1,2)*dN(2,i);
   Bm(4,j) = -(gamma(1,1)*dN(1,i)+gamma(1,2)*dN(2,i));
   Bm(4,j+1) = NN(i,1);
   Bm(5,j) = -(gamma(2,1)*dN(1,i)+gamma(2,2)*dN(2,i));
   Bm(5,j+2) = NN(i,1);
end

end

% Function returning number of GP for three specified schemes
function [kbGauss, ksGauss] = numGaussPoints(intMethod)
    if strcmp(intMethod,'Reduced')
        kbGauss = 1;
        ksGauss = 1;
    elseif strcmp(intMethod,'Selective')
        kbGauss = 2;
        ksGauss = 1;
    elseif strcmp(intMethod,'Full')
        kbGauss = 2;
        ksGauss = 2;
    end
end

% First three sets of GP and W for a quadrilateral
function [gaussVec,w] = gaussFunc(gaussPoints)
    Gauss = [0 0 0 ;
            -1/sqrt(3) 1/sqrt(3) 0;
            -sqrt(0.6) 0 sqrt(0.6)];

     ww = [2 0 0;
           1 1 0;
           5/9 8/9 5/9];

    gaussVec = Gauss(gaussPoints,1:gaussPoints);
    w = ww(gaussPoints,1:gaussPoints);
end
