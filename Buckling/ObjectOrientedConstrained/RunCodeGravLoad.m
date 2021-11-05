clear all
close all
foldername = 'data7';
filename = foldername;
%% Initialization
Nn = 100;
L = 1;
Si = linspace(0,L,Nn);
dS = Si(2) - Si(1);

xinit = Si;
yinit = zeros(size(Si));
zinit = zeros(size(Si));

x(Nn) = Node;

for i = 1:Nn
    x(i).X = [xinit(i);yinit(i);zinit(i)];
    x(i).NodeNumber = i;
    x(i).globalNodes = 4*(i-1) + (1:3);
end
x(1).FreeNode = false; x(2).FreeNode = false;

f(Nn-1) = Frame;
fc(Nn-1) = Inextensibility;
pdIdoteI = 1 + [0;1;0]'*[1;0;0] + [0;0;1]'*[0;1;0] + [1;0;0]'*[0;0;1];
epsR = cross([0;1;0],[1;0;0]) + cross([0;0;1],[0;1;0]) + cross([1;0;0],[0;0;1]);
Di = qgen((1/2) * sqrt(pdIdoteI),-(1/2) * epsR/sqrt(pdIdoteI));
for i = 1:(Nn-1)
    f(i).SetNodes(x(i),x(i+1));
    f(i).D = Di;
    f(i).FrameNumber = i;
    f(i).globalNodes = (1:7) + 4*(i-1);
    fc(i).f = f(i);
    fc(i).globalNodes = f(i).globalNodes;
    fc(i).dS = dS;
end
f(1) .FreePhi = false;

e(Nn) = Element;

for i = 1:(Nn-2)
    e(i).SetFrames(f(i),f(i+1));
end

m(Nn) = Material;

for i = 1:(Nn-2)
    m(i).BMod = 2*10^1*ones(3,1)/dS; m(i).BMod(2) = 10*m(i).BMod(1); m(i).BMod(3) = 100*m(i).BMod(1);
    m(i).element = e(i);
end
lc(Nn-2) = LocalConstructor;

for i = 1:(Nn-2)
    lc(i).element = e(i);
    lc(i).material = m(i);
    %lc(i).Update(2);
    lc(i).globalNodes = [(4*i-3):(4*i+7)];
end

r = Rod;
%r.ParHessian = true;
r.lc = lc;
r.Nn = Nn;
r.Nodes = x;
r.Frames = f;
r.AxMod = 0*100000*ones(Nn-1,1);
r.dS = dS;
r.Update(2);
r.calcFreeNodeMask;
Constraints = cell(max(size(fc))-1,1);
for i = 1:max(size(fc))
    Constraints{i} = fc(i);
end
r.Constraints = Constraints(2:end);
%%
rho = 8/Nn;

s = Solver;
s.r = r;
s.dS = dS;
s.mu = 1000/Nn;
s.mumin = s.mu;
s.mumax = 10^3*s.mumin;
s.Nc = max(size(r.Constraints));
s.lambda = zeros(s.Nc,1);
s.homog = 000*dS;
Nk = 1500;

multl = linspace(0,1,Nk);
for i = 1:Nk
    i
    r.rho = (1-multl(i))*10^-1*dS;
    tic
    
    for k = 3:Nn
        x(k).Fnode = (30)*multl(i)*[-1;0;0];
    end
    
    s.Solve
    s.LastNumIter
    
    
    xv = zeros(Nn,1);
    yv = zeros(Nn,1);
    zv = zeros(Nn,1);
    for k = 1:Nn
        xv(k) = x(k).X(1);
        yv(k) = x(k).X(2);
        zv(k) = x(k).X(3);
    end
    writevtkline([foldername '/' filename '.' num2str(i) '.vtk'],xv,yv,zv)
    toc
    r.setd2D;
end
%%
Fxdist = zeros(Nk,1);
Xall = zeros(3,Nn,Nk);
for i = 1:Nk
    i+Nk
    r.rho = 0;
    tic
    
    for k = 3:Nn
        x(k).Fnode = (30)*(1-multl(i))*[-1;0;0];
    end
    Fxdist(i) = x(k).Fnode(1);
    
    s.Solve
    s.LastNumIter
    
    
    xv = zeros(Nn,1);
    yv = zeros(Nn,1);
    zv = zeros(Nn,1);
    for k = 1:Nn
        xv(k) = x(k).X(1);
        yv(k) = x(k).X(2);
        zv(k) = x(k).X(3);
        Xall(:,k,i) = x(k).X;
    end
    
    writevtkline([foldername '/' filename '.' num2str(Nk+i) '.vtk'],xv,yv,zv)
    toc
    r.setd2D;
end

figure()
plot(-Fxdist/(dS^2*m(1).BMod(1)),squeeze(Xall(1,end,:)),'.')
hold all
%plot(-Fxdist/(dS^2*m(1).BMod(1)),squeeze(Xall(2,end,:)),'.')
plot(-Fxdist/(dS^2*m(1).BMod(1)),squeeze(Xall(3,end,:)),'.')
%legend('x','y','z')