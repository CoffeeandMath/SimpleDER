clear all
close all
%% Initialization
Nn = 50;
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
end

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

e(Nn) = Element;
ec(Nn-2) = ElementConstraintTemplate;
for i = 1:(Nn-2)
    e(i).SetFrames(f(i),f(i+1));
    e(i).globalNodes = (1:11) + 4*(i-1);
    ec(i).e = e(i);
    ec(i).globalNodes = e(i).globalNodes;
    ec(i).dS = dS;
end

m(Nn) = Material;

for i = 1:(Nn-2)
    m(i).BMod = 2*10^1*ones(3,1)/dS; m(i).BMod(3) = 10*m(i).BMod(1);
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
Constraints = cell(max(size(fc)) + max(size(ec)),1);
for i = 1:max(size(fc))
    Constraints{i} = fc(i);
end
for i = (1):(max(size(ec)))
    Constraints{i+max(size(fc))} = ec(i);
end
r.Constraints = Constraints;
%%
rho = 80/Nn;

s = Solver;
s.r = r;
s.dS = dS;
s.mu = 500/Nn;
s.mumin = s.mumin;
s.mumax = 10^4*s.mumin;
s.Nc = max(size(fc)) + max(size(ec));
s.lambda = zeros(s.Nc,1);

Nk = 200;

multl = linspace(0,1,Nk);
for i = 1:Nk
    i
    r.rho = multl(i)*rho;
    tic
    for k = 1:(Nn-2)
        m(k).kappa0 = multl(i)*[0;0;0]*dS;
    end
    x(end).Fnode = 0*multl(i)*[-300;200;300];
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
    writevtkline(['data2/data.' num2str(i) '.vtk'],xv,yv,zv)
    toc
    r.setd2D;
end
Fz = zeros(3,Nk);
Xz = zeros(3,Nk);
for i = 1:Nk
    i+Nk
    x(end).Fnode = multl(i)*600*Normalize([-1;1;1]);
    tic
    s.Solve
    s.LastNumIter
    toc

    xv = zeros(Nn,1);
    yv = zeros(Nn,1);
    zv = zeros(Nn,1);
    for k = 1:Nn
        xv(k) = x(k).X(1);
        yv(k) = x(k).X(2);
        zv(k) = x(k).X(3);
    end
    writevtkline(['data2/data.' num2str(Nk+i) '.vtk'],xv,yv,zv)
   
    r.setd2D;
    Fz(:,i) = x(end).Fnode(3);
    Xz(:,i)= x(end).X(3);
    
end

figure()
dXnorm = zeros(Nn,1);
Fnorm = zeros(Nn,1);
for i = 1:Nk
    dXnorm(i) = norm(Xz(:,i)-Xz(:,1));
    Fnorm(i) = norm(Fz(:,i));
end
plot(dXnorm,Fnorm)
xlabel('|Xt - X1|')
ylabel('F')