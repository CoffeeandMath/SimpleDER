clear all

%%

Nn = 100;

x(Nn) = Node;

xvals = linspace(0,1,Nn).^2;
yvals = zeros(Nn,1);
zvals = linspace(0,1,Nn);

for i = 1:Nn
    x(i).X = [xvals(i);yvals(i);zvals(i)];
end

f(Nn-1) = Frame;

for i = 1:(Nn-1)
    f(i).SetNodes(x(i),x(i+1));
end

e(Nn) = Element;

for i = 2:(Nn-1)
    e(i).SetFrames(f(i-1),f(i));
end

m(Nn) = Material;

for i = 2:(Nn-1)
    m(i).BMod = 10^-2*ones(3,1);
    m(i).kappa0 = zeros(3,1);
    m(i).Eloc = 0;
    m(i).dEdkappaloc = zeros(3,1);
    m(i).ddEddkappaloc = zeros(3,3);
    m(i).element = e(i);
end
lc(Nn) = LocalConstructor;

for i = 2:(Nn-1)
    lc(i).element = e(i);
    lc(i).material = m(i);
    lc(i).dEbend = zeros(11,1);
    lc(i).ddEbend = zeros(11,11);
    lc(i).Update(2);
end


r = Rod;
r.ParHessian = true;
r.lc = lc;
r.Nn = Nn;
r.Nodes = x;
r.Frames = f;
r.AxMod = ones(Nn-1,1);
r.dS = 0.01;
%tic;
xnew = [xvals;yvals';zvals];
for i = 1:80
    xnew = xnew + 0.1*rand(size(xnew));
    tic
    r.SetNodeValues(xnew);
    r.Update(2);
    H = r.ddEbend;
    toc
end
%toc

fprintf('Done \n')