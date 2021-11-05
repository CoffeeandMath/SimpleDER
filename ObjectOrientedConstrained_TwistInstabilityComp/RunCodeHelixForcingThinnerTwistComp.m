clear all
%close all
foldername = 'datathin';
filename = foldername;
%% Initialization
Nn = 350;
Nk = 500;
L = 1;
Si = linspace(0,L,Nn);
dS = Si(2) - Si(1);
h = L*.2/(pi*108);

w = L*8/(pi*108)
R = L/pi;
wstar = 2.2*sqrt(R*h)
modsc = 5;
nu = 0.3;
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
f(1).FreePhi = false;

e(Nn) = Element;

for i = 1:(Nn-2)
    e(i).SetFrames(f(i),f(i+1));
end

m(Nn) = Material;
for i = 1:(Nn-2)
    m(i).BMod = zeros(3,1);
    m(i).BMod(1) = (w^3 * h/12)*(10^modsc)/dS; m(i).BMod(2) = (h^3 * w/12)*(10^modsc)/dS;  m(i).BMod(3) = (h^3 * w/(6*(1+nu)))*(10^modsc)/dS;%m(i).BMod(3) = m(i).BMod(1) + m(i).BMod(2);
    xi = sqrt((1 - nu^2)*w^4/(60*h^2));
    m(i).xi = xi/dS;
    %m(i).BMod(1) = (a^3*b/12)*10^modsc/dS; m(i).BMod(2) = (a*b^3/12)*10^modsc/dS;  m(i).BMod(3) = (a*b^3/(6*(1+nu)))*10^modsc/dS;%m(i).BMod(3) = m(i).BMod(1) + m(i).BMod(2);
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
for i = 1:(max(size(fc))-1)
    Constraints{i} = fc(i+1);
end
r.Constraints = Constraints;
%%
rho = 0*8/Nn;

s = Solver;
s.r = r;
s.dS = dS;
s.mu = 100/Nn;
s.mumin = s.mu;
s.mumax = 10^3*s.mumin;
s.Nc = max(size(Constraints));
s.lambda = zeros(s.Nc,1);
s.homog = 0000*dS;
k0 = [0;pi*(1+0*1/Nn);0]/L;%k0 = [4;4;1];
multl = linspace(0,1,Nk);
for i = 1:Nk
    i
    %r.rho = sin(pi*multl(i))*rho;
    r.rho = multl(i)*rho;
    tic
    for k = 1:(Nn-2)
        m(k).kappa0 = multl(i)*k0*dS;
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
    d1loc = cell(Nn,1);
    d2loc = cell(Nn,1);
    for k = 1:(Nn-1)
        d1loc{k} = qvec(qmult(f(k).d,[0;1;0;0],qconj(f(k).d)));
        d2loc{k} = qvec(qmult(f(k).d,[0;0;1;0],qconj(f(k).d)));
    end
    d1loc{end} = d1loc{end-1};
    d2loc{end} = d2loc{end-1};
    
    writevtkline([foldername '/' filename '.' num2str(i) '.vtk'],xv,yv,zv,'a',h*ones(Nn,1),'b',w*ones(Nn,1),'d1',d1loc,'d2',d2loc)
    toc
    r.setd2D;
    
end
save([num2str(Nn) 'halfway.mat'])
%%
load([num2str(Nn) 'halfway.mat'])
foldername = 'datathin';
filename = foldername;
Nk = 1000;
h = L*.2/(pi*108);

w = L*8/(pi*108)
R = L/pi;
wstar = 2.2*sqrt(R*h)
modsc = 5;
nu = 0.3;

for i = 1:(Nn-2)
    m(i).BMod = zeros(3,1);
    m(i).BMod(1) = (w^3 * h/12)*(10^modsc)/dS; m(i).BMod(2) = (h^3 * w/12)*(10^modsc)/dS;  m(i).BMod(3) = (h^3 * w/(6*(1+nu)))*(10^modsc)/dS;%m(i).BMod(3) = m(i).BMod(1) + m(i).BMod(2);
    xi = sqrt((1 - nu^2)*w^4/(60*h^2));
    m(i).xi = (1)*xi/dS;
    %m(i).BMod(1) = (a^3*b/12)*10^modsc/dS; m(i).BMod(2) = (a*b^3/12)*10^modsc/dS;  m(i).BMod(3) = (a*b^3/(6*(1+nu)))*10^modsc/dS;%m(i).BMod(3) = m(i).BMod(1) + m(i).BMod(2);
    m(i).element = e(i);
end
xi

s.mu = .1*s.mu;
s.mumin = s.mu;
s.mumax = 1*s.mumax;
nc1 = NodeConstraint;
nc1.n = x(end);
nc1.X0 =  x(1).X + [0;2*L/pi;0];
nc1.globalNodes = x(end).globalNodes;



%r.Constraints = [r.Constraints;nc1;nc2];
r.Constraints{end+1} = nc1;

s.Nc = max(size(r.Constraints));
s.lambda = [s.lambda;0];
tic
s.Solve
toc


r.Constraints = r.Constraints(1:(end-2));
%s.homog = 1000*dS;
s.Nc = max(size(r.Constraints));
s.lambda = s.lambda(1:(end-2));
s.lambda
%s.mu = 100*s.mu;
%s.mumin = 100*s.mumin;
%s.mumax = 100*s.mumax;
x(end).FreeNode = false; x(end-1).FreeNode = false;f(end).FreePhi = false; f(end).phi = 0;
x(end-1).X = x(2).X + [0;2*L/pi;0];
x(end).X = x(1).X + [0;2*L/pi;0];

r.calcFreeNodeMask;
s.Solve
q0 = r.ExtractGenCoord;
Dold = cell(Nn-1,1);
for i = 1:max(size(f))
    Dold{i} = f(i).D;
end
% s.mu = 10000/Nn;
% s.mumin = s.mu;
% s.mumax = 10^3*s.mumin;
s.homog = 0*.03;
Xall = zeros(3,Nn,Nk);
dend = [0;0;0];
dphi = pi;
for k = 1:(Nn-2)
    m(k).kappa0 =0* 0.1*k0*dS;
end
y0 = zeros(Nk,1);
gamma = 0;
phivals = (dphi/Nk):(dphi/(Nk)):(dphi);
ksolvemax = 0;
StrEn = zeros(Nk,1);
minEval = zeros(Nk,1);

i = 0;
figure()
while i < Nk
    i = i+1;
    % Nk+i
    %r.rho = 0;
    tic
    
    
    
    f(1).phi = f(1).phi + dphi/Nk;
    f(end).phi = f(end).phi + dphi/Nk;
    %     f(1).phi = dphi/Nk;
    %     f(end).phi = dphi/Nk;
    s.Solve
    
    %s.LastNumIter
    StrEn(i) = r.E;
    H = s.ddEAL();
    evaltemp = eigs(H(r.FreeNodeMask,r.FreeNodeMask),1,'smallestabs')
    minEval(i) = evaltemp;
    %toc
    fprintf(['Step: ' num2str(Nk+i) ', Elapsed Time: ' num2str(toc) '\n'])
    
    xv = zeros(Nn,1);
    yv = zeros(Nn,1);
    zv = zeros(Nn,1);
    for k = 1:Nn
        xv(k) = x(k).X(1);
        yv(k) = x(k).X(2);
        zv(k) = x(k).X(3);
        Xall(:,k,i) = x(k).X;
    end
    y0(i) = zv(floor(Nn/2));
    d1loc = cell(Nn,1);
    d2loc = cell(Nn,1);
    for k = 1:(Nn-1)
        d1loc{k} = qvec(qmult(f(k).d,[0;1;0;0],qconj(f(k).d)));
        d2loc{k} = qvec(qmult(f(k).d,[0;0;1;0],qconj(f(k).d)));
    end
    d1loc{end} = d1loc{end-1};
    d2loc{end} = d2loc{end-1};
    
    writevtkline([foldername '/' filename '.' num2str(i) '.vtk'],xv,yv,zv,'a',h*ones(Nn,1),'b',w*ones(Nn,1),'d1',d1loc,'d2',d2loc)
    
    %r.setd2D;
    %     if mod(i,50)==0
    %         r.setd2D;
    %         f(1).phi = 0;
    %         f(end).phi = 0;
    %     end
    %f(end).phi
    
    %     if s.LastNumIter < 45
    %         ksolvemax = i;
    %         save([num2str(Nn) '_a_' num2str(h) '_b_' num2str(w) '_phivsy0.mat'],'phivals','y0')
    %     else
    %         for k = 1:max(size(f))
    %             f(end).phi = 0;
    %         end
    %         RunCodeHelixForcing5Neg;
    %         i = Nk+1;
    %
    %     end
    if i > 1
    evalslope = minEval(1)^(-1)*(evaltemp - minEval(i-1))/(dphi)
    if and(evalslope < -.003,evaltemp < .1*minEval(1))
        for k = 1:max(size(f))
            f(end).phi = 0;
        end
        ksolvemax = i;
        RunCodeHelixForcing5Neg;
        i = Nk+1;
    end
    end
    plot(phivals*180/pi,minEval,'.')
    drawnow
end
figure()


plot(phivals*180/pi,-y0*pi/L)
xticks([0 60 120 180])
xticklabels({'0','60','120','180'})
ylim([-0.9 0.5])
yticks([-0.8 -0.4 0 .4 0.8])


figure
yyaxis left
plot(phivals(abs(StrEn)<1)*180/pi,StrEn(abs(StrEn)<1),'.')
ylabel('Strain Energy')

yyaxis right
plot(phivals(abs(minEval)<10^-3)*180/pi,minEval(abs(minEval)<10^-4),'.')
ylabel('Smallest Eigenvalue')
xlabel('\phi')