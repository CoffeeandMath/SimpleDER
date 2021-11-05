
xmesh = linspace(0,1,50);

f0 = pi^2/4+1;
Niter = 1000;
fval1 = linspace(f0,0,Niter)';
xend1 = zeros(Niter,1); 
yend1 = zeros(Niter,1);
%figure()
opts = bvpset('RelTol',10^-7,'AbsTol',10^-8);
for i = 1:Niter
    i
    if i == 1
        solinit = bvpinit(xmesh, @(X) guess(X,fval1(i)));
    else
        guesslast = @(x) [spline(sol.x,sol.y(1,:),x) spline(sol.x,sol.y(2,:),x)];
        solinit = bvpinit(xmesh, guesslast);
    end
    bvpfun = @(x,y) bvpfcn(x,y,fval1(i));
    
    sol = bvp4c(bvpfun, @bcfun, solinit,opts);
%     plot(sol.x, sol.y(1,:), '-o')
%     title(num2str(fval1(i)))
%     ylim([-pi, pi])
%     drawnow
    
    yend1(i) = trapz(sol.x,sin(sol.y(1,:)));
    xend1(i) = trapz(sol.x,cos(sol.y(1,:)));
end

fval2 = linspace(f0,100,Niter)';
xend2 = zeros(Niter,1); 
yend2 = zeros(Niter,1);
for i = 1:Niter
    i
    if i == 1
        solinit = bvpinit(xmesh, @(X) guess(X,fval2(i)));
    else
        guesslast = @(x) [spline(sol.x,sol.y(1,:),x) spline(sol.x,sol.y(2,:),x)];
        solinit = bvpinit(xmesh, guesslast);
    end
    bvpfun = @(x,y) bvpfcn(x,y,fval2(i));
    
    sol = bvp4c(bvpfun, @bcfun, solinit,opts);
%     plot(sol.x, sol.y(1,:), '-o')
%     title(num2str(fval2(i)))
%     ylim([-pi, pi])
%     drawnow
    
    yend2(i) = trapz(sol.x,sin(sol.y(1,:)));
    xend2(i) = trapz(sol.x,cos(sol.y(1,:)));
end
%figure()
fvalsall = [flipud(fval1(2:end)); fval2];
xvalsall = [flipud(xend1(2:end));xend2];
yvalsall = -[flipud(yend1(2:end));yend2];
xhan = @(F) interp1(fvalsall,xvalsall,F);
yhan = @(F) interp1(fvalsall,yvalsall,F);
sp =200;
plot(linspace(min(fvalsall),max(fvalsall),sp),xhan(linspace(min(fvalsall),max(fvalsall),sp)),'x','MarkerSize',8)
hold all
plot(linspace(min(fvalsall),max(fvalsall),sp),yhan(linspace(min(fvalsall),max(fvalsall),sp)),'x','MarkerSize',8)
%xlim([-30,0])
%ylim([-1,1])
%legend('xend','yend')