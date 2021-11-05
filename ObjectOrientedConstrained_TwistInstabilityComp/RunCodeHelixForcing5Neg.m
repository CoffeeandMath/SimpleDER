%load(['a_' num2str(a) '_b_' num2str(b) '_phivsy0.mat'])
s.mu = 10000/Nn;
s.mumin = s.mu;
s.mumax = 10^3*s.mumin;
for i = 1:(Nn-1)
    f(i).D = Dold{i};
end
r.ApplyGenCoord(q0);
r.Update(2);
dphi = -dphi;
for i = Nk:-1:(ksolvemax+1)
    %Nk+i
    %r.rho = 0;
    tic
    
    
    x(end-1).X = x(end-1).X + dend/Nk;
    x(end).X = x(end).X + dend/Nk;
    f(1).phi = dphi/Nk;
    f(end).phi = dphi/Nk;
    s.Solve
    %s.LastNumIter
    fprintf(['Step: ' num2str(Nk+i) ', Elapsed Time: ' num2str(toc) '\n'])
    StrEn(i) = r.E;
    H = s.ddEAL();
    evaltemp = eigs(H(r.FreeNodeMask,r.FreeNodeMask),1,'smallestabs')
    minEval(i) = evaltemp;
    
    
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
    %toc
    r.setd2D;
    %f(end).phi
    save([num2str(Nn) '_a_' num2str(h) '_b_' num2str(w) '_phivsy0.mat'],'phivals','y0')
    plot(phivals*180/pi,minEval,'.')
    drawnow
end
figure()


plot(phivals*180/pi,-y0*pi)
xticks([0 60 120 180])
xticklabels({'0','60','120','180'})
yticks([-0.4 0 .4])

