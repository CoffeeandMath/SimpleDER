 function dkappadxout = dkappadxfun(dim,di,Bim,Bi)
dkappadxout = zeros(3,3,3);

dqidxim = [qmult(qconj(dim),[0;Bim(:,1)],di),qmult(qconj(dim),[0;Bim(:,2)],di),qmult(qconj(dim),[0;Bim(:,3)],di)];
dkappadxout(:,:,1) = 2*dqidxim(2:4,:);


dqidxip = [qmult(qconj(dim),[0;Bi(:,1)],di),qmult(qconj(dim),[0;Bi(:,2)],di),qmult(qconj(dim),[0;Bi(:,3)],di)];
dkappadxout(:,:,3) = 2*dqidxip(2:4,:);

dkappadxout(:,:,2) = -dkappadxout(:,:,1) - dkappadxout(:,:,3);


end

