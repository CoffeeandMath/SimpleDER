function ddkout = ddkappadxdphifun(dnm,dn,Bnm,Bn)

%%


e3 = [0;0;0;1];

dqndphinmdxnm = -[qmult(e3,qconj(dnm),[0;Bnm(:,1)],dn),qmult(e3,qconj(dnm),[0;Bnm(:,2)],dn),qmult(e3,qconj(dnm),[0;Bnm(:,3)],dn)];
ddkout(:,:,1,1) = qvec(dqndphinmdxnm);


dqndphinmdxnp = -[qmult(e3,qconj(dnm),[0;Bn(:,1)],dn),qmult(e3,qconj(dnm),[0;Bn(:,2)],dn),qmult(e3,qconj(dnm),[0;Bn(:,3)],dn)];
ddkout(:,:,3,1) = qvec(dqndphinmdxnp);

ddkout(:,:,2,1) = -ddkout(:,:,1,1) - ddkout(:,:,3,1);


dqndphindxnm = [qmult(qconj(dnm),[0;Bnm(:,1)],dn,e3),qmult(qconj(dnm),[0;Bnm(:,2)],dn,e3),qmult(qconj(dnm),[0;Bnm(:,3)],dn,e3)];
ddkout(:,:,1,2) = qvec(dqndphindxnm);

dqndphindxnp = [qmult(qconj(dnm),[0;Bn(:,1)],dn,e3),qmult(qconj(dnm),[0;Bn(:,2)],dn,e3),qmult(qconj(dnm),[0;Bn(:,3)],dn,e3)];
ddkout(:,:,3,2) = qvec(dqndphindxnp);

ddkout(:,:,2,2) = -ddkout(:,:,1,2) - ddkout(:,:,3,2);

end

