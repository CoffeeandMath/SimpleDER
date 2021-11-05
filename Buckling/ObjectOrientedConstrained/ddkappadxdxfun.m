function ddkout = ddkappadxdxfun(dim,di,Bim,Bi,DBim,DBi)


%% ddkappaIdxIndxKm
%upper case are components, lower case are nodes
% 1)component of kappa
% 2)component of x
% 3)component of x
% 4)node x
% 5)node x
ddkout = zeros(3,3,3,3,3);



 
dBimBimd = zeros(4,3,3);
dBimBid = zeros(4,3,3);
dBiBid = zeros(4,3,3);

cdim = qconj(dim);

for n = 1:3
    for q = 1:3
        
        dBimBimd(:,n,q) = qmult(cdim,qmult([0;Bim(:,q)],[0;Bim(:,n)]) - [0;DBim(:,n,q)],di);
        dBimBid(:,n,q) = qmult(cdim,[0;Bim(:,n)],[0;Bi(:,q)],di);
        dBiBid(:,n,q) = qmult(cdim,qmult([0;Bi(:,n)],[0;Bi(:,q)]) + [0;DBi(:,n,q)],di);
    end
end
ddkout(:,:,:,1,1) = 2*dBimBimd(2:4,:,:);
ddkout(:,:,:,1,3) = 2*dBimBid(2:4,:,:);
ddkout(:,:,:,3,3) = 2*dBiBid(2:4,:,:);
ddkout(:,:,:,3,1) = permute(ddkout(:,:,:,1,3),[1,3,2]);
ddkout(:,:,:,1,2) = -ddkout(:,:,:,1,1) - ddkout(:,:,:,1,3);
ddkout(:,:,:,2,1) = permute(ddkout(:,:,:,1,2),[1,3,2]);
ddkout(:,:,:,2,3) = -ddkout(:,:,:,1,3) - ddkout(:,:,:,3,3);
ddkout(:,:,:,3,2) = permute(ddkout(:,:,:,2,3),[1,3,2]);
ddkout(:,:,:,2,2) = -ddkout(:,:,:,2,3) - ddkout(:,:,:,2,1);



end

