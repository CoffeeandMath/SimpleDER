function mout = vtimes3D(v,M)
% multiples vector v times a 3 dimensional matrix and returns a 2D matrix
N = max(size(v));
sM = size(M);
mout = zeros(sM(2:end));

for i = 1:N
    mout = mout + squeeze(v(i)*M(i,:,:));
end
end

