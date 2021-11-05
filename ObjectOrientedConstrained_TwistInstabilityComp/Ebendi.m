function [Etot] = Ebendi(thetai,Sys)

D = Sys.D;
xi = Sys.xi;

if isfield(Sys,'kappa0')
    kappa0 = Sys.kappa0;
else
    kappa0 = zeros(size(thetai));
end

dkappasquared = (thetai - kappa0).^2;
erib = zeros(max(size(thetai)),1);
for i = 1:(max(size(thetai)))
    erib(i) = 0.5*D(2)*dkappasquared(i,3)^4/(1/(xi^2) + dkappasquared(i,2));
end

Etot = 0.5*sum(sum(D.*dkappasquared)) + sum(erib);

end

