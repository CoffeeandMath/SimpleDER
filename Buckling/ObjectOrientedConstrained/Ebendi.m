function [Etot] = Ebendi(thetai,Sys)

D = Sys.D;

if isfield(Sys,'kappa0')
    kappa0 = Sys.kappa0;
else
    kappa0 = zeros(size(thetai));
end

dkappasquared = (thetai - kappa0).^2;

Etot = 0.5*sum(sum(D.*dkappasquared));

end

