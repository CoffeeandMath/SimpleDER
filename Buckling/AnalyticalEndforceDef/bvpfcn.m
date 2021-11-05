function dydx = bvpfcn(s,y,fend)
dydx = zeros(2,1);
dydx = [y(2);-fend*sin(y(1))];
end

