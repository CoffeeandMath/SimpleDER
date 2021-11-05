function dydx = bvpfcn(s,y,fend)
dydx = zeros(2,1);
dydx = [y(2);fend*(s-1)*sin(y(1))];
end

