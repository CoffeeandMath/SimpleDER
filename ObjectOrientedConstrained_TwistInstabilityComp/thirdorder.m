function [d] = thirdorder(a,b,varargin)
as = size(a);
bs = size(b);

if as(2) > 1
    d = zeros([as max(bs)]);
    for i = 1:as(1)
        for j = 1:as(2)
            for k = 1:max(bs)
                d(i,j,k) = a(i,j)*b(k);
            end
        end
    end
elseif bs(2) > 1
    d = zeros([max(as) bs]);
    for i = 1:max(as)
        for j = 1:bs(1)
            for k = 1:bs(2)
                d(i,j,k) = a(i)*b(j,k);
            end
        end
    end
else
    c = varargin{1};
    cs = size(c);
    d = zeros([max(as) max(bs) max(cs)]);
    for i = 1:max(as)
        for j = 1:max(bs)
            for k = 1:max(cs)
                d(i,j,k) = a(i)*b(j)*c(k);
            end
        end
    end
    
end
end

