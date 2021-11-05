function c = thirdorderin(a,b)
as = size(a);
bs = size(b);

c = zeros([as(1) max(bs) as(2)]);

for i = 1:as(1)
    for j = 1:max(bs)
        for k = 1:as(2)
            c(i,j,k) = a(i,k)*b(j);
        end
    end
end
end

