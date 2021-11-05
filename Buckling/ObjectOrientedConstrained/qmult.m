function qout = qmult(q1,q2,varargin)

r1 = q1(1);
v1 = q1(2:end);

r2 = q2(1);
v2 = q2(2:end);

qout = zeros(4,1);

qout(1) = r1*r2 - dot(v1,v2);
qout(2:end) = r1*v2 + r2*v1 + cross(v1,v2);
if max(size(varargin)) > 0
    for i = 1:max(size(varargin))
        qout = qmult(qout,varargin{i});
    end
    
end

end

