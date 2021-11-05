function qout = qmult(q1,q2,varargin)



% qout = zeros(4,1);
% 
% 
% 
% qout(1) = q1(1)*q2(1) - q1(2:end)'*q2(2:end);
% qout(2:end) = q1(1)*q2(2:end) + q2(1)*q1(2:end) + [q1(3)*q2(4) - q1(4)*q2(3);q1(4)*q2(2) - q1(2)*q2(4);q1(2)*q2(3) - q1(3)*q2(2)];

qout = [q1(1)*q2(1) - q1(2:end)'*q2(2:end);q1(1)*q2(2:end) + q2(1)*q1(2:end) + [q1(3)*q2(4) - q1(4)*q2(3);q1(4)*q2(2) - q1(2)*q2(4);q1(2)*q2(3) - q1(3)*q2(2)]];
if max(size(varargin)) > 0
    for i = 1:max(size(varargin))
        qout = qmult(qout,varargin{i});
    end
    
end

end

