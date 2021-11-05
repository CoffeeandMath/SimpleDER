function qout = qconj(q)
% 
% qout = zeros(4,1);
% qout(1) = q(1);
% qout(2:end) = -q(2:end);
qout = [q(1);-q(2:end)];
end

