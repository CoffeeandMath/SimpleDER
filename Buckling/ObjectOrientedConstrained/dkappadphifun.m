function [dkappadphiout] = dkappadphifun(q_i)

dkappadphiout = zeros(3,2);
%%dkappadphiim
dkappai = -qmult([0;0;0;1],q_i);


dkappadphiout(:,1) = dkappai(2:end);

%% dkappadphii

dkappai = qmult(q_i,[0;0;0;1]);
dkappadphiout(:,2) = dkappai(2:end);

end

