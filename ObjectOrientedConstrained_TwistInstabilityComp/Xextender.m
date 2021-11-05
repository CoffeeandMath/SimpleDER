function [ X ] = Xextender( Xconst,dS )
%This function takes the input vector Xconst and appends 7 values to the
%beginning of it. It adds in the the position of the first node and fixes
%gamma1 to be zero then adds the position of the second node.

x1 = [0;0;0];
x2 = [dS;0;0];
X = [x1;0;x2;Xconst];

end

