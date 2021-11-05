 classdef Node < handle
     properties
         X = [0;0;0];
         Fnode = [0;0;0];
         NodeNumber = 0;
         globalNodes = zeros(1,3);
         FreeNode = true;
     end
 end