classdef FrameConstraintTemplate < Constraint
    properties
        f;
        globalNodes;       
        dS;
    end
    methods
        
        function Update(obj,Order)
            obj.Calc_c;
            if Order > 0
                obj.Calc_dcdq;
            end
            if Order > 1
                obj.Calc_ddcddq;
            end
        end
        
        function Calc_c(obj)
            obj.c = obj.f.eps/obj.dS - 1;
        end
        function Calc_dcdq(obj)
            obj.dc = (1/obj.dS)*[obj.f.depsdx(:,1);0;obj.f.depsdx(:,2)];
        end
        function Calc_ddcddq(obj)
            obj.ddc = (1/obj.dS)*[obj.f.ddepsdxdx(:,:,1,1),zeros(3,1),obj.f.ddepsdxdx(:,:,1,2);...
                zeros(1,7);...
                obj.f.ddepsdxdx(:,:,2,1),zeros(3,1),obj.f.ddepsdxdx(:,:,2,2)];
        end
        
        
    end
end