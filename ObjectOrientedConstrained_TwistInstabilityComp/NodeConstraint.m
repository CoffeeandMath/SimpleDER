classdef NodeConstraint < Constraint
    properties
        n;
        X0;
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
            obj.c = (1/2)*norm(obj.n.X - obj.X0)^2;
           
        end
        function Calc_dcdq(obj)
            obj.dc = obj.n.X-obj.X0;
        end
        function Calc_ddcddq(obj)
            obj.ddc = eye(3);
        end
        
        
    end
end