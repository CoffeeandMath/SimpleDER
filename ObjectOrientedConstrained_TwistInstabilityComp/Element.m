classdef Element < handle
    properties
        Frame1
        Frame2
        
        globalNodes = zeros(1,11);
        
        q = zeros(4,1);
        kappa = zeros(3,1);
        dkappadphi = zeros(3,2);
        dkappadx = zeros(3,3,3);
        ddkappadphidphi = zeros(3,2,2);
        ddkappadxdphi = zeros(3,3,3,2);
        ddkappadxdx = zeros(3,3,3,3,3)
    end
    methods
        function SetFrames(obj,f1,f2)
            obj.Frame1 = f1;
            obj.Frame2 = f2;
        end
        
        function Update(obj,Order)
            %obj.Frame1.Update(Order);
            %obj.Frame2.Update(Order);
            obj.q = qmult(qconj(obj.Frame1.d),obj.Frame2.d);
            obj.kappa = 2*qvec(obj.q);
            if Order > 0
                obj.dkappadx = dkappadxfun(obj.Frame1.d,obj.Frame2.d,obj.Frame1.B,obj.Frame2.B);
                obj.dkappadphi = dkappadphifun(obj.q);
            end
            if Order > 1
                obj.ddkappadphidphi = ddkappadphidphifun(obj.q);
                obj.ddkappadxdphi = ddkappadxdphifun(obj.Frame1.d,obj.Frame2.d,obj.Frame1.B,obj.Frame2.B);
                obj.ddkappadxdx = ddkappadxdxfun(obj.Frame1.d,obj.Frame2.d,obj.Frame1.B,obj.Frame2.B,obj.Frame1.DB,obj.Frame2.DB);
            end
        end
    end
end