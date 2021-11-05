classdef ElementConstraintTemplate < Constraint
    properties
        e;
        globalNodes;
        dS
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
            obj.c = obj.e.kappa(2)/obj.dS;
        end
        function Calc_dcdq(obj)
            obj.dc = (1/obj.dS)*[obj.e.dkappadx(2,:,1)';obj.e.dkappadphi(2,1);obj.e.dkappadx(2,:,2)';obj.e.dkappadphi(2,2);obj.e.dkappadx(2,:,3)'];
            
        end
        function Calc_ddcddq(obj)
            dSinv = (1/obj.dS);
            obj.ddc = zeros(11,11);
            obj.ddc(1:3,1:3) = dSinv*squeeze(obj.e.ddkappadxdx(2,:,:,1,1));
            obj.ddc(1:3,5:7) = dSinv*squeeze(obj.e.ddkappadxdx(2,:,:,1,2));
            obj.ddc(5:7,1:3) = obj.ddc(1:3,5:7)';
            obj.ddc(1:3,9:11) = dSinv*squeeze(obj.e.ddkappadxdx(2,:,:,1,3));
            obj.ddc(9:11,1:3) = obj.ddc(1:3,9:11)';
            obj.ddc(5:7,5:7) = dSinv*squeeze(obj.e.ddkappadxdx(2,:,:,2,2));
            obj.ddc(5:7,9:11) = dSinv*squeeze(obj.e.ddkappadxdx(2,:,:,2,3));
            obj.ddc(9:11,5:7) = obj.ddc(5:7,9:11)';
            obj.ddc(9:11,9:11) = dSinv*squeeze(obj.e.ddkappadxdx(2,:,:,3,3));
            
            
            obj.ddc(1:3,[4,8]) = dSinv*squeeze(obj.e.ddkappadxdphi(2,:,1,:));
            obj.ddc([4,8],1:3) = obj.ddc(1:3,[4,8])';
            obj.ddc(5:7,[4,8]) = dSinv*squeeze(obj.e.ddkappadxdphi(2,:,2,:));
            obj.ddc([4,8],5:7) = obj.ddc(5:7,[4,8])';
            obj.ddc(9:11,[4,8]) = dSinv*squeeze(obj.e.ddkappadxdphi(2,:,3,:));
            obj.ddc([4,8],9:11) = obj.ddc(9:11,[4,8])';
            
            
            obj.ddc([4,8],[4,8]) = dSinv*squeeze(obj.e.ddkappadphidphi(2,:,:));
            
        end
    end
end