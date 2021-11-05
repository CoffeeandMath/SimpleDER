classdef Material < handle
    properties
        element
        
        BMod
        kappa0 = [0;0;0];
        Eloc = 0;
        dEdkappaloc = zeros(3,1);
        ddEddkappaloc = zeros(3,3);
        xi = 0;
    end
    methods
        function Update(obj,Order)
            obj.Energy;
            if Order > 0
                obj.dEnergy;
            end
            if Order > 1
                obj.ddEnergy;
            end
        end
        function Energy(obj)
            
            dkappa = obj.element.kappa - obj.kappa0;
            obj.Eloc = 0.5*obj.BMod'*(dkappa.^2) + 0.5*obj.BMod(2)*(dkappa(3)^4/((1/obj.xi^2) + dkappa(2)^2));
            
        end
        
        function dEnergy(obj)
            dkappa = obj.element.kappa - obj.kappa0;
            obj.dEdkappaloc = obj.BMod.*(dkappa) + obj.BMod(2)*[0;-dkappa(2)*dkappa(3)^4/(((1/obj.xi^2) + dkappa(2)^2)^2);2*dkappa(3)^3/((1/obj.xi^2)+dkappa(2)^2)];
            
        end
        
        function ddEnergy(obj)
            dkappa = obj.element.kappa - obj.kappa0;
            DDEstr22 = obj.BMod(2)*obj.xi^4*(3*obj.xi^2 * dkappa(2)^2 - 1)*dkappa(3)^4 / ((1 + obj.xi^2 * dkappa(2)^2)^3);
            DDEstr23 = -4*obj.BMod(2)*dkappa(2)*dkappa(3)^3/((1/obj.xi^2 + dkappa(2)^2)^2);
            DDEstr33 = 6*obj.BMod(2)*dkappa(3)^2/((1/obj.xi)^2 + dkappa(2)^2);
            obj.ddEddkappaloc = diag(obj.BMod) + [0,0,0;...
                0, DDEstr22, DDEstr23;...
                0, DDEstr23, DDEstr33];
        end
    end
end