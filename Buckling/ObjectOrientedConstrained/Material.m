classdef Material < handle
    properties
        element
        
        BMod
        kappa0 = [0;0;0];
        Eloc = 0;
        dEdkappaloc = zeros(3,1);
        ddEddkappaloc = zeros(3,3);
        
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
                obj.Eloc = 0.5*obj.BMod'*(dkappa.^2);

        end
        
        function dEnergy(obj)

                obj.dEdkappaloc = obj.BMod.*(obj.element.kappa - obj.kappa0);

        end
        
        function ddEnergy(obj)
                obj.ddEddkappaloc = diag(obj.BMod);
        end
    end
end