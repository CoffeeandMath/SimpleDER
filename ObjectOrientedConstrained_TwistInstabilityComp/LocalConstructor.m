classdef LocalConstructor < handle
    properties
        element
        material
        E = 0;
        dEbend = zeros(11,1);
        ddEbend = zeros(11,11);
        
        globalNodes;
    end
    
    methods
         
        function Update(obj,Order)
            %obj.element.Update(Order);
            %obj.material.Update(Order);
            
            
            
            
            obj.Energy;
            if Order > 0
                obj.dEnergy;
            end
            if Order > 1
                obj.ddEnergy;
            end
        end
        
        
        function Energy(obj)
            obj.E = obj.material.Eloc;
        end
        
        function dEnergy(obj)
            
            dEbdphi = obj.element.dkappadphi'*obj.material.dEdkappaloc;
            dEbdxim = obj.element.dkappadx(:,:,1)'*obj.material.dEdkappaloc;
            dEbdxi = obj.element.dkappadx(:,:,2)'*obj.material.dEdkappaloc;
            dEbdxip = obj.element.dkappadx(:,:,3)'*obj.material.dEdkappaloc;
            
            obj.dEbend = [dEbdxim;dEbdphi(1);dEbdxi;dEbdphi(2);dEbdxip];
        end
        
        function ddEnergy(obj)
            
            %dximdxim
            obj.ddEbend(1:3,1:3) = obj.element.dkappadx(:,:,1)'*obj.material.ddEddkappaloc*obj.element.dkappadx(:,:,1) + vtimes3D(obj.material.dEdkappaloc,obj.element.ddkappadxdx(:,:,:,1,1));
            %dximdxi
            obj.ddEbend(1:3,5:7) = obj.element.dkappadx(:,:,1)'*obj.material.ddEddkappaloc*obj.element.dkappadx(:,:,2) + vtimes3D(obj.material.dEdkappaloc,obj.element.ddkappadxdx(:,:,:,1,2));
            %dxidxim
            obj.ddEbend(5:7,1:3) = obj.ddEbend(1:3,5:7)';
            %dximdxip
            obj.ddEbend(1:3,9:11) = obj.element.dkappadx(:,:,1)'*obj.material.ddEddkappaloc*obj.element.dkappadx(:,:,3) + vtimes3D(obj.material.dEdkappaloc,obj.element.ddkappadxdx(:,:,:,1,3));
            %dxipdxim
            obj.ddEbend(9:11,1:3) = obj.ddEbend(1:3,9:11)';
            %dxidxi
            obj.ddEbend(5:7,5:7) = obj.element.dkappadx(:,:,2)'*obj.material.ddEddkappaloc*obj.element.dkappadx(:,:,2) + vtimes3D(obj.material.dEdkappaloc,obj.element.ddkappadxdx(:,:,:,2,2));
            %dxidxip
            obj.ddEbend(5:7,9:11) = obj.element.dkappadx(:,:,2)'*obj.material.ddEddkappaloc*obj.element.dkappadx(:,:,3) + vtimes3D(obj.material.dEdkappaloc,obj.element.ddkappadxdx(:,:,:,2,3));
            %dxipdxi
            obj.ddEbend(9:11,5:7) = obj.ddEbend(5:7,9:11)';
            %dxipdxip
            obj.ddEbend(9:11,9:11) = obj.element.dkappadx(:,:,3)'*obj.material.ddEddkappaloc*obj.element.dkappadx(:,:,3) + vtimes3D(obj.material.dEdkappaloc,obj.element.ddkappadxdx(:,:,:,3,3));
            
            %dximdphi
            obj.ddEbend(1:3,[4,8]) = obj.element.dkappadx(:,:,1)'*obj.material.ddEddkappaloc*obj.element.dkappadphi + vtimes3D(obj.material.dEdkappaloc,squeeze(obj.element.ddkappadxdphi(:,:,1,:)));
            obj.ddEbend([4,8],1:3) = obj.ddEbend(1:3,[4,8])';
            %dxidphi
            obj.ddEbend(5:7,[4,8]) = obj.element.dkappadx(:,:,2)'*obj.material.ddEddkappaloc*obj.element.dkappadphi + vtimes3D(obj.material.dEdkappaloc,squeeze(obj.element.ddkappadxdphi(:,:,2,:)));
            obj.ddEbend([4,8],5:7) = obj.ddEbend(5:7,[4,8])';
            %dxipdphi
            obj.ddEbend(9:11,[4,8]) =  obj.element.dkappadx(:,:,3)'*obj.material.ddEddkappaloc*obj.element.dkappadphi + vtimes3D(obj.material.dEdkappaloc,squeeze(obj.element.ddkappadxdphi(:,:,3,:)));
            obj.ddEbend([4,8],9:11) = obj.ddEbend(9:11,[4,8])';
            
            
            obj.ddEbend([4,8],[4,8]) = obj.element.dkappadphi'*obj.material.ddEddkappaloc*obj.element.dkappadphi + vtimes3D(obj.material.dEdkappaloc,obj.element.ddkappadphidphi);
            
        end
        
    end
end