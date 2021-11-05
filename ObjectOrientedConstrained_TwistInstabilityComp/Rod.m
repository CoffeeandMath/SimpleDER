classdef Rod < handle
    properties
        lc
        
        ParHessian = false;
        Nodes
        Frames
        Constraints
        
        E
        Ebend
        Egrav
        Estrain
        Enode
        
        dE
        dEbend
        dEgrav
        dEstrain
        dEnode
        
        ddE
        ddEbend
        ddEstrain
        
        Nn
        Nelements
        dS
        rho = 0;
        AxMod
        
        FreeNodeMask
    end
    methods
        function Update(obj,Order)
            obj.Nelements = max(size(obj.lc));
            if and(Order>1,obj.ParHessian)
                Floc = obj.Frames;
                lcloc = obj.lc;
                parfor i = 1:(max(size(Floc)))
                    Floc(i).Update(Order);
                end
                parfor i = 1:max(size(lcloc))
                    lcloc(i).element.Update(Order);
                end
                for i = 1:max(size(lcloc))
                    lcloc(i).material.Update(Order);
                end
                for i = 1:max(size(lcloc))
                    lcloc(i).Update(Order);
                end
                
                for i = 1:max(size(obj.Constraints))
                    ci = obj.Constraints{i};
                    ci.Update(Order);
                end
            else
                for i = 1:(max(size(obj.Frames)))
                    obj.Frames(i).Update(Order);
                end
                for i = 1:(max(size(obj.lc)))
                    obj.lc(i).element.Update(Order);
                    obj.lc(i).material.Update(Order);
                end
                for i = 1:(max(size(obj.lc)))
                    obj.lc(i).Update(Order);
                end
                
                for i = 1:max(size(obj.Constraints))
                    ci = obj.Constraints{i};
                    ci.Update(Order);
                end
                %                 for i = 2:(obj.Nn-1)
                %                     obj.lc(i).Update(Order);
                %                 end
            end
            obj.Ebendcalc;
            obj.Egravcalc;
            obj.Estrcalc;
            obj.Enodecalc;
            
            obj.E = obj.Ebend + obj.Egrav + obj.Estrain - obj.Enode;
            if Order > 0
                obj.dEbendcalc;
                obj.dEstrcalc;
                obj.dEgravcalc;
                obj.dEnodecalc;
                
                obj.dE = obj.dEbend + obj.dEgrav + obj.dEstrain - obj.dEnode;
            end
            
            if Order > 1
                obj.ddEbendcalc;
                obj.ddEstrcalc;
                
                obj.ddE = obj.ddEbend + obj.ddEstrain;
            end
        end
        
        function SetNodeValues(obj,x)
            sz = size(x);
            if sz(2) == 1
                for i = 1:obj.Nn
                    obj.Nodes(i).X = x((3*i-2):(3*i));
                end
            else
                for i = 1:obj.Nn
                    obj.Nodes(i).X = x(:,i);
                end
            end
        end
        
        function SetPhiValues(obj,phiv)
            for i = 1:(obj.Nn - 1)
                obj.Frames(i).phi = phiv(i);
            end
        end
        
        function ApplyGenCoord(obj,Y)
            Yall = obj.ExtractGenCoord('All');
            Yall(obj.FreeNodeMask) = Y;
            xvals = zeros(3*obj.Nn,1);
            for i = 1:obj.Nn
                xvals((3*i-2):(3*i)) = Yall((4*i-3):(4*i-1));
            end
            obj.SetNodeValues(xvals);
            obj.SetPhiValues(Yall(4:4:end));
        end
        
        function Yc = ExtractGenCoord(obj,varargin)
            if max(size(varargin))==0
                Y = zeros(4*obj.Nn-1,1);
                for i = 1:obj.Nn
                    Y((4*i-3):(4*i-1)) = obj.Nodes(i).X;
                    if i < obj.Nn
                        Y(4*i) = obj.Frames(i).phi;
                    end
                end
                Yc = Y(obj.FreeNodeMask);
            elseif strcmp(varargin{1},'All')
                Y = zeros(4*obj.Nn-1,1);
                for i = 1:obj.Nn
                    Y((4*i-3):(4*i-1)) = obj.Nodes(i).X;
                    if i < obj.Nn
                        Y(4*i) = obj.Frames(i).phi;
                    end
                end
                Yc = Y;
            else
                fprintf('error\n')
                Yc = 0;
                
            end
        end
        
        function setd2D(obj)
            for i = 1:(obj.Nn - 1)
                obj.Frames(i).D = obj.Frames(i).d;
            end
        end
        
        %% Bending Energy
        function Ebendcalc(obj)
            obj.Ebend = sum([obj.lc(:).E]);
        end
        
        function dEbendcalc(obj)
            obj.dEbend = zeros(4*obj.Nn-1,1);
            
            for i = 1:(obj.Nelements)
                obj.dEbend(obj.lc(i).globalNodes) = obj.dEbend(obj.lc(i).globalNodes) + obj.lc(i).dEbend;
            end
        end
        
        function ddEbendcalc(obj)
            obj.ddEbend = spalloc(4*obj.Nn-1,4*obj.Nn-1,20*obj.Nn);
            
            for i = 1:(obj.Nelements)
                obj.ddEbend(obj.lc(i).globalNodes,obj.lc(i).globalNodes) = obj.ddEbend(obj.lc(i).globalNodes,obj.lc(i).globalNodes) + obj.lc(i).ddEbend;
            end
            
        end
        
        function Egravcalc(obj)
            Eg = 0;
            for i = 1:obj.Nn
                Eg = Eg + obj.rho*obj.Nodes(i).X(3);
            end
            obj.Egrav = Eg;
        end
        
        function dEgravcalc(obj)
            dEg = zeros(4*obj.Nn - 1,1);
            dEg(3:4:end) = obj.rho;
            obj.dEgrav = dEg;
        end
        
        function Estrcalc(obj)
            Estr = 0;
            for i = 1:(obj.Nn-1)
                Estr = Estr + 0.5*obj.AxMod(i)*(obj.Frames(i).eps - obj.dS)^2;
            end
            obj.Estrain = Estr;
        end
        
        function dEstrcalc(obj)
            obj.dEstrain = zeros(4*obj.Nn-1,1);
            for i = 1:(obj.Nn-1)
                obj.dEstrain((4*i-3):(4*i-1)) = obj.dEstrain((4*i-3):(4*i-1)) + obj.AxMod(i)*(obj.Frames(i).eps - obj.dS)*obj.Frames(i).depsdx(:,1);
                obj.dEstrain((4*i+1):(4*i+3)) = obj.dEstrain((4*i+1):(4*i+3)) + obj.AxMod(i)*(obj.Frames(i).eps - obj.dS)*obj.Frames(i).depsdx(:,2);
            end
        end
        
        function ddEstrcalc(obj)
            obj.ddEstrain = spalloc(4*obj.Nn-1,4*obj.Nn-1,20*obj.Nn);
            for i = 1:(obj.Nn-1)
                obj.ddEstrain((4*i-3):(4*i-1),(4*i-3):(4*i-1)) = obj.ddEstrain((4*i-3):(4*i-1),(4*i-3):(4*i-1)) + obj.AxMod(i)*obj.Frames(i).depsdx(:,1)*obj.Frames(i).depsdx(:,1)' + obj.AxMod(i)*(obj.Frames(i).eps - obj.dS)*obj.Frames(i).ddepsdxdx(:,:,1,1);
                obj.ddEstrain((4*i-3):(4*i-1),(4*i+1):(4*i+3)) = obj.ddEstrain((4*i-3):(4*i-1),(4*i+1):(4*i+3)) + obj.AxMod(i)*obj.Frames(i).depsdx(:,1)*obj.Frames(i).depsdx(:,2)' + obj.AxMod(i)*(obj.Frames(i).eps - obj.dS)*obj.Frames(i).ddepsdxdx(:,:,1,2);
                obj.ddEstrain((4*i+1):(4*i+3),(4*i-3):(4*i-1)) = obj.ddEstrain((4*i+1):(4*i+3),(4*i-3):(4*i-1)) + obj.AxMod(i)*obj.Frames(i).depsdx(:,2)*obj.Frames(i).depsdx(:,1)' + obj.AxMod(i)*(obj.Frames(i).eps - obj.dS)*obj.Frames(i).ddepsdxdx(:,:,2,1);
                obj.ddEstrain((4*i+1):(4*i+3),(4*i+1):(4*i+3)) = obj.ddEstrain((4*i+1):(4*i+3),(4*i+1):(4*i+3)) + obj.AxMod(i)*obj.Frames(i).depsdx(:,2)*obj.Frames(i).depsdx(:,2)' + obj.AxMod(i)*(obj.Frames(i).eps - obj.dS)*obj.Frames(i).ddepsdxdx(:,:,2,2);
            end
        end
        
        function Enodecalc(obj)
            obj.Enode = 0;
            for i = 1:obj.Nn
                obj.Enode = obj.Enode + obj.Nodes(i).X'*obj.Nodes(i).Fnode;
            end
        end
        function dEnodecalc(obj)
            obj.dEnode = zeros(4*obj.Nn - 1,1);
            for i = 1:obj.Nn
                obj.dEnode((4*i-3):(4*i-1)) = obj.Nodes(i).Fnode;
            end
        end
        
        function calcFreeNodeMask(obj)
            obj.FreeNodeMask = [];
            for i = 1:max(size(obj.Nodes))
                if obj.Nodes(i).FreeNode
                    obj.FreeNodeMask = [obj.FreeNodeMask, obj.Nodes(i).globalNodes];
                end
            end
            
            for i = 1:max(size(obj.Frames))
                if obj.Frames(i).FreePhi
                    obj.FreeNodeMask = [obj.FreeNodeMask,obj.Frames(i).globalNodes(4)];
                end
            end
            obj.FreeNodeMask = sort(obj.FreeNodeMask);
        end
    end
end
