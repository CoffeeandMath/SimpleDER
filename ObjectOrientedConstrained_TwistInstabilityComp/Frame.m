classdef Frame < handle
    properties
        Node1
        Node2
        phi = 0;
        
        FrameNumber = 0;
        globalNodes = zeros(1,7);
        
        t = zeros(3,1);
        T = zeros(3,1);
        d = zeros(4,1);
        D = [1;0;0;0];
        B = zeros(3,3);
        DB = zeros(3,3,3);
        ptw = zeros(4,1);
        ppar = zeros(4,1);
        eps = 0;
        depsdx = zeros(3,2);
        ddepsdxdx = zeros(3,3,2,2);
        FreePhi = true;
    end
    methods
        function SetNodes(obj,n1,n2)
            obj.Node1 = n1;
            obj.Node2 = n2;
            obj.Update(0);
            
        end
        function Update(obj,Order)
            obj.eps = norm(obj.Node2.X - obj.Node1.X);
            obj.t = (obj.Node2.X - obj.Node1.X)/obj.eps;
            obj.T = Normalize(qvec(qmult(obj.D,[0;0;0;1],qconj(obj.D))));
            obj.ptw = qgen(cos(obj.phi/2),sin(obj.phi/2)*obj.T);
            obj.ppar = (sqrt(2)^(-1))*qgen(sqrt(1 + obj.T'*obj.t),cross(obj.T,obj.t)/sqrt(1 + obj.T'*obj.t));
            obj.d = Normalize(qmult(obj.ppar,obj.ptw,obj.D));
            if Order > 0
                obj.B = (0.5/obj.eps)*((1/(1+obj.T'*obj.t))*(obj.T*cross(obj.T,obj.t)' - cross(obj.T,obj.t)*((obj.t+obj.T)')) + skw(obj.T));
                obj.Calcdepsdx;
            end
            if Order > 1
                
                f2i = 1 + obj.T'*obj.t;
                dBideps = thirdorder(obj.B,-obj.t/obj.eps);
                Ticti = cross(obj.T,obj.t);
                dBidf2 = thirdorder(obj.T*Ticti' - Ticti*(obj.t'+obj.T'),(0.5/(obj.eps^2))*(-1/(f2i^2))*(obj.T-obj.t*(f2i-1)));
                Tiskw = skw(obj.T);
                tiprojeps = eye(3) - obj.t*obj.t';
                dBidti = (0.5/(obj.eps^2*f2i))*(thirdorder(obj.T,Tiskw - Ticti*obj.t') - thirdorderin(Tiskw - Ticti*obj.t',obj.t+obj.T) - thirdorder(Ticti,tiprojeps));
                
                obj.DB = dBidti + dBidf2 + dBideps;
                
                obj.Calcddepsdxdx;
            end
        end
        
        function Calcdepsdx(obj)
            obj.depsdx = [-obj.t,obj.t];
        end
        
        function Calcddepsdxdx(obj)
            dt = (obj.eps)^(-1)*(eye(3) - obj.t*obj.t');
            obj.ddepsdxdx(:,:,1,1) = dt;
            obj.ddepsdxdx(:,:,1,2) = -dt;
            obj.ddepsdxdx(:,:,2,1) = -dt;
            obj.ddepsdxdx(:,:,2,2) = dt;
        end
    end
end