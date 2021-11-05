classdef Solver < handle
    properties
        r
        Nc = 0;
        dS
        Xmin
        LastNumIter
        lambda
        mu
        mumin = 100;
        mumax = 10^6;
        homog = 10;
        q0;
        q0full;
    end
    methods
        function Solve(obj)
            if obj.Nc == 0
                obj.NR
            else
                obj.AugmentedLagrangian
            end
        end
        
        function NR(obj)
            dfnorm = 10;
            tol = sqrt(obj.r.Nn)*10^-6;
            kmax = 100;
            obj.q0 = obj.r.ExtractGenCoord;
            qn = obj.q0;
            k = 0;
            obj.r.Update(2);
            while and(dfnorm > tol, k < kmax)
                k = k+1;
                
                
                Amat = obj.r.ddE(obj.r.FreeNodeMask,obj.r.FreeNodeMask);
                
                Eold = obj.r.E;
                stepdir = - (Amat\obj.r.dE(obj.r.FreeNodeMask));
                qntemp = qn + stepdir;
                obj.r.ApplyGenCoord(qntemp);
                obj.r.Update(0);
                while obj.r.E > Eold
                    stepdir = 0.95*stepdir;
                    qntemp = qn + stepdir;
                    obj.r.ApplyGenCoord(qntemp);
                    obj.r.Update(0);
                end
                qn = qntemp;
                obj.r.ApplyGenCoord(qn);
                
                obj.r.Update(2);
                dfnorm = norm(obj.r.dE(obj.r.FreeNodeMask));
            end
            obj.LastNumIter = k;
        end
        
        
        function NRLS(obj)
            
            dfnorm = 10;
            tol = sqrt(obj.r.Nn)*10^-6;
            kmax = 100;
            obj.q0 = Xconstrainer(obj.r.ExtractGenCoord);
            qn = obj.q0;
            k = 0;
            Nls = 10;
            alphai = linspace(-0.1,1.5,Nls);
            %alphai = 1;
            while and(dfnorm > tol, k < kmax)
                k = k+1;
                if k > 1
                    obj.r.ApplyGenCoord(Xextender(qn,obj.dS));
                end
                obj.r.Update(2);
                Amat = Hconstrainer(obj.r.ddE);
                
                stepdir = - (Amat\Xconstrainer(obj.r.dE));
                
                if or(norm(1.0*isnan(stepdir))>0,norm(stepdir)==Inf)
                    stepdir = - Xconstrainer(obj.r.dE);
                end
                fval = zeros(Nls,1);
                
                for i = 1:Nls
                    obj.r.ApplyGenCoord(Xextender(qn+alphai(i)*stepdir,obj.dS));
                    obj.r.Update(0);
                    fval(i) = obj.r.E;
                    if or(fval(i)> 10^8, isnan(fval(i)))
                        fval(i) = Inf;
                    end
                    
                end
                [~,I] = min(fval);
                
                step = alphai(I)*stepdir;
                
                stepn = norm(step);
                qn = qn+step;
                obj.r.ApplyGenCoord(Xextender(qn+step,obj.dS));
                obj.r.Update(1);
                while (norm(1.0*isnan(obj.r.dE)) > 1)
                    step = step/2;
                    obj.r.ApplyGenCoord(Xextender(qn+step,obj.dS));
                    obj.r.Update(1);
                end
                qn = qn + step;
                obj.r.ApplyGenCoord(Xextender(qn,obj.dS));
                obj.r.Update(1);
                dfnorm = norm(Xconstrainer(obj.r.dE));
            end
            obj.LastNumIter = k;
        end
        
        function AugmentedLagrangian(obj)
            tau0 = 10^-7;
            ALnorm = 10;
            obj.mu = max(obj.mu - (1/5)*obj.mumin,obj.mumin);
            %obj.mu = obj.mumin;
            %obj.lambda = zeros(size(obj.lambda));
            q0init = obj.r.ExtractGenCoord;
            kal = 0;
            while ALnorm > tau0
                if ALnorm > 10
                    fprintf('something went wrong\n');
                    obj.mu = 100*obj.mu;
                    obj.mumin = 100*obj.mumin;
                    obj.mumax = 100*obj.mumax;
                    obj.lambda = zeros(size(obj.lambda));
                    obj.r.ApplyGenCoord(q0init);
                end
                obj.q0 = obj.r.ExtractGenCoord;
                obj.q0full = obj.r.ExtractGenCoord('All');
                kal = kal+1;
                
                dfnorm = 10;
                tol = sqrt(obj.r.Nn)*10^-8;
                kmax = 50;
                
                qn = obj.r.ExtractGenCoord;
                k = 0;
                obj.r.Update(2);
                
                
                while and(dfnorm > tol, k < kmax)
                    k = k+1;
                    
                    ddEfull = obj.ddEAL();
                    Amat = ddEfull(obj.r.FreeNodeMask,obj.r.FreeNodeMask);
                    Amat = 0.5*(Amat + transpose(Amat));
                    Eold = obj.EAL;
                    dEfull = obj.dEAL();
                    stepdir = - (Amat\dEfull(obj.r.FreeNodeMask));
%                     if and(mod(k,10)==0,k>20)
%                         [~,lambdamag] = eigs(Amat,1,'largestabs');
%                         for ii = 1:10
%                          
%                          dE = dEfull(obj.r.FreeNodeMask);
%                          obj.r.ApplyGenCoord(obj.r.ExtractGenCoord - lambdamag^-1 * dE);
%                          obj.r.Update(1);
%                         end
%                         
%                     end
                    qntemp = qn + stepdir;
                    obj.r.ApplyGenCoord(qntemp);
                    obj.r.Update(0);
                    %                     kls = 0;
                    %                     while and(and(obj.EAL > Eold,abs(obj.EAL-Eold) > 10^-15),norm(obj.dEAL)>10^-5)
                    %                         kls = kls+1;
                    %                         stepdir = -0.95*stepdir;
                    %                         qntemp = qn + stepdir;
                    %                         obj.r.ApplyGenCoord(Xextender(qntemp,obj.dS));
                    %                         obj.r.Update(1);
                    %                         obj.EAL
                    %                         kls
                    %                     end
                    qn = qntemp;
                    obj.r.ApplyGenCoord(qn);
                    
                    %obj.r.Update(1);
                    %obj.r.Update(2);
                    obj.r.Update(2);
                    dEfull = obj.dEAL();
                    dfnorm = norm(dEfull(obj.r.FreeNodeMask));
                end
                obj.LastNumIter = k;
                
                %update constraint error
                ALnorm = obj.sum_c();
                %update lagrange multipliers
                for i = 1:obj.Nc
                    rc = obj.r.Constraints{i};
                    obj.lambda(i) = obj.lambda(i) - obj.mu*rc.c;
                end
                obj.mu = min(obj.mu + obj.mumin/10,obj.mumax);
                
                
            end
            %obj.mu
            
            
        end
        
        
        function Eout = EAL(obj,varargin)
            if max(size(varargin)) == 0
                
            else
                obj.r.ApplyGenCoord(varargin{1});
                obj.r.Update(0);
            end
            Eout = obj.r.E + 0.5*obj.homog*(norm(obj.r.ExtractGenCoord-obj.q0)^2);
            for i = 1:obj.Nc
                rc = obj.r.Constraints{i};
                Eout = Eout - obj.lambda(i)* rc.c + 0.5*obj.mu*(rc.c)^2;
            end
        end
        
        function dEout = dEAL(obj,varargin) 
            if max(size(varargin)) == 0
                
            else
                obj.r.ApplyGenCoord(varargin{1});
                obj.r.Update(1);
            end
            dEout = obj.r.dE + obj.homog*(obj.r.ExtractGenCoord('All')-obj.q0full);
            for i = 1:obj.Nc
                rc = obj.r.Constraints{i};
                dEout(rc.globalNodes) = dEout(rc.globalNodes)  - obj.lambda(i)* rc.dc + obj.mu*rc.c*rc.dc;
            end
        end
        
        function ddEout = ddEAL(obj,varargin)
            if max(size(varargin)) == 0
                
            else
                obj.r.ApplyGenCoord(varargin{1});
                obj.r.Update(2);
            end
            ddEout = obj.r.ddE + obj.homog*eye(size(obj.r.ddE));
            for i = 1:obj.Nc
                rc = obj.r.Constraints{i};
                ddEout(rc.globalNodes,rc.globalNodes) = ddEout(rc.globalNodes,rc.globalNodes)  - obj.lambda(i)* rc.ddc + obj.mu*(rc.dc*rc.dc' + rc.c*rc.ddc);
            end
            ddEout = sparse(ddEout);
        end
        
        function cout = sum_c(obj)
            cout = 0;
            
            for i = 1:obj.Nc
                rc = obj.r.Constraints{i};
                cout = cout + (rc.c)^2;
            end
        end
        
        
        
    end
    
end