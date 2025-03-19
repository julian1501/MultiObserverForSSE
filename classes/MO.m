classdef MO
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        sys % the system that this mo observes
        attack % the attack that is used on the system
        numOutputs % the number of sensors sensing the system
        numObservers % number of observers in this multi observer
        numOutputsObserver % number of outputs that each observer has
        eigenvalues % eigenvalues of the observer
        CSet % Output matrix of each observer i
        CSetIndices % Matrix containing the origin of each Ci
        A % Matrix containing (possibly different) A matrices for each observer
        L % Matrix containing the L matrices for each observer
        K % Matrix containing the K matrices for each observer
        mua
    end

    methods
        function obj = MO(sys,attack,numOutputsObserver,LMIconsts)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.sys = sys;
            obj.attack = attack;
            obj.numOutputs = obj.sys.nx;
            obj.A = obj.sys.A;

            if 2*obj.attack.numAttacks > obj.numOutputs
                warning("More than half of the sensors are attacked 2*numAttacks (%d) > numOutputs (%d)",...
                        obj.attack.numAttacks,obj.numOutputs);
            end
            
            obj.numOutputsObserver  = numOutputsObserver;
            obj.numObservers        = nchoosek(obj.numOutputs,obj.numOutputsObserver);
            
            [obj.CSet,obj.CSetIndices] = obj.CSetSetup();

            obj.eigenvalues = -1:-1:-obj.numOutputs;
            [obj.L, obj.K] = obj.defineObservers(LMIconsts,true);

            
        end

        function [CSet,CSetIndices] = CSetSetup(obj)
            
            CSetIndices = nchoosek(1:obj.numOutputs,obj.numOutputsObserver);
            CSet = zeros(obj.numOutputsObserver,...
                         obj.sys.nx,...
                         obj.numObservers);
            for i = 1:obj.numObservers
                CSet(:,:,i) = obj.sys.C(CSetIndices(i,:),:);
            end
        end

        function [Li, Ki] = defineObservers(obj,LMIconsts,diagnosticCheck)
            % This function defines the observers based on the LMI as in
            % section 5 of (Chong, 2025, Observer design for nonlinear 
            % systems with asynchronously sampled and maliciously corrupted
            % measurements)
            Li = zeros(obj.sys.nx,obj.numOutputsObserver,obj.numObservers);
            Ki = zeros(obj.sys.nx,obj.numOutputsObserver,obj.numObservers);
            ebar = 1; % max slope nonlinearity
            ebarM = diag(repmat(ebar^-1,1,obj.numOutputs));

            for i = 1:1:obj.numObservers
                % Generate observer construction matrix (OCM) that has to
                % satisfy LMI: OCM =< 0
                Ci = obj.CSet(:,:,i);
                setlmis([])
                L = lmivar(2,[obj.sys.nx obj.numOutputsObserver]);
                K = lmivar(2,[obj.sys.nx obj.numOutputsObserver]);

                Rs = lmivar(1,[obj.sys.nx 1]);
                Nu = lmivar(1,[1 0]);
                Ds = lmivar(1,repmat([1 0],obj.numOutputs,1));
                Mua = lmivar(1,[1 0]);
                Mud = lmivar(1,[1 0]);

                MuaLMI = newlmi; % 0 < Mua
                lmiterm([-MuaLMI 1 1 Mua],eye(obj.sys.nx),1);
                lmiterm([MuaLMI 1 1 0],LMIconsts.mua*eye(obj.sys.nx));

                MudLMIlo = newlmi; % 0 < Mud
                lmiterm([-MudLMIlo 1 1 Mud],eye(obj.sys.nx),1);

                MudLMIhi = newlmi; % Mud < mud
                lmiterm([MudLMIhi 1 1 Mud],eye(obj.sys.nx),1);
                lmiterm([-MudLMIhi 1 1 0],LMIconsts.mud*eye(obj.sys.nx));

                DsLMI = newlmi; % 0 < Ds
                lmiterm([-DsLMI 1 1 Ds],1,1);

                RsLMI = newlmi; % 0 < Rs
                lmiterm([-RsLMI 1 1 Rs],1,1);
                
                NuLMI = newlmi; % 0 < Mu
                lmiterm([-NuLMI 1 1 Nu],eye(obj.sys.nx),1);
                
                LKLMI = newlmi; % Whole big thing
                lmiterm([LKLMI 1 1 Rs],1,obj.A,'s');
                lmiterm([LKLMI 1 1 L],Rs,Ci,'s');
                lmiterm([LKLMI 1 1 Nu],1,1);

                lmiterm([LKLMI 1 2 Rs],obj.sys.B,1);
                lmiterm([LKLMI 1 2 -Ds],obj.sys.C',1);
                lmiterm([LKLMI 1 2 -K],Ci',Ds');

                lmiterm([LKLMI 1 3 Rs],-1,1);
                lmiterm([LKLMI 1 4 Rs],1,1);

                lmiterm([LKLMI 2 2 Ds],-2,ebarM);
                lmiterm([LKLMI 3 3 Mua],-1,1);
                lmiterm([LKLMI 4 4 Mud],-1,1);
                
                % construct lmi and solve
                LMISYS = getlmis;
                [tmin,xfeas] = feasp(LMISYS,[0,0,0,0,1]);

                assert(tmin <= 0,"Observer %d is not observable. tmin = %2.4f\n",i,tmin);
                
                L = dec2mat(LMISYS,xfeas,L);
                K = dec2mat(LMISYS,xfeas,K);
                Li(:,:,i) = L;
                Ki(:,:,i) = K;

                if diagnosticCheck
                    Rs  = dec2mat(LMISYS,xfeas,Rs);
                    assert(any(eig(Rs)>0),"Observer %d, Rs eigenvalues not larger than zero",i)
                    assert(any(Rs == Rs',"all"),"Observer %d, Rs is not symmetric",i)
                    Nu  = dec2mat(LMISYS,xfeas,Nu);
                    assert(Nu>=0,"Observer %d, Nu (%.2f) is not larger than zero",i,Nu)
                    Ds  = dec2mat(LMISYS,xfeas,Ds);
                    assert(any(eig(Ds)>0),"Observer %d, Ds is not positive definite",i)
                    Mua = dec2mat(LMISYS,xfeas,Mua);
                    assert(Mua>=LMIconsts.mua,"Observer %d, Mu_a (%.2f) is not larger than the desired size mua (%.2f)",i,Mua,LMIconsts.mua)
                    Mud = dec2mat(LMISYS,xfeas,Mud);
                    assert(Mud>=0,"Observer %d, Mu_d (%2.f) is not larger than zero",i,Mud)
                    assert(Mud<=LMIconsts.mud,"Observer %d, Mu_d (%2.f) is not smaller than the desired maximum mud (%2.f)",i,Mud,LMIconsts.mud)
                    

                    OCM = [(Rs*(obj.sys.A + L*Ci) + (obj.sys.A + L*Ci)'*Rs + Nu*eye(obj.sys.nx)) (Rs*obj.sys.B + (obj.sys.C + K*Ci)'*Ds) -Rs Rs;
                           (Rs*obj.sys.B + (obj.sys.C + K*Ci)'*Ds) -2*Ds*ebarM zeros(obj.sys.nx) zeros(obj.sys.nx);
                           -Rs zeros(obj.sys.nx) -Mua*eye(obj.sys.nx) zeros(obj.sys.nx);
                           Rs zeros(obj.sys.nx) zeros(obj.sys.nx) -Mud*eye(obj.sys.nx)];
                    assert(any(eig(OCM)<0,"all"),"Observer %d does not satisfy the main LMI",i)
                end
            end
        end
    end
end