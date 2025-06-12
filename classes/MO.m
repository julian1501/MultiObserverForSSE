classdef MO
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        sys % the system that this mo observes
        attack % the attack that is used on the system
        numOutputs % the number of sensors sensing the system
        numObservers % number of observers in this multi observer
        numOutputsObserver % number of outputs that each observer has
        C % Output matrix of each observer i
        CIndices % Matrix containing the origin of each Ci
        A % Matrix containing (possibly different) A matrices for each observer
        L % Matrix containing the L matrices for each observer
        K % Matrix containing the K matrices for each observer
        mua
        ai
        bi
        ci
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
            
            [obj.C,obj.CIndices] = obj.CSetup();

            [obj.L, obj.K, obj.ai, obj.bi, obj.ci] = obj.defineObservers(LMIconsts);

            
        end

        function [C,CIndices] = CSetup(obj)
            
            CIndices = nchoosek(1:obj.numOutputs,obj.numOutputsObserver);
            C = zeros(obj.numOutputsObserver,...
                         obj.sys.nx,...
                         obj.numObservers);
            for i = 1:obj.numObservers
                C(:,:,i) = obj.sys.C(CIndices(i,:),:);
            end
        end

        function [Li, Ki, ai, bi, ci] = defineObservers(obj,LMIconsts)
            % This function defines the observers based on the LMI as in
            % section 5 of (Chong, 2025, Observer design for nonlinear 
            % systems with asynchronously sampled and maliciously corrupted
            % measurements)
            Li = zeros(obj.sys.nx,obj.numOutputsObserver,obj.numObservers);
            Ki = zeros(obj.sys.nx,obj.numOutputsObserver,obj.numObservers);
            ai = zeros(1,obj.numObservers);
            bi = zeros(1,obj.numObservers);
            ci = zeros(1,obj.numObservers);
            ebar = 1; % max slope nonlinearity
            ebarM = diag(repmat(ebar^-1,1,obj.numOutputs));

            for i = 1:1:obj.numObservers
                % Generate observer construction matrix (OCM) that has to
                % satisfy LMI: OCM =< 0
                Ci = obj.C(:,:,i);
                setlmis([])
                % in order to keep the norm of L small i define two
                % components a positive and a negative component both will
                % have an objective to remain (somewhat) small and only the
                % sum will be used
                Lh = lmivar(2,[obj.sys.nx obj.numOutputsObserver]);
                % Lnegh = lmivar(2,[obj.sys.nx obj.numOutputsObserver]);
                % Lposh = lmivar(2,[obj.sys.nx obj.numOutputsObserver]);
                maxL = lmivar(1, [1 0]);
                Kh = lmivar(2,[obj.sys.nx obj.numOutputsObserver]);

                Rsh = lmivar(1,[obj.sys.nx 1]);
                Nuh = lmivar(1,[1 0]);
                Dsh = lmivar(1,repmat([1 0],obj.numOutputs,1));
                Muah = lmivar(1,[1 0]);
                Mudh = lmivar(1,[1 0]);

                MuaLMI = newlmi; % LMIconsts.mua < Mua
                lmiterm([-MuaLMI 1 1 Muah],eye(obj.sys.nx),1);
                lmiterm([MuaLMI 1 1 0],LMIconsts.mua*eye(obj.sys.nx));

                MudLMIlo = newlmi; % 0 < Mud
                lmiterm([-MudLMIlo 1 1 Mudh],eye(obj.sys.nx),1);

                MudLMIhi = newlmi; % Mud < mud
                lmiterm([MudLMIhi 1 1 Mudh],eye(obj.sys.nx),1);
                lmiterm([-MudLMIhi 1 1 0],LMIconsts.mud*eye(obj.sys.nx));

                DsLMI = newlmi; % 0 < Ds
                lmiterm([-DsLMI 1 1 Dsh],1,1);

                RsLMI = newlmi; % 0 < Rs
                lmiterm([-RsLMI 1 1 Rsh],1,1);

                LLMI = newlmi; % L^T L < maxL I
                lmiterm([-LLMI 1 1 0],eye(obj.numOutputsObserver));
                lmiterm([-LLMI,2,1,Lh],1,1);
                lmiterm([-LLMI,2,2,maxL],eye(obj.sys.nx),1)

                % LposLMI = newlmi; % 0 < Lpos
                % lmiterm([-LposLMI,1,1,Lposh],1,1);
                % LnegLMI = newlmi; % Lneg < 0
                % lmiterm([LnegLMI 1 1 Lnegh],1,1);
                
                NuLMI = newlmi; % 0 < Nu
                lmiterm([-NuLMI 1 1 Nuh],eye(obj.sys.nx),1);
                lmiterm([NuLMI 1 1 0],LMIconsts.nu*eye(obj.sys.nx));
                
                LKLMI = newlmi; % Whole big thing
                lmiterm([LKLMI 1 1 Rsh],1,obj.A,'s');
                lmiterm([LKLMI 1 1 Lh],Rsh,Ci,'s');
                lmiterm([LKLMI 1 1 Nuh],eye(obj.sys.nx),1);

                lmiterm([LKLMI 1 2 Rsh],1,obj.sys.B);
                lmiterm([LKLMI 1 2 Dsh],obj.sys.C',1);
                lmiterm([LKLMI 1 2 -Kh],Ci',Dsh); % use the previously defined Kht

                lmiterm([LKLMI 1 3 Rsh],-1,1);
                lmiterm([LKLMI 1 4 Rsh],1,1);

                lmiterm([LKLMI 2 2 Dsh],-2,ebarM);
                lmiterm([LKLMI 3 3 Muah],-1,eye(obj.sys.nx));
                lmiterm([LKLMI 4 4 Mudh],-1,eye(obj.sys.nx));
                
                % construct lmi and solve
                LMISYS = getlmis;

                NuInfo = decinfo(LMISYS, Nuh);  % Get index info for Nu
                MuaInfo = decinfo(LMISYS, Muah);
                numDecVars = decnbr(LMISYS);
                c = zeros(1, numDecVars);
                c(NuInfo(1):NuInfo(end)) = -10;  % Maximize Nu
                c(MuaInfo(1):MuaInfo(end)) = 1; % Minimize Mua

                Rs_info = decinfo(LMISYS, Rsh);
                c(Rs_info(1):Rs_info(end)) = 0.1; %  minimize the trace of Rs

                maxLInfo = decinfo(LMISYS, maxL);
                c(maxLInfo(1):maxLInfo(end)) = 1; % minimize the L

                % [tmin,xfeas] = feasp(LMISYS,[0,0,0,0,1]);
                [xopt, fopt] = mincx(LMISYS, c);


                assert(~isempty(xopt), "Observer %d is not observable. Optimization failed.", i);
                
                LVal = dec2mat(LMISYS,fopt,Lh);
                Kval = dec2mat(LMISYS,fopt,Kh);
                Li(:,:,i) = LVal;
                Ki(:,:,i) = Kval;

                RsVal  = dec2mat(LMISYS,fopt,Rsh);
                assert(any(eig(RsVal)>0),"Observer %d, Rs eigenvalues not larger than zero",i)
                assert(any(RsVal == RsVal',"all"),"Observer %d, Rs is not symmetric",i)
                NuVal  = dec2mat(LMISYS,fopt,Nuh);
                assert(NuVal>=LMIconsts.nu,"Observer %d, Nu (%.2f) is not larger than the desired size nu (%.2f)",i,NuVal,LMIconsts.nu)
                DsVal  = dec2mat(LMISYS,fopt,Dsh);
                assert(any(eig(DsVal)>0),"Observer %d, Ds is not positive definite",i)
                MuaVal = dec2mat(LMISYS,fopt,Muah);
                assert(MuaVal>=LMIconsts.mua,"Observer %d, Mu_a (%.2f) is not larger than the desired size mua (%.2f)",i,MuaVal,LMIconsts.mua)
                MudVal = dec2mat(LMISYS,fopt,Mudh);
                assert(MudVal>=0,"Observer %d, Mu_d (%2.f) is not larger than zero",i,MudVal)
                assert(MudVal<=LMIconsts.mud,"Observer %d, Mu_d (%2.f) is not smaller than the desired maximum mud (%2.f)",i,MudVal,LMIconsts.mud)
                

                OCM = [(RsVal*(obj.sys.A + LVal*Ci) + (obj.sys.A + LVal*Ci)'*RsVal + NuVal*eye(obj.sys.nx)) (RsVal*obj.sys.B + (obj.sys.C + Kval*Ci)'*DsVal) -RsVal RsVal;
                       (RsVal*obj.sys.B + (obj.sys.C + Kval*Ci)'*DsVal)' -2*DsVal*ebarM zeros(obj.sys.nx) zeros(obj.sys.nx);
                       -RsVal zeros(obj.sys.nx) -MuaVal*eye(obj.sys.nx) zeros(obj.sys.nx);
                       RsVal zeros(obj.sys.nx) zeros(obj.sys.nx) -MudVal*eye(obj.sys.nx)];
                assert(any(eig(OCM)<0,"all"),"Observer %d does not satisfy the main LMI",i)
                
                eigRs = eig(RsVal);
                lambdaMinRs = min(eigRs);
                ai(i) = lambdaMinRs;
                bi(i) = NuVal/lambdaMinRs;
                ci(i) = MuaVal*(norm(LVal)^2 + norm(obj.sys.B)^2*ebar^2);

            end
        end
    end
end