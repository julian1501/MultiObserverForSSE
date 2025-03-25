classdef PowerSystem
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        A % State matrix
        B % Input matrix
        C % Output matrix
        nx
        ny
        nu
        Vref % reference voltage
        actImp
        reaImp
        charInverters
    end

    methods
        function obj = PowerSystem(numCustomers,sysConsts)
            % PowerSystem Construct an instance of this class
            obj.A = diag(-1.*sysConsts.inverters);
            obj.B = diag(sysConsts.inverters);
            obj.actImp = sysConsts.actImp;
            obj.reaImp = sysConsts.reaImp;
            obj.Vref = sysConsts.vref;
            obj.charInverters = sysConsts.charInverters;

            reaImpMainLine   =   obj.reaImp(:,1);
            reaImpBranches   =   obj.reaImp(:,2);
            
            % Generate C matrix template
            obj.C = zeros(numCustomers,numCustomers);
            for i = 1:1:numCustomers
                obj.C(i:end,i:end) = obj.C(i:end,i:end) + repmat(reaImpMainLine(i),numCustomers+1-i,numCustomers+1-i);
            end
            obj.C = obj.C + diag(reaImpBranches);
            obj.C = -2.*obj.C;

            obj.nx = size(obj.A,1);
            obj.ny = size(obj.C,1);
            obj.nu = size(obj.B,2);

            if obj.isStable == 0
                fprintf("The system is not stable. \n")
            else
                fprintf("The system is stable. \n")
            end
            
            if obj.isObsv == 0
                error("The system (A,C) is not observable. \n")
            else
                fprintf("The system (A,C) is observable. \n")
            end

            if obj.isCtrb == 0
                error("The system (A,B) is not controlable. \n")
            else
                fprintf("The system (A,C) is controlable. \n")
            end
        end

        function stableBool = isStable(obj)
            % isStable
            % Check if system is stable
            eigenvalues = eig(obj.A);
            if any(eigenvalues > 0)
                stableBool = false;
            else
                stableBool = true;
            end
        end

        function obsvBool = isObsv(obj)
            % isObsv
            % Check if system is observable through PBH test
            ev = eig(obj.A);
            I = eye(obj.nx);
            % PBH test
            obsvBool = true;
            for i = 1:1:obj.nx
                if rank([ev(i)*I - obj.A; obj.C]) ~= obj.nx
                    obsvBool = false;
                    break
                end
            end
        end

        function ctrbBool = isCtrb(obj)
            % isObsv
            % Check if system is observable through PBH test
            ev = eig(obj.A);
            I = eye(obj.nx);
            % PBH test
            ctrbBool = true;
            for i = 1:1:obj.nx
                if rank([ev(i)*I - obj.A, obj.B]) ~= obj.nx
                    ctrbBool = false;
                    break
                end
            end
        end
    end
end